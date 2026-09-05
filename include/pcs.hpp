#ifndef MITM_PCS
#define MITM_PCS

#include <cassert>

#include "tools.hpp"
#include "problem.hpp"
#include "parameters.hpp"
#include "engine.hpp"
#include "pcs_common.hpp"
#include "walker.hpp"
#include "inserter.hpp"

/*
 * The PCS scheme: van Oorschot-Wiener parallel collision search over the distributed dictionary
 * (PROTOCOL.md §7).  The problem wrappers, the Scheme's functions and the two entry points,
 * claw_search / collision_search.  A driver includes this one, plus its own problem definition.
 */

namespace mitm::pcs {

/*
 * A collision problem as one random function: x -> f(i ^ x).  Needs |Range| >= |Domain|.
 * INCOMPLETE (CLAUDE.md, Gotchas): for some seeds the search plateaus and never retires the golden
 * pair.  Suspected: mix() is a plain xor, so every version i has the same collision set, and
 * mix_good_pair() does not normalise the pair's order the way swapmix() does.
 */
template <class Problem>
class CollisionWrapper {
public:
    const Problem &pb;
    const int n;                    /* domain bits */
    const int m;                    /* range bits */
    const u64 in_mask;              /* the low n bits */
    const u64 out_mask;             /* the low m bits */
    static constexpr int vlen = Problem::vlen;


    CollisionWrapper(const Problem &pb) : pb(pb), n(pb.n), m(pb.m), in_mask(make_mask(pb.n)), out_mask(make_mask(pb.m))
    {
        static_assert(std::is_base_of<AbstractCollisionProblem, Problem>::value,
            "problem not derived from mitm::AbstractCollisionProblem");
        assert(m <= 64);

        /* self-test: vmixf agrees with mixf */
        PRNG vprng;
        u64 i = vprng.rand() & out_mask;
        u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        for (int j = 0; j < vlen; j++)
            x[j] = vprng.rand() & out_mask;
        vmixf(i, x, y);
        for (int j = 0; j < vlen; j++)
            assert(y[j] == mixf(i, x[j]));
    }

    /* randomization by a family of permutations of {0, 1}^n */
    u64 mix(u64 i, u64 x) const   /* return σ_i(x) */
    {
        return i ^ x;
    }

    /* evaluates f o σ_i(x) */
    u64 mixf(u64 i, u64 x) const
    {
        return pb.f(mix(i, x));
    }

    void vmixf(u64 i, u64 x[], u64 r[]) const
    {
        /* careful: vlen can be more than one SIMD vector */
        u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        for (int j = 0; j < vlen; j++)
            y[j] = mix(i, x[j]);
        /* a scalar problem provides no vf(): call f() straight away */
        if constexpr (vlen == 1)
            r[0] = pb.f(y[0]);
        else
            pb.vf(y, r);
    }

    bool mix_good_pair(u64 i, u64 x0, u64 x1) const
    {
        return pb.is_good_pair(mix(i, x0), mix(i, x1));
    }

    /* one value every rank must agree on (PROTOCOL.md §4.1): the function has no uninitialised state */
    u64 self_test(u64 a, u64 b) const
    {
        return mixf(a & out_mask, b & out_mask);
    }
};


/****************************************************************************************/

template <class Problem>
class EqualSizeClawWrapper {
public:
    const Problem &pb;
    const int n;                    /* domain bits */
    const int m;                    /* range bits */
    const u64 in_mask;              /* the low n bits */
    const u64 out_mask;             /* the low m bits */
    const u64 choice_mask;          /* the bit of the mixed value that picks f or g */
    static constexpr int vlen = Problem::vlen;

    EqualSizeClawWrapper(const Problem& pb)
        : pb(pb), n(pb.n), m(pb.m), in_mask(make_mask(pb.n)), out_mask(make_mask(pb.m)), choice_mask(1ull << (pb.m - 1))
    {
        static_assert(std::is_base_of<AbstractClawProblem, Problem>::value,
            "problem not derived from mitm::AbstractClawProblem");
        assert(m <= 64);
        assert(pb.n == pb.m);

        /* self-test: vmixf agrees with mixf */
        PRNG vprng;
        u64 i = vprng.rand() & out_mask;
        u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        for (int j = 0; j < vlen; j++)
            x[j] = vprng.rand() & out_mask;
        vmixf(i, x, y);
        for (int j = 0; j < vlen; j++)
            assert(y[j] == mixf(i, x[j]));
    }

    /* pick either f() or g() */
    bool choose(u64 i, u64 x) const
    {
        return (x * (i | 1)) & choice_mask;
    }

    u64 mix(u64 i, u64 x) const
    {
        return i ^ x;
    }

    u64 mixf(u64 i, u64 x) const
    {
        u64 y = mix(i, x);
        if (choose(i, x))
            return pb.f(y);
        else
            return pb.g(y);
    }

    void vmixf(u64 i, u64 x[], u64 r[]) const
    {
        /* careful: vlen can be more than one SIMD vector */
        u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        bool choices[vlen];
        for (int j = 0; j < vlen; j++) {
            y[j] = mix(i, x[j]);
            choices[j] = choose(i, x[j]);
        }
        /* a scalar problem provides no vfg(): call f() / g() straight away */
        if constexpr (vlen == 1)
            r[0] = choices[0] ? pb.f(y[0]) : pb.g(y[0]);
        else
            pb.vfg(y, choices, r);
    }

    pair<u64, u64> swapmix(u64 i, u64 a, u64 b) const
    {
        u64 x0 = choose(i, a) ? a : b;
        u64 x1 = choose(i, b) ? a : b;
        assert(choose(i, x0));
        assert(not choose(i, x1));
        return pair(mix(i, x0), mix(i, x1));
    }

    bool mix_good_pair(u64 i, u64 a, u64 b) const
    {
        if (choose(i, a) == choose(i, b))
            return false;
        auto [x0, x1] = swapmix(i, a, b);
        return pb.is_good_pair(x0, x1);
    }

    /* one value every rank must agree on (PROTOCOL.md §4.1): the functions have no uninitialised state */
    u64 self_test(u64 a, u64 b) const
    {
        return mixf(a & out_mask, b & out_mask);
    }
};

template <class Problem>
class LargerRangeClawWrapper {
public:
    const Problem &pb;
    const int n;                    /* domain bits: pb.n + 1, the choice bit included */
    const int m;                    /* range bits */
    const u64 in_mask;              /* the low pb.n bits */
    const u64 out_mask;             /* the low m bits */
    static constexpr int vlen = Problem::vlen;
    u64 choice_mask;                /* the bit of the mixed value that picks f or g */

    LargerRangeClawWrapper(const Problem& pb)
        : pb(pb), n(pb.n + 1), m(pb.m), in_mask(make_mask(pb.n)), out_mask(make_mask(pb.m))
    {
        static_assert(std::is_base_of<AbstractClawProblem, Problem>::value,
            "problem not derived from mitm::AbstractClawProblem");
        assert(m <= 64);
        assert(n <= m);
        choice_mask = 1ull << n;

        /* self-test: vmixf agrees with mixf */
        PRNG vprng;
        u64 i = vprng.rand() & out_mask;
        u64 x[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        for (int j = 0; j < vlen; j++)
            x[j] = vprng.rand() & out_mask;
        vmixf(i, x, y);
        for (int j = 0; j < vlen; j++)
            assert(y[j] == mixf(i, x[j]));
    }

    inline u64 full_mix(u64 i, u64 x) const
    {
        u64 y = x * (i | 1);
        return (y ^ (y >> n));
    }

    /* pick either f() or g() */
    bool choose(u64 i, u64 x) const
    {
        return full_mix(i, x) & choice_mask;
    }

    u64 mix(u64 i, u64 x) const   // {0, 1}^m  x  {0, 1}^m ---> {0, 1}^n
    {
        return full_mix(i, x) & in_mask;
    }

    u64 mixf(u64 i, u64 x) const  // {0, 1}^m  x  {0, 1}^m ---> {0, 1}^m
    {
        u64 z = full_mix(i, x);
        if (z & choice_mask)
            return pb.f(z & in_mask);
        else
            return pb.g(z & in_mask);
    }

    void vmixf(u64 i, u64 x[], u64 r[]) const
    {
        /* careful: vlen can be more than one SIMD vector */
        u64 y[vlen] __attribute__ ((aligned(sizeof(u64) * vlen)));
        bool choice[vlen];
        for (int j = 0; j < vlen; j++) {
            y[j] = mix(i, x[j]);
            choice[j] = choose(i, x[j]);
        }
        /* a scalar problem provides no vfg(): call f() / g() straight away */
        if constexpr (vlen == 1)
            r[0] = choice[0] ? pb.f(y[0]) : pb.g(y[0]);
        else
            pb.vfg(y, choice, r);
    }

    pair<u64, u64> swapmix(u64 i, u64 a, u64 b) const
    {
        u64 x0 = choose(i, a) ? a : b;
        u64 x1 = choose(i, b) ? a : b;
        assert(choose(i, x0));
        assert(not choose(i, x1));
        return pair(mix(i, x0), mix(i, x1));
    }

    bool mix_good_pair(u64 i, u64 a, u64 b) const
    {
        if (choose(i, a) == choose(i, b))
            return false;
        auto [x0, x1] = swapmix(i, a, b);
        return pb.is_good_pair(x0, x1);
    }

    /* one value every rank must agree on (PROTOCOL.md §4.1): the functions have no uninitialised state */
    u64 self_test(u64 a, u64 b) const
    {
        return mixf(a & out_mask, b & out_mask);
    }
};


/******************************** the scheme's functions ********************************/

/* the version and the chains' root are drawn afresh; `stop` stays the controller's (PROTOCOL.md §4.2) */
template <class Wrapper>
void Scheme::next_header(const Params &, const Wrapper &wrapper, PRNG &prng, Header &h, u64 &)
{
    h.i = prng.rand() & wrapper.out_mask;
    h.root_seed = prng.rand();
}

inline void Scheme::build_dict(SharedContext<Scheme> &shared, const Params &params, int d)
{
    shared.shards[d] = std::make_unique<PcsDict>(params.jbits, params.w_shard);
    shared.scheme.coll_q[d] = std::make_unique<CollisionQueue>(params.coll_queue_capacity);
}

/* every round starts from an empty dictionary (PROTOCOL.md §5); a collective flush is future work */
inline void Scheme::after_round(SharedContext<Scheme> &shared, const Params &, int d, const Header &)
{
    shared.shards[d]->flush();
}

/* a version of the function is used for beta*w distinguished points, then dropped (PROTOCOL.md §2.2) */
inline bool Scheme::round_complete(const Params &params, const u64 reported[])
{
    return reported[N_DP] >= params.points_per_version;
}

/* the startup report, printed by run() before the team exists: the plan is on record even if the allocation fails */
inline void Scheme::banner(const Params &params, u64 seed)
{
    printf("Starting MPI+OpenMP parallel collision search (PCS) with seed=%016" PRIx64 "\n", seed);
    layout_banner(params);
    printf("Generating %.1f*w = %" PRIu64 " = 2^%0.2f distinguished points / version\n",
        params.beta, params.points_per_version, std::log2((double) params.points_per_version));
    printf("DP == %d words: endpoint + (%d-bit length | %d-bit chain index).  ",
        POINT_WORDS, params.lenbits, params.jbits);
    if (params.dp_max_it > params.len_sat)
        printf("Lengths >= %" PRIu64 " re-walked\n", params.len_sat);
    else
        printf("No length can saturate (dp_max_it == %" PRIu64 ")\n", params.dp_max_it);
    if (params.theta_auto)
        printf("AUTO-TUNING: setting 1/theta == %.2f\n", 1 / params.theta);
    else
        printf("NOTICE: using 1/theta == %.2f vs ``optimal'' 1/theta == %.2f\n",
            1 / params.theta, 1 / params.auto_theta);
    if (params.theta == 1) {
        printf("***** WARNING *****\n***** WARNING *****\n***** WARNING *****\n");
        printf("---> zero difficulty (use the direct scheme!)\n");
        printf("***** WARNING *****\n***** WARNING *****\n***** WARNING *****\n");
    }
    fflush(stdout);
}

/* the live one-line display, from the reports so far */
inline void Scheme::display(const Params &params, const u64 reported[], double delta, u64 nround)
{
    u64 ndp = reported[N_DP];
    double dp_rate = ndp / delta;
    double completion = (double) ndp / params.points_per_version;
    char hrate[8], hnrate[8], hprobe[8];
    human_format(dp_rate / params.theta / params.n_producers, hrate);
    human_format(ndp * POINT_WORDS * sizeof(u64) / params.n_nodes / delta, hnrate);
    human_format((double) reported[N_PROBE] / params.n_dicts / delta, hprobe);
    printf("\rRound %" PRId64 ":  %.1fs (%.1f%%, ETA %.1fs).  %.2f*w #DP.  %s #f/s per producer.  "
           "%s probe/s per dict thread.  node-->%sB/s   ",
        nround, delta, 100. * completion,
        (completion > 0) ? delta * (1 - completion) / completion : 0.,
        (double) ndp / params.w, hrate, hprobe, hnrate);
    fflush(stdout);
}

/*
 * The round report: `r` is the exact SUM of the counters over every thread of every node, `total` the
 * all-time sums; `round` and `all` the HyperLogLogs, this round's and every round's merged.
 */
inline void Scheme::round_report(const Params &params, const u64 r[], const u64 total[], const RoundStats &round,
                                 const RoundStats &all, double delta, u64 nround)
{
    u64 ndp = r[N_DP];
    char hrate[8], hnrate[8];
    human_format((double) r[N_EVAL] / params.n_producers / delta, hrate);
    human_format((double) ndp * POINT_WORDS * sizeof(u64) / params.n_nodes / delta, hnrate);

    printf("\n");
    printf("Round %" PRId64 ".  %.1fs.  #DP %.2f*w (total 2^%.2f).  #coll %.2f*w (total 2^%.2f).  "
           "Total #f=2^%.3f.  %s #f/s per producer.  node-->%sB/s\n",
        nround, delta,
        (double) ndp / params.w, std::log2((double) total[N_DP] ? (double) total[N_DP] : 1.),
        (double) r[N_COLLISIONS] / params.w,
        std::log2((double) total[N_COLLISIONS] ? (double) total[N_COLLISIONS] : 1.),
        std::log2((double) total[N_EVAL] ? (double) total[N_EVAL] : 1.), hrate, hnrate);

    if (ndp > 0) {
        double avglen = (double) r[N_POINTS_TRAILS] / ndp;
        printf("            %.2f avg trail length", avglen);
        if (r[N_COLLISIONS] > 0 && avglen > 0)
            printf(" (x%.2f & x%.2f colliding)",
                (double) r[COLLIDING_LEN_MIN] / r[N_COLLISIONS] / avglen,
                (double) r[COLLIDING_LEN_MAX] / r[N_COLLISIONS] / avglen);
        /* BAD_PROBE is a dict-thread counter: its denominator is the probes retired, not the DPs
           found.  The other four are producer-side, where ndp is right. */
        printf(".  %.2f%% probe failure.  %.2f%% walk-robinhood.  %.2f%% walk-noncolliding.  "
               "%.2f%% same-value.  %.2f%% DP failure.  %.2f re-walked/collision\n",
            r[N_PROBE] ? 100. * r[BAD_PROBE] / r[N_PROBE] : 0., 100. * r[BAD_WALK_ROBINHOOD] / ndp,
            100. * r[BAD_WALK_NONCOLLIDING] / ndp, 100. * r[BAD_COLLISION] / ndp,
            100. * r[BAD_DP] / ndp,
            r[N_COLLISIONS] ? (double) r[N_MEASURE] / r[N_COLLISIONS] : 0.);

        /*
         * What the round actually routed.  A round closes on points *found* (against
         * points_per_version), so a run whose queues overflow burns rounds against a nearly empty
         * dictionary and still reports them complete: the DPs that reached a shard, and the share of
         * those found, is the number to read.  PROBLEM.md §2: it is the attack's speed, one for one.
         */
        char hoff[8], hins[8], hprobe[8];
        human_format((double) ndp / delta, hoff);
        human_format((double) r[N_PROBE] / delta, hins);
        human_format((double) r[N_PROBE] / params.n_dicts / delta, hprobe);
        printf("            ROUTED  %s DP/s found --> %s DP/s inserted (%.2f%% reached a shard,"
               " %s probe/s per dict thread).  dict load %.2f/slot\n",
            hoff, hins, 100. * r[N_PROBE] / ndp, hprobe,
            (double) r[N_PROBE] / params.w);
    }

    if (r[DROP_PRODUCERQ] | r[DROP_OUT] | r[DROP_DICTQ] | r[DROP_COLL])
        printf("            DROPPED  %" PRId64 " producer-queue / %" PRId64 " output-buffer / "
               "%" PRId64 " dict-queue / %" PRId64 " collision-queue\n",
            r[DROP_PRODUCERQ], r[DROP_OUT], r[DROP_DICTQ], r[DROP_COLL]);

    u64 E_i = distinct_collisions_estimation(round.hll);
    u64 E = distinct_collisions_estimation(all.hll);
    printf("            #distinct coll (this i / total) %.02f*w / 2^%.2f\n",
        (double) E_i / params.w, std::log2((double) E ? (double) E : 1.));
    printf("\n");
    fflush(stdout);
}

/* the last line: found, or gave up after nround versions */
inline void Scheme::done(const Params &, u64 nround, bool found, double seconds)
{
    if (found)
        printf("Completed in %.2fs\n", seconds);
    else
        printf("Gave up after %" PRId64 " versions of the function (%.2fs)\n", nround, seconds);
    fflush(stdout);
}


/******************************** the entry points ********************************/

/* find x0 != x1 with f(x0) == f(x1) and is_good_pair(x0, x1) */
template <class Problem>
optional<pair<u64, u64>> collision_search(const Problem& pb, u64 nbytes_memory, const Options &opts, PRNG &prng)
{
    int rank;
    MPI_Comm_rank(opts.mpi_comm, &rank);
    if (opts.verbose && rank == 0)
        printf("Starting collision search with f : {0,1}^%d --> {0, 1}^%d (vlen=%d)\n",
            pb.n, pb.m, Problem::vlen);

    CollisionWrapper<Problem> wrapper(pb);
    auto collision = run<Scheme>(wrapper, nbytes_memory, opts, prng);
    if (not collision)
        return nullopt;

    auto [i, x, y] = *collision;
    u64 x0 = wrapper.mix(i, x);
    u64 x1 = wrapper.mix(i, y);

    assert(x0 != x1);
    assert(pb.f(x0) == pb.f(x1));
    assert(pb.is_good_pair(x0, x1));
    return optional(pair(x0, x1));
}


/* find x0, x1 with f(x0) == g(x1) and is_good_pair(x0, x1) */
template <class Problem>
optional<pair<u64, u64>> claw_search(const Problem& pb, u64 nbytes_memory, const Options &opts, PRNG &prng)
{
    int rank;
    MPI_Comm_rank(opts.mpi_comm, &rank);
    bool verbose = opts.verbose && (rank == 0);
    if (verbose)
        printf("Starting claw search with f, g : {0,1}^%d --> {0, 1}^%d (vlen=%d)\n",
            pb.n, pb.m, Problem::vlen);

    optional<tuple<u64,u64,u64>> claw;
    u64 x0, x1;

    /* distinct wrapper types, so each branch carries the whole search */
    if (pb.n == pb.m) {
        if (verbose)
            printf("  - using |Domain| == |Range| mode.  Expecting 1.8*n/w rounds.\n");
        EqualSizeClawWrapper<Problem> wrapper(pb);
        claw = run<Scheme>(wrapper, nbytes_memory, opts, prng);
        if (claw) {
            auto [i, a, b] = *claw;
            std::tie(x0, x1) = wrapper.swapmix(i, a, b);
        }
    } else if (pb.n < pb.m) {
        if (verbose)
            printf("  - using |Domain| << |Range| mode.  Expecting 0.9*n/w rounds.\n");
        LargerRangeClawWrapper<Problem> wrapper(pb);
        claw = run<Scheme>(wrapper, nbytes_memory, opts, prng);
        if (claw) {
            auto [i, a, b] = *claw;
            std::tie(x0, x1) = wrapper.swapmix(i, a, b);
        }
    } else {
        errx(1, "Larger domain not yet supported...");
    }

    if (not claw)
        return nullopt;

    assert((x0 & make_mask(pb.n)) == x0);
    assert((x1 & make_mask(pb.n)) == x1);
    assert(pb.f(x0) == pb.g(x1));
    assert(pb.is_good_pair(x0, x1));
    return optional(pair(x0, x1));
}

}
#endif
