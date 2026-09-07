#ifndef MITM_ROUTER_DRIVER
#define MITM_ROUTER_DRIVER

#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <err.h>

#include "router/router.hpp"

/*
 * The command line shared by router_test and router_bench: the team, the traffic, every Router_Opts
 * field.  Part of the examples, not of the library.
 */
namespace mitm {

struct RouterArgs {
	int senders = 4;                 /* per node */
	int receivers = 2;               /* per node */
	u64 points = 100000;             /* per sender and per round (test) */
	int rounds = 3;
	bool lossy = false;
	Router_Opts opts;
	int pop_max = 0;                 /* points per Router_Pop; 0 == the test's own schedule */
	int skew = 0;                    /* hot destinations; 0 == uniform */
	bool stagger = false;            /* sender s pushes (s+1)/S of the points */
	bool partial = false;            /* 1..7 points per destination, round robin */
	bool slow_recv = false;          /* a receiver dawdles after each pop */
	u64 seed = 1;                    /* bench: the senders' PRNG streams */
	double seconds = 3;              /* bench: duration of a round */
	bool local_only = false;         /* bench: every destination on the sender's node */
	std::string test = "all";
};

static void router_usage()
{
	printf("Router drivers.  Options:\n");
	printf("  --senders N --receivers N        threads per node (default 4 and 2)\n");
	printf("  --points N --rounds N            per sender per round; rounds\n");
	printf("  --lossy                          drop instead of waiting\n");
	printf("  --block N --swc N --n-recv N --inbox BLOCKS --credit N\n");
	printf("  --sweep N                        sealed blocks handled per turn\n");
	printf("  --dests N                        destinations per node: a virtual fan-out (default: the receivers)\n");
	printf("  --pop-max N --skew K --stagger --partial --slow-recv --seed S\n");
	printf("  --no-bind --cache-level L --group N     pinning (on by default); cores per group (default 16)\n");
	printf("  --seconds X --local-only               (bench)\n");
	printf("  --test NAME                      (test) connect_only basic tiny rounds zero skew partial_lines\n");
	printf("                                   pop_sizes swc_is_block staggered_close slow_recv groups colors\n");
	printf("                                   asymmetric all\n");
	printf("  --quiet\n");
	exit(EXIT_SUCCESS);
}

static void router_parse(int argc, char **argv, RouterArgs &a)
{
	enum {
		O_SENDERS = 1000, O_RECEIVERS, O_POINTS, O_ROUNDS, O_LOSSY, O_BLOCK, O_SWC, O_NRECV, O_INBOX,
		O_SWEEP, O_DESTS, O_CREDIT, O_POPMAX, O_SKEW, O_STAGGER, O_PARTIAL, O_SLOW,
		O_SEED, O_SECONDS, O_NOBIND, O_CACHE, O_GROUP, O_LOCAL, O_TEST, O_QUIET, O_HELP
	};
	struct option longopts[] = {
		{"senders", required_argument, NULL, O_SENDERS}, {"receivers", required_argument, NULL, O_RECEIVERS},
		{"points", required_argument, NULL, O_POINTS}, {"rounds", required_argument, NULL, O_ROUNDS},
		{"lossy", no_argument, NULL, O_LOSSY}, {"block", required_argument, NULL, O_BLOCK},
		{"swc", required_argument, NULL, O_SWC},
		{"n-recv", required_argument, NULL, O_NRECV}, {"inbox", required_argument, NULL, O_INBOX},
		{"sweep", required_argument, NULL, O_SWEEP},
		{"dests", required_argument, NULL, O_DESTS}, {"credit", required_argument, NULL, O_CREDIT},
		{"pop-max", required_argument, NULL, O_POPMAX}, {"skew", required_argument, NULL, O_SKEW},
		{"stagger", no_argument, NULL, O_STAGGER}, {"partial", no_argument, NULL, O_PARTIAL},
		{"slow-recv", no_argument, NULL, O_SLOW}, {"seed", required_argument, NULL, O_SEED},
		{"seconds", required_argument, NULL, O_SECONDS},
		{"no-bind", no_argument, NULL, O_NOBIND}, {"cache-level", required_argument, NULL, O_CACHE},
		{"group", required_argument, NULL, O_GROUP},
		{"local-only", no_argument, NULL, O_LOCAL}, {"test", required_argument, NULL, O_TEST},
		{"quiet", no_argument, NULL, O_QUIET}, {"help", no_argument, NULL, O_HELP},
		{NULL, 0, NULL, 0}
	};
	for (;;) {
		int ch = getopt_long(argc, argv, "", longopts, NULL);
		switch (ch) {
		case -1: return;
		case O_SENDERS: a.senders = atoi(optarg); break;
		case O_RECEIVERS: a.receivers = atoi(optarg); break;
		case O_POINTS: a.points = human_parse(optarg); break;
		case O_ROUNDS: a.rounds = atoi(optarg); break;
		case O_LOSSY: a.lossy = true; break;
		case O_BLOCK: a.opts.block_points = human_parse(optarg); break;
		case O_SWC: a.opts.swc_linesize = human_parse(optarg); break;
		case O_NRECV: a.opts.n_recv = atoi(optarg); break;
		case O_INBOX: a.opts.inbox_blocks = atoi(optarg); break;
		case O_SWEEP: a.opts.sweep_blocks = atoi(optarg); break;
		case O_DESTS: a.opts.dests_per_node = atoi(optarg); break;
		case O_CREDIT: a.opts.credit = atoi(optarg); break;
		case O_POPMAX: a.pop_max = atoi(optarg); break;
		case O_SKEW: a.skew = atoi(optarg); break;
		case O_STAGGER: a.stagger = true; break;
		case O_PARTIAL: a.partial = true; break;
		case O_SLOW: a.slow_recv = true; break;
		case O_SEED: a.seed = strtoull(optarg, NULL, 0); break;
		case O_SECONDS: a.seconds = atof(optarg); break;
		case O_NOBIND: a.opts.pin = false; break;
		case O_CACHE: a.opts.cache_level = atoi(optarg); break;
		case O_GROUP: a.opts.group_size = atoi(optarg); break;
		case O_LOCAL: a.local_only = true; break;
		case O_TEST: a.test = optarg; break;
		case O_QUIET: a.opts.verbose = false; break;
		case O_HELP: router_usage(); break;
		default: errx(1, "Unknown option (try --help)");
		}
	}
}

}
#endif
