#ifndef MITM_DIRECT
#define MITM_DIRECT

/*
 * The direct engine: the exhaustive meet-in-the-middle over a distributed dictionary, on the Router.  One
 * include for a driver: the engine pulls in its parameters, its shared state, the problem wrappers, the
 * dictionary, the producer and the round/epilogue machinery.  Entry points: mitm::direct::claw_search and
 * mitm::direct::collision_search.
 */

#include "direct/params.hpp"
#include "direct/shared.hpp"
#include "direct/wrappers.hpp"
#include "direct/dict.hpp"
#include "direct/producer.hpp"
#include "direct/engine.hpp"

#endif
