#ifndef MITM_PCS
#define MITM_PCS

/*
 * PCS, the parallel collision search, on the Router.  One include for a driver: the engine pulls in its
 * parameters, its shared state, the trails, the dictionary, the two worker threads and the control
 * channel.  Entry points: mitm::pcs::claw_search and mitm::pcs::collision_search.
 */

#include "pcs/params.hpp"
#include "pcs/shared.hpp"
#include "pcs/wrappers.hpp"
#include "pcs/trail.hpp"
#include "pcs/dict.hpp"
#include "pcs/walker.hpp"
#include "pcs/control.hpp"
#include "pcs/engine.hpp"

#endif
