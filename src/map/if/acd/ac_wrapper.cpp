/**C++File**************************************************************

  FileName    [ac_wrapper.cpp]

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [Ashenhurst-Curtis decomposition.]

  Synopsis    [Interface with the FPGA mapping package.]

  Author      [Alessandro Tempia Calvino]
  
  Affiliation [EPFL]

  Date        [Ver. 1.0. Started - November 20, 2023.]

***********************************************************************/

#include "ac_wrapper.h"
#include "ac_decomposition.hpp"
#include "acd66.hpp"
#include "acdXX.hpp"
#include "ac_dc_decomposition.hpp"

ABC_NAMESPACE_IMPL_START

int acd_evaluate( word * pTruth, unsigned nVars, int lutSize, unsigned *pdelay, unsigned *cost, int try_no_late_arrival )
{
  using namespace acd;

  ac_decomposition_params ps;
  ps.lut_size = lutSize;
  ps.use_first = false;
  ps.try_no_late_arrival = static_cast<bool>( try_no_late_arrival );
  ac_decomposition_stats st;

  ac_decomposition_impl acd( nVars, ps, &st );
  int val = acd.run( pTruth, *pdelay );

  if ( val < 0 )
  {
    *pdelay = 0;
    return -1;
  }

  *pdelay = acd.get_profile();
  *cost = st.num_luts;

  return val;
}

int acd_decompose( word * pTruth, unsigned nVars, int lutSize, unsigned *pdelay, unsigned char *decomposition )
{
  using namespace acd;

  ac_decomposition_params ps;
  ps.lut_size = lutSize;
  ps.use_first = true;
  ac_decomposition_stats st;

  ac_decomposition_impl acd( nVars, ps, &st );
  acd.run( pTruth, *pdelay );
  int val = acd.compute_decomposition();

  if ( val < 0 )
  {
    *pdelay = 0;
    return -1;
  }

  *pdelay = acd.get_profile();
  acd.get_decomposition( decomposition );
  return 0;
}

int acd2_evaluate( word * pTruth, unsigned nVars, int lutSize, unsigned *pdelay, unsigned *cost, int try_no_late_arrival )
{
  using namespace acd;

  acdXX_params ps;
  ps.lut_size = lutSize;
  ps.max_shared_vars = lutSize - 2;
  acdXX_impl acd( nVars, ps );
  int val = acd.run( pTruth, *pdelay );

  if ( val == 0 )
  {
    *pdelay = 0;
    return -1;
  }

  acd.compute_decomposition();
  *pdelay = acd.get_profile();
  *cost = 2;

  return val;
}

int acd2_decompose( word * pTruth, unsigned nVars, int lutSize, unsigned *pdelay, unsigned char *decomposition )
{
  using namespace acd;

  acdXX_params ps;
  ps.lut_size = lutSize;
  ps.max_shared_vars = lutSize - 2;
  acdXX_impl acd( nVars, ps );
  acd.run( pTruth, *pdelay );
  int val = acd.compute_decomposition();

  if ( val != 0 )
  {
    *pdelay = 0;
    return -1;
  }

  *pdelay = acd.get_profile();

  acd.get_decomposition( decomposition );
  return 0;
}

inline int acd66_evaluate( word * pTruth, unsigned nVars )
{
  using namespace acd;

  acd66_impl acd( nVars, true, false );

  if ( acd.run( pTruth ) == 0 )
    return 0;

  return 1;
}

int acdXX_evaluate( word * pTruth, unsigned lutSize, unsigned nVars )
{
  using namespace acd;

  if ( lutSize == 6 )
  {
    return acd66_evaluate( pTruth, nVars );
  }
  
  acdXX_params ps;
  ps.lut_size = lutSize;
  ps.max_shared_vars = lutSize - 2;
  acdXX_impl acd( nVars, ps );

  if ( acd.run( pTruth ) == 0 )
    return 0;

  return 1;
}

int acdXX_decompose( word * pTruth, unsigned lutSize, unsigned nVars, unsigned char *decomposition )
{
  using namespace acd;

  acdXX_params ps;
  ps.lut_size = lutSize;

  for ( int i = 0; i <= lutSize - 2; ++i )
  {
    ps.max_shared_vars = i;
    ps.min_shared_vars = i;
    acdXX_impl acd( nVars, ps );

    if ( acd.run( pTruth ) == 0 )
      continue;
    acd.compute_decomposition();
    acd.get_decomposition( decomposition );
    return 0;
  }

  return 1;
}

int acd_dc_evaluate( word * pTruth, word * pCareSet, unsigned nVars, int lutSize, unsigned *pdelay, unsigned *cost, int try_no_late_arrival )
{
  using namespace acd;

  ac_dc_decomposition_params ps;
  ps.lut_size = lutSize;
  ps.max_free_set_vars = 4; // same limit as ac_decomposition_params, so -B differs from -Z only by the don't cares
  ps.use_first = false;
  ps.try_no_late_arrival = static_cast<bool>( try_no_late_arrival );
  ac_dc_decomposition_stats st;

  ac_dc_decomposition_impl acd( nVars, ps, &st );
  int val = acd.run( pTruth, pCareSet, *pdelay );

  if ( val < 0 )
  {
    *pdelay = 0;
    return -1;
  }

  *pdelay = acd.get_profile();
  *cost = st.num_luts;
  return val;
}

/* Replays acd_dc_evaluate: *pdelay is the late-arriving mask and try_no_late_arrival the flag of
   the evaluation, not the free set found by it. The column multiplicity with don't cares is a
   heuristic whose result depends on the variable order, so running the engine with the free set
   as profile (as acd_decompose does for the exact plain engine) can fail to find the
   decomposition the evaluation found; the replay is deterministic and returns the same one. */
int acd_dc_decompose( word * pTruth, word * pCareSet, unsigned nVars, int lutSize, unsigned *pdelay, int try_no_late_arrival, unsigned char *decomposition )
{
  using namespace acd;

  ac_dc_decomposition_params ps;
  ps.lut_size = lutSize;
  ps.max_free_set_vars = 4;
  ps.use_first = false;
  ps.try_no_late_arrival = static_cast<bool>( try_no_late_arrival );
  ac_dc_decomposition_stats st;

  ac_dc_decomposition_impl acd( nVars, ps, &st );
  int val = acd.run( pTruth, pCareSet, *pdelay );
  if ( val >= 0 )
    val = acd.compute_decomposition();

  if ( val < 0 )
  {
    *pdelay = 0;
    return -1;
  }

  *pdelay = acd.get_profile();
  acd.get_decomposition( decomposition );
  return 0;
}

ABC_NAMESPACE_IMPL_END
