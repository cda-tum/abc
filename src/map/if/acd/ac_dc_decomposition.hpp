/**C++File**************************************************************

  FileName    [ac_dc_decomposition.hpp]

  SystemName  [ABC: Logic synthesis and verification system.]

  PackageName [Ashenhurst-Curtis decomposition using don't cares.]

  Synopsis    [Interface with the FPGA mapping package.]

  Author      [Benjamin Hien]

  Affiliation [TUM]

  Date        [Ver. 1.0. Started - November 20, 2024.]

***********************************************************************/
/*!
  \file ac_dc_decomposition.hpp
  \brief Ashenhurst-Curtis decomposition using don't cares

  Generated from mockturtle include/mockturtle/algorithms/acd_dc.hpp (commit 0903240) by
  port_acd_dc_to_abc.py; do not edit by hand, change the mockturtle header and regenerate.

  \author Benjamin Hien
*/

#ifndef _ACD_DC_H_
#define _ACD_DC_H_
#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "kitty_constants.hpp"
#include "kitty_constructors.hpp"
#include "kitty_dynamic_tt.hpp"
#include "kitty_operations.hpp"
#include "kitty_operators.hpp"
#include "kitty_static_tt.hpp"

#ifdef _MSC_VER
#  include <intrin.h>
#  define __builtin_popcount __popcnt
#endif

ABC_NAMESPACE_CXX_HEADER_START

namespace acd
{

struct ac_dc_decomposition_params
{
 /*! \brief LUT size for decomposition (3 < num < 7). */
 uint32_t lut_size{ 6 };

 /*! \brief Maximum size of the free set (1 < num < 6). */
 uint32_t max_free_set_vars{ 5 };

 /*! \brief Maximum multiplicity value allowed (2 <= num < 16). */
 uint32_t max_multiplicity{ 16 };

 /*! \brief Perform only support reducing (2-level) decompositions. */
 bool support_reducing_only{ true };

 /*! \brief Use the first feasible decomposition found. */
 bool use_first{ false };

 /*! \brief If decomposition with delay profile fails, try without. */
 bool try_no_late_arrival{ false };
};

struct ac_dc_decomposition_stats
{
 uint32_t num_luts{ 0 };
 uint32_t num_edges{ 0 };
 uint32_t num_levels{ 0 };
};

struct ac_dc_decomposition_result
{
 kitty::dynamic_truth_table tt;
 std::vector<uint32_t> support;
};

class ac_dc_decomposition_impl
{
private:
  struct encoding_column
  {
    uint64_t column[2];
    uint32_t cost;
    uint32_t index;
    float sort_cost;
  };

private:
  static constexpr uint32_t max_num_vars = 11;
  using STT = kitty::static_truth_table<max_num_vars>;
  using word = uint64_t;

public:
  explicit ac_dc_decomposition_impl( uint32_t num_vars, ac_dc_decomposition_params const& ps, ac_dc_decomposition_stats* pst = nullptr )
      : num_vars( num_vars ), ps( ps ), pst( pst )
  {
    std::iota( permutations.begin(), permutations.end(), 0 );
  }

  /*! \brief Runs ACD using late arriving variables */
  int run( word* ptt, word* pcs, unsigned delay_profile )
  {
    /* truth table is too large for the settings */
    if ( num_vars > max_num_vars )
    {
      return -1;
    }

    uint32_t late_arriving = __builtin_popcount( delay_profile );

    /* relax maximum number of free set variables if a function has more variables */
    if ( num_vars > ps.max_free_set_vars + ps.lut_size )
    {
      ps.max_free_set_vars = num_vars - ps.lut_size;
    }
    if ( late_arriving > ps.max_free_set_vars )
    {
      return -1;
      // ps.max_free_set_vars = late_arriving;
    }

    /* return a high cost if too many late arriving variables */
    if ( late_arriving > ps.lut_size - 1 || late_arriving > ps.max_free_set_vars )
    {
      return -1;
    }

    /* convert to static TT */
    init_truth_table( ptt );
    init_care_set( pcs );

    /* permute late arriving variables to be the least significant */
    reposition_late_arriving_variables( delay_profile, late_arriving );

    /* best_cs follows every permutation of best_tt (find_decomposition overwrites best_tt with the
     * permuted table of a feasible free set before trying the next free-set size) */
    best_cs = best_tt_cs;

    /* run ACD trying different bound sets and free sets */
    if ( !find_decomposition( delay_profile, late_arriving ) )
    {
      return -1;
    }

    /* return number of levels (one LUT only if the multiplicity is 1) */
    return ( delay_profile == 0 && best_multiplicity > 1 ) ? 2 : 1;
  }

  /* Completes best_tt on its don't cares according to the column cover found for best_free_set:
   * every partially specified column becomes the base function it was matched with, columns
   * without care minterms become the first base function. Returns the multiplicity of the
   * completed function (equal to best_multiplicity) or UINT32_MAX on failure. */
  uint32_t complete_truth_table_on_dc()
  {
    STT completed = best_tt;
    uint32_t multiplicity = UINT32_MAX;
    switch ( best_free_set )
    {
    case 1:
      multiplicity = column_multiplicity_dc1<1u, true>( best_tt, best_cs, UINT32_MAX, &completed );
      break;
    case 2:
      multiplicity = column_multiplicity_dc2<2u, true>( best_tt, best_cs, UINT32_MAX, &completed );
      break;
    case 3:
      multiplicity = column_multiplicity_dc5<3u, true>( best_tt, best_cs, UINT32_MAX, &completed );
      break;
    case 4:
      multiplicity = column_multiplicity_dc5<4u, true>( best_tt, best_cs, UINT32_MAX, &completed );
      break;
    case 5:
      multiplicity = column_multiplicity_dc5<5u, true>( best_tt, best_cs, UINT32_MAX, &completed );
      break;
    default:
      return UINT32_MAX;
    }

    if ( multiplicity == UINT32_MAX || multiplicity > best_multiplicity )
      return UINT32_MAX;

    best_tt = completed;
    best_multiplicity = multiplicity;
    return multiplicity;
  }

  int compute_decomposition()
  {
    if ( best_multiplicity == UINT32_MAX )
      return -1;

    /* complete the truth table on don't cares, so that the column count equals best_multiplicity */
    if ( complete_truth_table_on_dc() == UINT32_MAX )
    {
      assert( false && "don't-care completion does not reproduce the column multiplicity" );
      return -1;
    }

    /* compute isets */
    std::vector<STT> isets = compute_isets();
    if ( isets.size() != best_multiplicity )
    {
      assert( false && "number of distinct columns differs from the column multiplicity" );
      return -1;
    }

    /* multiplicity 1 (possible only with don't cares): the function depends on the free set only,
     * a single LUT over the free-set variables implements it; no bound-set LUT, no encoding */
    if ( best_multiplicity == 1 )
    {
      best_bound_sets.clear();
      best_care_sets.clear();
      best_iset_onset.clear();
      best_iset_offset.clear();
      return 0;
    }

    generate_support_minimization_encodings();

    /* solves exactly only for small multiplicities */
    if ( best_multiplicity <= 4u )
      solve_min_support_exact( isets );
    else
      solve_min_support_heuristic( isets );

    /* unfeasible decomposition */
    assert( !best_bound_sets.empty() );

    return 0;
  }

  unsigned get_profile()
  {
    unsigned profile = 0;

    if ( best_free_set > num_vars )
      return -1;

    for ( uint32_t i = 0; i < best_free_set; ++i )
    {
      profile |= 1 << permutations[i];
    }

    return profile;
  }

  void get_decomposition( unsigned char* decompArray )
  {
    if ( best_free_set > num_vars )
      return;

    generate_decomposition();
    get_decomposition_abc( decompArray );
  }

  /*std::vector<ac_dc_decomposition_result> const& get_acd_result()
  {
    if ( best_free_set > num_vars )
      return;

    generate_decomposition();
    return dec_result;
  }*/

private:
  bool find_decomposition( unsigned& delay_profile, uint32_t late_arriving )
  {
    best_multiplicity = UINT32_MAX;
    best_free_set = UINT32_MAX;
    uint32_t best_cost = UINT32_MAX;
    uint32_t offset = static_cast<uint32_t>( late_arriving );
    uint32_t start = std::max( offset, 1u );

    /* perform only support reducing decomposition */
    if ( ps.support_reducing_only )
    {
      start = std::max( start, num_vars - ps.lut_size );
    }

    /* array of functions to compute the column multiplicity */
    std::function<uint32_t( STT const& tt )> column_multiplicity_fn[5] = {
        [this]( STT const& tt ) { return column_multiplicity<1u>( tt ); },
        [this]( STT const& tt ) { return column_multiplicity<2u>( tt ); },
        [this]( STT const& tt ) { return column_multiplicity<3u>( tt ); },
        [this]( STT const& tt ) { return column_multiplicity5<4u>( tt ); },
        [this]( STT const& tt ) { return column_multiplicity5<5u>( tt ); } };

    std::function<uint32_t( STT const& tt, STT const& cs, uint32_t loc_cost )> column_multiplicity_fn_dc[5] = {
        [this]( STT const& tt, STT const& cs, uint32_t loc_cost ) { return column_multiplicity_dc1<1u>( tt, cs, loc_cost ); },
        [this]( STT const& tt, STT const& cs, uint32_t loc_cost ) { return column_multiplicity_dc2<2u>( tt, cs, loc_cost ); },
        [this]( STT const& tt, STT const& cs, uint32_t loc_cost ) { return column_multiplicity_dc5<3u>( tt, cs, loc_cost ); },
        [this]( STT const& tt, STT const& cs, uint32_t loc_cost ) { return column_multiplicity_dc5<4u>( tt, cs, loc_cost ); },
        [this]( STT const& tt, STT const& cs, uint32_t loc_cost ) { return column_multiplicity_dc5<5u>( tt, cs, loc_cost ); } };

    /* find a feasible AC decomposition */
    for ( uint32_t i = start; i <= ps.lut_size - 1 && i <= ps.max_free_set_vars; ++i )
    {
      auto [tt_p, cs_p, perm, multiplicity] = enumerate_iset_combinations_offset( i, offset, column_multiplicity_fn_dc[i - 1] );

      /* additional cost if not support reducing */
      uint32_t additional_cost = ( num_vars - i > ps.lut_size ) ? 128 : 0;

      /* check for feasible solution that improves the cost */
      if ( multiplicity <= ( 1 << ( ps.lut_size - i ) ) && multiplicity + additional_cost < best_cost && ( multiplicity <= ps.max_multiplicity || multiplicity == 16 ) )
      {
        best_tt = tt_p;
        best_cs = cs_p;
        permutations = perm;
        best_multiplicity = multiplicity;
        best_cost = multiplicity + additional_cost;
        best_free_set = i;

        if ( !ps.use_first && multiplicity > 2 )
        {
          continue;
        }

        break;
      }
    }

    if ( best_multiplicity == UINT32_MAX && ( !ps.try_no_late_arrival || late_arriving == 0 ) )
      return false;

    /* try without the delay profile */
    if ( best_multiplicity == UINT32_MAX )
    {
      delay_profile = 0;
      if ( ps.support_reducing_only )
      {
        start = std::max( 1u, num_vars - ps.lut_size );
      }

      for ( uint32_t i = start; i <= ps.lut_size - 1 && i <= ps.max_free_set_vars; ++i )
      {
        auto [tt_p, cs_p, perm, multiplicity] = enumerate_iset_combinations_offset( i, 0, column_multiplicity_fn_dc[i - 1] );

        /* additional cost if not support reducing */
        uint32_t additional_cost = ( num_vars - i > ps.lut_size ) ? 128 : 0;

        /* check for feasible solution that improves the cost */
        if ( multiplicity <= ( 1 << ( ps.lut_size - i ) ) && multiplicity + additional_cost < best_cost && ( multiplicity <= ps.max_multiplicity || multiplicity == 16 ) )
        {
          best_tt = tt_p;
          best_cs = cs_p;
          permutations = perm;
          best_multiplicity = multiplicity;
          best_cost = multiplicity + additional_cost;
          best_free_set = i;

          if ( !ps.use_first )
          {
            continue;
          }

          break;
        }
      }
    }

    if ( best_multiplicity == UINT32_MAX )
      return false;

    /* estimation on number of LUTs (multiplicity 1: a single LUT over the free set) */
    if ( pst )
    {
      pst->num_luts = best_multiplicity <= 1 ? 1 : best_multiplicity <= 2 ? 2 : best_multiplicity <= 4 ? 3
                                                                          : best_multiplicity <= 8   ? 4
                                                                                                     : 5;
      pst->num_edges = ( pst->num_luts - 1 ) * ( num_vars - best_free_set ) + ( pst->num_luts - 1 ) + best_free_set;
    }

    return true;
  }

  void init_truth_table( word* ptt )
  {
    uint32_t const num_blocks = ( num_vars <= 6 ) ? 1 : ( 1 << ( num_vars - 6 ) );

    for ( uint32_t i = 0; i < num_blocks; ++i )
    {
      best_tt._bits[i] = ptt[i];
    }

    local_extend_to( best_tt, num_vars );
  }

  void init_care_set( const word* pcs )
  {
    uint32_t const num_blocks = ( num_vars <= 6 ) ? 1 : ( 1 << ( num_vars - 6 ) );

    for ( uint32_t i = 0; i < num_blocks; ++i )
    {
      best_tt_cs._bits[i] = pcs[i];
    }

    local_extend_to( best_tt_cs, num_vars );
  }

  template<uint32_t free_set_size>
  uint32_t column_multiplicity( STT tt )
  {
    uint64_t multiplicity_set[4] = { 0u, 0u, 0u, 0u };
    uint32_t multiplicity = 0;
    uint32_t const num_blocks = ( num_vars > 6 ) ? ( 1u << ( num_vars - 6 ) ) : 1;
    uint64_t constexpr masks_bits[] = { 0x0, 0x3, 0xF, 0x3F };
    uint64_t constexpr masks_idx[] = { 0x0, 0x0, 0x0, 0x3 };

    /* supports up to 64 values of free set (256 for |FS| == 3)*/
    static_assert( free_set_size <= 3 );

    /* extract iset functions */
    auto it = std::begin( tt );
    for ( auto i = 0u; i < num_blocks; ++i )
    {
      for ( auto j = 0; j < ( 64 >> free_set_size ); ++j )
      {
        multiplicity_set[( *it >> 6 ) & masks_idx[free_set_size]] |= UINT64_C( 1 ) << ( *it & masks_bits[free_set_size] );
        *it >>= ( 1u << free_set_size );
      }
      ++it;
    }

    multiplicity = __builtin_popcountl( multiplicity_set[0] );

    if constexpr ( free_set_size == 3 )
    {
      multiplicity += __builtin_popcountl( multiplicity_set[1] );
      multiplicity += __builtin_popcountl( multiplicity_set[2] );
      multiplicity += __builtin_popcountl( multiplicity_set[3] );
    }

    return multiplicity;
  }

  // Helper: Assign partial functions to selected representative values
  static void encode_dc1( uint64_t* mapping, uint64_t mask, uint64_t cof_masked )
  {
    while ( mask )
    {
      uint32_t pos = __builtin_ctzll( mask );
      mask &= mask - 1;

      if ( mapping[pos] == 0xF )
      {
        mapping[pos] = cof_masked;
      }
    }
  }

  // Greedy approximation of the minimum hitting set problem
  static void min_hitting_set1( uint64_t& selected_set, uint64_t uncovered_mask )
  {
    constexpr uint64_t coverage_masks[] = { 0x5, 0x6, 0x9, 0xA };

    while ( uncovered_mask )
    {
      uint32_t best_coverage = 0;
      uint8_t best_index = 0xFF;

      for ( uint8_t i = 0; i < 4; ++i )
      {
        if ( selected_set & ( 1ULL << i ) )
          continue;

        uint64_t coverage = coverage_masks[i] & uncovered_mask;
        uint32_t count = __builtin_popcountll( coverage );

        if ( count > best_coverage )
        {
          best_coverage = count;
          best_index = i;
        }
      }

      if ( best_index == 0xFF )
        break;

      uncovered_mask &= ~coverage_masks[best_index];
      selected_set |= ( 1ULL << best_index );
    }
  }

  template<uint32_t free_set_size, bool manipulate = false>
  uint32_t column_multiplicity_dc1( STT const& tt, STT const& cs, uint32_t loc_cost, STT* out_tt = nullptr )
  {
    static_assert( free_set_size == 1, "Expected free_set_size to be 1 for DC1 optimization." );

    constexpr uint64_t coverage_masks[] = { 0x5, 0x6, 0x9, 0xA };

    uint64_t multiplicity_set = 0;
    uint64_t encoding_mask = 0;
    uint32_t multiplicity = 0;

    uint32_t const num_blocks = ( num_vars > 6 ) ? ( 1u << ( num_vars - 6 ) ) : 1;

    // Analyze TT + CS to determine base functions and coverage
    for ( uint32_t i = 0; i < num_blocks; ++i )
    {
      uint64_t cof = tt._bits[i];
      uint64_t care = cs._bits[i];

      for ( uint32_t j = 0; j < ( 64 >> free_set_size ); ++j )
      {
        uint64_t cs_bits = care & 0x3;
        uint64_t tt_bits = cof & 0x3;

        if ( cs_bits == 0x3 ) // fully specified
        {
          multiplicity_set |= 1ULL << tt_bits;
          encoding_mask |= coverage_masks[tt_bits];
        }
        else if ( cs_bits != 0 ) // partial care
        {
          uint32_t index = 2u + ( cs_bits << 1u ) + __builtin_popcountl( cs_bits & tt_bits );
          multiplicity_set |= 1ULL << index;
        }

        cof >>= 2;
        care >>= 2;
      }
    }

    // Invert encoding mask to find uncovered elements
    encoding_mask = ~encoding_mask & ( multiplicity_set >> 4 );
    min_hitting_set1( multiplicity_set, encoding_mask );

    // Count number of selected base functions
    multiplicity_set &= 0xF;
    multiplicity = __builtin_popcountl( multiplicity_set );

    /* a function without any care minterm has one (arbitrary) column */
    if ( multiplicity == 0 )
    {
      multiplicity_set = 0x1;
      multiplicity = 1;
    }

    // If better than current best, optionally manipulate the TT

    if constexpr ( manipulate )
    {
      STT new_tt = tt;
      uint64_t mapping[4] = { 0xF, 0xF, 0xF, 0xF };

      // Assign unassigned partials to selected base values
      uint64_t default_dc_val = 0xF;
      uint64_t temp_set = multiplicity_set;

      while ( temp_set )
      {
        uint32_t pos = __builtin_ctzl( temp_set );
        temp_set &= temp_set - 1;

        encode_dc1( mapping, coverage_masks[pos], pos );

        if ( default_dc_val == 0xF )
        {
          default_dc_val = pos;
        }
      }

      for ( uint32_t i = 0; i < num_blocks; ++i )
      {
        uint64_t old_tt = tt._bits[i];
        uint64_t care = cs._bits[i];
        uint64_t& new_tt_block = new_tt._bits[i];

        for ( uint32_t j = 0; j < ( 64 >> free_set_size ); ++j )
        {
          uint64_t cs_bits = care & 0x3;
          uint64_t tt_bits = old_tt & 0x3;

          uint64_t new_val = 0;

          if ( cs_bits != 0 )
          {
            if ( cs_bits != 3 )
            {
              uint32_t idx = ( cs_bits << 1 ) + __builtin_popcountl( cs_bits & tt_bits ) - 2;
              new_val = mapping[idx] & 0x3;
            }
            else
            {
              new_val = tt_bits;
            }
          }
          else
          {
            new_val = default_dc_val;
          }

          uint64_t shift = j * 2;
          new_tt_block = ( new_tt_block & ~( 0x3ULL << shift ) ) | ( new_val << shift );

          old_tt >>= 2;
          care >>= 2;
        }

        // Ensure new TT agrees with original where care set is 1
        assert( ( ( new_tt_block ^ tt._bits[i] ) & cs._bits[i] ) == 0 );
      }

      if ( out_tt )
        *out_tt = new_tt;
    }
    return multiplicity;
  }

  static inline uint32_t extract_relevant_bits( uint64_t ccs_mask, uint64_t cof_mask, uint32_t pop )
  {
    uint32_t compacted_value = 0;
    uint32_t bit_position = 0;

    for ( uint32_t i = 0; i < pop; ++i )
    {                                                                        // Loop runs at most 4 times
      uint32_t lowest_bit = __builtin_ctzll( ccs_mask );                     // Find rightmost set bit in ccs_mask
      compacted_value |= ( ( cof_mask >> lowest_bit ) & 1 ) << bit_position; // Extract and shift
      ccs_mask ^= ( UINT64_C( 1 ) << lowest_bit );                           // Remove lowest set bit efficiently
      ++bit_position;
    }

    return compacted_value;
  }

  // Helper: Map uncovered DC terms to a selected representative
  static void encode_dc2( uint64_t* mapping, uint64_t mask, uint64_t cof_masked )
  {
    while ( mask )
    {
      uint32_t pos = __builtin_ctzll( mask );
      mask &= mask - 1;
      if ( mapping[pos] == 0xFF )
        mapping[pos] = cof_masked;
    }
  }

  // Greedy Minimum Hitting Set solver for 4-variable functions
  static void min_hitting_set2( uint64_t& selected_set, uint64_t uncovered_mask )
  {
    constexpr uint64_t coverage_masks[] = {
        0x0101010111111155, 0x0102020211212256, 0x0201040412121459, 0x020208081222285A,
        0x0404011021144165, 0x0408022021248266, 0x0804044022184469, 0x080808802228886A,
        0x1010100144411195, 0x1020200244812296, 0x2010400448421499, 0x202080084882289A,
        0x40401010844441A5, 0x40802020848482A6, 0x80404040884844A9, 0x80808080888888AA };

    std::vector<uint8_t> active_indices;
    for ( uint8_t i = 0; i < 16; ++i )
    {
      if ( coverage_masks[i] & uncovered_mask )
      {
        active_indices.push_back( i );
      }
    }

    while ( uncovered_mask )
    {
      uint32_t best_coverage = 0;
      uint8_t best_index = 0xFF;

      for ( size_t j = 0; j < active_indices.size(); )
      {
        uint8_t i = active_indices[j];

        if ( selected_set & ( UINT64_C( 1 ) << i ) )
        {
          active_indices.erase( active_indices.begin() + j );
          continue;
        }

        uint64_t coverage = coverage_masks[i] & uncovered_mask;
        uint32_t count = __builtin_popcountll( coverage );

        if ( count > best_coverage )
        {
          best_coverage = count;
          best_index = i;
        }

        if ( count == 0 )
        {
          active_indices.erase( active_indices.begin() + j );
        }
        else
        {
          ++j;
        }
      }

      if ( best_index == 0xFF )
        break;

      uncovered_mask &= ~coverage_masks[best_index];
      selected_set |= ( UINT64_C( 1 ) << best_index );
    }
  }

  // Compute the minimal column multiplicity with optional TT manipulation for free_set_size = 2
  template<uint32_t free_set_size, bool manipulate = false>
  uint32_t column_multiplicity_dc2( const STT& tt, const STT& cs, uint32_t loc_cost, STT* out_tt = nullptr )
  {
    static_assert( free_set_size == 2, "Wrong free set size for method used, expected 2" );

    constexpr uint64_t encode[16] = { 255, 0, 1, 4, 2, 5, 6, 10, 3, 7, 8, 11, 9, 12, 13, 14 };
    uint64_t constexpr coverage_masks[] = {
        0x0101010111111155, 0x0102020211212256, 0x0201040412121459, 0x020208081222285A,
        0x0404011021144165, 0x0408022021248266, 0x0804044022184469, 0x080808802228886A,
        0x1010100144411195, 0x1020200244812296, 0x2010400448421499, 0x202080084882289A,
        0x40401010844441A5, 0x40802020848482A6, 0x80404040884844A9, 0x80808080888888AA };

    const uint32_t num_blocks = ( num_vars > 6 ) ? ( 1u << ( num_vars - 6 ) ) : 1;
    uint64_t selected_set = 0;
    uint64_t dc_set = 0;
    uint64_t uncovered_mask = 0;
    uint64_t mapping[64];
    std::fill( std::begin( mapping ), std::end( mapping ), 0xFF );

    for ( uint32_t i = 0; i < num_blocks; ++i )
    {
      uint64_t tt_block = tt._bits[i];
      uint64_t cs_block = cs._bits[i];

      for ( uint32_t j = 0; j < ( 64 >> free_set_size ); ++j )
      {
        uint64_t care = cs_block & 0xF;
        uint64_t value = tt_block & 0xF;

        if ( care == 0xF )
        {
          selected_set |= UINT64_C( 1 ) << value;
          uncovered_mask |= coverage_masks[value];
        }
        else if ( care )
        {
          uint32_t count = __builtin_popcountll( care );
          uint32_t idx = extract_relevant_bits( care, value, count );

          if ( count == 1 )
          {
            dc_set |= UINT64_C( 1 ) << ( encode[care] * 2u + idx );
          }
          else if ( count == 2 )
          {
            dc_set |= UINT64_C( 1 ) << ( 8u + ( encode[care] - 4u ) * 4u + idx );
          }
          else if ( count == 3 )
          {
            dc_set |= UINT64_C( 1 ) << ( 32u + ( encode[care] - 10u ) * 8u + idx );
          }
        }

        tt_block >>= 4;
        cs_block >>= 4;
      }
    }

    if ( __builtin_popcountll( selected_set ) >= loc_cost )
    {
      return UINT32_MAX;
    }

    uncovered_mask = ~uncovered_mask & dc_set;
    min_hitting_set2( selected_set, uncovered_mask );

    uint32_t multiplicity = __builtin_popcountll( selected_set );
    assert( multiplicity <= 16 && "Bug" );

    /* a function without any care minterm has one (arbitrary) column */
    if ( multiplicity == 0 )
    {
      selected_set = 0x1;
      multiplicity = 1;
    }

    if constexpr ( manipulate )
    {
      STT new_tt = tt;
      uint64_t default_value = 0xFF;

      uint64_t tmp_set = selected_set;
      while ( tmp_set )
      {
        uint32_t index = __builtin_ctzll( tmp_set );
        tmp_set &= tmp_set - 1;
        encode_dc2( mapping, coverage_masks[index], index );
        if ( default_value == 0xFF )
          default_value = index;
      }

      for ( uint32_t i = 0; i < num_blocks; ++i )
      {
        uint64_t tt_block = tt._bits[i];
        uint64_t cs_block = cs._bits[i];
        uint64_t& new_block = new_tt._bits[i];

        for ( uint32_t j = 0; j < ( 64 >> free_set_size ); ++j )
        {
          uint64_t care = cs_block & 0xF;
          uint64_t value = tt_block & 0xF;

          if ( care )
          {
            if ( care != 0xF )
            {
              uint32_t count = __builtin_popcountll( care );
              uint32_t idx = extract_relevant_bits( care, value, count );
              uint32_t map_idx = 0;

              if ( count == 1 )
                map_idx = encode[care] * 2u + idx;
              else if ( count == 2 )
                map_idx = 8u + ( encode[care] - 4u ) * 4u + idx;
              else if ( count == 3 )
                map_idx = 32u + ( encode[care] - 10u ) * 8u + idx;

              new_block = ( new_block & ~( UINT64_C( 0xF ) << ( j * 4 ) ) ) |
                          ( mapping[map_idx] << ( j * 4 ) );
            }
          }
          else
          {
            new_block = ( new_block & ~( UINT64_C( 0xF ) << ( j * 4 ) ) ) |
                        ( default_value << ( j * 4 ) );
          }

          tt_block >>= 4;
          cs_block >>= 4;
        }
        // Ensure new TT agrees with original where care set is 1
        assert( ( ( new_block ^ tt._bits[i] ) & cs._bits[i] ) == 0 );
      }

      if ( out_tt )
        *out_tt = new_tt;
    }

    return multiplicity;
  }

  template<uint32_t free_set_size>
  uint32_t column_multiplicity5( STT tt )
  {
    uint32_t const num_blocks = ( num_vars > 6 ) ? ( 1u << ( num_vars - 6 ) ) : 1;
    uint64_t constexpr masks[] = { 0x0, 0x3, 0xF, 0xFF, 0xFFFF, 0xFFFFFFFF };

    static_assert( free_set_size == 5 || free_set_size == 4 );

    uint32_t size = 0;
    uint64_t prev = -1;
    std::array<uint32_t, 64> multiplicity_set;

    /* extract iset functions */
    auto it = std::begin( tt );
    for ( auto i = 0u; i < num_blocks; ++i )
    {
      for ( auto j = 0; j < ( 64 >> free_set_size ); ++j )
      {
        uint32_t fs_fn = static_cast<uint32_t>( *it & masks[free_set_size] );
        if ( fs_fn != prev )
        {
          multiplicity_set[size++] = fs_fn;
          prev = fs_fn;
        }
        *it >>= ( 1u << free_set_size );
      }
      ++it;
    }

    std::sort( multiplicity_set.begin(), multiplicity_set.begin() + size );

    /* count unique */
    uint32_t multiplicity = 1;
    for ( auto i = 1u; i < size; ++i )
    {
      multiplicity += multiplicity_set[i] != multiplicity_set[i - 1] ? 1 : 0;
    }

    return multiplicity;
  }

  template<uint32_t free_set_size, bool manipulate = false>
  uint32_t column_multiplicity_dc5( STT const& tt, STT const& cs, uint32_t loc_cost, STT* out_tt = nullptr )
  {
    uint32_t const num_blocks = ( num_vars > 6 ) ? ( 1u << ( num_vars - 6 ) ) : 1;
    uint64_t constexpr masks[] = { 0x0, 0x3, 0xF, 0xFF, 0xFFFF, 0xFFFFFFFF };

    uint32_t size = 0;
    uint32_t partial_size = 0;
    uint64_t prev = -1;

    std::array<uint32_t, 64> base_set{};
    std::array<uint32_t, 64> partial_fn{};
    std::array<uint32_t, 64> partial_cs{};
    std::array<uint32_t, 64> block_index{};
    std::array<uint32_t, 64> entry_index{};
    std::array<uint32_t, 64> dc_block_index{};
    std::array<uint32_t, 64> dc_entry_index{};
    uint32_t dc_size = 0;

    STT new_tt = tt;

    for ( auto i = 0u; i < num_blocks; ++i )
    {
      uint64_t cof = tt._bits[i];
      uint64_t ccs = cs._bits[i];

      for ( auto j = 0; j < ( 64 >> free_set_size ); ++j )
      {
        uint32_t fs_fn = static_cast<uint32_t>( cof & masks[free_set_size] );
        uint32_t fs_cs = static_cast<uint32_t>( ccs & masks[free_set_size] );

        if ( fs_cs == masks[free_set_size] )
        {
          if ( fs_fn != prev )
          {
            base_set[size++] = fs_fn;
            prev = fs_fn;
          }
        }
        else if ( fs_cs )
        {
          partial_fn[partial_size] = fs_fn;
          partial_cs[partial_size] = fs_cs;
          block_index[partial_size] = i;
          entry_index[partial_size] = j;
          ++partial_size;
        }
        else if constexpr ( manipulate )
        {
          /* column without care minterms: completed later with the first base function */
          dc_block_index[dc_size] = i;
          dc_entry_index[dc_size] = j;
          ++dc_size;
        }

        cof >>= ( 1u << free_set_size );
        ccs >>= ( 1u << free_set_size );
      }
    }

    std::sort( base_set.begin(), base_set.begin() + size );
    uint32_t unique_size = ( size == 0 ) ? 0 : 1;

    for ( size_t i = 1; i < size; ++i )
    {
      if ( base_set[i] != base_set[unique_size - 1] )
      {
        base_set[unique_size++] = base_set[i];
      }
    }
    if ( unique_size >= loc_cost )
    {
      return UINT32_MAX;
    }
    // uint32_t unique_size = std::unique( base_set.begin(), base_set.begin() + size ) - base_set.begin();

    for ( uint32_t i = 0; i < partial_size; ++i )
    {
      bool matched = false;

      for ( uint32_t j = 0; j < unique_size; ++j )
      {
        if ( ( base_set[j] & partial_cs[i] ) == ( partial_fn[i] & partial_cs[i] ) )
        {
          matched = true;

          if constexpr ( manipulate )
          {
            /* the partially specified column becomes the matched base function */
            uint32_t block = block_index[i];
            uint32_t pos = entry_index[i] * ( 1u << free_set_size );
            uint64_t mask = masks[free_set_size] << pos;
            uint64_t value = static_cast<uint64_t>( base_set[j] ) << pos;

            new_tt._bits[block] = ( new_tt._bits[block] & ~mask ) | value;
          }

          break;
        }
      }

      if ( !matched )
      {
        base_set[unique_size++] = partial_fn[i];
        if ( unique_size >= loc_cost )
        {
          return UINT32_MAX;
        }
      }
    }

    /* a function without any care minterm has one (arbitrary) column */
    if ( unique_size == 0 )
    {
      base_set[unique_size++] = 0;
    }

    if constexpr ( manipulate )
    {
      /* columns without care minterms become the first base function */
      for ( uint32_t i = 0; i < dc_size; ++i )
      {
        uint32_t pos = dc_entry_index[i] * ( 1u << free_set_size );
        uint64_t mask = masks[free_set_size] << pos;
        new_tt._bits[dc_block_index[i]] = ( new_tt._bits[dc_block_index[i]] & ~mask ) | ( static_cast<uint64_t>( base_set[0] ) << pos );
      }

      for ( size_t i = 0; i < num_blocks; ++i )
      {
        assert( ( new_tt._bits[i] & cs._bits[i] ) == ( tt._bits[i] & cs._bits[i] ) &&
                "Modified truth table is not equivalent under care set!" );
      }

      if ( out_tt )
        *out_tt = new_tt;
    }

    return unique_size;
  }

  uint32_t column_multiplicity2( STT const& tt, uint32_t free_set_size )
  {
    assert( free_set_size <= 5 );

    uint32_t const num_blocks = ( num_vars > 6 ) ? ( 1u << ( num_vars - 6 ) ) : 1;
    uint64_t const shift = UINT64_C( 1 ) << free_set_size;
    uint64_t const mask = ( UINT64_C( 1 ) << shift ) - 1;
    uint32_t cofactors[4];
    uint32_t size = 0;

    /* extract iset functions */
    for ( auto i = 0u; i < num_blocks; ++i )
    {
      uint64_t sub = tt._bits[i];
      for ( auto j = 0; j < ( 64 >> free_set_size ); ++j )
      {
        uint32_t fs_fn = static_cast<uint32_t>( sub & mask );
        uint32_t k;
        for ( k = 0; k < size; ++k )
        {
          if ( fs_fn == cofactors[k] )
            break;
        }
        if ( k == 2 )
          return 3;
        if ( k == size )
          cofactors[size++] = fs_fn;
        sub >>= shift;
      }
    }

    return size;
  }

  inline bool combinations_offset_next( uint32_t k, uint32_t offset, uint32_t* pComb, uint32_t* pInvPerm, STT& tt, STT& cs )
  {
    uint32_t i;

    for ( i = k - 1; pComb[i] == num_vars - k + i; --i )
    {
      if ( i == offset )
        return false;
    }

    /* move vars */
    uint32_t var_old = pComb[i];
    uint32_t pos_new = pInvPerm[var_old + 1];
    std::swap( pInvPerm[var_old + 1], pInvPerm[var_old] );
    std::swap( pComb[i], pComb[pos_new] );
    kitty::swap_inplace( tt, i, pos_new );
    kitty::swap_inplace( cs, i, pos_new );

    for ( uint32_t j = i + 1; j < k; j++ )
    {
      var_old = pComb[j];
      pos_new = pInvPerm[pComb[j - 1] + 1];
      std::swap( pInvPerm[pComb[j - 1] + 1], pInvPerm[var_old] );
      std::swap( pComb[j], pComb[pos_new] );
      kitty::swap_inplace( tt, j, pos_new );
      kitty::swap_inplace( cs, j, pos_new );
    }

    return true;
  }

  template<typename Fn>
  std::tuple<STT, STT, std::array<uint32_t, max_num_vars>, uint32_t> enumerate_iset_combinations_offset( uint32_t free_set_size, uint32_t offset, Fn&& fn )
  {
    STT tt = best_tt;
    STT cs = best_cs;

    /* TT with best cost */
    STT local_best_tt = tt;
    STT local_best_cs = cs;
    uint32_t best_cost = UINT32_MAX;

    assert( free_set_size >= offset );

    /* special case */
    if ( free_set_size == offset )
    {
      best_cost = fn( tt, cs, best_cost );
      return { tt, cs, permutations, best_cost };
    }

    /* works up to 16 input truth tables */
    assert( num_vars <= 16 );

    /* Search for column multiplicity of 2 */
    if ( free_set_size == ps.lut_size - 1 )
    {
      return enumerate_iset_combinations2( free_set_size, offset, fn );
    }

    /* bail early combination */
    uint32_t bail_early = 2;
    if ( best_multiplicity < UINT32_MAX )
      bail_early = ( best_multiplicity >> 1 ) + ( best_multiplicity & 1 );

    /* init combinations */
    uint32_t pComb[16], pInvPerm[16], bestPerm[16];
    for ( uint32_t i = 0; i < num_vars; ++i )
    {
      pComb[i] = pInvPerm[i] = i;
    }

    /* enumerate combinations */
    do
    {
      uint32_t cost = fn( tt, cs, best_cost );
      if ( cost < best_cost )
      {
        local_best_tt = tt;
        local_best_cs = cs;
        best_cost = cost;
        for ( uint32_t i = 0; i < num_vars; ++i )
        {
          bestPerm[i] = pComb[i];
        }

        if ( best_cost <= bail_early )
          break;
      }
    } while ( combinations_offset_next( free_set_size, offset, pComb, pInvPerm, tt, cs ) );

    std::array<uint32_t, max_num_vars> res_perm;
    for ( uint32_t i = 0; i < num_vars; ++i )
    {
      res_perm[i] = permutations[bestPerm[i]];
    }

    return std::make_tuple( local_best_tt, local_best_cs, res_perm, best_cost );
  }

  /* free set of lut_size - 1 variables: only a column multiplicity of at most 2 is feasible; the
   * don't-care-aware multiplicity function bails out (UINT32_MAX) as soon as 3 columns are needed */
  template<typename Fn>
  inline std::tuple<STT, STT, std::array<uint32_t, max_num_vars>, uint32_t> enumerate_iset_combinations2( uint32_t free_set_size, uint32_t offset, Fn&& fn )
  {
    STT tt = best_tt;
    STT cs = best_cs;

    /* TT with best cost */
    STT local_best_tt = tt;
    STT local_best_cs = cs;
    uint32_t best_cost = ( 1 << ( ps.lut_size - free_set_size ) ) + 1;

    assert( free_set_size >= offset );

    /* init combinations */
    uint32_t pComb[16], pInvPerm[16];
    for ( uint32_t i = 0; i < num_vars; ++i )
    {
      pComb[i] = pInvPerm[i] = i;
    }

    /* enumerate combinations */
    std::array<uint32_t, max_num_vars> res_perm;

    do
    {
      uint32_t cost = fn( tt, cs, 3u );
      if ( cost <= 2 )
      {
        local_best_tt = tt;
        local_best_cs = cs;
        best_cost = cost;
        for ( uint32_t i = 0; i < num_vars; ++i )
        {
          res_perm[i] = permutations[pComb[i]];
        }
        return std::make_tuple( local_best_tt, local_best_cs, res_perm, best_cost );
      }
    } while ( combinations_offset_next( free_set_size, offset, pComb, pInvPerm, tt, cs ) );

    return std::make_tuple( local_best_tt, local_best_cs, res_perm, UINT32_MAX );
  }

  std::vector<STT> compute_isets( bool verbose = false )
  {
    /* construct isets involved in multiplicity */
    uint32_t isets_support = num_vars - best_free_set;
    std::vector<STT> isets( best_multiplicity );

    /* construct isets */
    std::unordered_map<uint64_t, uint32_t> column_to_iset;
    STT tt = best_tt;
    uint32_t offset = 0;
    uint32_t num_blocks = ( num_vars > 6 ) ? ( 1u << ( num_vars - 6 ) ) : 1;
    uint64_t constexpr masks[] = { 0x0, 0x3, 0xF, 0xFF, 0xFFFF, 0xFFFFFFFF };

    auto it = std::begin( tt );
    for ( auto i = 0u; i < num_blocks; ++i )
    {
      for ( auto j = 0; j < ( 64 >> best_free_set ); ++j )
      {
        uint64_t val = *it & masks[best_free_set];

        if ( auto el = column_to_iset.find( val ); el != column_to_iset.end() )
        {
          isets[el->second]._bits[i / ( 1u << best_free_set )] |= UINT64_C( 1 ) << ( j + offset );
        }
        else
        {
          if ( column_to_iset.size() >= isets.size() )
          {
            /* more distinct columns than the multiplicity: caller reports failure */
            isets.emplace_back();
          }
          isets[column_to_iset.size()]._bits[i / ( 1u << best_free_set )] |= UINT64_C( 1 ) << ( j + offset );
          column_to_iset[val] = column_to_iset.size();
        }

        *it >>= ( 1u << best_free_set );
      }

      offset = ( offset + ( 64 >> best_free_set ) ) % 64;
      ++it;
    }

    /* extend isets to cover the whole truth table */
    for ( STT& iset : isets )
    {
      local_extend_to( iset, isets_support );
    }

    /* save free_set functions */
    std::vector<STT> free_set_tts( best_multiplicity );

    for ( auto const& pair : column_to_iset )
    {
      free_set_tts[pair.second]._bits[0] = pair.first;
      local_extend_to( free_set_tts[pair.second], best_free_set );
    }

    /* print isets  and free set*/
    if ( verbose )
    {
      std::cout << "iSets\n";
      uint32_t i = 0;
      for ( auto iset : isets )
      {
        kitty::print_hex( iset );
        std::cout << " of func ";
        kitty::print_hex( free_set_tts[i++] );
        std::cout << "\n";
      }
    }

    best_free_set_tts = std::move( free_set_tts );

    return isets;
  }

  void generate_decomposition()
  {
    dec_result.clear();

    uint32_t num_edges = 0;
    for ( uint32_t i = 0; i < best_bound_sets.size(); ++i )
    {
      ac_dc_decomposition_result dec;
      auto tt = best_bound_sets[i];
      auto care = best_care_sets[i];

      /* compute and minimize support for bound set variables */
      uint32_t k = 0;
      for ( uint32_t j = 0; j < num_vars - best_free_set; ++j )
      {
        if ( !kitty::has_var( tt, care, j ) )
        {
          /* fix truth table */
          adjust_truth_table_on_dc( tt, care, j );
          continue;
        }

        if ( k < j )
        {
          kitty::swap_inplace( tt, k, j );
          kitty::swap_inplace( care, k, j );
        }
        dec.support.push_back( permutations[best_free_set + j] );
        ++k;
      }

      dec.tt = kitty::shrink_to( tt, dec.support.size() );
      dec_result.push_back( dec );
      num_edges += dec.support.size() > 1 ? dec.support.size() : 0;
    }

    /* compute the decomposition for the top-level LUT */
    compute_top_lut_decomposition();

    if ( pst )
    {
      pst->num_luts = dec_result.size();
      pst->num_edges = num_edges + dec_result.back().support.size();
    }
  }

  void compute_top_lut_decomposition()
  {
    uint32_t top_vars = best_bound_sets.size() + best_free_set;
    assert( top_vars <= ps.lut_size );

    /* extend bound set functions with free_set_size LSB vars */
    kitty::dynamic_truth_table tt( top_vars );

    /* compute support */
    dec_result.emplace_back();
    for ( uint32_t i = 0; i < best_free_set; ++i )
    {
      dec_result.back().support.push_back( permutations[i] );
    }

    /* create functions for bound set */
    std::vector<kitty::dynamic_truth_table> bound_set_vars;
    auto res_it = dec_result.begin();
    uint32_t offset = 0;
    for ( uint32_t i = 0; i < best_bound_sets.size(); ++i )
    {
      bound_set_vars.emplace_back( top_vars );
      kitty::create_nth_var( bound_set_vars[i], best_free_set + i );

      /* add bound-set variables to the support, remove buffers (shared set) */
      if ( res_it->support.size() == 1 )
      {
        dec_result.back().support.push_back( res_it->support.front() );
        /* it is a NOT */
        if ( ( res_it->tt._bits[0] & 1 ) == 1 )
        {
          bound_set_vars[i] = ~bound_set_vars[i];
        }
        dec_result.erase( res_it );
        ++offset;
      }
      else
      {
        dec_result.back().support.push_back( num_vars + i - offset );
        ++res_it;
      }
    }

    /* create composition function */
    for ( uint32_t i = 0; i < best_free_set_tts.size(); ++i )
    {
      kitty::dynamic_truth_table free_set_tt = kitty::shrink_to( best_free_set_tts[i], top_vars );

      /* find MUX assignments */
      for ( uint32_t j = 0; j < bound_set_vars.size(); ++j )
      {
        /* AND with ONSET or OFFSET */
        if ( ( ( best_iset_onset[j] >> i ) & 1 ) )
        {
          free_set_tt &= bound_set_vars[j];
        }
        else if ( ( ( best_iset_offset[j] >> i ) & 1 ) )
        {
          free_set_tt &= ~bound_set_vars[j];
        }
      }

      tt |= free_set_tt;
    }

    /* add top-level LUT to result */
    dec_result.back().tt = tt;
  }

  inline void reposition_late_arriving_variables( unsigned delay_profile, uint32_t late_arriving )
  {
    uint32_t k = 0;
    for ( uint32_t i = 0; i < late_arriving; ++i )
    {
      while ( ( ( delay_profile >> k ) & 1 ) == 0 )
        ++k;

      if ( permutations[i] == k )
      {
        ++k;
        continue;
      }

      std::swap( permutations[i], permutations[k] );
      kitty::swap_inplace( best_tt, i, k );
      kitty::swap_inplace( best_tt_cs, i, k );
      ++k;
    }
  }

  template<class Iterator>
  void print_perm( Iterator begin, uint32_t free_set )
  {
    std::cout << "[";
    for ( uint32_t i = 0; i < num_vars; ++i )
    {
      if ( i == free_set )
      {
        std::cout << ", ";
      }
      std::cout << *begin << " ";
      ++begin;
    }
    std::cout << "]\n";
  }

  void generate_support_minimization_encodings()
  {
    uint32_t count = 0;

    /* enable don't cares only if not a power of 2 */
    uint32_t num_combs = 2;
    if ( __builtin_popcount( best_multiplicity ) == 1 )
    {
      uint32_t num_combs_exact[4] = { 1, 3, 35, 6435 };
      for ( uint32_t i = 0; i < 4; ++i )
      {
        if ( ( best_multiplicity >> i ) == 2u )
        {
          num_combs = num_combs_exact[i];
        }
      }
      support_minimization_encodings = std::vector<std::array<uint32_t, 2>>( num_combs );
      generate_support_minimization_encodings_rec<false, true>( 0, 0, 0, count, best_multiplicity >> 1, true );
      assert( count == num_combs );
      return;
    }

    /* constraint the size of a ON-set OFF-set class for a strict encoding */
    uint32_t min_set_size = 1;
    if ( best_multiplicity <= 4 )
      min_set_size = 2;
    else if ( best_multiplicity <= 8 )
      min_set_size = 4;
    else
      min_set_size = 8;
    min_set_size = best_multiplicity - min_set_size;

    if ( best_multiplicity > 8 )
    {
      uint32_t class_sizes[13] = { 3, 3, 15, 25, 35, 35, 255, 501, 957, 1749, 3003, 4719, 6435 };
      num_combs = class_sizes[best_multiplicity - 3];
      support_minimization_encodings = std::vector<std::array<uint32_t, 2>>( num_combs );
      generate_support_minimization_encodings_rec<false, false>( 0, 0, 0, count, min_set_size, true );
    }
    else
    {
      uint32_t class_sizes[13] = { 6, 3, 90, 130, 105, 35, 9330, 23436, 48708, 78474, 91377, 70785, 32175 };
      num_combs = class_sizes[best_multiplicity - 3];
      support_minimization_encodings = std::vector<std::array<uint32_t, 2>>( num_combs );
      generate_support_minimization_encodings_rec<true, false>( 0, 0, 0, count, min_set_size, true );
    }

    assert( count == num_combs );
  }

  template<bool enable_dcset, bool equal_size_partition>
  void generate_support_minimization_encodings_rec( uint32_t onset, uint32_t offset, uint32_t var, uint32_t& count, uint32_t min_set_size, bool first )
  {
    if ( var == best_multiplicity )
    {
      if ( equal_size_partition )
      {
        /* sets must be equally populated */
        if ( __builtin_popcount( onset ) != __builtin_popcount( offset ) )
        {
          return;
        }
      }
      else if ( __builtin_popcount( onset ) < min_set_size || __builtin_popcount( offset ) < min_set_size )
      {
        /* ON-set and OFF-set must be populated with at least min_set_size elements */
        return;
      }

      support_minimization_encodings[count][0] = onset;
      support_minimization_encodings[count][1] = offset;
      ++count;
      return;
    }

    /* var in DCSET */
    if ( enable_dcset )
    {
      generate_support_minimization_encodings_rec<enable_dcset, equal_size_partition>( onset, offset, var + 1, count, min_set_size, first );
    }

    /* move var in ONSET */
    onset |= 1 << var;
    generate_support_minimization_encodings_rec<enable_dcset, equal_size_partition>( onset, offset, var + 1, count, min_set_size, false );
    onset &= ~( 1 << var );

    /* remove symmetries */
    if ( first )
    {
      return;
    }

    /* move var in OFFSET */
    offset |= 1 << var;
    generate_support_minimization_encodings_rec<enable_dcset, equal_size_partition>( onset, offset, var + 1, count, min_set_size, false );
    offset &= ~( 1 << var );
  }

  void solve_min_support_exact( std::vector<STT> const& isets )
  {
    std::vector<encoding_column> matrix;
    matrix.reserve( support_minimization_encodings.size() );
    best_bound_sets.clear();

    /* create covering matrix */
    if ( !create_covering_matrix<false>( isets, matrix, false ) )
    {
      return;
    }

    /* solve the covering problem */
    std::array<uint32_t, 6> solution = covering_solve_exact( matrix );

    /* check for failed decomposition */
    if ( solution[0] == UINT32_MAX )
    {
      return;
    }

    /* compute best bound sets */
    uint32_t num_luts = 1 + solution[5];
    uint32_t num_levels = 2;
    uint32_t num_edges = best_free_set + solution[5];
    uint32_t isets_support = num_vars - best_free_set;
    best_care_sets.clear();
    best_iset_onset.clear();
    best_iset_offset.clear();
    for ( uint32_t i = 0; i < solution[5]; ++i )
    {
      STT tt;
      STT care;

      const uint32_t onset = support_minimization_encodings[matrix[solution[i]].index][0];
      const uint32_t offset = support_minimization_encodings[matrix[solution[i]].index][1];
      for ( uint32_t j = 0; j < best_multiplicity; ++j )
      {
        if ( ( ( onset >> j ) & 1 ) )
        {
          tt |= isets[j];
        }
        if ( ( ( offset >> j ) & 1 ) )
        {
          care |= isets[j];
        }
      }

      care |= tt;
      num_edges += matrix[solution[i]].cost & ( ( 1 << isets_support ) - 1 );

      best_bound_sets.push_back( tt );
      best_care_sets.push_back( care );
      best_iset_onset.push_back( onset );
      best_iset_offset.push_back( offset );
    }

    if ( pst )
    {
      pst->num_luts = num_luts;
      pst->num_levels = num_levels;
      pst->num_edges = num_edges;
    }
  }

  void solve_min_support_heuristic( std::vector<STT> const& isets )
  {
    std::vector<encoding_column> matrix;
    matrix.reserve( support_minimization_encodings.size() );
    best_bound_sets.clear();

    /* create covering matrix */
    if ( !create_covering_matrix<true>( isets, matrix, true ) )
    {
      return;
    }

    /* solve the covering problem: heuristic pass + local search */
    std::array<uint32_t, 6> solution = covering_solve_heuristic( matrix );

    /* check for failed decomposition */
    if ( solution[0] == UINT32_MAX )
    {
      return;
    }

    /* improve solution with local search */
    while ( covering_improve( matrix, solution ) )
      ;

    /* try to increase the encoding length */
    // if ( covering_increase_length( matrix, solution ) )
    //   std::cout << "[i] stat\n";

    /* compute best bound sets */
    uint32_t num_luts = 1 + solution[5];
    uint32_t num_levels = 2;
    uint32_t num_edges = best_free_set + solution[5];
    uint32_t isets_support = num_vars - best_free_set;
    best_care_sets.clear();
    best_iset_onset.clear();
    best_iset_offset.clear();
    for ( uint32_t i = 0; i < solution[5]; ++i )
    {
      STT tt;
      STT care;

      const uint32_t onset = support_minimization_encodings[matrix[solution[i]].index][0];
      const uint32_t offset = support_minimization_encodings[matrix[solution[i]].index][1];
      for ( uint32_t j = 0; j < best_multiplicity; ++j )
      {
        if ( ( ( onset >> j ) & 1 ) )
        {
          tt |= isets[j];
        }
        if ( ( ( offset >> j ) & 1 ) )
        {
          care |= isets[j];
        }
      }

      care |= tt;
      num_edges += matrix[solution[i]].cost & ( ( 1 << isets_support ) - 1 );

      best_bound_sets.push_back( tt );
      best_care_sets.push_back( care );
      best_iset_onset.push_back( onset );
      best_iset_offset.push_back( offset );
    }

    if ( pst )
    {
      pst->num_luts = num_luts;
      pst->num_levels = num_levels;
      pst->num_edges = num_edges;
    }
  }

  template<bool UseHeuristic>
  bool create_covering_matrix( std::vector<STT> const& isets, std::vector<encoding_column>& matrix, bool sort )
  {
    assert( best_multiplicity <= 16 );
    uint32_t combinations = ( best_multiplicity * ( best_multiplicity - 1 ) ) / 2;
    uint32_t iset_support = num_vars - best_free_set;

    /* insert dichotomies */
    for ( uint32_t i = 0; i < support_minimization_encodings.size(); ++i )
    {
      uint32_t const onset = support_minimization_encodings[i][0];
      uint32_t const offset = support_minimization_encodings[i][1];

      /* compute function and distinguishable seed dichotomies */
      uint64_t column[2] = { 0, 0 };
      STT tt;
      STT care;
      uint32_t pair_pointer = 0;
      for ( uint32_t j = 0; j < best_multiplicity; ++j )
      {
        auto onset_shift = ( onset >> j );
        auto offset_shift = ( offset >> j );
        if ( ( onset_shift & 1 ) )
        {
          tt |= isets[j];
        }

        if ( ( offset_shift & 1 ) )
        {
          care |= isets[j];
        }

        /* compute included seed dichotomies */
        for ( uint32_t k = j + 1; k < best_multiplicity; ++k )
        {
          /* if are in diffent sets */
          if ( ( ( ( onset_shift & ( offset >> k ) ) | ( ( onset >> k ) & offset_shift ) ) & 1 ) )
          {
            column[pair_pointer >> 6u] |= UINT64_C( 1 ) << ( pair_pointer & 0x3F );
          }

          ++pair_pointer;
        }
      }

      care |= tt;

      /* compute cost */
      uint32_t cost = 0;
      for ( uint32_t j = 0; j < iset_support; ++j )
      {
        cost += has_var_support( tt, care, iset_support, j ) ? 1 : 0;
        // if ( !has_var_support( tt, care, iset_support, j ) )
        // {
        //   /* adjust truth table and care set */
        //   adjust_truth_table_on_dc_support( tt, care, iset_support, j );
        //   continue;
        // }
        // ++cost;
      }

      /* discard solutions with support over LUT size */
      if ( cost > ps.lut_size )
        continue;

      /* buffers have zero cost */
      if ( cost == 1 )
        cost = 0;

      float sort_cost = 0;
      if constexpr ( UseHeuristic )
      {
        sort_cost = 1.0f / ( __builtin_popcountl( column[0] ) + __builtin_popcountl( column[1] ) );
      }
      else
      {
        sort_cost = cost + ( ( combinations - __builtin_popcountl( column[0] + __builtin_popcountl( column[1] ) ) ) << num_vars );
      }

      /* insert */
      matrix.emplace_back( encoding_column{ { column[0], column[1] }, cost, i, sort_cost } );
    }

    if ( !sort )
    {
      return true;
    }

    if constexpr ( UseHeuristic )
    {
      std::sort( matrix.begin(), matrix.end(), [&]( auto const& a, auto const& b ) {
        return a.cost < b.cost;
      } );
    }
    else
    {
      std::sort( matrix.begin(), matrix.end(), [&]( auto const& a, auto const& b ) {
        return a.sort_cost < b.sort_cost;
      } );
    }

    return true;
  }

  std::array<uint32_t, 6> covering_solve_exact( std::vector<encoding_column>& matrix )
  {
    /* last value of res contains the size of the bound set */
    std::array<uint32_t, 6> res = { UINT32_MAX };
    uint32_t best_cost = UINT32_MAX;
    uint32_t combinations = ( best_multiplicity * ( best_multiplicity - 1 ) ) / 2;

    assert( best_multiplicity <= 4 );

    /* determine the number of needed loops*/
    if ( best_multiplicity <= 2 )
    {
      res[5] = 1;
      res[0] = 0;
    }
    else if ( best_multiplicity <= 4 )
    {
      res[5] = 2;
      for ( uint32_t i = 0; i < matrix.size() - 1; ++i )
      {
        for ( uint32_t j = 1; j < matrix.size(); ++j )
        {
          /* filter by cost */
          if ( matrix[i].cost + matrix[j].cost >= best_cost )
            continue;

          /* check validity */
          if ( __builtin_popcountl( matrix[i].column[0] | matrix[j].column[0] ) + __builtin_popcountl( matrix[i].column[1] | matrix[j].column[1] ) == combinations )
          {
            res[0] = i;
            res[1] = j;
            best_cost = matrix[i].cost + matrix[j].cost;
          }
        }
      }
    }

    return res;
  }

  std::array<uint32_t, 6> covering_solve_heuristic( std::vector<encoding_column>& matrix )
  {
    /* last value of res contains the size of the bound set */
    std::array<uint32_t, 6> res = { UINT32_MAX };
    uint32_t combinations = ( best_multiplicity * ( best_multiplicity - 1 ) ) / 2;
    uint64_t column0 = 0, column1 = 0;

    uint32_t best = 0;
    float best_cost = std::numeric_limits<float>::max();
    for ( uint32_t i = 0; i < matrix.size(); ++i )
    {
      if ( matrix[i].sort_cost < best_cost )
      {
        best = i;
        best_cost = matrix[i].sort_cost;
      }
    }

    /* select */
    column0 = matrix[best].column[0];
    column1 = matrix[best].column[1];
    std::swap( matrix[0], matrix[best] );

    /* get max number of BS's */
    uint32_t iter = 1;

    while ( iter < ps.lut_size - best_free_set && __builtin_popcountl( column0 ) + __builtin_popcountl( column1 ) != combinations )
    {
      /* select column that minimizes the cost */
      best = 0;
      best_cost = std::numeric_limits<float>::max();
      for ( uint32_t i = iter; i < matrix.size(); ++i )
      {
        float local_cost = 1.0f / ( __builtin_popcountl( matrix[i].column[0] & ~column0 ) + __builtin_popcountl( matrix[i].column[1] & ~column1 ) );
        if ( local_cost < best_cost )
        {
          best = i;
          best_cost = local_cost;
        }
      }

      column0 |= matrix[best].column[0];
      column1 |= matrix[best].column[1];
      std::swap( matrix[iter], matrix[best] );
      ++iter;
    }

    if ( __builtin_popcountl( column0 ) + __builtin_popcountl( column1 ) == combinations )
    {
      for ( uint32_t i = 0; i < iter; ++i )
      {
        res[i] = i;
      }
      res[5] = iter;
    }

    return res;
  }

  bool covering_improve( std::vector<encoding_column> const& matrix, std::array<uint32_t, 6>& solution )
  {
    /* performs one iteration of local search */
    uint32_t best_cost = 0, local_cost = 0;
    uint32_t num_elements = solution[5];
    uint32_t combinations = ( best_multiplicity * ( best_multiplicity - 1 ) ) / 2;
    bool improved = false;

    /* compute current cost */
    for ( uint32_t i = 0; i < num_elements; ++i )
    {
      best_cost += matrix[solution[i]].cost;
    }

    uint64_t column0, column1;
    for ( uint32_t i = 0; i < num_elements; ++i )
    {
      /* remove element i */
      local_cost = 0;
      column0 = 0;
      column1 = 0;
      for ( uint32_t j = 0; j < num_elements; ++j )
      {
        if ( j == i )
          continue;
        local_cost += matrix[solution[j]].cost;
        column0 |= matrix[solution[j]].column[0];
        column1 |= matrix[solution[j]].column[1];
      }

      /* search for a better replecemnts */
      for ( uint32_t j = 0; j < matrix.size(); ++j )
      {
        if ( __builtin_popcount( column0 | matrix[j].column[0] ) + __builtin_popcount( column1 | matrix[j].column[1] ) != combinations )
          continue;
        if ( local_cost + matrix[j].cost < best_cost )
        {
          solution[i] = j;
          best_cost = local_cost + matrix[j].cost;
          improved = true;
        }
      }
    }

    return improved;
  }

  bool covering_increase_length( std::vector<encoding_column> const& matrix, std::array<uint32_t, 6>& solution )
  {
    /* count the number of solutions with cost > 0 */
    uint32_t best_cost = 0, local_cost = 0;
    uint32_t replaceable_columns = 0;
    uint32_t num_elements = solution[5];

    for ( uint32_t i = 0; i < num_elements; ++i )
    {
      replaceable_columns += !!matrix[solution[i]].cost;
    }

    /* not possible to increase the encoding length */
    if ( replaceable_columns == 1 || ps.lut_size - best_free_set - num_elements == 0 )
      return false;

    /* extract shared set BS candidates */
    std::vector<uint32_t> ss_candidates;
    for ( uint32_t i = 0; i < matrix.size(); ++i )
    {
      if ( matrix[i].cost == 0 )
        ss_candidates.push_back( i );
    }

    // if ( ss_candidates.size() < 2 )
    //   return false;

    /* find a column and try to replace it using two shared set variables */
    uint64_t column0, column1;
    uint32_t combinations = ( best_multiplicity * ( best_multiplicity - 1 ) ) / 2;
    bool improved = false;

    uint32_t i = 0;
    while ( replaceable_columns > 1 && ps.lut_size - best_free_set - num_elements > 0 )
    {
      /* find a column */
      while ( i < num_elements && matrix[solution[i]].cost == 0 )
        ++i;

      if ( i == num_elements )
        break;

      /* remove the column */
      column0 = 0;
      column1 = 0;
      for ( uint32_t j = 0; j < num_elements; ++j )
      {
        if ( j == i )
          continue;

        column0 |= matrix[solution[j]].column[0];
        column1 |= matrix[solution[j]].column[1];
      }

      /* find two buffers */
      if ( covering_increase_length_find_ss( matrix, solution, ss_candidates, column0, column1, i ) )
      {
        ++num_elements;
        --replaceable_columns;
        improved = true;
      }

      ++i;
    }

    return improved;
  }

  inline bool covering_increase_length_find_ss( std::vector<encoding_column> const& matrix, std::array<uint32_t, 6>& solution, std::vector<uint32_t> const& ss_candidates, uint64_t column0, uint64_t column1, uint32_t i )
  {
    uint32_t combinations = ( best_multiplicity * ( best_multiplicity - 1 ) ) / 2;

    for ( uint32_t j = 0; j < ss_candidates.size() - 1; ++j )
    {
      for ( uint32_t k = j + 1; k < ss_candidates.size(); ++k )
      {
        uint64_t test_column0 = matrix[ss_candidates[j]].column[0] | matrix[ss_candidates[k]].column[0];
        uint64_t test_column1 = matrix[ss_candidates[j]].column[1] | matrix[ss_candidates[k]].column[1];

        if ( __builtin_popcount( column0 | test_column0 ) + __builtin_popcount( column1 | test_column1 ) != combinations )
          continue;

        /* perform the substitution */
        solution[i] = ss_candidates[j];
        solution[solution[5]] = ss_candidates[k];
        ++solution[5];
        return true;
      }
    }

    return false;
  }

  void adjust_truth_table_on_dc( STT& tt, STT& care, uint32_t var_index )
  {
    assert( var_index < tt.num_vars() );
    assert( tt.num_vars() == care.num_vars() );

    if ( tt.num_vars() <= 6 || var_index < 6 )
    {
      auto it_tt = std::begin( tt._bits );
      auto it_care = std::begin( care._bits );
      while ( it_tt != std::end( tt._bits ) )
      {
        uint64_t new_bits = *it_tt & *it_care;
        *it_tt = ( ( new_bits | ( new_bits >> ( uint64_t( 1 ) << var_index ) ) ) & kitty::detail::projections_neg[var_index] ) |
                 ( ( new_bits | ( new_bits << ( uint64_t( 1 ) << var_index ) ) ) & kitty::detail::projections[var_index] );
        *it_care = ( *it_care | ( *it_care >> ( uint64_t( 1 ) << var_index ) ) ) & kitty::detail::projections_neg[var_index];
        *it_care = *it_care | ( *it_care << ( uint64_t( 1 ) << var_index ) );

        ++it_tt;
        ++it_care;
      }
      return;
    }

    const auto step = 1 << ( var_index - 6 );
    for ( auto i = 0u; i < static_cast<uint32_t>( tt.num_blocks() ); i += 2 * step )
    {
      for ( auto j = 0; j < step; ++j )
      {
        tt._bits[i + j] = ( tt._bits[i + j] & care._bits[i + j] ) | ( tt._bits[i + j + step] & care._bits[i + j + step] );
        tt._bits[i + j + step] = tt._bits[i + j];
        care._bits[i + j] = care._bits[i + j] | care._bits[i + j + step];
        care._bits[i + j + step] = care._bits[i + j];
      }
    }
  }

  void local_extend_to( STT& tt, uint32_t real_num_vars )
  {
    if ( real_num_vars < 6 )
    {
      auto mask = *tt.begin();

      for ( auto i = real_num_vars; i < num_vars; ++i )
      {
        mask |= ( mask << ( 1 << i ) );
      }

      std::fill( tt.begin(), tt.end(), mask );
    }
    else
    {
      uint32_t num_blocks = ( 1u << ( real_num_vars - 6 ) );
      auto it = tt.begin();
      while ( it != tt.end() )
      {
        it = std::copy( tt.cbegin(), tt.cbegin() + num_blocks, it );
      }
    }
  }

  bool has_var_support( const STT& tt, const STT& care, uint32_t real_num_vars, uint8_t var_index )
  {
    assert( var_index < real_num_vars );
    assert( real_num_vars <= tt.num_vars() );
    assert( tt.num_vars() == care.num_vars() );

    const uint32_t num_blocks = real_num_vars <= 6 ? 1 : ( 1 << ( real_num_vars - 6 ) );
    if ( real_num_vars <= 6 || var_index < 6 )
    {
      auto it_tt = std::begin( tt._bits );
      auto it_care = std::begin( care._bits );
      while ( it_tt != std::begin( tt._bits ) + num_blocks )
      {
        if ( ( ( ( *it_tt >> ( uint64_t( 1 ) << var_index ) ) ^ *it_tt ) & kitty::detail::projections_neg[var_index] & ( *it_care >> ( uint64_t( 1 ) << var_index ) ) & *it_care ) != 0 )
        {
          return true;
        }
        ++it_tt;
        ++it_care;
      }

      return false;
    }

    const auto step = 1 << ( var_index - 6 );
    for ( auto i = 0u; i < num_blocks; i += 2 * step )
    {
      for ( auto j = 0; j < step; ++j )
      {
        if ( ( ( tt._bits[i + j] ^ tt._bits[i + j + step] ) & care._bits[i + j] & care._bits[i + j + step] ) != 0 )
        {
          return true;
        }
      }
    }

    return false;
  }

  void adjust_truth_table_on_dc_support( STT& tt, STT& care, uint32_t real_num_vars, uint32_t var_index )
  {
    assert( var_index < real_num_vars );
    assert( tt.num_vars() == care.num_vars() );

    const uint32_t num_blocks = real_num_vars <= 6 ? 1 : ( 1 << ( real_num_vars - 6 ) );
    if ( real_num_vars <= 6 || var_index < 6 )
    {
      auto it_tt = std::begin( tt._bits );
      auto it_care = std::begin( care._bits );
      while ( it_tt != std::begin( tt._bits ) + num_blocks )
      {
        uint64_t new_bits = *it_tt & *it_care;
        *it_tt = ( ( new_bits | ( new_bits >> ( uint64_t( 1 ) << var_index ) ) ) & kitty::detail::projections_neg[var_index] ) |
                 ( ( new_bits | ( new_bits << ( uint64_t( 1 ) << var_index ) ) ) & kitty::detail::projections[var_index] );
        *it_care = ( *it_care | ( *it_care >> ( uint64_t( 1 ) << var_index ) ) ) & kitty::detail::projections_neg[var_index];
        *it_care = *it_care | ( *it_care << ( uint64_t( 1 ) << var_index ) );

        ++it_tt;
        ++it_care;
      }
      return;
    }

    const auto step = 1 << ( var_index - 6 );
    for ( auto i = 0u; i < static_cast<uint32_t>( num_blocks ); i += 2 * step )
    {
      for ( auto j = 0; j < step; ++j )
      {
        tt._bits[i + j] = ( tt._bits[i + j] & care._bits[i + j] ) | ( tt._bits[i + j + step] & care._bits[i + j + step] );
        tt._bits[i + j + step] = tt._bits[i + j];
        care._bits[i + j] = care._bits[i + j] | care._bits[i + j + step];
        care._bits[i + j + step] = care._bits[i + j];
      }
    }
  }

  /* Decomposition format for ABC
   *
   * The record is an array of unsigned chars where:
   *   - the first unsigned char entry stores the number of unsigned chars in the record
   *   - the second entry stores the number of LUTs
   * After this, several sub-records follow, each representing one LUT as follows:
   *   - an unsigned char entry listing the number of fanins
   *   - a list of fanins, from the LSB to the MSB of the truth table. The N inputs of the original function
   *     have indexes from 0 to N-1, followed by the internal signals in a topological order
   *   - the LUT truth table occupying 2^(M-3) bytes, where M is the fanin count of the LUT, from the LSB to the MSB.
   *     A 2-input LUT, which takes 4 bits, should be stretched to occupy 8 bits (one unsigned char)
   *     A 0- or 1-input LUT can be represented similarly but it is not expected that such LUTs will be represented
   */
  void get_decomposition_abc( unsigned char* decompArray )
  {
    unsigned char* pArray = decompArray;
    unsigned char bytes = 2;

    /* write number of LUTs */
    pArray++;
    *pArray = dec_result.size();
    pArray++;

    /* write LUTs */
    for ( ac_dc_decomposition_result const& lut : dec_result )
    {
      /* write fanin size*/
      *pArray = lut.support.size();
      pArray++;
      ++bytes;

      /* write support */
      for ( uint32_t i : lut.support )
      {
        *pArray = (unsigned char)i;
        pArray++;
        ++bytes;
      }

      /* write truth table */
      uint32_t tt_num_bytes = ( lut.tt.num_vars() <= 3 ) ? 1 : ( 1 << ( lut.tt.num_vars() - 3 ) );
      tt_num_bytes = std::min( tt_num_bytes, 8u );
      for ( uint32_t i = 0; i < lut.tt.num_blocks(); ++i )
      {
        for ( uint32_t j = 0; j < tt_num_bytes; ++j )
        {
          *pArray = (unsigned char)( ( lut.tt._bits[i] >> ( 8 * j ) ) & 0xFF );
          pArray++;
          ++bytes;
        }
      }
    }

    /* write numBytes */
    *decompArray = bytes;
  }

private:
  uint32_t best_multiplicity{ UINT32_MAX };
  uint32_t best_free_set{ UINT32_MAX };
  STT best_tt;
  STT best_tt_cs;
  STT best_cs; /* care set permuted like best_tt (set by find_decomposition) */
  std::vector<STT> best_bound_sets;
  std::vector<STT> best_care_sets;
  std::vector<STT> best_free_set_tts;
  std::vector<uint64_t> best_iset_onset;
  std::vector<uint64_t> best_iset_offset;
  std::vector<ac_dc_decomposition_result> dec_result;

  std::vector<std::array<uint32_t, 2>> support_minimization_encodings;

  uint32_t num_vars;
  ac_dc_decomposition_params ps;
  ac_dc_decomposition_stats* pst;
  std::array<uint32_t, max_num_vars> permutations;
};

} // namespace acd

ABC_NAMESPACE_CXX_HEADER_END

#endif // _ACD_DC_H_
