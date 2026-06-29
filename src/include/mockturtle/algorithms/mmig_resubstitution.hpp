/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2022  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file mmig_resubstitution.hpp
  \brief mMIG-oriented resubstitution wrapper with selective candidate CEC
*/

#pragma once

#include "cleanup.hpp"
#include "equivalence_checking.hpp"
#include "mig_resub.hpp"
#include "miter.hpp"

#include "../traits.hpp"
#include "../utils/stopwatch.hpp"
#include "../views/depth_view.hpp"
#include "../views/fanout_view.hpp"

#include <cstdint>
#include <optional>
#include <type_traits>

namespace mockturtle
{

struct mmig_resubstitution_params
{
  mmig_resubstitution_params()
  {
    resub_ps.max_pis = 8u;
    resub_ps.max_divisors = 128u;
    resub_ps.max_inserts = 2u;
    resub_ps.use_dont_cares = false;
    resub_ps.preserve_depth = true;
    resub_ps.window_size = 8u;
    resub_ps.verbose = false;
  }

  resubstitution_params resub_ps{};
  bool enable_candidate_cec{ false };
  bool verify_mixed_risk_only{ true };
  uint32_t max_candidate_cec_checks{ 0u };
  bool skip_candidate_when_cec_budget_exhausted{ true };
  bool verbose{ false };
};

struct mmig_resubstitution_stats
{
  stopwatch<>::duration time_total{ 0 };
  resubstitution_stats resub{};
  uint32_t num_candidate_checks{ 0 };
  uint32_t num_candidate_accepted{ 0 };
  uint32_t num_candidate_rejected{ 0 };
  uint32_t num_candidate_inconclusive{ 0 };
  uint32_t num_candidate_skipped_non_mixed{ 0 };
  uint32_t num_candidate_skipped_budget{ 0 };
  stopwatch<>::duration time_candidate_cec{ 0 };
};

namespace detail
{

template<class Ntk>
bool is_mixed_risk_root( Ntk const& ntk, node<Ntk> const& n )
{
  if ( ntk.is_min( n ) )
  {
    return true;
  }

  bool mixed = false;

  ntk.foreach_fanin( n, [&]( auto const& f ) {
    auto const child = ntk.get_node( f );
    if ( ntk.is_constant( child ) || ntk.is_pi( child ) )
    {
      return;
    }
    if ( ntk.is_min( child ) )
    {
      mixed = true;
    }
  } );

  if ( mixed )
  {
    return true;
  }

  ntk.foreach_fanout( n, [&]( auto const& parent ) {
    if ( ntk.is_min( parent ) )
    {
      mixed = true;
      return false;
    }
    return true;
  } );

  return mixed;
}

template<class Ntk>
bool try_candidate_substitution( Ntk& ntk,
                                 node<Ntk> const& n,
                                 signal<Ntk> const& g,
                                 mmig_resubstitution_params const& ps,
                                 mmig_resubstitution_stats& st )
{
  if ( !ps.enable_candidate_cec )
  {
    ntk.substitute_node( n, g );
    return true;
  }

  if ( ps.verify_mixed_risk_only && !is_mixed_risk_root( ntk, n ) )
  {
    ++st.num_candidate_skipped_non_mixed;
    ntk.substitute_node( n, g );
    return true;
  }

  if ( ps.max_candidate_cec_checks > 0u && st.num_candidate_checks >= ps.max_candidate_cec_checks )
  {
    ++st.num_candidate_skipped_budget;
    if ( ps.skip_candidate_when_cec_budget_exhausted )
    {
      return false;
    }
    ntk.substitute_node( n, g );
    return true;
  }

  ++st.num_candidate_checks;

  auto const current = ntk.clone();
  auto candidate = current.clone();
  candidate.substitute_node( n, g );
  candidate = cleanup_dangling<mig_network>( candidate );

  auto const maybe_miter = miter<mig_network>( current, candidate );
  if ( !maybe_miter.has_value() )
  {
    ++st.num_candidate_inconclusive;
    ++st.num_candidate_rejected;
    return false;
  }

  std::optional<bool> equivalent;
  {
    stopwatch t( st.time_candidate_cec );
    equivalent = equivalence_checking( *maybe_miter );
  }

  if ( !equivalent.has_value() )
  {
    ++st.num_candidate_inconclusive;
    ++st.num_candidate_rejected;
    return false;
  }

  if ( *equivalent )
  {
    ntk.substitute_node( n, g );
    ++st.num_candidate_accepted;
    return true;
  }

  ++st.num_candidate_rejected;
  return false;
}

template<class Ntk>
void run_mmig_resub_impl( Ntk& ntk,
                          resubstitution_params const& resub_ps,
                          mmig_resubstitution_params const& ps,
                          mmig_resubstitution_stats& st )
{
  auto callback = [&]( Ntk& net, node<Ntk> const& n, signal<Ntk> const& g ) {
    return try_candidate_substitution( net, n, g, ps, st );
  };

  if ( resub_ps.max_pis == 8u )
  {
    using truthtable_t = kitty::static_truth_table<8u>;
    using truthtable_dc_t = kitty::dynamic_truth_table;
    using functor_t = mig_enumerative_resub_functor<Ntk, detail::window_simulator<Ntk, truthtable_t>, truthtable_dc_t>;
    using resub_impl_t = detail::resubstitution_impl<Ntk, detail::window_based_resub_engine<Ntk, truthtable_t, truthtable_dc_t, functor_t>>;

    typename resub_impl_t::engine_st_t engine_st;
    typename resub_impl_t::collector_st_t collector_st;
    resub_impl_t p( ntk, resub_ps, st.resub, engine_st, collector_st );
    p.run( callback );

    if ( resub_ps.verbose )
    {
      st.resub.report();
      collector_st.report();
      engine_st.report();
    }
    return;
  }

  using truthtable_t = kitty::dynamic_truth_table;
  using truthtable_dc_t = kitty::dynamic_truth_table;
  using functor_t = mig_enumerative_resub_functor<Ntk, detail::window_simulator<Ntk, truthtable_t>, truthtable_dc_t>;
  using resub_impl_t = detail::resubstitution_impl<Ntk, detail::window_based_resub_engine<Ntk, truthtable_t, truthtable_dc_t, functor_t>>;

  typename resub_impl_t::engine_st_t engine_st;
  typename resub_impl_t::collector_st_t collector_st;
  resub_impl_t p( ntk, resub_ps, st.resub, engine_st, collector_st );
  p.run( callback );

  if ( resub_ps.verbose )
  {
    st.resub.report();
    collector_st.report();
    engine_st.report();
  }
}

} // namespace detail

template<class Ntk>
void mmig_resubstitution( Ntk& ntk,
                          mmig_resubstitution_params const& ps = {},
                          mmig_resubstitution_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( std::is_same_v<typename Ntk::base_type, mig_network>, "Network type is not MIG-based" );
  static_assert( has_is_min_v<Ntk>, "Ntk does not implement is_min" );
  static_assert( has_create_min_v<Ntk>, "Ntk does not implement create_min" );

  mmig_resubstitution_stats st{};
  {
    stopwatch t( st.time_total );
    auto resub_ps = ps.resub_ps;
    if ( ps.verbose )
    {
      resub_ps.verbose = true;
    }
    resub_ps.enable_min_resubstitution = true;

    /* Dynamic-table resubstitution with don't-cares at large cuts is unstable in
     * some networks. Keep don't-cares on the proven-safe 8-input setting. */
    if ( resub_ps.use_dont_cares && resub_ps.max_pis > 8u )
    {
      resub_ps.use_dont_cares = false;
    }

    fanout_view<Ntk> fntk{ ntk };
    depth_view<fanout_view<Ntk>> dfntk{ fntk };
    detail::run_mmig_resub_impl( dfntk, resub_ps, ps, st );
    ntk = cleanup_dangling( ntk );
  }

  if ( pst != nullptr )
  {
    *pst = st;
  }
}

} // namespace mockturtle
