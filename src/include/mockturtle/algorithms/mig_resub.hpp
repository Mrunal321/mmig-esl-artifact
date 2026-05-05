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
  \file mig_resub.hpp
  \brief Majority-specific resustitution rules

  \author Eleonora Testa
  \author Heinz Riener
  \author Mathias Soeken
  \author Siang-Yun (Sonia) Lee
*/

#pragma once

#include "cleanup.hpp"
#include "../networks/mig.hpp"
#include "../utils/index_list/index_list.hpp"
#include "../utils/truth_table_utils.hpp"
#include "resubstitution.hpp"
#include "resyn_engines/mig_resyn.hpp"

#include <kitty/kitty.hpp>

namespace mockturtle
{

struct mig_enumerative_resub_stats
{
  /*! \brief Accumulated runtime for const-resub */
  stopwatch<>::duration time_resubC{ 0 };

  /*! \brief Accumulated runtime for zero-resub */
  stopwatch<>::duration time_resub0{ 0 };

  /*! \brief Accumulated runtime for collecting unate divisors. */
  stopwatch<>::duration time_collect_unate_divisors{ 0 };

  /*! \brief Accumulated runtime for one-resub */
  stopwatch<>::duration time_resub1{ 0 };

  /*! \brief Accumulated runtime for relevance resub */
  stopwatch<>::duration time_resubR{ 0 };

  /*! \brief Accumulated runtime for collecting unate divisors. */
  stopwatch<>::duration time_collect_binate_divisors{ 0 };

  /*! \brief Accumulated runtime for two-resub. */
  stopwatch<>::duration time_resub2{ 0 };

  /*! \brief Number of accepted constant resubsitutions */
  uint32_t num_const_accepts{ 0 };

  /*! \brief Number of accepted zero resubsitutions */
  uint32_t num_div0_accepts{ 0 };

  /*! \brief Number of accepted one resubsitutions */
  uint64_t num_div1_accepts{ 0 };

  /*! \brief Number of accepted relevance resubsitutions */
  uint32_t num_divR_accepts{ 0 };

  /*! \brief Number of accepted two resubsitutions */
  uint64_t num_div2_accepts{ 0 };

  void report() const
  {
    std::cout << "[i] kernel: mig_enumerative_resub_functor\n";
    std::cout << fmt::format( "[i]     constant-resub {:6d}                                   ({:>5.2f} secs)\n",
                              num_const_accepts, to_seconds( time_resubC ) );
    std::cout << fmt::format( "[i]            0-resub {:6d}                                   ({:>5.2f} secs)\n",
                              num_div0_accepts, to_seconds( time_resub0 ) );
    std::cout << fmt::format( "[i]            R-resub {:6d}                                   ({:>5.2f} secs)\n",
                              num_divR_accepts, to_seconds( time_resubR ) );
    std::cout << fmt::format( "[i]            collect unate divisors                           ({:>5.2f} secs)\n", to_seconds( time_collect_unate_divisors ) );
    std::cout << fmt::format( "[i]            1-resub {:6d} = {:6d} MAJ/MIN                  ({:>5.2f} secs)\n",
                              num_div1_accepts, num_div1_accepts, to_seconds( time_resub1 ) );
    std::cout << fmt::format( "[i]            collect binate divisors                          ({:>5.2f} secs)\n", to_seconds( time_collect_binate_divisors ) );
    std::cout << fmt::format( "[i]            2-resub {:6d} = {:6d} 2MAJ/MIN                 ({:>5.2f} secs)\n",
                              num_div2_accepts, num_div2_accepts, to_seconds( time_resub2 ) );
    std::cout << fmt::format( "[i]            total   {:6d}\n",
                              ( num_const_accepts + num_div0_accepts + num_divR_accepts + num_div1_accepts + num_div2_accepts ) );
  }
}; /* mig_enumerative_resub_stats */

template<typename Ntk, typename Simulator, typename TT, bool use_constant = true>
struct mig_enumerative_resub_functor
{
public:
  using node = mig_network::node;
  using signal = mig_network::signal;
  using stats = mig_enumerative_resub_stats;

  struct unate_divisors
  {
    std::vector<signal> u0;
    std::vector<signal> u1;
    std::vector<signal> next_candidates;

    void clear()
    {
      u0.clear();
      u1.clear();
      next_candidates.clear();
    }
  };

  struct binate_divisors
  {
    std::vector<signal> b0;
    std::vector<signal> b1;
    std::vector<signal> b2;

    void clear()
    {
      b0.clear();
      b1.clear();
      b2.clear();
    }
  };

  enum class candidate_stage : uint8_t
  {
    constant = 0u,
    div0,
    divR,
    div1,
    div2
  };

  enum class candidate_kind : uint8_t
  {
    direct = 0u,
    maj3,
    maj5
  };

  struct resub_candidate
  {
    candidate_stage stage{ candidate_stage::constant };
    candidate_kind kind{ candidate_kind::direct };
    signal direct{};
    std::array<signal, 3u> outer{};
    std::array<signal, 3u> inner{};
    bool outer_is_min{ false };
    bool inner_is_min{ false };
    bool complement_inner_in_outer{ false };
    bool complement_output{ false };
    uint32_t gain{ 0u };
    uint32_t inserts{ 0u };
    uint32_t discovery_index{ 0u };
  };

  struct cheap_candidate_score
  {
    int32_t level_gain{ std::numeric_limits<int32_t>::min() };
    int32_t estimated_level{ std::numeric_limits<int32_t>::max() };
    int32_t inv_cost{ std::numeric_limits<int32_t>::max() };
    int32_t duplication_penalty{ std::numeric_limits<int32_t>::max() };
    int32_t share_bonus{ std::numeric_limits<int32_t>::min() };
    uint32_t gain{ 0u };
    uint32_t discovery_index{ 0u };
  };

  struct ranked_candidate
  {
    resub_candidate candidate;
    cheap_candidate_score cheap;
  };

  struct full_candidate_score
  {
    ranked_candidate ranked;
    int32_t depth_gain{ std::numeric_limits<int32_t>::min() };
    int32_t gate_gain{ std::numeric_limits<int32_t>::min() };
    int32_t inv_gain{ std::numeric_limits<int32_t>::min() };
  };

  struct network_cost
  {
    int32_t depth{ 0 };
    int32_t gates{ 0 };
    int32_t inverted_edges{ 0 };
  };

public:
  explicit mig_enumerative_resub_functor( Ntk& ntk, Simulator const& sim, std::vector<node> const& divs, uint32_t num_divs, stats& st )
      : ntk( ntk ), sim( sim ), divs( divs ), num_divs( num_divs ), st( st ), ps( default_ps )
  {
  }

  explicit mig_enumerative_resub_functor( Ntk& ntk, Simulator const& sim, std::vector<node> const& divs, uint32_t num_divs, resubstitution_params const& ps, stats& st )
      : ntk( ntk ), sim( sim ), divs( divs ), num_divs( num_divs ), st( st ), ps( ps )
  {
  }

  std::optional<signal> operator()( node const& root, TT care, uint32_t required, uint32_t max_inserts, uint32_t num_mffc, uint32_t& last_gain )
  {
    if ( !ps.rank_mig_candidates )
    {
      return run_legacy( root, care, required, max_inserts, num_mffc, last_gain );
    }
    return run_ranked( root, care, required, max_inserts, num_mffc, last_gain );
  }

  std::optional<signal> run_legacy( node const& root, TT care, uint32_t required, uint32_t max_inserts, uint32_t num_mffc, uint32_t& last_gain )
  {
    (void)care;
    assert( is_const0( ~care ) );

    /* consider constants */
    auto g = call_with_stopwatch( st.time_resubC, [&]() {
      return resub_const( root, required );
    } );
    if ( g )
    {
      ++st.num_const_accepts;
      last_gain = num_mffc;
      return g; /* accepted resub */
    }

    /* consider equal nodes */
    g = call_with_stopwatch( st.time_resub0, [&]() {
      return resub_div0( root, required );
    } );
    if ( g )
    {
      ++st.num_div0_accepts;
      last_gain = num_mffc;
      return g; /* accepted resub */
    }

    /* consider relevance optimization */
    g = call_with_stopwatch( st.time_resubR, [&]() {
      return resub_divR( root, required );
    } );
    if ( g )
    {
      ++st.num_divR_accepts;
      last_gain = num_mffc;
      return g; /* accepted resub */
    }

    if ( max_inserts == 0 || num_mffc == 1 )
      return std::nullopt;

    /* collect level one divisors */
    call_with_stopwatch( st.time_collect_unate_divisors, [&]() {
      collect_unate_divisors( root, required );
    } );

    /* consider equal nodes */
    g = call_with_stopwatch( st.time_resub1, [&]() {
      return resub_div1( root, required );
    } );
    if ( g )
    {
      ++st.num_div1_accepts;
      last_gain = num_mffc - 1;
      return g; /* accepted resub */
    }

    if ( max_inserts == 1 || num_mffc == 2 )
      return std::nullopt;

    /* collect level two divisors */
    call_with_stopwatch( st.time_collect_binate_divisors, [&]() {
      collect_binate_divisors( root, required );
    } );

    /* consider two nodes */
    g = call_with_stopwatch( st.time_resub2, [&]() { return resub_div2( root, required ); } );
    if ( g )
    {
      ++st.num_div2_accepts;
      last_gain = num_mffc - 2;
      return g; /* accepted resub */
    }

    return std::nullopt;
  }

  std::optional<signal> run_ranked( node const& root, TT care, uint32_t required, uint32_t max_inserts, uint32_t num_mffc, uint32_t& last_gain )
  {
    (void)care;
    assert( is_const0( ~care ) );

    std::vector<ranked_candidate> pool;
    pool.reserve( std::max( 1u, ps.mig_candidate_pool_size ) );
    uint32_t discovery_index = 0u;

    auto consider = [&]( resub_candidate cand ) {
      cand.discovery_index = discovery_index++;
      ranked_candidate ranked{ cand, cheap_score_candidate( root, cand ) };
      if ( pool.size() < std::max( 1u, ps.mig_candidate_pool_size ) )
      {
        pool.push_back( ranked );
        return;
      }

      auto worst = std::min_element( pool.begin(), pool.end(), [&]( auto const& a, auto const& b ) {
        return better_cheap_candidate( b, a, required );
      } );
      if ( worst != pool.end() && better_cheap_candidate( ranked, *worst, required ) )
      {
        *worst = ranked;
      }
    };

    call_with_stopwatch( st.time_resubC, [&]() {
      collect_resub_const( root, required, num_mffc, consider );
    } );
    call_with_stopwatch( st.time_resub0, [&]() {
      collect_resub_div0( root, required, num_mffc, consider );
    } );
    call_with_stopwatch( st.time_resubR, [&]() {
      collect_resub_divR( root, required, consider, num_mffc );
    } );

    if ( max_inserts > 0 && num_mffc > 1 )
    {
      call_with_stopwatch( st.time_collect_unate_divisors, [&]() {
        collect_unate_divisors( root, required );
      } );
      call_with_stopwatch( st.time_resub1, [&]() {
        collect_resub_div1( root, required, consider, num_mffc );
      } );

      if ( max_inserts > 1 && num_mffc > 2 )
      {
        call_with_stopwatch( st.time_collect_binate_divisors, [&]() {
          collect_binate_divisors( root, required );
        } );
        call_with_stopwatch( st.time_resub2, [&]() {
          collect_resub_div2( root, required, consider, num_mffc );
        } );
      }
    }

    if ( pool.empty() )
    {
      return std::nullopt;
    }

    std::sort( pool.begin(), pool.end(), [&]( auto const& a, auto const& b ) {
      return better_cheap_candidate( a, b, required );
    } );

    auto const full_limit = std::min<uint32_t>( static_cast<uint32_t>( pool.size() ), std::max( 1u, ps.mig_full_score_limit ) );
    auto const before = compute_network_cost();

    std::optional<full_candidate_score> best;
    for ( uint32_t i = 0u; i < full_limit; ++i )
    {
      if ( auto scored = full_score_candidate( root, pool[i], before, required ) )
      {
        if ( !best.has_value() || better_full_candidate( *scored, *best, required ) )
        {
          best = *scored;
        }
      }
    }

    if ( !best.has_value() )
    {
      auto const& fallback = pool.front();
      last_gain = fallback.candidate.gain;
      update_accept_stats( fallback.candidate.stage );
      return materialize_candidate( ntk, fallback.candidate );
    }

    last_gain = best->ranked.candidate.gain;
    update_accept_stats( best->ranked.candidate.stage );
    return materialize_candidate( ntk, best->ranked.candidate );
  }

  std::optional<signal> resub_const( node const& root, uint32_t required ) const
  {
    (void)required;
    auto const tt = sim.get_tt( ntk.make_signal( root ) );
    if ( tt == sim.get_tt( ntk.get_constant( false ) ) )
    {
      return sim.get_phase( root ) ? ntk.get_constant( true ) : ntk.get_constant( false );
    }
    return std::nullopt;
  }

  std::optional<signal> resub_div0( node const& root, uint32_t required ) const
  {
    (void)required;
    auto const tt = sim.get_tt( ntk.make_signal( root ) );
    for ( auto i = 0u; i < num_divs; ++i )
    {
      auto const d = divs.at( i );

      if ( tt != sim.get_tt( ntk.make_signal( d ) ) )
        continue; /* next */

      return ( sim.get_phase( d ) ^ sim.get_phase( root ) ) ? !ntk.make_signal( d ) : ntk.make_signal( d );
    }

    return std::nullopt;
  }

  std::optional<signal> resub_divR( node const& root, uint32_t required )
  {
    (void)required;

    std::vector<signal> fs;
    ntk.foreach_fanin( root, [&]( const auto& f ) {
      fs.emplace_back( f );
    } );

    for ( auto i = 0u; i < divs.size(); ++i )
    {
      auto const& d0 = divs.at( i );
      auto const& s = ntk.make_signal( d0 );
      auto const& tt = sim.get_tt( s );

      if ( d0 == root )
        break;

      auto const tt0 = sim.get_tt( fs[0] );
      auto const tt1 = sim.get_tt( fs[1] );
      auto const tt2 = sim.get_tt( fs[2] );

      if ( ntk.get_node( fs[0] ) != d0 && ntk.fanout_size( ntk.get_node( fs[0] ) ) == 1 && can_replace_majority_fanin( tt0, tt1, tt2, tt ) )
      {
        auto const b = sim.get_phase( ntk.get_node( fs[1] ) ) ? !fs[1] : fs[1];
        auto const c = sim.get_phase( ntk.get_node( fs[2] ) ) ? !fs[2] : fs[2];

        return sim.get_phase( root ) ? !ntk.create_maj( sim.get_phase( d0 ) ? !s : s, b, c ) : ntk.create_maj( sim.get_phase( d0 ) ? !s : s, b, c );
      }
      else if ( ntk.get_node( fs[1] ) != d0 && ntk.fanout_size( ntk.get_node( fs[1] ) ) == 1 && can_replace_majority_fanin( tt1, tt0, tt2, tt ) )
      {
        auto const a = sim.get_phase( ntk.get_node( fs[0] ) ) ? !fs[0] : fs[0];
        auto const c = sim.get_phase( ntk.get_node( fs[2] ) ) ? !fs[2] : fs[2];

        return sim.get_phase( root ) ? !ntk.create_maj( sim.get_phase( d0 ) ? !s : s, a, c ) : ntk.create_maj( sim.get_phase( d0 ) ? !s : s, a, c );
      }
      else if ( ntk.get_node( fs[2] ) != d0 && ntk.fanout_size( ntk.get_node( fs[2] ) ) == 1 && can_replace_majority_fanin( tt2, tt0, tt1, tt ) )
      {
        auto const a = sim.get_phase( ntk.get_node( fs[0] ) ) ? !fs[0] : fs[0];
        auto const b = sim.get_phase( ntk.get_node( fs[1] ) ) ? !fs[1] : fs[1];

        return sim.get_phase( root ) ? !ntk.create_maj( sim.get_phase( d0 ) ? !s : s, a, b ) : ntk.create_maj( sim.get_phase( d0 ) ? !s : s, a, b );
      }
      else if ( ntk.get_node( fs[0] ) != d0 && ntk.fanout_size( ntk.get_node( fs[0] ) ) == 1 && can_replace_majority_fanin( ~tt0, tt1, tt2, tt ) )
      {
        auto const b = sim.get_phase( ntk.get_node( fs[1] ) ) ? !fs[1] : fs[1];
        auto const c = sim.get_phase( ntk.get_node( fs[2] ) ) ? !fs[2] : fs[2];

        return sim.get_phase( root ) ? !ntk.create_maj( sim.get_phase( d0 ) ? s : !s, b, c ) : ntk.create_maj( sim.get_phase( d0 ) ? s : !s, b, c );
      }
      else if ( ntk.get_node( fs[1] ) != d0 && ntk.fanout_size( ntk.get_node( fs[1] ) ) == 1 && can_replace_majority_fanin( ~tt1, tt0, tt2, tt ) )
      {
        auto const a = sim.get_phase( ntk.get_node( fs[0] ) ) ? !fs[0] : fs[0];
        auto const c = sim.get_phase( ntk.get_node( fs[2] ) ) ? !fs[2] : fs[2];

        return sim.get_phase( root ) ? !ntk.create_maj( sim.get_phase( d0 ) ? s : !s, a, c ) : ntk.create_maj( sim.get_phase( d0 ) ? s : !s, a, c );
      }
      else if ( ntk.get_node( fs[2] ) != d0 && ntk.fanout_size( ntk.get_node( fs[2] ) ) == 1 && can_replace_majority_fanin( ~tt2, tt0, tt1, tt ) )
      {
        auto const a = sim.get_phase( ntk.get_node( fs[0] ) ) ? !fs[0] : fs[0];
        auto const b = sim.get_phase( ntk.get_node( fs[1] ) ) ? !fs[1] : fs[1];

        return sim.get_phase( root ) ? !ntk.create_maj( sim.get_phase( d0 ) ? s : !s, a, b ) : ntk.create_maj( sim.get_phase( d0 ) ? s : !s, a, b );
      }
    }

    return std::nullopt;
  }

  void collect_unate_divisors( node const& root, uint32_t required )
  {
    udivs.clear();

    auto const& tt = sim.get_tt( ntk.make_signal( root ) );
    auto const& one = sim.get_tt( ntk.get_constant( true ) );
    for ( auto i = 0u; i < num_divs; ++i )
    {
      auto const d0 = divs.at( i );
      if ( ntk.level( d0 ) > required - 1 )
        continue;
      auto const& tt_s0 = sim.get_tt( ntk.make_signal( d0 ) );

      for ( auto j = i + 1; j < num_divs; ++j )
      {
        auto const d1 = divs.at( j );
        if ( ntk.level( d1 ) > required - 1 )
          continue;
        auto const& tt_s1 = sim.get_tt( ntk.make_signal( d1 ) );

        /* Boolean filtering rule for MAJ-3 */
        if ( kitty::ternary_majority( tt_s0, tt_s1, tt ) == tt )
        {
          udivs.u0.emplace_back( ntk.make_signal( d0 ) );
          udivs.u1.emplace_back( ntk.make_signal( d1 ) );
          continue;
        }

        if ( kitty::ternary_majority( ~tt_s0, tt_s1, tt ) == tt )
        {
          udivs.u0.emplace_back( !ntk.make_signal( d0 ) );
          udivs.u1.emplace_back( ntk.make_signal( d1 ) );
          continue;
        }

        if ( kitty::ternary_majority( tt_s0, ~tt_s1, tt ) == tt )
        {
          udivs.u0.emplace_back( ntk.make_signal( d0 ) );
          udivs.u1.emplace_back( !ntk.make_signal( d1 ) );
          continue;
        }

        if ( std::find( udivs.next_candidates.begin(), udivs.next_candidates.end(), ntk.make_signal( d1 ) ) == udivs.next_candidates.end() )
          udivs.next_candidates.emplace_back( ntk.make_signal( d1 ) );
      }

      if constexpr ( use_constant ) /* allowing "not real" MAJ gates (one fanin is constant) */
      {
        if ( kitty::ternary_majority( tt_s0, one, tt ) == tt )
        {
          udivs.u0.emplace_back( ntk.make_signal( d0 ) );
          udivs.u1.emplace_back( ntk.get_constant( true ) );
          continue;
        }

        if ( kitty::ternary_majority( ~tt_s0, one, tt ) == tt )
        {
          udivs.u0.emplace_back( !ntk.make_signal( d0 ) );
          udivs.u1.emplace_back( ntk.get_constant( true ) );
          continue;
        }

        if ( kitty::ternary_majority( tt_s0, ~one, tt ) == tt )
        {
          udivs.u0.emplace_back( ntk.make_signal( d0 ) );
          udivs.u1.emplace_back( ntk.get_constant( false ) );
          continue;
        }
      }

      if ( std::find( udivs.next_candidates.begin(), udivs.next_candidates.end(), ntk.make_signal( d0 ) ) == udivs.next_candidates.end() )
        udivs.next_candidates.emplace_back( ntk.make_signal( d0 ) );
    }

    if constexpr ( use_constant )
    {
      udivs.next_candidates.emplace_back( ntk.get_constant( true ) );
    }
  }

  std::optional<signal> resub_div1( node const& root, uint32_t required )
  {
    (void)required;
    auto const& tt = sim.get_tt( ntk.make_signal( root ) );

    for ( auto i = 0u; i < udivs.u0.size(); ++i )
    {
      auto const s0 = udivs.u0.at( i );
      auto const s1 = udivs.u1.at( i );
      auto const& tt_s0 = sim.get_tt( s0 );
      auto const& tt_s1 = sim.get_tt( s1 );

      for ( auto j = i + 1; j < udivs.u0.size(); ++j )
      {
        auto s2 = udivs.u0.at( j );
        auto tt_s2 = sim.get_tt( s2 );

        if ( kitty::ternary_majority( tt_s0, tt_s1, tt_s2 ) == tt )
        {
          auto const a = sim.get_phase( ntk.get_node( s0 ) ) ? !s0 : s0;
          auto const b = sim.get_phase( ntk.get_node( s1 ) ) ? !s1 : s1;
          auto const c = sim.get_phase( ntk.get_node( s2 ) ) ? !s2 : s2;
          return sim.get_phase( root ) ? !ntk.create_maj( a, b, c ) : ntk.create_maj( a, b, c );
        }

        s2 = udivs.u1.at( j );
        tt_s2 = sim.get_tt( s2 );

        if ( kitty::ternary_majority( tt_s0, tt_s1, tt_s2 ) == tt )
        {
          auto const a = sim.get_phase( ntk.get_node( s0 ) ) ? !s0 : s0;
          auto const b = sim.get_phase( ntk.get_node( s1 ) ) ? !s1 : s1;
          auto const c = sim.get_phase( ntk.get_node( s2 ) ) ? !s2 : s2;
          return sim.get_phase( root ) ? !ntk.create_maj( a, b, c ) : ntk.create_maj( a, b, c );
        }
      }
    }

    return std::nullopt;
  }

  void collect_binate_divisors( node const& root, uint32_t required )
  {
    bdivs.clear();

    auto const& tt = sim.get_tt( ntk.make_signal( root ) );
    for ( auto i = 0u; i < udivs.next_candidates.size(); ++i )
    {
      auto const& s0 = udivs.next_candidates.at( i );
      if ( ntk.level( ntk.get_node( s0 ) ) > required - 2 )
        continue;

      auto const& tt_s0 = sim.get_tt( s0 );

      for ( auto j = i + 1; j < udivs.next_candidates.size(); ++j )
      {
        auto const& s1 = udivs.next_candidates.at( j );
        if ( ntk.level( ntk.get_node( s1 ) ) > required - 2 )
          continue;

        auto const& tt_s1 = sim.get_tt( s1 );

        for ( auto k = j + 1; k < udivs.next_candidates.size(); ++k )
        {
          auto const& s2 = udivs.next_candidates.at( k );
          if ( ntk.level( ntk.get_node( s2 ) ) > required - 2 )
            continue;

          auto const& tt_s2 = sim.get_tt( s2 );

          /* Note: the implication relation is actually not necessary for majority; this is an over-filtering */
          if ( kitty::implies( kitty::ternary_majority( tt_s0, tt_s1, tt_s2 ), tt ) )
          {
            bdivs.b0.emplace_back( s0 );
            bdivs.b1.emplace_back( s1 );
            bdivs.b2.emplace_back( s2 );
            continue;
          }

          if ( kitty::implies( kitty::ternary_majority( ~tt_s0, tt_s1, tt_s2 ), tt ) )
          {
            bdivs.b0.emplace_back( !s0 );
            bdivs.b1.emplace_back( s1 );
            bdivs.b2.emplace_back( s2 );
            continue;
          }

          if ( kitty::implies( kitty::ternary_majority( tt_s0, ~tt_s1, tt_s2 ), tt ) )
          {
            bdivs.b0.emplace_back( s0 );
            bdivs.b1.emplace_back( !s1 );
            bdivs.b2.emplace_back( s2 );
            continue;
          }

          if ( kitty::implies( kitty::ternary_majority( tt_s0, tt_s1, ~tt_s2 ), tt ) )
          {
            bdivs.b0.emplace_back( s0 );
            bdivs.b1.emplace_back( s1 );
            bdivs.b2.emplace_back( !s2 );
            continue;
          }

          if ( kitty::implies( kitty::ternary_majority( ~tt_s0, ~tt_s1, tt_s2 ), tt ) )
          {
            bdivs.b0.emplace_back( !s0 );
            bdivs.b1.emplace_back( !s1 );
            bdivs.b2.emplace_back( s2 );
            continue;
          }

          if ( kitty::implies( kitty::ternary_majority( tt_s0, ~tt_s1, ~tt_s2 ), tt ) )
          {
            bdivs.b0.emplace_back( s0 );
            bdivs.b1.emplace_back( !s1 );
            bdivs.b2.emplace_back( !s2 );
            continue;
          }

          if ( kitty::implies( kitty::ternary_majority( ~tt_s0, tt_s1, ~tt_s2 ), tt ) )
          {
            bdivs.b0.emplace_back( !s0 );
            bdivs.b1.emplace_back( s1 );
            bdivs.b2.emplace_back( !s2 );
            continue;
          }

          if ( kitty::implies( kitty::ternary_majority( ~tt_s0, ~tt_s1, ~tt_s2 ), tt ) )
          {
            bdivs.b0.emplace_back( !s0 );
            bdivs.b1.emplace_back( !s1 );
            bdivs.b2.emplace_back( !s2 );
            continue;
          }
        }
      }
    }
  }

  std::optional<signal> resub_div2( node const& root, uint32_t required )
  {
    (void)required;
    auto const& tt = sim.get_tt( ntk.make_signal( root ) );

    for ( auto i = 0u; i < udivs.u0.size(); ++i )
    {
      auto const& s0 = udivs.u0.at( i );
      auto const& s1 = udivs.u1.at( i );

      for ( auto j = 0u; j < bdivs.b0.size(); ++j )
      {
        auto const& s2 = bdivs.b0.at( j );
        auto const& s3 = bdivs.b1.at( j );
        auto const& s4 = bdivs.b2.at( j );

        auto const a = sim.get_phase( ntk.get_node( s0 ) ) ? !s0 : s0;
        auto const b = sim.get_phase( ntk.get_node( s1 ) ) ? !s1 : s1;
        auto const c = sim.get_phase( ntk.get_node( s2 ) ) ? !s2 : s2;
        auto const d = sim.get_phase( ntk.get_node( s3 ) ) ? !s3 : s3;
        auto const e = sim.get_phase( ntk.get_node( s4 ) ) ? !s4 : s4;

        auto const& tt_s0 = sim.get_tt( s0 );
        auto const& tt_s1 = sim.get_tt( s1 );
        auto const& tt_s2 = sim.get_tt( s2 );
        auto const& tt_s3 = sim.get_tt( s3 );
        auto const& tt_s4 = sim.get_tt( s4 );

        if ( kitty::ternary_majority( tt_s0, tt_s1, kitty::ternary_majority( tt_s2, tt_s3, tt_s4 ) ) == tt )
        {
          return sim.get_phase( root ) ? !ntk.create_maj( a, b, ntk.create_maj( c, d, e ) ) : ntk.create_maj( a, b, ntk.create_maj( c, d, e ) );
        }
      }
    }

    return std::nullopt;
  }

  template<class Fn>
  void collect_resub_const( node const& root, uint32_t required, uint32_t num_mffc, Fn&& fn ) const
  {
    if ( auto const g = resub_const( root, required ) )
    {
      resub_candidate cand{};
      cand.stage = candidate_stage::constant;
      cand.kind = candidate_kind::direct;
      cand.direct = *g;
      cand.gain = num_mffc;
      fn( cand );
    }
  }

  template<class Fn>
  void collect_resub_div0( node const& root, uint32_t required, uint32_t num_mffc, Fn&& fn ) const
  {
    (void)required;
    auto const tt = sim.get_tt( ntk.make_signal( root ) );
    for ( auto i = 0u; i < num_divs; ++i )
    {
      auto const d = divs.at( i );
      if ( tt != sim.get_tt( ntk.make_signal( d ) ) )
      {
        continue;
      }

      resub_candidate cand{};
      cand.stage = candidate_stage::div0;
      cand.kind = candidate_kind::direct;
      cand.direct = ( sim.get_phase( d ) ^ sim.get_phase( root ) ) ? !ntk.make_signal( d ) : ntk.make_signal( d );
      cand.gain = num_mffc;
      fn( cand );
    }
  }

  template<class Fn>
  void collect_resub_divR( node const& root, uint32_t required, Fn&& fn, uint32_t num_mffc ) const
  {
    (void)required;

    std::vector<signal> fs;
    ntk.foreach_fanin( root, [&]( const auto& f ) {
      fs.emplace_back( f );
    } );

    auto const tt0 = sim.get_tt( fs[0] );
    auto const tt1 = sim.get_tt( fs[1] );
    auto const tt2 = sim.get_tt( fs[2] );

    auto emit = [&]( signal const& a, signal const& b, signal const& c ) {
      resub_candidate cand{};
      cand.stage = candidate_stage::divR;
      cand.kind = candidate_kind::maj3;
      cand.outer = { a, b, c };
      cand.complement_output = sim.get_phase( root );
      cand.gain = num_mffc;
      cand.inserts = 1u;
      fn( cand );

      emit_min_dual_maj3( fn, cand );
    };

    for ( auto i = 0u; i < divs.size(); ++i )
    {
      auto const& d0 = divs.at( i );
      auto const& s = ntk.make_signal( d0 );
      auto const& tt = sim.get_tt( s );

      if ( d0 == root )
      {
        break;
      }

      if ( ntk.get_node( fs[0] ) != d0 && ntk.fanout_size( ntk.get_node( fs[0] ) ) == 1 && can_replace_majority_fanin( tt0, tt1, tt2, tt ) )
      {
        auto const a = sim.get_phase( d0 ) ? !s : s;
        auto const b = sim.get_phase( ntk.get_node( fs[1] ) ) ? !fs[1] : fs[1];
        auto const c = sim.get_phase( ntk.get_node( fs[2] ) ) ? !fs[2] : fs[2];
        emit( a, b, c );
      }
      if ( ntk.get_node( fs[1] ) != d0 && ntk.fanout_size( ntk.get_node( fs[1] ) ) == 1 && can_replace_majority_fanin( tt1, tt0, tt2, tt ) )
      {
        auto const a = sim.get_phase( d0 ) ? !s : s;
        auto const b = sim.get_phase( ntk.get_node( fs[0] ) ) ? !fs[0] : fs[0];
        auto const c = sim.get_phase( ntk.get_node( fs[2] ) ) ? !fs[2] : fs[2];
        emit( a, b, c );
      }
      if ( ntk.get_node( fs[2] ) != d0 && ntk.fanout_size( ntk.get_node( fs[2] ) ) == 1 && can_replace_majority_fanin( tt2, tt0, tt1, tt ) )
      {
        auto const a = sim.get_phase( d0 ) ? !s : s;
        auto const b = sim.get_phase( ntk.get_node( fs[0] ) ) ? !fs[0] : fs[0];
        auto const c = sim.get_phase( ntk.get_node( fs[1] ) ) ? !fs[1] : fs[1];
        emit( a, b, c );
      }
      if ( ntk.get_node( fs[0] ) != d0 && ntk.fanout_size( ntk.get_node( fs[0] ) ) == 1 && can_replace_majority_fanin( ~tt0, tt1, tt2, tt ) )
      {
        auto const a = sim.get_phase( d0 ) ? s : !s;
        auto const b = sim.get_phase( ntk.get_node( fs[1] ) ) ? !fs[1] : fs[1];
        auto const c = sim.get_phase( ntk.get_node( fs[2] ) ) ? !fs[2] : fs[2];
        emit( a, b, c );
      }
      if ( ntk.get_node( fs[1] ) != d0 && ntk.fanout_size( ntk.get_node( fs[1] ) ) == 1 && can_replace_majority_fanin( ~tt1, tt0, tt2, tt ) )
      {
        auto const a = sim.get_phase( d0 ) ? s : !s;
        auto const b = sim.get_phase( ntk.get_node( fs[0] ) ) ? !fs[0] : fs[0];
        auto const c = sim.get_phase( ntk.get_node( fs[2] ) ) ? !fs[2] : fs[2];
        emit( a, b, c );
      }
      if ( ntk.get_node( fs[2] ) != d0 && ntk.fanout_size( ntk.get_node( fs[2] ) ) == 1 && can_replace_majority_fanin( ~tt2, tt0, tt1, tt ) )
      {
        auto const a = sim.get_phase( d0 ) ? s : !s;
        auto const b = sim.get_phase( ntk.get_node( fs[0] ) ) ? !fs[0] : fs[0];
        auto const c = sim.get_phase( ntk.get_node( fs[1] ) ) ? !fs[1] : fs[1];
        emit( a, b, c );
      }
    }
  }

  template<class Fn>
  void collect_resub_div1( node const& root, uint32_t required, Fn&& fn, uint32_t num_mffc ) const
  {
    (void)required;
    auto const& tt = sim.get_tt( ntk.make_signal( root ) );

    for ( auto i = 0u; i < udivs.u0.size(); ++i )
    {
      auto const s0 = udivs.u0.at( i );
      auto const s1 = udivs.u1.at( i );
      auto const& tt_s0 = sim.get_tt( s0 );
      auto const& tt_s1 = sim.get_tt( s1 );

      for ( auto j = i + 1; j < udivs.u0.size(); ++j )
      {
        auto emit = [&]( signal const& s2 ) {
          auto const a = sim.get_phase( ntk.get_node( s0 ) ) ? !s0 : s0;
          auto const b = sim.get_phase( ntk.get_node( s1 ) ) ? !s1 : s1;
          auto const c = sim.get_phase( ntk.get_node( s2 ) ) ? !s2 : s2;

          resub_candidate cand{};
          cand.stage = candidate_stage::div1;
          cand.kind = candidate_kind::maj3;
          cand.outer = { a, b, c };
          cand.complement_output = sim.get_phase( root );
          cand.gain = num_mffc - 1u;
          cand.inserts = 1u;
          fn( cand );

          emit_min_dual_maj3( fn, cand );
        };

        auto s2 = udivs.u0.at( j );
        auto tt_s2 = sim.get_tt( s2 );
        if ( kitty::ternary_majority( tt_s0, tt_s1, tt_s2 ) == tt )
        {
          emit( s2 );
        }

        s2 = udivs.u1.at( j );
        tt_s2 = sim.get_tt( s2 );
        if ( kitty::ternary_majority( tt_s0, tt_s1, tt_s2 ) == tt )
        {
          emit( s2 );
        }
      }
    }
  }

  template<class Fn>
  void collect_resub_div2( node const& root, uint32_t required, Fn&& fn, uint32_t num_mffc ) const
  {
    (void)required;
    auto const& tt = sim.get_tt( ntk.make_signal( root ) );

    for ( auto i = 0u; i < udivs.u0.size(); ++i )
    {
      auto const& s0 = udivs.u0.at( i );
      auto const& s1 = udivs.u1.at( i );
      auto const& tt_s0 = sim.get_tt( s0 );
      auto const& tt_s1 = sim.get_tt( s1 );

      for ( auto j = 0u; j < bdivs.b0.size(); ++j )
      {
        auto const& s2 = bdivs.b0.at( j );
        auto const& s3 = bdivs.b1.at( j );
        auto const& s4 = bdivs.b2.at( j );

        auto const& tt_s2 = sim.get_tt( s2 );
        auto const& tt_s3 = sim.get_tt( s3 );
        auto const& tt_s4 = sim.get_tt( s4 );

        if ( kitty::ternary_majority( tt_s0, tt_s1, kitty::ternary_majority( tt_s2, tt_s3, tt_s4 ) ) != tt )
        {
          continue;
        }

        auto const a = sim.get_phase( ntk.get_node( s0 ) ) ? !s0 : s0;
        auto const b = sim.get_phase( ntk.get_node( s1 ) ) ? !s1 : s1;
        auto const c = sim.get_phase( ntk.get_node( s2 ) ) ? !s2 : s2;
        auto const d = sim.get_phase( ntk.get_node( s3 ) ) ? !s3 : s3;
        auto const e = sim.get_phase( ntk.get_node( s4 ) ) ? !s4 : s4;

        resub_candidate cand{};
        cand.stage = candidate_stage::div2;
        cand.kind = candidate_kind::maj5;
        cand.outer = { a, b, signal{} };
        cand.inner = { c, d, e };
        cand.complement_output = sim.get_phase( root );
        cand.gain = num_mffc - 2u;
        cand.inserts = 2u;
        fn( cand );

        emit_min_dual_maj5( fn, cand );
      }
    }
  }

  bool enable_min_resubstitution() const
  {
    if constexpr ( has_create_min_v<Ntk> )
    {
      return ps.enable_min_resubstitution;
    }
    return false;
  }

  template<class Fn>
  void emit_min_dual_maj3( Fn&& fn, resub_candidate const& maj_cand ) const
  {
    if ( !enable_min_resubstitution() )
    {
      return;
    }

    resub_candidate min_cand = maj_cand;
    min_cand.outer = { !maj_cand.outer[0], !maj_cand.outer[1], !maj_cand.outer[2] };
    min_cand.outer_is_min = true;
    fn( min_cand );
  }

  template<class Fn>
  void emit_min_dual_maj5( Fn&& fn, resub_candidate const& maj_cand ) const
  {
    if ( !enable_min_resubstitution() )
    {
      return;
    }

    /* Replace the inner MAJ by its equivalent MIN dual. */
    resub_candidate inner_min = maj_cand;
    inner_min.inner = { !maj_cand.inner[0], !maj_cand.inner[1], !maj_cand.inner[2] };
    inner_min.inner_is_min = true;
    fn( inner_min );

    /* Replace the outer MAJ by its equivalent MIN dual. */
    resub_candidate outer_min = maj_cand;
    outer_min.outer = { !maj_cand.outer[0], !maj_cand.outer[1], signal{} };
    outer_min.outer_is_min = true;
    outer_min.complement_inner_in_outer = true;
    fn( outer_min );

    /* Keep both duals when both levels benefit from absorbed inversions. */
    resub_candidate both_min = outer_min;
    both_min.inner = inner_min.inner;
    both_min.inner_is_min = true;
    fn( both_min );
  }

  template<class Net>
  signal remap_signal( Net& net, signal const& s ) const
  {
    auto remapped = net.make_signal( ntk.get_node( s ) );
    return ntk.is_complemented( s ) ? !remapped : remapped;
  }

  template<class Net>
  signal materialize_candidate( Net& net, resub_candidate const& cand ) const
  {
    if ( cand.kind == candidate_kind::direct )
    {
      return remap_signal( net, cand.direct );
    }

    if ( cand.kind == candidate_kind::maj3 )
    {
      auto const a = remap_signal( net, cand.outer[0] );
      auto const b = remap_signal( net, cand.outer[1] );
      auto const c = remap_signal( net, cand.outer[2] );
      auto g = create_resub_gate( net, a, b, c, cand.outer_is_min );
      return cand.complement_output ? !g : g;
    }

    auto const inner = create_resub_gate( net,
                                          remap_signal( net, cand.inner[0] ),
                                          remap_signal( net, cand.inner[1] ),
                                          remap_signal( net, cand.inner[2] ),
                                          cand.inner_is_min );
    auto const outer_inner = cand.complement_inner_in_outer ? !inner : inner;
    auto g = create_resub_gate( net,
                                remap_signal( net, cand.outer[0] ),
                                remap_signal( net, cand.outer[1] ),
                                outer_inner,
                                cand.outer_is_min );
    return cand.complement_output ? !g : g;
  }

  template<class Net>
  signal create_resub_gate( Net& net, signal const& a, signal const& b, signal const& c, bool use_min ) const
  {
    if ( use_min )
    {
      if constexpr ( has_create_min_v<Net> )
      {
        return net.create_min( a, b, c );
      }
      else
      {
        return net.create_maj( !a, !b, !c );
      }
    }
    return net.create_maj( a, b, c );
  }

  cheap_candidate_score cheap_score_candidate( node const& root, resub_candidate const& cand ) const
  {
    cheap_candidate_score score{};
    score.estimated_level = estimate_candidate_level( cand );
    score.level_gain = static_cast<int32_t>( ntk.level( root ) ) - score.estimated_level;
    score.inv_cost = count_candidate_inversions( cand );
    score.duplication_penalty = candidate_duplication_penalty( cand );
    score.share_bonus = candidate_share_bonus( cand );
    score.gain = cand.gain;
    score.discovery_index = cand.discovery_index;
    return score;
  }

  bool better_cheap_candidate( ranked_candidate const& lhs, ranked_candidate const& rhs, uint32_t required ) const
  {
    auto const depth_mode = required != std::numeric_limits<uint32_t>::max();
    if ( depth_mode )
    {
      if ( lhs.cheap.level_gain != rhs.cheap.level_gain )
      {
        return lhs.cheap.level_gain > rhs.cheap.level_gain;
      }
      if ( lhs.cheap.share_bonus != rhs.cheap.share_bonus )
      {
        return lhs.cheap.share_bonus > rhs.cheap.share_bonus;
      }
      if ( lhs.cheap.gain != rhs.cheap.gain )
      {
        return lhs.cheap.gain > rhs.cheap.gain;
      }
    }
    else
    {
      if ( lhs.cheap.gain != rhs.cheap.gain )
      {
        return lhs.cheap.gain > rhs.cheap.gain;
      }
      if ( lhs.cheap.share_bonus != rhs.cheap.share_bonus )
      {
        return lhs.cheap.share_bonus > rhs.cheap.share_bonus;
      }
      if ( lhs.cheap.level_gain != rhs.cheap.level_gain )
      {
        return lhs.cheap.level_gain > rhs.cheap.level_gain;
      }
    }

    if ( lhs.cheap.inv_cost != rhs.cheap.inv_cost )
    {
      return lhs.cheap.inv_cost < rhs.cheap.inv_cost;
    }
    if ( lhs.cheap.duplication_penalty != rhs.cheap.duplication_penalty )
    {
      return lhs.cheap.duplication_penalty < rhs.cheap.duplication_penalty;
    }
    if ( lhs.candidate.inserts != rhs.candidate.inserts )
    {
      return lhs.candidate.inserts < rhs.candidate.inserts;
    }
    if ( lhs.cheap.estimated_level != rhs.cheap.estimated_level )
    {
      return lhs.cheap.estimated_level < rhs.cheap.estimated_level;
    }
    return lhs.cheap.discovery_index < rhs.cheap.discovery_index;
  }

  std::optional<full_candidate_score> full_score_candidate( node const& root, ranked_candidate const& ranked, network_cost const& before, uint32_t required ) const
  {
    (void)required;
    if constexpr ( has_clone_v<Ntk> && has_foreach_po_v<Ntk> )
    {
      auto candidate_ntk = ntk.clone();
      auto const g = materialize_candidate( candidate_ntk, ranked.candidate );
      candidate_ntk.substitute_node( root, g );
      candidate_ntk = cleanup_dangling( candidate_ntk );

      auto const after = compute_network_cost( candidate_ntk );

      full_candidate_score scored{};
      scored.ranked = ranked;
      scored.depth_gain = before.depth - after.depth;
      scored.gate_gain = before.gates - after.gates;
      scored.inv_gain = before.inverted_edges - after.inverted_edges;
      return scored;
    }

    return std::nullopt;
  }

  bool better_full_candidate( full_candidate_score const& lhs, full_candidate_score const& rhs, uint32_t required ) const
  {
    auto const depth_mode = required != std::numeric_limits<uint32_t>::max();
    if ( depth_mode )
    {
      if ( lhs.depth_gain != rhs.depth_gain )
      {
        return lhs.depth_gain > rhs.depth_gain;
      }
      if ( lhs.gate_gain != rhs.gate_gain )
      {
        return lhs.gate_gain > rhs.gate_gain;
      }
    }
    else
    {
      if ( lhs.gate_gain != rhs.gate_gain )
      {
        return lhs.gate_gain > rhs.gate_gain;
      }
      if ( lhs.depth_gain != rhs.depth_gain )
      {
        return lhs.depth_gain > rhs.depth_gain;
      }
    }

    if ( lhs.inv_gain != rhs.inv_gain )
    {
      return lhs.inv_gain > rhs.inv_gain;
    }
    if ( lhs.ranked.cheap.share_bonus != rhs.ranked.cheap.share_bonus )
    {
      return lhs.ranked.cheap.share_bonus > rhs.ranked.cheap.share_bonus;
    }
    if ( lhs.ranked.candidate.inserts != rhs.ranked.candidate.inserts )
    {
      return lhs.ranked.candidate.inserts < rhs.ranked.candidate.inserts;
    }
    return lhs.ranked.cheap.discovery_index < rhs.ranked.cheap.discovery_index;
  }

  int32_t estimate_candidate_level( resub_candidate const& cand ) const
  {
    if ( cand.kind == candidate_kind::direct )
    {
      return signal_level( cand.direct );
    }

    if ( cand.kind == candidate_kind::maj3 )
    {
      return 1 + std::max( { signal_level( cand.outer[0] ), signal_level( cand.outer[1] ), signal_level( cand.outer[2] ) } );
    }

    auto const inner_level = 1 + std::max( { signal_level( cand.inner[0] ), signal_level( cand.inner[1] ), signal_level( cand.inner[2] ) } );
    return 1 + std::max( { signal_level( cand.outer[0] ), signal_level( cand.outer[1] ), inner_level } );
  }

  int32_t signal_level( signal const& s ) const
  {
    auto const n = ntk.get_node( s );
    if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
    {
      return 0;
    }
    return static_cast<int32_t>( ntk.level( n ) );
  }

  int32_t count_candidate_inversions( resub_candidate const& cand ) const
  {
    int32_t inv = cand.complement_output ? 1 : 0;
    auto add_signal = [&]( signal const& s ) {
      if ( ntk.is_complemented( s ) )
      {
        ++inv;
      }
    };

    if ( cand.kind == candidate_kind::direct )
    {
      add_signal( cand.direct );
      return inv;
    }

    add_signal( cand.outer[0] );
    add_signal( cand.outer[1] );
    if ( cand.kind == candidate_kind::maj3 )
    {
      add_signal( cand.outer[2] );
      return inv;
    }

    add_signal( cand.inner[0] );
    add_signal( cand.inner[1] );
    add_signal( cand.inner[2] );
    if ( cand.complement_inner_in_outer )
    {
      ++inv;
    }
    return inv;
  }

  int32_t candidate_duplication_penalty( resub_candidate const& cand ) const
  {
    std::array<node, 5u> seen{};
    uint32_t seen_size = 0u;
    int32_t penalty = 0;

    auto visit_signal = [&]( signal const& s ) {
      auto const n = ntk.get_node( s );
      if ( ntk.is_constant( n ) )
      {
        return;
      }
      if ( std::find( seen.begin(), seen.begin() + seen_size, n ) != seen.begin() + seen_size )
      {
        return;
      }
      seen[seen_size++] = n;
      penalty += static_cast<int32_t>( ntk.fanout_size( n ) );
    };

    if ( cand.kind == candidate_kind::direct )
    {
      visit_signal( cand.direct );
      return penalty;
    }

    visit_signal( cand.outer[0] );
    visit_signal( cand.outer[1] );
    if ( cand.kind == candidate_kind::maj3 )
    {
      visit_signal( cand.outer[2] );
      return penalty;
    }

    visit_signal( cand.inner[0] );
    visit_signal( cand.inner[1] );
    visit_signal( cand.inner[2] );
    return penalty;
  }

  int32_t candidate_share_bonus( resub_candidate const& cand ) const
  {
    if constexpr ( has_has_maj_v<Ntk> )
    {
      auto existing_bonus = [&]( signal const& a, signal const& b, signal const& c ) {
        auto const existing = ntk.has_maj( a, b, c );
        if ( existing.has_value() && ntk.get_node( *existing ) != ntk.get_node( ntk.get_constant( false ) ) )
        {
          return 1;
        }
        return 0;
      };

      if ( cand.kind == candidate_kind::maj3 )
      {
        if ( cand.outer_is_min )
        {
          return 0;
        }
        return existing_bonus( cand.outer[0], cand.outer[1], cand.outer[2] );
      }

      if ( cand.kind == candidate_kind::maj5 )
      {
        if ( cand.outer_is_min || cand.inner_is_min || cand.complement_inner_in_outer )
        {
          return 0;
        }
        auto bonus = existing_bonus( cand.inner[0], cand.inner[1], cand.inner[2] );
        if ( auto const inner = ntk.has_maj( cand.inner[0], cand.inner[1], cand.inner[2] ) )
        {
          bonus += existing_bonus( cand.outer[0], cand.outer[1], *inner );
        }
        return bonus;
      }
    }

    return 0;
  }

  template<class Net>
  int32_t count_inverted_edges( Net const& net ) const
  {
    int32_t inv_edges = 0;
    net.foreach_gate( [&]( auto const& n ) {
      net.foreach_fanin( n, [&]( auto const& f ) {
        if ( net.is_complemented( f ) )
        {
          ++inv_edges;
        }
      } );
    } );
    net.foreach_po( [&]( auto const& f, auto ) {
      if ( net.is_complemented( f ) )
      {
        ++inv_edges;
      }
    } );
    return inv_edges;
  }

  template<class Net>
  network_cost compute_network_cost( Net const& net ) const
  {
    depth_view depth_net{ net };
    network_cost cost{};
    cost.depth = static_cast<int32_t>( depth_net.depth() );
    cost.gates = static_cast<int32_t>( net.num_gates() );
    cost.inverted_edges = count_inverted_edges( net );
    return cost;
  }

  network_cost compute_network_cost() const
  {
    return compute_network_cost( ntk );
  }

  void update_accept_stats( candidate_stage stage )
  {
    switch ( stage )
    {
    case candidate_stage::constant:
      ++st.num_const_accepts;
      break;
    case candidate_stage::div0:
      ++st.num_div0_accepts;
      break;
    case candidate_stage::divR:
      ++st.num_divR_accepts;
      break;
    case candidate_stage::div1:
      ++st.num_div1_accepts;
      break;
    case candidate_stage::div2:
      ++st.num_div2_accepts;
      break;
    }
  }

private:
  Ntk& ntk;
  Simulator const& sim;
  std::vector<node> const& divs;
  uint32_t const num_divs;
  stats& st;
  resubstitution_params default_ps{};
  resubstitution_params const& ps;

  unate_divisors udivs;
  binate_divisors bdivs;
}; /* mig_enumerative_resub_functor */

struct mig_resyn_resub_stats
{
  /*! \brief Time for finding dependency function. */
  stopwatch<>::duration time_compute_function{ 0 };

  /*! \brief Number of found solutions. */
  uint32_t num_success{ 0 };

  /*! \brief Number of times that no solution can be found. */
  uint32_t num_fail{ 0 };

  void report() const
  {
    fmt::print( "[i]     <ResubFn: mig_resyn_functor>\n" );
    fmt::print( "[i]         #solution = {:6d}\n", num_success );
    fmt::print( "[i]         #invoke   = {:6d}\n", num_success + num_fail );
    fmt::print( "[i]         engine time: {:>5.2f} secs\n", to_seconds( time_compute_function ) );
  }
}; /* mig_resyn_resub_stats */

/*! \brief Interfacing resubstitution functor with MIG resynthesis engines for `window_based_resub_engine`.
 */
template<typename Ntk, typename Simulator, typename TTcare, typename ResynEngine = mig_resyn_topdown<typename Simulator::truthtable_t>>
struct mig_resyn_functor
{
public:
  using node = mig_network::node;
  using signal = mig_network::signal;
  using stats = mig_resyn_resub_stats;
  using TT = typename ResynEngine::truth_table_t;

  static_assert( std::is_same_v<TT, typename Simulator::truthtable_t>, "truth table type of the simulator does not match" );

public:
  explicit mig_resyn_functor( Ntk& ntk, Simulator const& sim, std::vector<node> const& divs, uint32_t num_divs, stats& st )
      : ntk( ntk ), sim( sim ), tts( ntk ), divs( divs ), st( st )
  {
    assert( divs.size() == num_divs );
    (void)num_divs;
    div_signals.reserve( divs.size() );
  }

  std::optional<signal> operator()( node const& root, TTcare care, uint32_t required, uint32_t max_inserts, uint32_t potential_gain, uint32_t& real_gain )
  {
    (void)required;
    TT target = sim.get_tt( sim.get_phase( root ) ? !ntk.make_signal( root ) : ntk.make_signal( root ) );
    TT care_transformed = target.construct();
    care_transformed = care;

    typename ResynEngine::stats st_eng;
    ResynEngine engine( st_eng );
    for ( auto const& d : divs )
    {
      div_signals.emplace_back( sim.get_phase( d ) ? !ntk.make_signal( d ) : ntk.make_signal( d ) );
      tts[d] = sim.get_tt( div_signals.back() );
    }

    auto const res = call_with_stopwatch( st.time_compute_function, [&]() {
      return engine( target, care_transformed, divs.begin(), divs.end(), tts, std::min( potential_gain - 1, max_inserts ) );
    } );
    if ( res )
    {
      ++st.num_success;
      signal ret;
      real_gain = potential_gain - ( *res ).num_gates();
      insert( ntk, div_signals.begin(), div_signals.end(), *res, [&]( signal const& s ) { ret = s; } );
      return ret;
    }
    else
    {
      ++st.num_fail;
      return std::nullopt;
    }
  }

private:
  Ntk& ntk;
  Simulator const& sim;
  unordered_node_map<TT, Ntk> tts;
  std::vector<node> const& divs;
  std::vector<signal> div_signals;
  stats& st;
}; /* mig_resyn_functor */

/*! \brief MIG-specific resubstitution algorithm.
 *
 * This algorithms iterates over each node, creates a
 * reconvergence-driven cut, and attempts to re-express the node's
 * function using existing nodes from the cut.  Node which are no
 * longer used (including nodes in their transitive fanins) can then
 * be removed.  The objective is to reduce the size of the network as
 * much as possible while maintaining the global input-output
 * functionality.
 *
 * **Required network functions:**
 *
 * - `clear_values`
 * - `fanout_size`
 * - `foreach_fanin`
 * - `foreach_fanout`
 * - `foreach_gate`
 * - `foreach_node`
 * - `get_constant`
 * - `get_node`
 * - `is_complemented`
 * - `is_pi`
 * - `level`
 * - `make_signal`
 * - `set_value`
 * - `set_visited`
 * - `size`
 * - `substitute_node`
 * - `value`
 * - `visited`
 *
 * \param ntk A network type derived from mig_network
 * \param ps Resubstitution parameters
 * \param pst Resubstitution statistics
 */
template<class Ntk>
void mig_resubstitution( Ntk& ntk, resubstitution_params const& ps = {}, resubstitution_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( std::is_same_v<typename Ntk::base_type, mig_network>, "Network type is not mig_network" );

  static_assert( has_clear_values_v<Ntk>, "Ntk does not implement the clear_values method" );
  static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement the fanout_size method" );
  static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );
  static_assert( has_foreach_gate_v<Ntk>, "Ntk does not implement the foreach_gate method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_get_constant_v<Ntk>, "Ntk does not implement the get_constant method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_make_signal_v<Ntk>, "Ntk does not implement the make_signal method" );
  static_assert( has_num_gates_v<Ntk>, "Ntk does not implement the num_gates method" );
  static_assert( has_set_value_v<Ntk>, "Ntk does not implement the set_value method" );
  static_assert( has_set_visited_v<Ntk>, "Ntk does not implement the set_visited method" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the has_size method" );
  static_assert( has_substitute_node_v<Ntk>, "Ntk does not implement the has substitute_node method" );
  static_assert( has_value_v<Ntk>, "Ntk does not implement the has_value method" );
  static_assert( has_visited_v<Ntk>, "Ntk does not implement the has_visited method" );
  static_assert( has_level_v<Ntk>, "Ntk does not implement the level method" );
  static_assert( has_foreach_fanout_v<Ntk>, "Ntk does not implement the foreach_fanout method" );

  if ( ps.max_pis == 8 )
  {
    using truthtable_t = kitty::static_truth_table<8u>;
    using truthtable_dc_t = kitty::dynamic_truth_table;
    using functor_t = mig_enumerative_resub_functor<Ntk, detail::window_simulator<Ntk, truthtable_t>, truthtable_dc_t>;
    using resub_impl_t = detail::resubstitution_impl<Ntk, detail::window_based_resub_engine<Ntk, truthtable_t, truthtable_dc_t, functor_t>>;

    resubstitution_stats st;
    typename resub_impl_t::engine_st_t engine_st;
    typename resub_impl_t::collector_st_t collector_st;

    resub_impl_t p( ntk, ps, st, engine_st, collector_st );
    p.run();

    if ( ps.verbose )
    {
      st.report();
      collector_st.report();
      engine_st.report();
    }

    if ( pst )
    {
      *pst = st;
    }
  }
  else
  {
    using truthtable_t = kitty::dynamic_truth_table;
    using truthtable_dc_t = kitty::dynamic_truth_table;
    using functor_t = mig_enumerative_resub_functor<Ntk, detail::window_simulator<Ntk, truthtable_t>, truthtable_dc_t>;
    using resub_impl_t = detail::resubstitution_impl<Ntk, detail::window_based_resub_engine<Ntk, truthtable_t, truthtable_dc_t, functor_t>>;

    resubstitution_stats st;
    typename resub_impl_t::engine_st_t engine_st;
    typename resub_impl_t::collector_st_t collector_st;

    resub_impl_t p( ntk, ps, st, engine_st, collector_st );
    p.run();

    if ( ps.verbose )
    {
      st.report();
      collector_st.report();
      engine_st.report();
    }

    if ( pst )
    {
      *pst = st;
    }
  }
}

/*! \brief MIG-specific resubstitution algorithm.
 *
 * This algorithms iterates over each node, creates a
 * reconvergence-driven cut, and attempts to re-express the node's
 * function using existing nodes from the cut.  Node which are no
 * longer used (including nodes in their transitive fanins) can then
 * be removed.  The objective is to reduce the size of the network as
 * much as possible while maintaining the global input-output
 * functionality.
 *
 * **Required network functions:**
 *
 * - `clear_values`
 * - `fanout_size`
 * - `foreach_fanin`
 * - `foreach_fanout`
 * - `foreach_gate`
 * - `foreach_node`
 * - `get_constant`
 * - `get_node`
 * - `is_complemented`
 * - `is_pi`
 * - `level`
 * - `make_signal`
 * - `set_value`
 * - `set_visited`
 * - `size`
 * - `substitute_node`
 * - `value`
 * - `visited`
 *
 * \param ntk A network type derived from mig_network
 * \param ps Resubstitution parameters
 * \param pst Resubstitution statistics
 */
template<class Ntk>
void mig_resubstitution2( Ntk& ntk, resubstitution_params const& ps = {}, resubstitution_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( std::is_same_v<typename Ntk::base_type, mig_network>, "Network type is not mig_network" );

  static_assert( has_clear_values_v<Ntk>, "Ntk does not implement the clear_values method" );
  static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement the fanout_size method" );
  static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement the foreach_fanin method" );
  static_assert( has_foreach_gate_v<Ntk>, "Ntk does not implement the foreach_gate method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_get_constant_v<Ntk>, "Ntk does not implement the get_constant method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement the is_complemented method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_make_signal_v<Ntk>, "Ntk does not implement the make_signal method" );
  static_assert( has_set_value_v<Ntk>, "Ntk does not implement the set_value method" );
  static_assert( has_set_visited_v<Ntk>, "Ntk does not implement the set_visited method" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the has_size method" );
  static_assert( has_substitute_node_v<Ntk>, "Ntk does not implement the has substitute_node method" );
  static_assert( has_value_v<Ntk>, "Ntk does not implement the has_value method" );
  static_assert( has_visited_v<Ntk>, "Ntk does not implement the has_visited method" );
  static_assert( has_level_v<Ntk>, "Ntk does not implement the level method" );
  static_assert( has_foreach_fanout_v<Ntk>, "Ntk does not implement the foreach_fanout method" );

  using truthtable_t = kitty::dynamic_truth_table;
  using truthtable_dc_t = kitty::dynamic_truth_table;
  using functor_t = mig_resyn_functor<Ntk, detail::window_simulator<Ntk, truthtable_t>, truthtable_dc_t>;

  using resub_impl_t = detail::resubstitution_impl<Ntk, detail::window_based_resub_engine<Ntk, truthtable_t, truthtable_dc_t, functor_t>>;

  resubstitution_stats st;
  typename resub_impl_t::engine_st_t engine_st;
  typename resub_impl_t::collector_st_t collector_st;

  resub_impl_t p( ntk, ps, st, engine_st, collector_st );
  p.run();

  if ( ps.verbose )
  {
    st.report();
    collector_st.report();
    engine_st.report();
  }

  if ( pst )
  {
    *pst = st;
  }
}

} /* namespace mockturtle */
