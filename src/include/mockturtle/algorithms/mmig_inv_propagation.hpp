#pragma once

#include <cstdint>

#include "../utils/stopwatch.hpp"
#include "../views/fanout_view.hpp"
#include "../traits.hpp"

namespace mockturtle
{

struct mmig_inv_propagation_params
{
  bool enable_dual_inversion{ true };
  uint32_t dual_inv_min_fanout{ 2u };
};

struct mmig_inv_propagation_stats
{
  int32_t num_rewritten{ 0 };
  int32_t num_dual_rewritten{ 0 };
  int32_t estimated_gain{ 0 };
  stopwatch<>::duration time_total{ 0 };
};

namespace detail
{

template<class Ntk>
class mmig_inv_propagation_impl
{
public:
  mmig_inv_propagation_impl( Ntk& ntk, mmig_inv_propagation_params const& ps, mmig_inv_propagation_stats& st )
      : ntk( ntk ),
        ps( ps ),
        st( st )
  {
  }

  void run()
  {
    stopwatch t( st.time_total );
    ntk.foreach_gate( [&]( auto const& n ) {
      signal<Ntk> fanins[3];
      uint32_t inv_count = 0;
      ntk.foreach_fanin( n, [&]( auto const& f, auto i ) {
        fanins[i] = f;
        if ( ntk.is_complemented( f ) )
        {
          ++inv_count;
        }
      } );

      if ( inv_count == 3u )
      {
        rewrite_triple_inversion( n, fanins );
      }
      else if ( inv_count == 2u && ps.enable_dual_inversion )
      {
        rewrite_dual_inversion( n, fanins );
      }
    } );
  }

private:
  /* --- Triple inversion ---
   * MAJ(!a,!b,!c) = MIN(a,b,c)   [no output complement needed, saves 3]
   * MIN(!a,!b,!c) = MAJ(a,b,c)   [no output complement needed, saves 3]
   *
   * Fallback for pure-MAJ networks (no MIN gate): uses De Morgan
   *   MAJ(!a,!b,!c) = !MAJ(a,b,c), net saving = 2.
   */
  void rewrite_triple_inversion( node<Ntk> const& n, signal<Ntk> fanins[3] )
  {
    auto const a = !fanins[0];
    auto const b = !fanins[1];
    auto const c = !fanins[2];

    signal<Ntk> replacement{};
    if constexpr ( has_is_min_v<Ntk> && has_create_min_v<Ntk> )
    {
      if ( ntk.is_min( n ) )
      {
        /* MIN(!a,!b,!c) = MAJ(a,b,c) — flip gate type, no output inversion */
        replacement = ntk.create_maj( a, b, c );
      }
      else
      {
        /* MAJ(!a,!b,!c) = MIN(a,b,c) — flip gate type, no output inversion */
        replacement = ntk.create_min( a, b, c );
      }
    }
    else
    {
      /* Pure-MAJ fallback: MAJ(!a,!b,!c) = !MAJ(a,b,c), saves 2 */
      replacement = !ntk.create_maj( a, b, c );
    }

    ntk.substitute_node( n, replacement );
    ntk.replace_in_outputs( n, replacement );
    ++st.num_rewritten;
    st.estimated_gain += 3;
  }

  /* --- Dual inversion: MAJ(!a,!b,c) = !MIN(a,b,!c) ---
   *
   * De Morgan for majority:
   *   MAJ(!x, !y, z)  =  !MIN(x, y, !z)
   *   MIN(!x, !y, z)  =  !MAJ(x, y, !z)
   *
   * Net inversion change:
   *   Before: 2 inverted inputs
   *   After:  1 output inversion + 1 inverted input (on the non-inverted child)
   *   = saves (2 - 1 - 1) = 0 inversions locally
   *
   * BUT: the output inversion is shared across ALL fanouts.
   *   - Before: each consuming gate sees a non-inverted edge from n
   *   - After:  the replacement has 1 output inversion, but at the gate level
   *             this trades 2 input inversions + F non-inverted fanout edges
   *             for 1 input inversion + F inverted fanout edges + 1 new gate type
   *
   * When fanout F >= 2, AND downstream consumers already have complemented
   * connections to this node, the trade is net-positive. We conservatively
   * require F >= dual_inv_min_fanout.
   */
  void rewrite_dual_inversion( node<Ntk> const& n, signal<Ntk> fanins[3] )
  {
    if constexpr ( !has_is_min_v<Ntk> || !has_create_min_v<Ntk> )
    {
      return;
    }

    /* Count downstream references to estimate benefit */
    fanout_view<Ntk> fntk{ ntk };
    uint32_t total_refs = 0u;
    uint32_t complemented_refs = 0u;
    fntk.foreach_fanout( n, [&]( auto const& parent ) {
      fntk.foreach_fanin( parent, [&]( auto const& f ) {
        if ( fntk.get_node( f ) == n )
        {
          ++total_refs;
          if ( fntk.is_complemented( f ) )
          {
            ++complemented_refs;
          }
        }
      } );
    } );

    /* Include PO references */
    ntk.foreach_po( [&]( auto const& s, auto ) {
      if ( ntk.get_node( s ) == n )
      {
        ++total_refs;
        if ( ntk.is_complemented( s ) )
        {
          ++complemented_refs;
        }
      }
    } );

    if ( total_refs < ps.dual_inv_min_fanout )
    {
      return;
    }

    /* Gain calculation:
     *   Remove 2 input inversions, add 1 input inversion on the formerly
     *   non-inverted child => local inversion reduction = +1
     *   Adding one output inversion flips every downstream edge:
     *     non-complemented refs become complemented  -> +1 inversion each
     *     complemented refs become non-complemented  -> -1 inversion each
     *   Fanout contribution to reduction:
     *     complemented_refs - non_complemented_refs
     *   Net reduction:
     *     1 + ( complemented_refs - ( total_refs - complemented_refs ) )
     *   = 1 + 2*complemented_refs - total_refs
     * We want gain > 0:
     *   1 + 2*complemented_refs - total_refs > 0
     */
    int32_t estimated_gain = 1 + 2 * static_cast<int32_t>( complemented_refs ) - static_cast<int32_t>( total_refs );
    if ( estimated_gain <= 0 )
    {
      return;
    }

    /* Identify the two inverted fanins and one non-inverted */
    uint32_t inv_indices[2];
    uint32_t non_inv_index = 0;
    uint32_t ip = 0;
    for ( uint32_t i = 0u; i < 3u; ++i )
    {
      if ( ntk.is_complemented( fanins[i] ) )
      {
        inv_indices[ip++] = i;
      }
      else
      {
        non_inv_index = i;
      }
    }

    /* Strip inversions from the two inverted fanins, add inversion to the
     * non-inverted one */
    auto const a = !fanins[inv_indices[0]]; // strip inversion
    auto const b = !fanins[inv_indices[1]]; // strip inversion
    auto const c = !fanins[non_inv_index];  // add inversion

    /* MAJ(!a,!b,c) = !MIN(a,b,!c)
     * MIN(!a,!b,c) = !MAJ(a,b,!c) */
    signal<Ntk> replacement{};
    if constexpr ( has_is_min_v<Ntk> && has_create_min_v<Ntk> )
    {
      if ( ntk.is_min( n ) )
      {
        /* MIN(!a, !b, c) = !MAJ(a, b, !c) */
        replacement = !ntk.create_maj( a, b, c );
      }
      else
      {
        /* MAJ(!a, !b, c) = !MIN(a, b, !c) */
        replacement = !ntk.create_min( a, b, c );
      }
    }

    ntk.substitute_node( n, replacement );
    ntk.replace_in_outputs( n, replacement );
    ++st.num_dual_rewritten;
    ++st.num_rewritten;
    st.estimated_gain += estimated_gain;
  }

  Ntk& ntk;
  mmig_inv_propagation_params const& ps;
  mmig_inv_propagation_stats& st;
};

} // namespace detail

template<class Ntk>
void mmig_inv_propagation( Ntk& ntk,
                           mmig_inv_propagation_params const& ps = {},
                           mmig_inv_propagation_stats* pst = nullptr )
{
  static_assert( has_foreach_gate_v<Ntk>, "Ntk does not implement foreach_gate" );
  static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement foreach_fanin" );
  static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement is_complemented" );
  static_assert( has_create_maj_v<Ntk>, "Ntk does not implement create_maj" );
  static_assert( has_substitute_node_v<Ntk>, "Ntk does not implement substitute_node" );
  static_assert( has_replace_in_outputs_v<Ntk>, "Ntk does not implement replace_in_outputs" );

  mmig_inv_propagation_stats st;
  detail::mmig_inv_propagation_impl<Ntk> impl( ntk, ps, st );
  impl.run();
  if ( pst != nullptr )
  {
    *pst = st;
  }
}

template<class Ntk>
void mmig_inv_propagation( Ntk& ntk, mmig_inv_propagation_stats* pst )
{
  mmig_inv_propagation( ntk, mmig_inv_propagation_params{}, pst );
}

} // namespace mockturtle
