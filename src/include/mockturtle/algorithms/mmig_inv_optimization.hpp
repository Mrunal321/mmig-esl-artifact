#pragma once

#include <cstdint>

#include "../views/fanout_view.hpp"
#include "../utils/stopwatch.hpp"
#include "../traits.hpp"

namespace mockturtle
{

struct mmig_inv_optimization_stats
{
  int32_t num_inverted{ 0 };
  int32_t total_gain{ 0 };
  stopwatch<>::duration time_total{ 0 };
};

namespace detail
{

template<class Ntk>
class mmig_inv_optimization_impl
{
public:
  mmig_inv_optimization_impl( Ntk& ntk, mmig_inv_optimization_stats& st )
      : ntk( ntk ),
        st( st )
  {
  }

  void run()
  {
    stopwatch t( st.time_total );
    bool changed = true;
    while ( changed )
    {
      changed = false;
      ntk.foreach_gate( [&]( auto const& n ) {
        const auto gain = calculate_gain( n );
        if ( gain > 0 )
        {
          if ( invert_node( n ) )
          {
            st.num_inverted++;
            st.total_gain += gain;
            changed = true;
          }
        }
      } );
    }
  }

private:
  int32_t calculate_gain( node<Ntk> n )
  {
    int32_t gain = 0;

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      if ( ntk.is_constant( ntk.get_node( f ) ) )
      {
        return;
      }
      gain += ntk.is_complemented( f ) ? 1 : -1;
    } );

    ntk.foreach_fanout( n, [&]( auto const& parent ) {
      gain += is_complemented_parent( parent, n ) ? 1 : -1;
    } );

    ntk.foreach_po( [&]( auto const& s, auto ) {
      if ( ntk.get_node( s ) == n )
      {
        gain += ntk.is_complemented( s ) ? 1 : -1;
      }
    } );

    return gain;
  }

  bool is_complemented_parent( node<Ntk> parent, node<Ntk> child )
  {
    bool value = false;
    ntk.foreach_fanin( parent, [&]( auto const& f ) {
      if ( ntk.get_node( f ) == child )
      {
        value = ntk.is_complemented( f );
      }
    } );
    return value;
  }

  bool invert_node( node<Ntk> n )
  {
    signal<Ntk> fanins[3]{};
    ntk.foreach_fanin( n, [&]( auto const& f, auto i ) { fanins[i] = f; } );

    auto const a = !fanins[0];
    auto const b = !fanins[1];
    auto const c = !fanins[2];

    signal<Ntk> replacement{};
    if constexpr ( has_is_min_v<Ntk> && has_create_min_v<Ntk> )
    {
      if ( ntk.is_min( n ) )
      {
        replacement = ntk.create_min( a, b, c );
      }
      else
      {
        replacement = ntk.create_maj( a, b, c );
      }
    }
    else
    {
      replacement = ntk.create_maj( a, b, c );
    }

    auto const rewritten = !replacement;
    if ( rewritten == ntk.make_signal( n ) )
    {
      return false;
    }

    ntk.substitute_node( n, rewritten );
    ntk.replace_in_outputs( n, rewritten );
    return true;
  }

private:
  fanout_view<Ntk> ntk;
  mmig_inv_optimization_stats& st;
};

} // namespace detail

template<class Ntk>
void mmig_inv_optimization( Ntk& ntk, mmig_inv_optimization_stats* pst = nullptr )
{
  static_assert( has_foreach_gate_v<Ntk>, "Ntk does not implement foreach_gate" );
  static_assert( has_foreach_fanin_v<Ntk>, "Ntk does not implement foreach_fanin" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement foreach_po" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement is_constant" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement get_node" );
  static_assert( has_is_complemented_v<Ntk>, "Ntk does not implement is_complemented" );
  static_assert( has_make_signal_v<Ntk>, "Ntk does not implement make_signal" );
  static_assert( has_create_maj_v<Ntk>, "Ntk does not implement create_maj" );
  static_assert( !has_is_min_v<Ntk> || has_create_min_v<Ntk>, "Ntk implements is_min but not create_min" );
  static_assert( has_substitute_node_v<Ntk>, "Ntk does not implement substitute_node" );
  static_assert( has_replace_in_outputs_v<Ntk>, "Ntk does not implement replace_in_outputs" );

  mmig_inv_optimization_stats st;
  detail::mmig_inv_optimization_impl<Ntk> impl( ntk, st );
  impl.run();
  if ( pst != nullptr )
  {
    *pst = st;
  }
}

} // namespace mockturtle
