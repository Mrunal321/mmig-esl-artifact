#include <catch.hpp>

#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/equivalence_checking.hpp>
#include <mockturtle/algorithms/miter.hpp>
#include <mockturtle/algorithms/mmig_inv_propagation.hpp>
#include <mockturtle/networks/mig.hpp>

using namespace mockturtle;

namespace
{

void check_equivalent( mig_network const& a, mig_network const& b )
{
  auto const m = miter<mig_network>( a, b );
  REQUIRE( m.has_value() );
  auto const eq = equivalence_checking( *m );
  REQUIRE( eq.has_value() );
  CHECK( *eq );
}

} // namespace

/* Helper: count complement edges across the whole network */
static uint32_t count_inv_edges( mig_network const& ntk )
{
  uint32_t inv = 0;
  ntk.foreach_gate( [&]( auto const& n ) {
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      if ( ntk.is_complemented( f ) ) ++inv;
    } );
  } );
  ntk.foreach_po( [&]( auto const& f ) {
    if ( ntk.is_complemented( f ) ) ++inv;
  } );
  return inv;
}

/* Build a MAJ node with exactly 3 complement fanins at the edge level,
 * bypassing create_maj normalisation (which folds >=2 complemented inputs
 * into a single output complement).
 *
 * Strategy:
 *   1. Create three DISTINCT non-trivial nodes d, e, f from disjoint PI sets.
 *   2. Create the central node n = MAJ(d, e, f) — 0 complement fanins.
 *   3. Create distinct replacement nodes d2, e2, f2.
 *   4. substitute_node(d, !d2) flips the complement bit on n's d-edge.
 *      Repeat for e and f.  Result: n has 3 complement fanins without
 *      triggering create_maj's normalisation path.
 *
 * PIs used: a,b,c,g,h,k (6 inputs total).
 */
static mig_network::signal force_triple_inv_maj( mig_network& ntk,
                                                  mig_network::signal a,
                                                  mig_network::signal b,
                                                  mig_network::signal c,
                                                  mig_network::signal g,
                                                  mig_network::signal h,
                                                  mig_network::signal k )
{
  /* Three structurally distinct nodes so none collapse to a simpler signal */
  auto const d = ntk.create_maj( a, b, c );
  auto const e = ntk.create_maj( a, b, g );
  auto const f = ntk.create_maj( a, b, h );

  /* Central node: d, e, f all non-complemented fanins */
  auto const n = ntk.create_maj( d, e, f );

  /* Distinct replacements (different PIs so not hash-equal to d/e/f) */
  auto const d2 = ntk.create_maj( b, c, k );
  auto const e2 = ntk.create_maj( b, g, k );
  auto const f2 = ntk.create_maj( b, h, k );

  /* Each substitution flips the complement bit on n's corresponding fanin */
  ntk.substitute_node( ntk.get_node( d ), !d2 );
  ntk.substitute_node( ntk.get_node( e ), !e2 );
  ntk.substitute_node( ntk.get_node( f ), !f2 );

  return n;
}

TEST_CASE( "mmig triple inversion: MAJ rewrites to MIN saving 3 inverters (fanout-1)", "[mmig_inv_propagation]" )
{
  /* Force MAJ(!d,!e,!f) at the edge level (bypassing normalisation) then run
   * inv_propagation and verify: function preserved, net saving >= 3,
   * and the replacement is a MIN gate. */
  mig_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const g = ntk.create_pi();
  auto const h = ntk.create_pi();
  auto const k = ntk.create_pi();

  auto const n = force_triple_inv_maj( ntk, a, b, c, g, h, k );
  ntk.create_po( n );
  /* NOTE: no cleanup_dangling here — that would rebuild via create_maj and
   * re-normalise the forced triple-inversion pattern away. */

  /* Verify the central node really has 3 complement fanins */
  uint32_t pre_inv_on_n = 0;
  ntk.foreach_gate( [&]( auto const& nd ) {
    uint32_t cnt = 0;
    ntk.foreach_fanin( nd, [&]( auto const& f ) { if ( ntk.is_complemented( f ) ) ++cnt; } );
    if ( cnt == 3u ) ++pre_inv_on_n;
  } );
  REQUIRE( pre_inv_on_n >= 1u );   /* at least one triple-inv node exists */

  auto const original = ntk.clone();
  uint32_t inv_before = count_inv_edges( ntk );

  mmig_inv_propagation_stats st{};
  mmig_inv_propagation( ntk, &st );
  ntk = cleanup_dangling( ntk );

  uint32_t inv_after = count_inv_edges( ntk );

  check_equivalent( original, ntk );
  CHECK( st.num_rewritten >= 1 );

  /* New MIN-gate rewrite saves 3 per node; no output complement added.
   * Old !MAJ rewrite saved only 2 (added 1 output complement). */
  CHECK( (int32_t)inv_before - (int32_t)inv_after >= 3 );

  /* At least one MIN gate must exist after the rewrite */
  uint32_t min_count = 0;
  ntk.foreach_gate( [&]( auto const& nd ) { if ( ntk.is_min( nd ) ) ++min_count; } );
  CHECK( min_count >= 1u );
}

TEST_CASE( "mmig triple inversion: high-fanout MAJ saves 3 inverters (not 0 as old code)", "[mmig_inv_propagation]" )
{
  /* With the OLD !MAJ code, a triple-inv node with fanout=3 would add 3 output
   * complement edges, exactly cancelling the 3 removed input edges — net 0.
   * With the NEW MIN-gate code, no output edges are added — net +3 always. */
  mig_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const g = ntk.create_pi();
  auto const h = ntk.create_pi();
  auto const k = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const e = ntk.create_pi();

  auto const n  = force_triple_inv_maj( ntk, a, b, c, g, h, k );
  auto const p1 = ntk.create_maj( n, d, e );
  auto const p2 = ntk.create_maj( n, a, b );
  auto const p3 = ntk.create_maj( n, c, d );
  ntk.create_po( p1 );
  ntk.create_po( p2 );
  ntk.create_po( p3 );
  /* No cleanup_dangling before inv_propagation — it would re-normalise. */

  auto const original = ntk.clone();
  uint32_t inv_before = count_inv_edges( ntk );

  mmig_inv_propagation_stats st{};
  mmig_inv_propagation( ntk, &st );
  ntk = cleanup_dangling( ntk );

  uint32_t inv_after = count_inv_edges( ntk );

  check_equivalent( original, ntk );
  CHECK( st.num_rewritten >= 1 );
  /* MIN-gate rewrite: 3 input invs gone, 0 output invs added regardless of fanout */
  CHECK( (int32_t)inv_before - (int32_t)inv_after >= 3 );
}

TEST_CASE( "mmig inv propagation dual-inversion skips when fanout polarity would regress", "[mmig_inv_propagation]" )
{
  mig_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const e = ntk.create_pi();

  auto const n = ntk.create_maj( !a, !b, c );
  auto const p1 = ntk.create_maj( n, d, e );
  auto const p2 = ntk.create_maj( n, a, d );
  ntk.create_po( p1 );
  ntk.create_po( p2 );

  auto const original = ntk.clone();

  mmig_inv_propagation_params ps;
  ps.enable_dual_inversion = true;
  ps.dual_inv_min_fanout = 2u;

  mmig_inv_propagation_stats st{};
  mmig_inv_propagation( ntk, ps, &st );

  CHECK( st.num_dual_rewritten == 0 );
  CHECK( st.num_rewritten == 0 );
  check_equivalent( original, ntk );
}

TEST_CASE( "mmig inv propagation keeps legacy stats-only call signature", "[mmig_inv_propagation]" )
{
  mig_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const e = ntk.create_pi();

  auto const n = ntk.create_maj( !a, !b, c );
  auto const p1 = ntk.create_maj( !n, d, e );
  auto const p2 = ntk.create_min( !n, a, d );
  ntk.create_po( p1 );
  ntk.create_po( p2 );
  ntk.create_po( !n );

  auto const original = ntk.clone();

  mmig_inv_propagation_stats st{};
  mmig_inv_propagation( ntk, &st );

  CHECK( st.num_rewritten >= 0 );
  CHECK( st.num_dual_rewritten >= 0 );
  check_equivalent( original, ntk );
}
