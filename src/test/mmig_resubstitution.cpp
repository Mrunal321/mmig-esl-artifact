#include <catch.hpp>

#include <mockturtle/algorithms/equivalence_checking.hpp>
#include <mockturtle/algorithms/miter.hpp>
#include <mockturtle/algorithms/mmig_resubstitution.hpp>
#include <mockturtle/networks/mig.hpp>

using namespace mockturtle;

namespace
{

uint32_t count_min_nodes( mig_network const& ntk )
{
  uint32_t count = 0u;
  ntk.foreach_gate( [&]( auto const& n ) {
    if ( ntk.is_min( n ) )
    {
      ++count;
    }
  } );
  return count;
}

void check_equivalent( mig_network const& a, mig_network const& b )
{
  auto const m = miter<mig_network>( a, b );
  REQUIRE( m.has_value() );
  auto const eq = equivalence_checking( *m );
  REQUIRE( eq.has_value() );
  CHECK( *eq );
}

} // namespace

TEST_CASE( "mmig resubstitution materializes polarity-aware MIN candidates", "[mmig_resubstitution]" )
{
  mig_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();

  auto const n1 = ntk.create_maj( !a, !b, !c );
  auto const root = ntk.create_maj( !a, !b, n1 ); // functionally equal to MIN(a,b,c)
  ntk.create_po( root );

  auto const original = ntk.clone();
  auto const min_before = count_min_nodes( ntk );

  mmig_resubstitution_params ps;
  ps.enable_candidate_cec = true;
  ps.verify_mixed_risk_only = false;
  ps.resub_ps.max_pis = 4u;
  ps.resub_ps.max_divisors = 64u;
  ps.resub_ps.max_inserts = 1u;
  ps.resub_ps.preserve_depth = false;
  ps.resub_ps.rank_mig_candidates = true;

  mmig_resubstitution_stats st{};
  mmig_resubstitution( ntk, ps, &st );

  CHECK( st.num_candidate_accepted > 0u );
  CHECK( count_min_nodes( ntk ) > min_before );
  check_equivalent( original, ntk );
}

TEST_CASE( "mmig resubstitution preserves equivalence with mixed-risk candidate CEC", "[mmig_resubstitution]" )
{
  mig_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const e = ntk.create_pi();

  auto const x = ntk.create_maj( a, b, c );
  auto const y = ntk.create_maj( a, b, x ); // functionally equal to x
  auto const z = ntk.create_min( y, d, !c ); // mixed fanout for y
  auto const u = ntk.create_maj( a, d, b );
  auto const w = ntk.create_maj( a, d, u ); // functionally equal to u (non-mixed)
  auto const f = ntk.create_maj( w, z, e );

  ntk.create_po( f );
  ntk.create_po( y );
  ntk.create_po( w );

  auto const original = ntk.clone();

  mmig_resubstitution_params ps;
  ps.enable_candidate_cec = true;
  ps.verify_mixed_risk_only = true;
  ps.resub_ps.max_pis = 8u;
  ps.resub_ps.max_divisors = 128u;
  ps.resub_ps.max_inserts = 2u;

  mmig_resubstitution_stats st{};
  mmig_resubstitution( ntk, ps, &st );

  CHECK( st.num_candidate_checks > 0u );
  CHECK( st.num_candidate_accepted > 0u );
  CHECK( st.num_candidate_skipped_non_mixed > 0u );
  check_equivalent( original, ntk );
}

TEST_CASE( "mmig resubstitution selective mixed-only CEC skips non-mixed candidates", "[mmig_resubstitution]" )
{
  mig_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();

  auto const u = ntk.create_maj( a, b, c );
  auto const w = ntk.create_maj( a, b, u ); // functionally equal to u
  auto const f = ntk.create_maj( w, d, c );
  ntk.create_po( f );
  ntk.create_po( w );

  auto const original = ntk.clone();

  mmig_resubstitution_params ps;
  ps.enable_candidate_cec = true;
  ps.verify_mixed_risk_only = true;
  ps.resub_ps.max_pis = 8u;
  ps.resub_ps.max_divisors = 64u;

  mmig_resubstitution_stats st{};
  mmig_resubstitution( ntk, ps, &st );

  CHECK( st.num_candidate_checks == 0u );
  CHECK( st.num_candidate_skipped_non_mixed > 0u );
  check_equivalent( original, ntk );
}

TEST_CASE( "mmig resubstitution respects candidate CEC budget", "[mmig_resubstitution]" )
{
  mig_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const e = ntk.create_pi();

  auto const n1 = ntk.create_maj( a, b, c );
  auto const n2 = ntk.create_maj( a, b, n1 );
  auto const n3 = ntk.create_min( n2, d, !c );
  auto const n4 = ntk.create_maj( n3, e, b );
  auto const n5 = ntk.create_min( n4, a, !d );
  ntk.create_po( n2 );
  ntk.create_po( n5 );

  auto const original = ntk.clone();

  mmig_resubstitution_params ps;
  ps.enable_candidate_cec = true;
  ps.verify_mixed_risk_only = true;
  ps.max_candidate_cec_checks = 1u;
  ps.skip_candidate_when_cec_budget_exhausted = true;
  ps.resub_ps.max_pis = 8u;
  ps.resub_ps.max_divisors = 128u;
  ps.resub_ps.max_inserts = 2u;

  mmig_resubstitution_stats st{};
  mmig_resubstitution( ntk, ps, &st );

  CHECK( st.num_candidate_checks > 0u );
  CHECK( st.num_candidate_checks <= 1u );
  check_equivalent( original, ntk );
}

TEST_CASE( "mmig resubstitution keeps dont-cares safe for large cuts", "[mmig_resubstitution]" )
{
  mig_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const e = ntk.create_pi();
  auto const f = ntk.create_pi();

  auto const n1 = ntk.create_maj( a, b, c );
  auto const n2 = ntk.create_maj( n1, d, e );
  auto const n3 = ntk.create_min( n2, !c, f );
  auto const n4 = ntk.create_maj( n3, a, !d );
  auto const n5 = ntk.create_min( n4, b, e );
  ntk.create_po( n3 );
  ntk.create_po( n5 );

  auto const original = ntk.clone();

  mmig_resubstitution_params ps;
  ps.enable_candidate_cec = false;
  ps.verify_mixed_risk_only = true;
  ps.resub_ps.max_pis = 10u;
  ps.resub_ps.max_divisors = 128u;
  ps.resub_ps.max_inserts = 2u;
  ps.resub_ps.use_dont_cares = true;

  mmig_resubstitution_stats st{};
  mmig_resubstitution( ntk, ps, &st );

  CHECK( st.time_total.count() >= 0 );
  check_equivalent( original, ntk );
}
