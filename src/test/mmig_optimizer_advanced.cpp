#include <catch.hpp>

#include <mockturtle/algorithms/equivalence_checking.hpp>
#include <mockturtle/algorithms/miter.hpp>
#include <mockturtle/algorithms/mmig_optimizer.hpp>
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

mig_network build_mixed_network()
{
  mig_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();
  auto const d = ntk.create_pi();
  auto const e = ntk.create_pi();

  auto const n1 = ntk.create_maj( a, b, c );
  auto const n2 = ntk.create_min( n1, d, !c );
  auto const n3 = ntk.create_maj( n2, e, b );
  auto const n4 = ntk.create_min( n3, !a, d );
  auto const n5 = ntk.create_maj( n4, n2, a );
  ntk.create_po( n3 );
  ntk.create_po( n4 );
  ntk.create_po( n5 );
  return ntk;
}

mig_network build_simple_network()
{
  mig_network ntk;
  auto const a = ntk.create_pi();
  auto const b = ntk.create_pi();
  auto const c = ntk.create_pi();

  auto const n = ntk.create_maj( a, b, c );
  ntk.create_po( n );
  return ntk;
}

} // namespace

TEST_CASE( "mmig optimizer advanced disabled keeps baseline behavior", "[mmig_optimizer]" )
{
  auto input = build_mixed_network();
  auto ntk_default = input.clone();
  auto ntk_explicit = input.clone();

  mmig_optimizer_params ps_default;
  ps_default.max_iterations = 2u;
  ps_default.enable_minority_seeding = true;
  ps_default.run_cec = true;
  ps_default.enable_advanced = false;

  mmig_optimizer_stats st_default{};
  mmig_depth_optimization( ntk_default, ps_default, &st_default );

  auto ps_explicit = ps_default;
  ps_explicit.enable_advanced = false;

  mmig_optimizer_stats st_explicit{};
  mmig_depth_optimization( ntk_explicit, ps_explicit, &st_explicit );

  CHECK( ntk_default.num_gates() == ntk_explicit.num_gates() );
  CHECK( st_default.cec.rejected == st_explicit.cec.rejected );
  check_equivalent( ntk_default, ntk_explicit );
}

TEST_CASE( "mmig optimizer advanced mode runs with CEC and preserves equivalence", "[mmig_optimizer]" )
{
  auto original = build_mixed_network();
  auto ntk = original.clone();

  mmig_optimizer_params ps;
  ps.max_iterations = 2u;
  ps.enable_minority_seeding = true;
  ps.run_cec = true;
  ps.enable_advanced = true;
  ps.advanced_rounds = 1u;
  ps.advanced_run_resub = true;
  ps.advanced_run_cut = true;
  ps.advanced_run_refactor = true;

  mmig_optimizer_stats st{};
  mmig_depth_optimization( ntk, ps, &st );

  CHECK( st.cec.accepted > 0u );
  CHECK( st.advanced_majority_before + st.advanced_minority_before == st.advanced_gates_before );
  CHECK( st.advanced_majority_after + st.advanced_minority_after == st.advanced_gates_after );
  check_equivalent( original, ntk );
}

TEST_CASE( "mmig optimizer advanced stagnation stop avoids useless rounds", "[mmig_optimizer]" )
{
  auto original = build_simple_network();
  auto ntk = original.clone();

  mmig_optimizer_params ps;
  ps.max_iterations = 1u;
  ps.enable_minority_seeding = false;
  ps.run_cec = true;
  ps.enable_advanced = true;
  ps.advanced_rounds = 4u;
  ps.advanced_flow = mmig_advanced_flow::round_robin;
  ps.advanced_stop_on_stagnation = true;
  ps.advanced_max_stagnation_rounds = 1u;
  ps.advanced_rollback_on_objective_regression = true;

  mmig_optimizer_stats st{};
  mmig_depth_optimization( ntk, ps, &st );

  CHECK( st.advanced_rounds_executed == 1u );
  CHECK( st.advanced_rounds_improved == 0u );
  CHECK( st.advanced_stagnation_stops == 1u );
  check_equivalent( original, ntk );
}

TEST_CASE( "mmig optimizer exact stage gating skips expensive exact when configured", "[mmig_optimizer]" )
{
  auto original = build_mixed_network();
  auto ntk = original.clone();

  mmig_optimizer_params ps;
  ps.max_iterations = 1u;
  ps.enable_minority_seeding = false;
  ps.run_cec = true;
  ps.enable_advanced = true;
  ps.advanced_rounds = 1u;
  ps.advanced_flow = mmig_advanced_flow::round_robin;
  ps.advanced_run_resub = false;
  ps.advanced_run_cut = false;
  ps.advanced_run_refactor = false;
  ps.advanced_run_exact = true;
  ps.advanced_run_balance = false;
  ps.advanced_tuned_policy = true;
  ps.advanced_exact_gate_limit = 1u;

  mmig_optimizer_stats st{};
  mmig_depth_optimization( ntk, ps, &st );

  CHECK( st.advanced_exact.num_gate_skips > 0u );
  check_equivalent( original, ntk );
}

TEST_CASE( "mmig optimizer interleaved seeding handles zero interval safely", "[mmig_optimizer]" )
{
  auto original = build_mixed_network();
  auto ntk = original.clone();

  mmig_optimizer_params ps;
  ps.max_iterations = 1u;
  ps.enable_minority_seeding = false;
  ps.run_cec = true;
  ps.enable_advanced = true;
  ps.advanced_rounds = 2u;
  ps.advanced_flow = mmig_advanced_flow::round_robin;
  ps.enable_interleaved_seeding = true;
  ps.interleaved_seeding_interval = 0u;

  mmig_optimizer_stats st{};
  mmig_depth_optimization( ntk, ps, &st );

  CHECK( st.advanced_rounds_executed > 0u );
  check_equivalent( original, ntk );
}
