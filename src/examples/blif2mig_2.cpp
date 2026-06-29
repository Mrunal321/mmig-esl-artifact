#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <chrono>
#include <utility>
#include <queue>
#include <set>
#include <algorithm>
#include <iomanip>
#include <unordered_map>
#include <string_view>
#include <sstream>
#include <optional>
#include <cassert>
#include <fstream>
#include <cctype>
#include <kitty/print.hpp>
#include <filesystem>
#include <cstdlib>

#include <mockturtle/mockturtle.hpp>
#include <mockturtle/algorithms/mig_inv_propagation.hpp>
#include <mockturtle/algorithms/mig_inv_optimization.hpp>
#include <mockturtle/algorithms/mig_algebraic_rewriting.hpp>
#include <mockturtle/algorithms/mig_resub.hpp>
#include <mockturtle/algorithms/balancing.hpp>
#include <mockturtle/algorithms/mmig_optimizer.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/fanout_view.hpp>
#include <mockturtle/views/topo_view.hpp>
#include <mockturtle/io/write_verilog.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/algorithms/cut_rewriting.hpp>
#include <mockturtle/algorithms/cut_enumeration.hpp>
#include <mockturtle/algorithms/refactoring.hpp>
#include <mockturtle/algorithms/resubstitution.hpp>
#include <mockturtle/views/mapping_view.hpp>
#include <mockturtle/algorithms/lut_mapping.hpp>
#include <mockturtle/algorithms/collapse_mapped.hpp>
#include <mockturtle/algorithms/klut_to_graph.hpp>
#include <lorina/blif.hpp>

using namespace mockturtle;

#ifndef BLIF2MIG_NETWORK_METRICS_DEFINED
#define BLIF2MIG_NETWORK_METRICS_DEFINED
struct network_metrics
{
  uint32_t pis = 0;
  uint32_t pos = 0;
  uint32_t gates = 0;
  uint32_t majority_nodes = 0;
  uint32_t minority_nodes = 0;
  uint32_t depth = 0;
  uint64_t edges_total = 0;
  uint64_t edges_internal = 0;
  uint64_t edges_po = 0;
  uint64_t inv_edges_total = 0;
  uint64_t inv_edges_internal = 0;
  uint64_t inv_edges_po = 0;
  uint32_t nodes_with_inverted_fanout = 0;
};
#endif

struct lut_metrics
{
  uint32_t lut_count = 0;
  uint32_t depth = 0;
  uint32_t pis = 0;
  std::vector<uint32_t> fanin_histogram;
};

static uint32_t count_nodes_with_inverted_fanout( mig_network const& mig )
{
  std::vector<uint8_t> inverted( mig.size(), 0 );

  mig.foreach_gate( [&]( auto target ) {
    mig.foreach_fanin( target, [&]( auto const& fi ) {
      if ( mig.is_complemented( fi ) )
      {
        const auto source = mig.get_node( fi );
        inverted[mig.node_to_index( source )] = 1;
      }
    } );
  } );

  mig.foreach_po( [&]( auto const& s, auto ) {
    if ( mig.is_complemented( s ) )
    {
      const auto source = mig.get_node( s );
      inverted[mig.node_to_index( source )] = 1;
    }
  } );

  uint32_t count = 0;
  mig.foreach_node( [&]( auto n ) {
    if ( mig.is_constant( n ) )
    {
      return;
    }
    if ( inverted[mig.node_to_index( n )] )
    {
      ++count;
    }
  } );

  return count;
}

static network_metrics compute_network_metrics( mig_network const& mig )
{
  network_metrics stats{};
  stats.pis = mig.num_pis();
  stats.pos = mig.num_pos();
  stats.gates = mig.num_gates();
  mig.foreach_gate( [&]( auto n ) {
    if ( mig.is_min( n ) )
    {
      ++stats.minority_nodes;
    }
    else
    {
      ++stats.majority_nodes;
    }
  } );

  depth_view dv{ mig };
  stats.depth = dv.depth();

  uint64_t edges_internal = 0;
  uint64_t inv_edges_internal = 0;
  mig.foreach_gate( [&]( auto n ) {
    mig.foreach_fanin( n, [&]( auto const& fi ) {
      ++edges_internal;
      if ( mig.is_complemented( fi ) )
      {
        ++inv_edges_internal;
      }
    } );
  } );

  uint64_t edges_po = 0;
  uint64_t inv_edges_po = 0;
  mig.foreach_po( [&]( auto const& s, auto ) {
    ++edges_po;
    if ( mig.is_complemented( s ) )
    {
      ++inv_edges_po;
    }
  } );

  stats.edges_internal = edges_internal;
  stats.inv_edges_internal = inv_edges_internal;
  stats.edges_po = edges_po;
  stats.inv_edges_po = inv_edges_po;
  stats.edges_total = edges_internal + edges_po;
  stats.inv_edges_total = inv_edges_internal + inv_edges_po;
  stats.nodes_with_inverted_fanout = count_nodes_with_inverted_fanout( mig );

  return stats;
}

static bool g_use_exact_lut_mapping = false;
static int g_exact_lut_window = 0;
static bool g_resub_use_dont_cares = true;
static int g_resub_dc_window = 12;

static void print_stats( const mig_network& mig, const std::string& label )
{
  depth_view d{ mig };
  fanout_view f{ mig };

  auto metrics = compute_network_metrics( mig );
  const double inv_ratio = metrics.edges_total ? static_cast<double>( metrics.inv_edges_total ) / static_cast<double>( metrics.edges_total ) : 0.0;

  uint32_t max_fo = 0, ge2 = 0, ge3 = 0;
  uint64_t sum_fo = 0;
  f.foreach_gate( [&]( auto n ) {
    int fo = 0;
    f.foreach_fanout( n, [&]( auto const& ) { ++fo; } );
    sum_fo += fo;
    max_fo = std::max<uint32_t>( max_fo, fo );
    if ( fo >= 2 )
      ++ge2;
    if ( fo >= 3 )
      ++ge3;
  } );
  const double avg_fo = metrics.gates ? static_cast<double>( sum_fo ) / static_cast<double>( metrics.gates ) : 0.0;

  std::vector<uint32_t> lvl_counts( metrics.depth + 1, 0 );
  uint64_t lvl_sum = 0;
  d.foreach_gate( [&]( auto n ) {
    const auto L = d.level( n );
    lvl_counts[L]++;
    lvl_sum += L;
  } );
  const double avg_lvl = metrics.gates ? static_cast<double>( lvl_sum ) / static_cast<double>( metrics.gates ) : 0.0;
  const uint32_t cpw = lvl_counts.empty() ? 0u : lvl_counts.back();

  std::cout << "\n=== " << label << " ===\n";
  std::cout << "PIs=" << metrics.pis << "  POs=" << metrics.pos
            << "  Gates=" << metrics.gates << "  Depth=" << metrics.depth << "\n";
  std::cout << "GateTypes: MAJ=" << metrics.majority_nodes << "  MIN=" << metrics.minority_nodes << "\n";
  std::cout << "Edges(total)=" << metrics.edges_total
            << " (internal=" << metrics.edges_internal << ", PO=" << metrics.edges_po << ")  "
            << "InvEdges(total)=" << metrics.inv_edges_total
            << " (internal=" << metrics.inv_edges_internal << ", PO=" << metrics.inv_edges_po << ")"
            << "  InvEdgeRatio=" << inv_ratio << "\n";
  std::cout << "Nodes with inverted fanout=" << metrics.nodes_with_inverted_fanout << "\n";
  std::cout << "Fanout: max=" << max_fo << "  avg=" << avg_fo
            << "  nodes(fo>=2)=" << ge2 << "  nodes(fo>=3)=" << ge3 << "\n";
  std::cout << "Levels: avg=" << avg_lvl << "  crit_path_width=" << cpw << "\n";
}

// ---------- small helpers (GLOBAL SCOPE!) ----------
#ifndef BLIF2MIG_SCOPED_TIMER_DEFINED
#define BLIF2MIG_SCOPED_TIMER_DEFINED
struct scoped_timer {
  std::string name;
  std::chrono::steady_clock::time_point t0;
  bool enabled;
  explicit scoped_timer( std::string n, bool enabled_ = true )
      : name( std::move( n ) ),
        t0( std::chrono::steady_clock::now() ),
        enabled( enabled_ ) {
    if ( enabled ) {
      std::cout << "[pass] " << name << " ...\n";
    }
  }
  ~scoped_timer() {
    if ( enabled ) {
      using namespace std::chrono;
      auto ms = duration_cast<milliseconds>( steady_clock::now() - t0 ).count();
      std::cout << "[pass] " << name << " done in " << ms << " ms\n";
    }
  }
};
#endif

static bool has_flag(int argc, char** argv, const char* f) {
  for (int i = 3; i < argc; ++i)
    if (std::string(argv[i]) == f) return true;
  return false;
}

static std::optional<int> parse_int_option( int argc, char** argv, std::string_view prefix )
{
  for ( int i = 3; i < argc; ++i )
  {
    std::string_view arg{ argv[i] };
    if ( arg.substr( 0, prefix.size() ) == prefix )
    {
      auto value_str = arg.substr( prefix.size() );
      if ( value_str.empty() )
        return std::nullopt;
      try
      {
        size_t idx = 0;
        const int value = std::stoi( std::string( value_str ), &idx );
        if ( idx == value_str.size() )
        {
          return value;
        }
      }
      catch ( ... )
      {
        return std::nullopt;
      }
    }
  }
  return std::nullopt;
}

static std::optional<std::string> parse_string_option( int argc, char** argv, std::string_view prefix )
{
  for ( int i = 3; i < argc; ++i )
  {
    std::string_view arg{ argv[i] };
    if ( arg.substr( 0, prefix.size() ) == prefix )
    {
      return std::string( arg.substr( prefix.size() ) );
    }
  }
  return std::nullopt;
}

static std::string format_ms_value( double ms );

static mmig_stage_mode parse_mmig_stage_mode( std::string const& value )
{
  if ( value == "pre" )
  {
    return mmig_stage_mode::pre;
  }
  if ( value == "post" )
  {
    return mmig_stage_mode::post;
  }
  return mmig_stage_mode::both;
}

static mmig_seed_mode parse_mmig_seed_mode( std::string const& value )
{
  if ( value == "l1" )
    return mmig_seed_mode::level1;
  if ( value == "l2" )
    return mmig_seed_mode::level2;
  if ( value == "l3" )
    return mmig_seed_mode::level3;
  return mmig_seed_mode::both;
}

enum class majority_flow_kind : uint8_t
{
  identity,
  standard,
  dac19_default,
  dac19_area,
  dac19_compat,
  compress2rs,
  legacy
};

static majority_flow_kind parse_majority_flow_kind( std::string const& value )
{
  if ( value == "identity" || value == "none" )
  {
    return majority_flow_kind::identity;
  }
  if ( value == "dac19_default" )
  {
    return majority_flow_kind::dac19_default;
  }
  if ( value == "dac19_area" )
  {
    return majority_flow_kind::dac19_area;
  }
  if ( value == "dac19_compat" )
  {
    return majority_flow_kind::dac19_compat;
  }
  if ( value == "compress2rs" )
  {
    return majority_flow_kind::compress2rs;
  }
  if ( value == "legacy" )
  {
    return majority_flow_kind::legacy;
  }
  return majority_flow_kind::standard;
}

static mmig_advanced_flow parse_mmig_advanced_flow( std::string const& value )
{
  if ( value == "dac19_default" )
  {
    return mmig_advanced_flow::dac19_default;
  }
  if ( value == "dac19_area" )
  {
    return mmig_advanced_flow::dac19_area;
  }
  if ( value == "compress2rs" )
  {
    return mmig_advanced_flow::compress2rs;
  }
  if ( value == "paper2019" || value == "dac19_compat" )
  {
    return mmig_advanced_flow::paper2019;
  }
  if ( value == "epfl" )
  {
    return mmig_advanced_flow::epfl;
  }
  if ( value == "legacy" )
  {
    return mmig_advanced_flow::legacy;
  }
  return mmig_advanced_flow::round_robin;
}

static char const* to_string( mmig_stage_mode mode )
{
  switch ( mode )
  {
  case mmig_stage_mode::pre:
    return "pre";
  case mmig_stage_mode::post:
    return "post";
  case mmig_stage_mode::both:
  default:
    return "both";
  }
}

static char const* to_string( mmig_seed_mode mode )
{
  switch ( mode )
  {
  case mmig_seed_mode::level1:
    return "l1";
  case mmig_seed_mode::level2:
    return "l2";
  case mmig_seed_mode::level3:
    return "l3";
  case mmig_seed_mode::both:
  default:
    return "both";
  }
}

static char const* to_string( majority_flow_kind flow )
{
  switch ( flow )
  {
  case majority_flow_kind::identity:
    return "identity";
  case majority_flow_kind::dac19_default:
    return "dac19_default";
  case majority_flow_kind::dac19_area:
    return "dac19_area";
  case majority_flow_kind::dac19_compat:
    return "dac19_compat";
  case majority_flow_kind::compress2rs:
    return "compress2rs";
  case majority_flow_kind::legacy:
    return "legacy";
  case majority_flow_kind::standard:
  default:
    return "standard";
  }
}

static char const* to_string( mmig_advanced_flow flow )
{
  switch ( flow )
  {
  case mmig_advanced_flow::dac19_default:
    return "dac19_default";
  case mmig_advanced_flow::dac19_area:
    return "dac19_area";
  case mmig_advanced_flow::compress2rs:
    return "compress2rs";
  case mmig_advanced_flow::paper2019:
    return "paper2019";
  case mmig_advanced_flow::epfl:
    return "epfl";
  case mmig_advanced_flow::legacy:
    return "legacy";
  case mmig_advanced_flow::round_robin:
  default:
    return "round_robin";
  }
}

static uint32_t count_minority_nodes( mig_network const& mig )
{
  uint32_t count = 0u;
  mig.foreach_gate( [&]( auto const& n ) {
    if ( mig.is_min( n ) )
    {
      ++count;
    }
  } );
  return count;
}

static void print_mmig_stage_summary( mmig_optimizer_stats const& st )
{
  auto const ms = [&]( auto const& d ) {
    return std::chrono::duration_cast<std::chrono::duration<double, std::milli>>( d ).count();
  };

  std::cout << "\n=== mMIG stage summary ===\n";
  std::cout << "CEC accepted=" << st.cec.accepted
            << " rejected=" << st.cec.rejected
            << " inconclusive=" << st.cec.inconclusive
            << " cec_time_ms=" << format_ms_value( ms( st.cec.time_cec ) ) << "\n";

  std::cout << "Seeding(pre): scored=" << st.seeding_pre.num_scored
            << " selected=" << st.seeding_pre.num_selected
            << " evaluated=" << st.seeding_pre.num_evaluated
            << " applied=" << st.seeding_pre.num_applied
            << " eq_checks=" << st.seeding_pre.num_equivalence_checks
            << " eq_rej=" << st.seeding_pre.num_equivalence_rejected
            << " low_l1_filtered=" << st.seeding_pre.num_filtered_low_l1
            << " depth_gain=" << st.seeding_pre.total_depth_gain
            << " area_gain=" << st.seeding_pre.total_area_gain
            << " inv_gain=" << st.seeding_pre.total_inv_gain
            << " time_ms=" << format_ms_value( ms( st.seeding_pre.time_total ) ) << "\n";

  std::cout << "Rewrite(pre): applied=" << st.rewrite_pre.num_applied
            << " rejected=" << st.rewrite_pre.num_rejected
            << " iterations=" << st.rewrite_pre.num_iterations
            << " time_ms=" << format_ms_value( ms( st.rewrite_pre.time_total ) ) << "\n";
  std::cout << "Rewrite(post): applied=" << st.rewrite_post.num_applied
            << " rejected=" << st.rewrite_post.num_rejected
            << " iterations=" << st.rewrite_post.num_iterations
            << " time_ms=" << format_ms_value( ms( st.rewrite_post.time_total ) ) << "\n";

  std::cout << "Inv-prop(pre): rewritten=" << st.inv_propagation_pre.num_rewritten
            << " est_gain=" << st.inv_propagation_pre.estimated_gain
            << " time_ms=" << format_ms_value( ms( st.inv_propagation_pre.time_total ) ) << "\n";
  std::cout << "Inv-prop(post): rewritten=" << st.inv_propagation_post.num_rewritten
            << " est_gain=" << st.inv_propagation_post.estimated_gain
            << " time_ms=" << format_ms_value( ms( st.inv_propagation_post.time_total ) ) << "\n";

  std::cout << "Inv-opt(pre): inverted=" << st.inv_optimization_pre.num_inverted
            << " gain=" << st.inv_optimization_pre.total_gain
            << " time_ms=" << format_ms_value( ms( st.inv_optimization_pre.time_total ) ) << "\n";
  std::cout << "Inv-opt(post): inverted=" << st.inv_optimization_post.num_inverted
            << " gain=" << st.inv_optimization_post.total_gain
            << " time_ms=" << format_ms_value( ms( st.inv_optimization_post.time_total ) ) << "\n";

  std::cout << "Adv-mighty: applied=" << st.advanced_mighty.num_applied
            << " rejected=" << st.advanced_mighty.num_rejected
            << " iterations=" << st.advanced_mighty.num_iterations
            << " time_ms=" << format_ms_value( ms( st.advanced_mighty.time_total ) ) << "\n";

  std::cout << "Adv-resub: checks=" << st.advanced_resub.num_candidate_checks
            << " accepted=" << st.advanced_resub.num_candidate_accepted
            << " rejected=" << st.advanced_resub.num_candidate_rejected
            << " inconclusive=" << st.advanced_resub.num_candidate_inconclusive
            << " skipped_non_mixed=" << st.advanced_resub.num_candidate_skipped_non_mixed
            << " skipped_budget=" << st.advanced_resub.num_candidate_skipped_budget
            << " est_gain=" << st.advanced_resub.resub.estimated_gain
            << " time_ms=" << format_ms_value( ms( st.advanced_resub.time_total ) ) << "\n";

  std::cout << "Adv-cut: depth_before=" << st.advanced_cut.depth_before
            << " depth_after=" << st.advanced_cut.depth_after
            << " rollbacks=" << st.advanced_cut.num_depth_rollbacks
            << " time_ms=" << format_ms_value( ms( st.advanced_cut.time_total ) ) << "\n";

  std::cout << "Adv-refac: depth_before=" << st.advanced_refactor.depth_before
            << " depth_after=" << st.advanced_refactor.depth_after
            << " rollbacks=" << st.advanced_refactor.num_depth_rollbacks
            << " time_ms=" << format_ms_value( ms( st.advanced_refactor.time_total ) ) << "\n";

  std::cout << "Adv-exact: candidates=" << st.advanced_exact.rewrite.candidates
            << " est_gain=" << st.advanced_exact.rewrite.estimated_gain
            << " depth_before=" << st.advanced_exact.depth_before
            << " depth_after=" << st.advanced_exact.depth_after
            << " rollbacks=" << st.advanced_exact.num_depth_rollbacks
            << " gated_skips=" << st.advanced_exact.num_gate_skips
            << " obs=(g=" << st.advanced_exact.observed_gates
            << ",min=" << st.advanced_exact.observed_minority
            << ",inv=" << st.advanced_exact.observed_inverted_edges
            << ")"
            << " time_ms=" << format_ms_value( ms( st.advanced_exact.time_total ) ) << "\n";

  std::cout << "Adv-balance: depth_before=" << st.advanced_balance.depth_before
            << " depth_after=" << st.advanced_balance.depth_after
            << " rollbacks=" << st.advanced_balance.num_depth_rollbacks
            << " time_ms=" << format_ms_value( ms( st.advanced_balance.time_total ) ) << "\n";

  std::cout << "Adv-rounds: executed=" << st.advanced_rounds_executed
            << " improved=" << st.advanced_rounds_improved
            << " rolled_back=" << st.advanced_rounds_rolled_back
            << " stagnation_stops=" << st.advanced_stagnation_stops
            << " obj_before=(d=" << st.advanced_depth_before
            << ",g=" << st.advanced_gates_before
            << ",maj=" << st.advanced_majority_before
            << ",min=" << st.advanced_minority_before
            << ",i=" << st.advanced_inverted_edges_before
            << ",ei=" << st.advanced_effective_inverted_edges_before
            << ",ea=" << st.advanced_effective_area_score_before
            << ") obj_after=(d=" << st.advanced_depth_after
            << ",g=" << st.advanced_gates_after
            << ",maj=" << st.advanced_majority_after
            << ",min=" << st.advanced_minority_after
            << ",i=" << st.advanced_inverted_edges_after
            << ",ei=" << st.advanced_effective_inverted_edges_after
            << ",ea=" << st.advanced_effective_area_score_after
            << ")\n";
}

// ---------- inverted edge counter ----------
uint64_t count_inverted_edges(const mig_network& mig) {
  uint64_t inv_edges = 0;
  
  // Count inverted edges on gate inputs
  mig.foreach_gate([&](auto n){
    mig.foreach_fanin(n, [&](auto const& fi){
      if (mig.is_complemented(fi)) ++inv_edges;
    });
  });
  
  // Count inverted edges on primary outputs
  mig.foreach_po([&](auto const& s, auto) {
    if (mig.is_complemented(s)) ++inv_edges;
  });
  
  return inv_edges;
}

struct optimization_result {
  mig_network network;
  uint64_t common_inv_edges = 0;
  uint64_t final_inv_edges = 0;
  double common_time_ms = 0.0;
  double post_time_ms = 0.0;
  double total_time_ms = 0.0;
};

enum class pipeline_mode
{
  depth,
  area
};

struct pipeline_params {
  bool no_cut = false;
  bool no_resub = false;
  int cut_k = 2;
  int cut_limit = 6;
  int rs_max_pis = 8;
  int rs_divs = 128;
  int rs_window_pre = 8;
  int rs_window_post = 10;
  bool run_inv_optimization = false;
  int lut_report_k = 0;
  bool write_lut_netlists = false;
  bool resub_use_dont_cares = false;
  int resub_dc_window = 12;
  pipeline_mode mode = pipeline_mode::depth;
  majority_flow_kind flow = majority_flow_kind::standard;
  bool flow_allow_zero_gain = false;
  bool flow_enable_balance = false;
  bool flow_enable_exact = false;
  int rs_max_inserts = 2;
};

static std::string stage_label(std::string_view prefix, std::string pass)
{
  if (prefix.empty())
  {
    return pass;
  }
  return "[" + std::string(prefix) + "] " + pass;
}

static void mig_algebraic_rewriting( mig_network& mig )
{
  depth_view depth_mig{ mig };
  mig_algebraic_depth_rewriting_params ps;
  ps.allow_area_increase = false;
  mig_algebraic_depth_rewriting( depth_mig, ps );
  mig = cleanup_dangling<mig_network>( depth_mig );
}

template<typename Fn>
static void run_pass_with_guard( mig_network& mig, pipeline_params const& params, const std::string& stage, Fn&& fn )
{
  if ( params.mode == pipeline_mode::area )
  {
    auto backup = mig;
    const auto before = backup.num_gates();
    {
      scoped_timer T( stage );
      fn();
    }
    const auto after = mig.num_gates();
    if ( after > before )
    {
      std::cout << "[info] " << stage << " reverted (gates " << before << " -> " << after << ")\n";
      mig = std::move( backup );
    }
  }
  else
  {
    scoped_timer T( stage );
    fn();
  }
}

static void run_resubstitution_stage( mig_network& mig, pipeline_params const& params, int max_pis, int max_divs, int max_inserts, int window_size, std::string_view prefix, std::string pass_label )
{
  std::string stage = stage_label( prefix, std::move( pass_label ) );
  run_pass_with_guard( mig, params, stage, [&](){
    fanout_view<mig_network> fview{ mig };
    depth_view<fanout_view<mig_network>> dfview{ fview };

    resubstitution_params rs;
    rs.max_pis = max_pis;
    rs.max_divisors = max_divs;
    rs.max_inserts = static_cast<uint32_t>( std::max( 0, max_inserts ) );
    rs.window_size = g_resub_use_dont_cares ? static_cast<uint32_t>( g_resub_dc_window ) : window_size;
    rs.preserve_depth = true;
    rs.use_dont_cares = g_resub_use_dont_cares;
    rs.verbose = false;

    mig_resubstitution( dfview, rs );
    mig = cleanup_dangling( mig );
  } );
}

static double to_milliseconds( std::chrono::steady_clock::time_point start, std::chrono::steady_clock::time_point end )
{
  return std::chrono::duration_cast<std::chrono::duration<double, std::milli>>( end - start ).count();
}

static std::string format_ms_value( double ms )
{
  std::ostringstream oss;
  oss << std::fixed << std::setprecision( 1 ) << ms;
  return oss.str();
}

static std::string format_seconds_value( double ms )
{
  std::ostringstream oss;
  oss << std::fixed << std::setprecision( 3 ) << ( ms / 1000.0 );
  return oss.str();
}

struct lut_artifact
{
  klut_network network;
  lut_metrics metrics;
};

static std::string shell_quote_string( std::string const& input )
{
  std::string quoted;
  quoted.reserve( input.size() + 2 );
  quoted.push_back( '"' );
  for ( char ch : input )
  {
    if ( ch == '\\' || ch == '"' )
    {
      quoted.push_back( '\\' );
    }
    quoted.push_back( ch );
  }
  quoted.push_back( '"' );
  return quoted;
}

static std::string shell_quote_path( std::filesystem::path const& path )
{
  return shell_quote_string( path.string() );
}

static std::filesystem::path make_temp_path( std::string_view suffix )
{
  static uint64_t counter = 0;
  const auto base = std::chrono::steady_clock::now().time_since_epoch().count();
  auto filename = std::string( "blif2mig_" ) + std::to_string( base ) + "_" + std::to_string( counter++ ) + std::string( suffix );
  return std::filesystem::temp_directory_path() / filename;
}

static bool is_executable_file( std::filesystem::path const& path )
{
  std::error_code ec;
  auto status = std::filesystem::status( path, ec );
  if ( ec || !std::filesystem::is_regular_file( status ) )
  {
    return false;
  }
#if defined( _WIN32 )
  (void)status;
  return true;
#else
  const auto perms = status.permissions();
  const auto exec_bits = std::filesystem::perms::owner_exec | std::filesystem::perms::group_exec | std::filesystem::perms::others_exec;
  return ( perms & exec_bits ) != std::filesystem::perms::none;
#endif
}

static std::optional<std::string> find_on_path( std::string const& name )
{
  if ( const char* path_env = std::getenv( "PATH" ); path_env && *path_env )
  {
    std::string paths = path_env;
    std::stringstream ss( paths );
    std::string entry;
    while ( std::getline( ss, entry, ':' ) )
    {
      if ( entry.empty() )
      {
        continue;
      }
      std::filesystem::path candidate = std::filesystem::path( entry ) / name;
      if ( is_executable_file( candidate ) )
      {
        return candidate.string();
      }
    }
  }
  return std::nullopt;
}

static std::string find_abc_executable()
{
  if ( const char* env = std::getenv( "ABC_BIN" ); env && *env )
  {
    const std::filesystem::path user_path = env;
    if ( is_executable_file( user_path ) )
    {
      return user_path.string();
    }
    return env;
  }

  if ( auto found = find_on_path( "abc" ) )
  {
    return *found;
  }

  auto dir = std::filesystem::current_path();
  for ( int depth = 0; depth < 8 && !dir.empty(); ++depth )
  {
    for ( auto const& rel : { std::filesystem::path( "abc" ), std::filesystem::path( "abc/abc" ), std::filesystem::path( "third_party/abc/abc" ) } )
    {
      const auto candidate = dir / rel;
      if ( is_executable_file( candidate ) )
      {
        return candidate.string();
      }
    }
    dir = dir.parent_path();
  }

  return "abc";
}

static klut_network mig_to_klut_subject( mig_network const& mig )
{
  klut_network klut;
  node_map<klut_network::signal, mig_network> node_to_signal( mig );

  const auto n_false = mig.get_node( mig.get_constant( false ) );
  node_to_signal[n_false] = klut.get_constant( false );
  const auto n_true = mig.get_node( mig.get_constant( true ) );
  if ( n_true != n_false )
  {
    node_to_signal[n_true] = klut.get_constant( true );
  }

  mig.foreach_pi( [&]( auto const& pi ) {
    node_to_signal[pi] = klut.create_pi();
  } );

  mig.foreach_node( [&]( auto const& n ) {
    if ( mig.is_constant( n ) || mig.is_pi( n ) )
    {
      return;
    }

    std::vector<klut_network::signal> fanins;
    fanins.reserve( mig.fanin_size( n ) );
    mig.foreach_fanin( n, [&]( auto const& f ) {
      auto sig = node_to_signal[mig.get_node( f )];
      if ( mig.is_complemented( f ) )
      {
        sig = klut.create_not( sig );
      }
      fanins.push_back( sig );
    } );

    assert( fanins.size() == 3u );
    auto maj_sig = klut.create_maj( fanins[0], fanins[1], fanins[2] );
    if ( mig.is_min( n ) )
    {
      maj_sig = klut.create_not( maj_sig );
    }
    node_to_signal[n] = maj_sig;
  } );

  mig.foreach_po( [&]( auto const& f ) {
    auto sig = node_to_signal[mig.get_node( f )];
    if ( mig.is_complemented( f ) )
    {
      sig = klut.create_not( sig );
    }
    klut.create_po( sig );
  } );

  return klut;
}

template<int K>
static std::optional<lut_artifact> map_to_lut_artifact_via_abc( mig_network const& mig )
{
  const auto abc_exec = find_abc_executable();
  const auto tmp_in = make_temp_path( "_abc_in.blif" );
  const auto tmp_out = make_temp_path( "_abc_out.blif" );
  const auto tmp_script = make_temp_path( "_abc_script.txt" );

  write_blif( mig, tmp_in.string() );

  std::ofstream script( tmp_script );
  if ( !script.is_open() )
  {
    std::cout << "[warn] Unable to create temporary ABC script at " << tmp_script << "\n";
    std::error_code ec;
    std::filesystem::remove( tmp_in, ec );
    return std::nullopt;
  }

  script << "read_blif \"" << tmp_in.string() << "\"\n";
  script << "strash\n";
  script << "if -K " << K << " -a\n";
  script << "write_blif \"" << tmp_out.string() << "\"\n";
  script << "quit\n";
  script.close();

  const std::string command = shell_quote_string( abc_exec ) + " -s -f " + shell_quote_path( tmp_script );
  const int result = std::system( command.c_str() );

  auto cleanup = [&]() {
    std::error_code ec;
    std::filesystem::remove( tmp_in, ec );
    std::filesystem::remove( tmp_script, ec );
    std::filesystem::remove( tmp_out, ec );
  };

  if ( result != 0 )
  {
    cleanup();
    std::cout << "[warn] ABC LUT mapping command failed (code " << result << "); falling back to internal LUT mapper\n";
    return std::nullopt;
  }

  klut_network klut;
  if ( lorina::read_blif( tmp_out.string(), blif_reader( klut ) ) != lorina::return_code::success )
  {
    cleanup();
    std::cout << "[warn] Failed to parse ABC output BLIF '" << tmp_out.string() << "'; falling back to internal LUT mapper\n";
    return std::nullopt;
  }

  cleanup();

  lut_artifact artifact{ std::move( klut ), {} };
  depth_view depth_v{ artifact.network };
  lut_metrics metrics;
  metrics.lut_count = artifact.network.num_gates();
  metrics.depth = depth_v.depth();
  metrics.pis = artifact.network.num_pis();
  metrics.fanin_histogram.assign( K + 1, 0 );
  artifact.network.foreach_gate( [&]( auto const& n ) {
    const auto fi = artifact.network.fanin_size( n );
    if ( fi < metrics.fanin_histogram.size() )
    {
      metrics.fanin_histogram[fi]++;
    }
  } );
  artifact.metrics = metrics;
  return artifact;
}

template<int K>
static std::optional<lut_artifact> map_to_lut_artifact( mig_network const& mig )
{
  if ( !g_use_exact_lut_mapping )
  {
    if ( auto abc_artifact = map_to_lut_artifact_via_abc<K>( mig ) )
    {
      return abc_artifact;
    }
  }

  auto subject = mig_to_klut_subject( mig );
  mapping_view<klut_network, true> mapped{ subject };
  lut_mapping_params ps;
  ps.cut_enumeration_ps.cut_size = K;
  if ( g_use_exact_lut_mapping )
  {
    satlut_mapping_params sps;
    sps.cut_enumeration_ps.cut_size = K;
    if ( g_exact_lut_window > 0 )
    {
      satlut_mapping<decltype( mapped ), true>( mapped, static_cast<uint32_t>( g_exact_lut_window ), sps );
    }
    else
    {
      satlut_mapping<decltype( mapped ), true>( mapped, sps );
    }
  }
  else
  {
    lut_mapping<decltype( mapped ), true>( mapped, ps );
  }

  auto klut_opt = collapse_mapped_network<klut_network>( mapped );
  if ( !klut_opt.has_value() )
  {
    return std::nullopt;
  }

  lut_artifact artifact{ std::move( *klut_opt ), {} };
  depth_view dv{ artifact.network };

  lut_metrics metrics;
  metrics.lut_count = artifact.network.num_gates();
  metrics.depth = dv.depth();
  metrics.pis = artifact.network.num_pis();
  metrics.fanin_histogram.assign( K + 1, 0 );

  artifact.network.foreach_gate( [&]( auto n ) {
    const auto fi = artifact.network.fanin_size( n );
    if ( fi <= K )
    {
      metrics.fanin_histogram[fi]++;
    }
  } );

  artifact.metrics = std::move( metrics );
  return artifact;
}

static std::string lut_sanitize_identifier( std::string name )
{
  std::string result;
  result.reserve( name.size() );
  for ( char ch : name )
  {
    if ( std::isalnum( static_cast<unsigned char>( ch ) ) || ch == '_' )
    {
      result.push_back( ch );
    }
    else
    {
      result.push_back( '_' );
    }
  }
  if ( result.empty() )
  {
    result = "lut_net";
  }
  if ( std::isdigit( static_cast<unsigned char>( result.front() ) ) )
  {
    result.insert( result.begin(), '_' );
  }
  return result;
}

static std::string lut_module_name_from_path( std::string path )
{
  const auto sep_pos = path.find_last_of( "/\\" );
  std::string base = ( sep_pos == std::string::npos ) ? path : path.substr( sep_pos + 1 );
  const auto dot_pos = base.find_last_of( '.' );
  if ( dot_pos != std::string::npos )
  {
    base = base.substr( 0, dot_pos );
  }
  if ( base.empty() )
  {
    base = "lut_net";
  }
  return lut_sanitize_identifier( base );
}

static std::string lut_apply_complement_expr( std::string expr, bool complemented )
{
  if ( !complemented )
  {
    return expr;
  }
  const bool needs_paren = expr.find_first_of( " []()" ) != std::string::npos;
  if ( needs_paren )
  {
    return "~(" + expr + ")";
  }
  return "~" + expr;
}

static void write_lut_verilog_simple( klut_network const& klut, std::string const& path )
{
  std::ofstream ofs( path );
  if ( !ofs.is_open() )
  {
    std::cerr << "[warn] cannot open " << path << " for LUT Verilog output\n";
    return;
  }

  const auto module_name = lut_module_name_from_path( path );
  const uint32_t num_pis = klut.num_pis();
  const uint32_t num_pos = klut.num_pos();

  ofs << "`timescale 1ns/1ps\n";
  ofs << "module " << module_name << "(";
  bool first_port = true;
  if ( num_pis > 0 )
  {
    ofs << "in";
    first_port = false;
  }
  if ( num_pos > 0 )
  {
    if ( !first_port )
    {
      ofs << ", ";
    }
    ofs << "out";
  }
  ofs << ");\n";

  if ( num_pis > 0 )
  {
    ofs << "  input [" << ( num_pis - 1 ) << ":0] in;\n";
  }
  else
  {
    ofs << "  // No primary inputs\n";
  }

  if ( num_pos > 0 )
  {
    ofs << "  output [" << ( num_pos - 1 ) << ":0] out;\n";
  }
  else
  {
    ofs << "  // No primary outputs\n";
  }

  node_map<std::string, klut_network> node_names( klut );
  node_names[klut.get_constant( false )] = "1'b0";
  if ( klut.get_node( klut.get_constant( false ) ) != klut.get_node( klut.get_constant( true ) ) )
  {
    node_names[klut.get_constant( true )] = "1'b1";
  }

  klut.foreach_pi( [&]( auto const& pi, uint32_t index ) {
    node_names[pi] = num_pis > 0 ? "in[" + std::to_string( index ) + "]" : "1'b0";
  } );

  std::vector<std::string> wire_names;
  klut.foreach_gate( [&]( auto const& n ) {
    std::string wire_name = "n" + std::to_string( klut.node_to_index( n ) );
    wire_names.push_back( wire_name );
    node_names[n] = wire_name;
  } );

  if ( !wire_names.empty() )
  {
    ofs << "  wire ";
    for ( size_t i = 0; i < wire_names.size(); ++i )
    {
      if ( i != 0 )
      {
        ofs << ", ";
      }
      ofs << wire_names[i];
    }
    ofs << ";\n";
  }

  klut.foreach_gate( [&]( auto const& n ) {
    std::vector<std::string> fanins;
    fanins.reserve( klut.fanin_size( n ) );
    klut.foreach_fanin( n, [&]( auto const& f ) {
      auto expr = node_names[klut.get_node( f )];
      fanins.push_back( lut_apply_complement_expr( expr, klut.is_complemented( f ) ) );
    } );

    const auto tt = klut.node_function( n );
    const uint32_t fanin_count = static_cast<uint32_t>( fanins.size() );
    if ( fanin_count == 0 )
    {
      const bool value = kitty::get_bit( tt, 0 );
      ofs << "  assign " << node_names[n] << " = 1'b" << ( value ? '1' : '0' ) << ";\n";
      return;
    }

    const uint32_t table_size = 1u << fanin_count;
    auto hex = kitty::to_hex( tt );
    const auto expected_digits = ( table_size + 3u ) / 4u;
    if ( hex.size() < expected_digits )
    {
      hex = std::string( expected_digits - hex.size(), '0' ) + hex;
    }

    const std::string param_name = "LUT_INIT_" + std::to_string( klut.node_to_index( n ) );
    ofs << "  localparam [" << ( table_size - 1 ) << ":0] " << param_name << " = " << table_size << "'h" << hex << ";\n";
    ofs << "  assign " << node_names[n] << " = " << param_name << "[";
    if ( fanin_count == 1 )
    {
      ofs << fanins.front();
    }
    else
    {
      ofs << "{";
      for ( auto it = fanins.rbegin(); it != fanins.rend(); ++it )
      {
        if ( it != fanins.rbegin() )
        {
          ofs << ", ";
        }
        ofs << *it;
      }
      ofs << "}";
    }
    ofs << "];\n";
  } );

  if ( num_pos > 0 )
  {
    klut.foreach_po( [&]( auto const& f, uint32_t index ) {
      const auto expr = lut_apply_complement_expr( node_names[klut.get_node( f )], klut.is_complemented( f ) );
      ofs << "  assign out[" << index << "] = " << expr << ";\n";
    } );
  }

  ofs << "endmodule\n";
}

template<int K>
static std::optional<lut_metrics> compute_k_lut_metrics( mig_network const& mig )
{
  if ( auto artifact = map_to_lut_artifact<K>( mig ) )
  {
    return artifact->metrics;
  }
  return std::nullopt;
}

template<int K>
static bool print_k_lut_stats( const mig_network& mig, const std::string& label, const std::string* output_base = nullptr )
{
  auto artifact_opt = map_to_lut_artifact<K>( mig );
  if ( !artifact_opt.has_value() )
  {
    std::cout << "[warn] " << label << ": failed to derive " << K << "-LUT mapping\n";
    return false;
  }

  const auto& artifact = *artifact_opt;
  const auto& metrics = artifact.metrics;

  std::cout << "\n=== " << label << " : " << K << "-LUT mapping ===\n";
  std::cout << "PIs=" << metrics.pis << "  POs=" << mig.num_pos()
            << "  LUTs=" << metrics.lut_count << "  LUT-Depth=" << metrics.depth << "\n";
  std::cout << "LUT fanin histogram:\n";
  for ( int i = 0; i <= K; ++i )
  {
    if ( metrics.fanin_histogram[i] )
    {
      std::cout << "  " << i << "-LUTs: " << metrics.fanin_histogram[i] << "\n";
    }
  }

  if ( output_base )
  {
    const std::string file_base = *output_base + "_lut" + std::to_string( K );
    write_lut_verilog_simple( artifact.network, file_base + ".v" );
    write_blif( artifact.network, file_base + ".blif" );
    std::cout << "[info] Wrote " << file_base << ".{v,blif}\n";
  }
  return true;
}

[[maybe_unused]] static void print_k_lut_stats_runtime( const mig_network& mig, int K, const std::string& label, const std::string* output_base = nullptr )
{
  switch ( K )
  {
    case 2: print_k_lut_stats<2>( mig, label, output_base ); break;
    case 3: print_k_lut_stats<3>( mig, label, output_base ); break;
    case 4: print_k_lut_stats<4>( mig, label, output_base ); break;
    case 5: print_k_lut_stats<5>( mig, label, output_base ); break;
    case 6: print_k_lut_stats<6>( mig, label, output_base ); break;
    case 7: print_k_lut_stats<7>( mig, label, output_base ); break;
    case 8: print_k_lut_stats<8>( mig, label, output_base ); break;
    default:
      std::cout << "[warn] Unsupported --lut-k value " << K << "; no LUT stats computed\n";
      break;
  }
}

enum class majority_flow_stage : uint8_t
{
  mighty_area,
  balancing,
  resubstitution,
  rewriting,
  rewriting_zero,
  refactoring,
  refactoring_zero,
  exact_rewriting
};

struct majority_flow_step
{
  majority_flow_stage stage;
  std::optional<int> rs_cut_override;
  std::optional<int> rs_insert_override;
};

static std::vector<majority_flow_step> build_majority_flow( pipeline_params const& params )
{
  std::vector<majority_flow_step> steps;
  auto add_stage = [&]( majority_flow_stage stage ) {
    steps.push_back( majority_flow_step{ stage, std::nullopt, std::nullopt } );
  };
  auto add_resub = [&]( int cut, int inserts ) {
    steps.push_back( majority_flow_step{ majority_flow_stage::resubstitution, cut, inserts } );
  };
  auto add_resub_with_depth_cap = [&]( int cut, int inserts ) {
    const int capped_inserts = ( cut > 8 && inserts > 1 ) ? 1 : inserts;
    add_resub( cut, capped_inserts );
  };

  switch ( params.flow )
  {
  case majority_flow_kind::compress2rs:
    add_stage( majority_flow_stage::mighty_area );
    add_resub_with_depth_cap( 6, 1 );
    add_stage( majority_flow_stage::rewriting );
    add_resub_with_depth_cap( 6, 2 );
    add_stage( majority_flow_stage::refactoring );
    add_resub_with_depth_cap( 8, 1 );
    add_stage( majority_flow_stage::mighty_area );
    add_resub_with_depth_cap( 8, 2 );
    add_stage( majority_flow_stage::rewriting );
    add_resub_with_depth_cap( 10, 1 );
    add_stage( majority_flow_stage::rewriting_zero );
    add_resub_with_depth_cap( 10, 2 );
    add_stage( majority_flow_stage::mighty_area );
    add_resub_with_depth_cap( 12, 1 );
    add_stage( majority_flow_stage::refactoring_zero );
    add_resub_with_depth_cap( 12, 2 );
    add_stage( majority_flow_stage::rewriting_zero );
    add_stage( majority_flow_stage::mighty_area );
    if ( params.flow_enable_exact )
    {
      add_stage( majority_flow_stage::exact_rewriting );
    }
    if ( params.flow_enable_balance )
    {
      add_stage( majority_flow_stage::balancing );
    }
    break;
  case majority_flow_kind::dac19_compat:
    add_stage( majority_flow_stage::mighty_area );
    add_resub_with_depth_cap( 6, 1 );
    add_stage( majority_flow_stage::rewriting );
    add_resub_with_depth_cap( 6, 2 );
    add_stage( majority_flow_stage::refactoring );
    add_resub_with_depth_cap( 8, 1 );
    add_stage( majority_flow_stage::mighty_area );
    add_resub_with_depth_cap( 8, 2 );
    add_stage( majority_flow_stage::rewriting );
    add_resub_with_depth_cap( 10, 1 );
    add_stage( majority_flow_stage::rewriting_zero );
    add_resub_with_depth_cap( 10, 2 );
    add_stage( majority_flow_stage::mighty_area );
    add_resub_with_depth_cap( 12, 1 );
    add_stage( majority_flow_stage::refactoring_zero );
    add_resub_with_depth_cap( 12, 2 );
    add_stage( majority_flow_stage::rewriting_zero );
    add_stage( majority_flow_stage::mighty_area );
    break;
  case majority_flow_kind::dac19_default:
  case majority_flow_kind::dac19_area:
    add_resub( params.rs_max_pis, params.rs_max_inserts );
    add_stage( majority_flow_stage::rewriting );
    add_resub( params.rs_max_pis, params.rs_max_inserts );
    add_stage( majority_flow_stage::refactoring );
    add_stage( majority_flow_stage::rewriting );
    if ( params.flow_allow_zero_gain )
    {
      add_stage( majority_flow_stage::rewriting_zero );
    }
    add_resub( params.rs_max_pis, params.rs_max_inserts );
    add_stage( params.flow_allow_zero_gain ? majority_flow_stage::refactoring_zero : majority_flow_stage::refactoring );
    add_stage( majority_flow_stage::rewriting );
    if ( params.flow_enable_exact )
    {
      add_stage( majority_flow_stage::exact_rewriting );
    }
    if ( params.flow_enable_balance )
    {
      add_stage( majority_flow_stage::balancing );
    }
    break;
  case majority_flow_kind::legacy:
    if ( params.flow_enable_balance )
    {
      add_stage( majority_flow_stage::balancing );
    }
    add_resub( params.rs_max_pis, params.rs_max_inserts );
    add_stage( majority_flow_stage::rewriting );
    add_resub( params.rs_max_pis, params.rs_max_inserts );
    add_stage( majority_flow_stage::refactoring );
    add_stage( majority_flow_stage::rewriting );
    if ( params.flow_allow_zero_gain )
    {
      add_stage( majority_flow_stage::rewriting_zero );
    }
    if ( params.flow_enable_balance )
    {
      add_stage( majority_flow_stage::balancing );
    }
    add_resub( params.rs_max_pis, params.rs_max_inserts );
    add_stage( majority_flow_stage::refactoring );
    if ( params.flow_allow_zero_gain )
    {
      add_stage( majority_flow_stage::refactoring_zero );
    }
    add_stage( majority_flow_stage::rewriting );
    if ( params.flow_enable_exact )
    {
      add_stage( majority_flow_stage::exact_rewriting );
    }
    if ( params.flow_enable_balance )
    {
      add_stage( majority_flow_stage::balancing );
    }
    break;
  case majority_flow_kind::standard:
  default:
    break;
  }

  return steps;
}

static void run_majority_mighty_area_stage( mig_network& mig, pipeline_params const& params, std::string_view prefix )
{
  std::string stage = stage_label( prefix, "mighty(area-aware)" );
  run_pass_with_guard( mig, params, stage, [&]() {
    depth_view depth_mig{ mig };
    mig_algebraic_depth_rewriting_params ps;
    ps.allow_area_increase = false;
    mig_algebraic_depth_rewriting( depth_mig, ps );
    mig = cleanup_dangling( mig );
  } );
}

static void run_majority_rewriting_stage( mig_network& mig,
                                          pipeline_params const& params,
                                          mig_npn_resynthesis& npn_resyn,
                                          bool allow_zero_gain,
                                          std::string_view prefix )
{
  std::string stage = stage_label( prefix, allow_zero_gain ? "rewriting(zero-gain)" : "rewriting" );
  run_pass_with_guard( mig, params, stage, [&]() {
    cut_rewriting_params crp;
    crp.cut_enumeration_ps.cut_size = static_cast<uint32_t>( std::max( 2, params.cut_k ) );
    crp.cut_enumeration_ps.cut_limit = static_cast<uint32_t>( std::max( 1, params.cut_limit ) );
    crp.cut_enumeration_ps.minimize_truth_table = true;
    crp.allow_zero_gain = allow_zero_gain;
    crp.preserve_depth = false;
    cut_rewriting( mig, npn_resyn, crp );
    mig = cleanup_dangling( mig );
  } );
}

static void run_majority_refactoring_stage( mig_network& mig,
                                            pipeline_params const& params,
                                            mig_npn_resynthesis& npn_resyn,
                                            bool allow_zero_gain,
                                            std::string_view prefix )
{
  std::string stage = stage_label( prefix, allow_zero_gain ? "refactoring(zero-gain)" : "refactoring" );
  run_pass_with_guard( mig, params, stage, [&]() {
    fanout_view<mig_network> fview{ mig };
    depth_view<fanout_view<mig_network>> dfview{ fview };
    refactoring_params rps;
    rps.allow_zero_gain = allow_zero_gain;
    rps.use_dont_cares = false;
    rps.verbose = false;
    refactoring( dfview, npn_resyn, rps );
    mig = cleanup_dangling( mig );
  } );
}

static void run_majority_balancing_stage( mig_network& mig, pipeline_params const& params, std::string_view prefix )
{
  std::string stage = stage_label( prefix, "balancing" );
  run_pass_with_guard( mig, params, stage, [&]() {
    sop_rebalancing<mig_network> balance_fn;
    balancing_params bps;
    bps.progress = false;
    auto balanced = balancing( mig, { balance_fn }, bps );
    mig = cleanup_dangling( balanced );
  } );
}

static void run_majority_exact_rewriting_stage( mig_network& mig,
                                                pipeline_params const& params,
                                                exact_library<mig_network>& exact_lib,
                                                bool allow_zero_gain,
                                                std::string_view prefix )
{
  std::string stage = stage_label( prefix, allow_zero_gain ? "exact_rewriting(zero-gain)" : "exact_rewriting" );
  run_pass_with_guard( mig, params, stage, [&]() {
    rewrite_params rps;
    rps.cut_enumeration_ps.cut_size = 4u;
    rps.cut_enumeration_ps.cut_limit = static_cast<uint32_t>( std::max( 8, params.cut_limit ) );
    rps.cut_enumeration_ps.minimize_truth_table = true;
    rps.preserve_depth = ( params.mode == pipeline_mode::depth );
    rps.allow_zero_gain = allow_zero_gain;
    rps.use_dont_cares = false;
    rewrite( mig, exact_lib, rps );
    mig = cleanup_dangling( mig );
  } );
}

static bool run_majority_flow_dac19_compat( mig_network& mig, pipeline_params const& params, std::string_view prefix, mig_npn_resynthesis& npn_resyn )
{
  auto const steps = build_majority_flow( params );
  for ( auto const& step : steps )
  {
    switch ( step.stage )
    {
    case majority_flow_stage::mighty_area:
    {
      std::string stage = stage_label( prefix, "mighty(area-aware)" );
      scoped_timer T( stage );
      depth_view depth_mig{ mig };
      mig_algebraic_depth_rewriting_params ps;
      ps.allow_area_increase = false;
      mig_algebraic_depth_rewriting( depth_mig, ps );
      mig = cleanup_dangling( mig );
      break;
    }
    case majority_flow_stage::resubstitution:
      if ( params.no_resub )
      {
        break;
      }
      {
        std::string stage = stage_label( prefix, "resubstitution" );
        scoped_timer T( stage );
        fanout_view<mig_network> fview{ mig };
        depth_view depth_fview{ fview };
        resubstitution_params rs;
        rs.max_pis = static_cast<uint32_t>( std::max( 1, step.rs_cut_override.value_or( params.rs_max_pis ) ) );
        rs.max_inserts = static_cast<uint32_t>( std::max( 0, step.rs_insert_override.value_or( params.rs_max_inserts ) ) );
        rs.preserve_depth = true;
        rs.use_dont_cares = false;
        rs.progress = false;
        mig_resubstitution( depth_fview, rs );
        mig = cleanup_dangling( mig );
      }
      break;
    case majority_flow_stage::rewriting:
    case majority_flow_stage::rewriting_zero:
      if ( params.no_cut )
      {
        break;
      }
      {
        std::string stage = stage_label( prefix, step.stage == majority_flow_stage::rewriting_zero ? "rewriting(zero-gain)" : "rewriting" );
        scoped_timer T( stage );
        cut_rewriting_params crp;
        crp.cut_enumeration_ps.cut_size = static_cast<uint32_t>( std::max( 2, params.cut_k ) );
        crp.cut_enumeration_ps.cut_limit = static_cast<uint32_t>( std::max( 1, params.cut_limit ) );
        crp.cut_enumeration_ps.minimize_truth_table = true;
        crp.min_cand_cut_size = 3;
        crp.allow_zero_gain = ( step.stage == majority_flow_stage::rewriting_zero );
        crp.preserve_depth = false;
        crp.progress = false;
        crp.verbose = false;

        fanout_view_params fvps;
        fvps.update_on_delete = false;
        fanout_view<mig_network> f{ mig, fvps };
        cut_rewriting_with_compatibility_graph( f, npn_resyn, crp );
        mig = cleanup_dangling( mig );
      }
      break;
    case majority_flow_stage::refactoring:
    case majority_flow_stage::refactoring_zero:
    {
      std::string stage = stage_label( prefix, step.stage == majority_flow_stage::refactoring_zero ? "refactoring(zero-gain)" : "refactoring" );
      scoped_timer T( stage );
      fanout_view<mig_network> f{ mig };
      depth_view depth_mig{ f };
      refactoring_params rps;
      rps.allow_zero_gain = ( step.stage == majority_flow_stage::refactoring_zero );
      rps.progress = false;
      refactoring( depth_mig, npn_resyn, rps );
      mig = cleanup_dangling( mig );
      break;
    }
    case majority_flow_stage::balancing:
    case majority_flow_stage::exact_rewriting:
      break;
    }
  }

  return true;
}

static bool run_majority_flow_preset( mig_network& mig, pipeline_params const& params, std::string_view prefix, mig_npn_resynthesis& npn_resyn )
{
  if ( params.flow == majority_flow_kind::standard )
  {
    return false;
  }

  if ( params.flow == majority_flow_kind::dac19_compat )
  {
    return run_majority_flow_dac19_compat( mig, params, prefix, npn_resyn );
  }

  auto const steps = build_majority_flow( params );
  if ( steps.empty() )
  {
    return true;
  }

  std::optional<exact_library<mig_network>> exact_lib{};
  if ( params.flow_enable_exact )
  {
    exact_lib.emplace( npn_resyn );
  }

  for ( auto const& step : steps )
  {
    switch ( step.stage )
    {
    case majority_flow_stage::mighty_area:
      run_majority_mighty_area_stage( mig, params, prefix );
      break;
    case majority_flow_stage::balancing:
      run_majority_balancing_stage( mig, params, prefix );
      break;
    case majority_flow_stage::resubstitution:
      if ( params.no_resub )
      {
        break;
      }
      run_resubstitution_stage(
          mig, params,
          step.rs_cut_override.value_or( params.rs_max_pis ),
          params.rs_divs,
          step.rs_insert_override.value_or( params.rs_max_inserts ),
          params.rs_window_pre,
          prefix,
          "resubstitution" );
      break;
    case majority_flow_stage::rewriting:
      if ( params.no_cut )
      {
        break;
      }
      run_majority_rewriting_stage( mig, params, npn_resyn, false, prefix );
      break;
    case majority_flow_stage::rewriting_zero:
      if ( params.no_cut )
      {
        break;
      }
      run_majority_rewriting_stage( mig, params, npn_resyn, true, prefix );
      break;
    case majority_flow_stage::refactoring:
      run_majority_refactoring_stage( mig, params, npn_resyn, false, prefix );
      break;
    case majority_flow_stage::refactoring_zero:
      run_majority_refactoring_stage( mig, params, npn_resyn, true, prefix );
      break;
    case majority_flow_stage::exact_rewriting:
      if ( params.no_cut )
      {
        break;
      }
      if ( exact_lib.has_value() )
      {
        run_majority_exact_rewriting_stage( mig, params, *exact_lib, params.flow_allow_zero_gain, prefix );
      }
      break;
    }
  }

  return true;
}

static void apply_common_passes( mig_network& mig, pipeline_params const& params, std::string_view prefix, mig_npn_resynthesis& npn_resyn )
{
  if ( params.mode == pipeline_mode::depth )
  {
    {
      std::string stage = stage_label( prefix, "mig_algebraic_depth_rewriting (pre)" );
      scoped_timer T( stage );
      depth_view depth_mig{ mig };
      mig_algebraic_depth_rewriting_params ps;
      ps.allow_area_increase = true;
      mig_algebraic_depth_rewriting( depth_mig, ps );
      mig = cleanup_dangling<mig_network>( depth_mig );
    }

    if ( !params.no_resub )
    {
      run_resubstitution_stage(
          mig, params, params.rs_max_pis, params.rs_divs, params.rs_max_inserts, params.rs_window_pre, prefix,
          "resubstitution (pre)" );
    }

    return;
  }

  // Area-oriented flow
  if ( !params.no_cut )
  {
    std::string stage = stage_label( prefix, "cut_rewriting" );
    run_pass_with_guard( mig, params, stage, [&](){
      cut_rewriting_params crp;
      crp.cut_enumeration_ps.cut_size = params.cut_k;
      crp.cut_enumeration_ps.cut_limit = params.cut_limit;
      cut_rewriting( mig, npn_resyn, crp );
      mig = cleanup_dangling( mig );
    } );
  }

  if ( !params.no_resub )
  {
    run_resubstitution_stage(
        mig, params, params.rs_max_pis, params.rs_divs, params.rs_max_inserts, params.rs_window_pre, prefix,
        "resubstitution (pre)" );
  }

  {
    std::string stage = stage_label( prefix, "mig_algebraic_rewriting" );
    run_pass_with_guard( mig, params, stage, [&](){
      mig_algebraic_rewriting( mig );
      mig = cleanup_dangling( mig );
    } );
  }
}

static void apply_post_cleanup( mig_network& mig, pipeline_params const& params, std::string_view prefix )
{
  if ( params.mode == pipeline_mode::depth )
  {
    {
      std::string stage = stage_label( prefix, "mig_algebraic_depth_rewriting (post)" );
      scoped_timer T( stage );
      depth_view depth_mig{ mig };
      mig_algebraic_depth_rewriting_params ps;
      ps.allow_area_increase = false;
      mig_algebraic_depth_rewriting( depth_mig, ps );
      mig = cleanup_dangling<mig_network>( depth_mig );
    }

    if ( params.run_inv_optimization )
    {
      std::string stage = stage_label( prefix, "final mig_inv_optimization" );
      scoped_timer T( stage );
      mig_inv_optimization( mig );
      mig = cleanup_dangling( mig );
    }

    return;
  }

  // Area-oriented cleanup
  if ( !params.no_resub )
  {
    run_resubstitution_stage(
        mig, params, params.rs_max_pis, params.rs_divs, params.rs_max_inserts, params.rs_window_post, prefix,
        "resubstitution (post)" );
  }

  if ( params.run_inv_optimization )
  {
    std::string stage = stage_label( prefix, "final mig_inv_optimization" );
    run_pass_with_guard( mig, params, stage, [&](){
      mig_inv_optimization( mig );
      mig = cleanup_dangling( mig );
    } );
  }
}

static optimization_result run_pipeline( mig_network mig, pipeline_params const& params, std::string_view pipeline_tag, std::string_view final_label, mig_npn_resynthesis& npn_resyn )
{
  optimization_result result;

  print_stats( mig, stage_label( pipeline_tag, "Original network" ) );

  if ( params.flow == majority_flow_kind::standard )
  {
    const auto common_begin = std::chrono::steady_clock::now();
    apply_common_passes( mig, params, pipeline_tag, npn_resyn );
    const auto common_end = std::chrono::steady_clock::now();
    result.common_time_ms = to_milliseconds( common_begin, common_end );

    result.common_inv_edges = count_inverted_edges( mig );
    print_stats( mig, stage_label( pipeline_tag, "After common majority passes" ) );

    const auto post_begin = std::chrono::steady_clock::now();
    apply_post_cleanup( mig, params, pipeline_tag );
    const auto post_end = std::chrono::steady_clock::now();
    result.post_time_ms = to_milliseconds( post_begin, post_end );
    result.total_time_ms = result.common_time_ms + result.post_time_ms;
  }
  else
  {
    const auto flow_begin = std::chrono::steady_clock::now();
    auto const flow_ok = run_majority_flow_preset( mig, params, pipeline_tag, npn_resyn );
    const auto flow_end = std::chrono::steady_clock::now();
    if ( !flow_ok )
    {
      std::cout << "[warn] selected majority flow did not run; falling back to standard pipeline\n";
      const auto common_begin = std::chrono::steady_clock::now();
      apply_common_passes( mig, params, pipeline_tag, npn_resyn );
      const auto common_end = std::chrono::steady_clock::now();
      result.common_time_ms = to_milliseconds( common_begin, common_end );
      const auto post_begin = std::chrono::steady_clock::now();
      apply_post_cleanup( mig, params, pipeline_tag );
      const auto post_end = std::chrono::steady_clock::now();
      result.post_time_ms = to_milliseconds( post_begin, post_end );
      result.total_time_ms = result.common_time_ms + result.post_time_ms;
    }
    else
    {
      result.common_time_ms = to_milliseconds( flow_begin, flow_end );
      result.post_time_ms = 0.0;
      result.total_time_ms = result.common_time_ms;
    }
    result.common_inv_edges = count_inverted_edges( mig );
    print_stats( mig, stage_label( pipeline_tag, "After selected majority flow" ) );
  }

  result.final_inv_edges = count_inverted_edges( mig );

  print_stats( mig, std::string( final_label ) );
  result.network = std::move( mig );

  return result;
}

[[maybe_unused]] static void print_simple_summary( const optimization_result& result, uint64_t original_inv_edges )
{
  std::cout << "\n" << std::string( 60, '=' ) << "\n";
  std::cout << "INVERTED EDGE SUMMARY\n";
  std::cout << std::string( 60, '=' ) << "\n";
  std::cout << "Original network:             " << original_inv_edges << " inverted edges\n";
  std::cout << "After selected flow:          " << result.common_inv_edges << " inverted edges\n";
  std::cout << "Final output:                 " << result.final_inv_edges << " inverted edges\n";
  std::cout << std::string( 60, '=' ) << "\n";

  std::cout << "\nTIMING (ms)\n";
  std::cout << "Common majority passes:       " << format_ms_value( result.common_time_ms ) << "\n";
  std::cout << "Post-cleanup passes:          " << format_ms_value( result.post_time_ms ) << "\n";
  std::cout << "Total optimisation time:      " << format_ms_value( result.total_time_ms ) << "\n";
}

// ---------- main ----------=444
#ifndef BLIF2MIG_BUILD_LIBRARY
int main(int argc, char** argv)
{
  if (argc < 3) {
    std::cerr << "Usage: blif2mig <in.blif> <out_base> [--fast|--no-cut|--no-resub|--no-refac"
             << "|--window=<n>|--window-pre=<n>|--window-post=<n>"
             << "|--mode=area|--mode=depth|--mig-flow=identity|standard|dac19_default|dac19_area|dac19_compat|compress2rs|legacy"
             << "|--mig-allow-zero-gain|--mig-enable-balance|--mig-enable-exact|--inv-opt|--no-inv-opt|--lut-k=<n>|--write-lut-netlists"
             << "|--rs-cut=<n>|--rs-insert=<n>|--rw-cut=<n>|--rw-cut-limit=<n>"
             << "|--resub-dc|--resub-dc-window=<n>|--exact-lut|--exact-lut-window=<n>"
             << "|--enable-mmig|--mmig-max-iters=<n>|--mmig-stage=pre|post|both|--mmig-cec|--no-mmig-cec|--mmig-enable-sr5"
             << "|--mmig-seed=off|l1|l2|both|--mmig-seed-budget=<n>|--mmig-seed-rounds=<n>"
             << "|--mmig-advanced|--mmig-flow=round_robin|dac19_default|dac19_area|compress2rs|paper2019|epfl|legacy"
             << "|--mmig-allow-zero-gain|--mmig-advanced-rounds=<n>|--mmig-adv-no-resub|--mmig-adv-no-cut|--mmig-adv-no-refac|--mmig-adv-no-exact|--mmig-adv-no-balance"
             << "|--mmig-adv-no-objective-guard|--mmig-adv-no-stagnation-stop|--mmig-adv-stagnation=<n>"
             << "|--mmig-dont-cares|--mmig-dc-threshold=<n>|--mmig-interleaved-seeding|--mmig-interleaved-seed-interval=<n>"
             << "|--no-mmig-dual-inv|--no-mmig-auto-sr5|--no-mmig-rewrite-ranking|--no-mmig-normalize-inner"
             << "|--no-mmig-resub-ranking|--mmig-adv-no-tuned-policy]" << "\n";
    return 1;
  }

  const bool fast = has_flag(argc, argv, "--fast");
  const bool no_cut = has_flag(argc, argv, "--no-cut");
  const bool no_resub = has_flag(argc, argv, "--no-resub");
  const bool mig_allow_zero_gain = has_flag( argc, argv, "--mig-allow-zero-gain" );
  const bool mig_enable_balance = has_flag( argc, argv, "--mig-enable-balance" );
  const bool mig_enable_exact = has_flag( argc, argv, "--mig-enable-exact" );
  const bool inv_opt = has_flag( argc, argv, "--inv-opt" );
  const bool no_inv_opt = has_flag( argc, argv, "--no-inv-opt" );
  const bool enable_mmig = has_flag( argc, argv, "--enable-mmig" );
  const bool mmig_cec_requested = has_flag( argc, argv, "--mmig-cec" );
  const bool no_mmig_cec = has_flag( argc, argv, "--no-mmig-cec" );
  const bool mmig_enable_sr5 = has_flag( argc, argv, "--mmig-enable-sr5" );
  const bool mmig_advanced = has_flag( argc, argv, "--mmig-advanced" );
  const bool mmig_allow_zero_gain = has_flag( argc, argv, "--mmig-allow-zero-gain" );
  const bool mmig_adv_no_resub = has_flag( argc, argv, "--mmig-adv-no-resub" );
  const bool mmig_adv_no_cut = has_flag( argc, argv, "--mmig-adv-no-cut" );
  const bool mmig_adv_no_refac = has_flag( argc, argv, "--mmig-adv-no-refac" );
  const bool mmig_adv_no_exact = has_flag( argc, argv, "--mmig-adv-no-exact" );
  const bool mmig_adv_no_balance = has_flag( argc, argv, "--mmig-adv-no-balance" );
  const bool mmig_adv_no_objective_guard = has_flag( argc, argv, "--mmig-adv-no-objective-guard" );
  const bool mmig_adv_no_stagnation_stop = has_flag( argc, argv, "--mmig-adv-no-stagnation-stop" );
  const bool mmig_use_dont_cares = has_flag( argc, argv, "--mmig-dont-cares" ); /* P2 */
  const bool mmig_interleaved_seeding = has_flag( argc, argv, "--mmig-interleaved-seeding" ); /* P3 */
  const bool mmig_no_dual_inv = has_flag( argc, argv, "--no-mmig-dual-inv" ); /* P1 */
  const bool mmig_no_auto_sr5 = has_flag( argc, argv, "--no-mmig-auto-sr5" ); /* P7 */
  const bool mmig_no_rewrite_ranking = has_flag( argc, argv, "--no-mmig-rewrite-ranking" );
  const bool mmig_no_normalize_inner = has_flag( argc, argv, "--no-mmig-normalize-inner" );
  const bool mmig_no_resub_ranking = has_flag( argc, argv, "--no-mmig-resub-ranking" );
  const bool mmig_adv_no_tuned_policy = has_flag( argc, argv, "--mmig-adv-no-tuned-policy" );
  majority_flow_kind mig_flow = majority_flow_kind::standard;
  mmig_advanced_flow mmig_flow = mmig_advanced_flow::round_robin;
  bool mmig_seed_enabled = true;
  mmig_seed_mode mmig_seed_mode_ps = mmig_seed_mode::both;
  int mmig_seed_budget = 24;
  int mmig_seed_rounds = 1;
  bool mmig_seed_incremental = false;
  int mmig_advanced_rounds = 1;
  int mmig_adv_stagnation = 1;
  int mmig_dc_threshold = 10000;
  int mmig_interleaved_seed_interval = 2;

  int mmig_max_iters = 16;
  if ( auto mmig_iters = parse_int_option( argc, argv, "--mmig-max-iters=" ); mmig_iters )
  {
    mmig_max_iters = std::clamp( *mmig_iters, 1, 200 );
  }

  mmig_stage_mode mmig_stage = mmig_stage_mode::both;
  if ( auto mmig_stage_opt = parse_string_option( argc, argv, "--mmig-stage=" ); mmig_stage_opt )
  {
    if ( *mmig_stage_opt == "pre" || *mmig_stage_opt == "post" || *mmig_stage_opt == "both" )
    {
      mmig_stage = parse_mmig_stage_mode( *mmig_stage_opt );
    }
    else
    {
      std::cout << "[warn] Unknown --mmig-stage value '" << *mmig_stage_opt << "', defaulting to both\n";
    }
  }

  if ( auto mig_flow_opt = parse_string_option( argc, argv, "--mig-flow=" ); mig_flow_opt )
  {
    if ( *mig_flow_opt == "identity" || *mig_flow_opt == "none" || *mig_flow_opt == "standard" || *mig_flow_opt == "dac19_default" || *mig_flow_opt == "dac19_area" || *mig_flow_opt == "dac19_compat" || *mig_flow_opt == "compress2rs" || *mig_flow_opt == "legacy" )
    {
      mig_flow = parse_majority_flow_kind( *mig_flow_opt );
    }
    else
    {
      std::cout << "[warn] Unknown --mig-flow value '" << *mig_flow_opt << "', defaulting to standard\n";
    }
  }

  if ( auto mmig_flow_opt = parse_string_option( argc, argv, "--mmig-flow=" ); mmig_flow_opt )
  {
    if ( *mmig_flow_opt == "round_robin" || *mmig_flow_opt == "standard" || *mmig_flow_opt == "dac19_default" || *mmig_flow_opt == "dac19_area" || *mmig_flow_opt == "compress2rs" || *mmig_flow_opt == "paper2019" || *mmig_flow_opt == "dac19_compat" || *mmig_flow_opt == "epfl" || *mmig_flow_opt == "legacy" )
    {
      mmig_flow = parse_mmig_advanced_flow( *mmig_flow_opt );
    }
    else
    {
      std::cout << "[warn] Unknown --mmig-flow value '" << *mmig_flow_opt << "', defaulting to round_robin\n";
    }
  }

  if ( auto mmig_seed_opt = parse_string_option( argc, argv, "--mmig-seed=" ); mmig_seed_opt )
  {
    if ( *mmig_seed_opt == "off" || *mmig_seed_opt == "none" )
    {
      mmig_seed_enabled = false;
    }
    else if ( *mmig_seed_opt == "l1" || *mmig_seed_opt == "l2" || *mmig_seed_opt == "l3" || *mmig_seed_opt == "both" )
    {
      mmig_seed_enabled = true;
      mmig_seed_mode_ps = parse_mmig_seed_mode( *mmig_seed_opt );
    }
    else
    {
      std::cout << "[warn] Unknown --mmig-seed value '" << *mmig_seed_opt << "', defaulting to both\n";
    }
  }

  if ( auto mmig_seed_budget_opt = parse_int_option( argc, argv, "--mmig-seed-budget=" ); mmig_seed_budget_opt )
  {
    mmig_seed_budget = std::clamp( *mmig_seed_budget_opt, 1, 10000 );
  }

  if ( auto mmig_seed_rounds_opt = parse_int_option( argc, argv, "--mmig-seed-rounds=" ); mmig_seed_rounds_opt )
  {
    mmig_seed_rounds = std::clamp( *mmig_seed_rounds_opt, 1, 16 );
  }

  if ( has_flag( argc, argv, "--mmig-seed-incremental" ) )
  {
    mmig_seed_incremental = true;
  }

  if ( auto mmig_dc_threshold_opt = parse_int_option( argc, argv, "--mmig-dc-threshold=" ); mmig_dc_threshold_opt )
  {
    mmig_dc_threshold = std::clamp( *mmig_dc_threshold_opt, 0, 100000000 );
  }

  if ( auto mmig_interleaved_seed_interval_opt = parse_int_option( argc, argv, "--mmig-interleaved-seed-interval=" ); mmig_interleaved_seed_interval_opt )
  {
    mmig_interleaved_seed_interval = std::clamp( *mmig_interleaved_seed_interval_opt, 1, 64 );
  }

  if ( auto mmig_advanced_rounds_opt = parse_int_option( argc, argv, "--mmig-advanced-rounds=" ); mmig_advanced_rounds_opt )
  {
    mmig_advanced_rounds = std::clamp( *mmig_advanced_rounds_opt, 1, 8 );
  }
  if ( auto mmig_adv_stagnation_opt = parse_int_option( argc, argv, "--mmig-adv-stagnation=" ); mmig_adv_stagnation_opt )
  {
    mmig_adv_stagnation = std::clamp( *mmig_adv_stagnation_opt, 0, 16 );
  }

  bool mmig_cec_enabled = true;
  if ( no_mmig_cec )
  {
    mmig_cec_enabled = false;
  }
  if ( mmig_cec_requested )
  {
    mmig_cec_enabled = true;
  }

  g_resub_use_dont_cares = has_flag( argc, argv, "--resub-dc" );
  if ( auto dc_win = parse_int_option( argc, argv, "--resub-dc-window=" ); dc_win )
  {
    g_resub_dc_window = std::clamp( *dc_win, 2, 32 );
  }
  if ( auto dc_win_alt = parse_int_option( argc, argv, "--resub-window=" ); dc_win_alt )
  {
    g_resub_dc_window = std::clamp( *dc_win_alt, 2, 32 );
  }
  g_use_exact_lut_mapping = has_flag( argc, argv, "--exact-lut" );
  if ( auto exact_win = parse_int_option( argc, argv, "--exact-lut-window=" ); exact_win )
  {
    g_exact_lut_window = std::max( 0, *exact_win );
  }

  int rs_max_pis = 8;
  int rs_max_inserts = 2;
  const int rs_divs = fast ? 64 : 128;
  int rs_window = fast ? 6 : 8;

  if ( auto override_window = parse_int_option( argc, argv, "--window=" ); override_window )
  {
    rs_window = std::clamp( *override_window, 2, 32 );
  }
  if ( auto rs_insert_opt = parse_int_option( argc, argv, "--rs-insert=" ); rs_insert_opt )
  {
    rs_max_inserts = std::clamp( *rs_insert_opt, 0, 4 );
  }
  if ( auto rs_cut_opt = parse_int_option( argc, argv, "--rs-cut=" ); rs_cut_opt )
  {
    rs_max_pis = std::clamp( *rs_cut_opt, 4, 16 );
  }

  pipeline_params params;
  params.no_cut = no_cut;
  params.no_resub = no_resub;
  params.cut_k = fast ? 2 : 4;
  params.cut_limit = fast ? 6 : 10;
  params.rs_max_pis = rs_max_pis;
  params.rs_max_inserts = rs_max_inserts;
  params.rs_divs = rs_divs;
  params.rs_window_pre = rs_window;
  params.run_inv_optimization = inv_opt && !no_inv_opt;
  params.resub_use_dont_cares = g_resub_use_dont_cares;
  params.resub_dc_window = g_resub_dc_window;
  params.flow = mig_flow;
  params.flow_allow_zero_gain = mig_allow_zero_gain;
  params.flow_enable_balance = mig_enable_balance;
  params.flow_enable_exact = mig_enable_exact;
  if ( params.flow == majority_flow_kind::dac19_compat )
  {
    params.cut_k = 4;
    params.cut_limit = 12;
  }
  params.lut_report_k = 0;
  if ( auto mode_opt = parse_string_option( argc, argv, "--mode=" ); mode_opt )
  {
    if ( *mode_opt == "area" )
    {
      params.mode = pipeline_mode::area;
    }
    else if ( *mode_opt == "depth" )
    {
      params.mode = pipeline_mode::depth;
    }
    else
    {
      std::cout << "[warn] Unknown --mode value '" << *mode_opt << "', defaulting to depth\n";
    }
  }

  if ( auto lut_k_override = parse_int_option( argc, argv, "--lut-k=" ); lut_k_override )
  {
    params.lut_report_k = std::max( 0, *lut_k_override );
  }
  if ( auto rw_cut_override = parse_int_option( argc, argv, "--rw-cut=" ); rw_cut_override )
  {
    params.cut_k = std::clamp( *rw_cut_override, 2, 8 );
  }
  if ( auto rw_limit_override = parse_int_option( argc, argv, "--rw-cut-limit=" ); rw_limit_override )
  {
    params.cut_limit = std::clamp( *rw_limit_override, 1, 64 );
  }
  if ( params.lut_report_k > 0 && ( params.lut_report_k < 2 || params.lut_report_k > 8 ) )
  {
    std::cout << "[warn] --lut-k value " << params.lut_report_k << " not supported (valid range 2..8); skipping LUT stats\n";
    params.lut_report_k = 0;
  }
  params.write_lut_netlists = has_flag( argc, argv, "--write-lut-netlists" );
  if ( params.write_lut_netlists && params.lut_report_k == 0 )
  {
    std::cout << "[warn] --write-lut-netlists requires --lut-k; ignoring LUT netlist request\n";
    params.write_lut_netlists = false;
  }

  // 1) Read BLIF -> KLUT
  klut_network klut;
  if (lorina::read_blif(argv[1], blif_reader(klut)) != lorina::return_code::success) {
    std::cerr << "Error: could not read BLIF '" << argv[1] << "'\n";
    return 2;
  }

  // 2) KLUT -> MIG
  const bool import_lut = has_flag( argc, argv, "--import-lut" );
  mig_npn_resynthesis npn_resyn;
  mig_network original_mig = import_lut ? convert_klut_to_graph<mig_network>( klut ) : node_resynthesis<mig_network>( klut, npn_resyn );
  print_stats( original_mig, "Original MIG" );
  if ( params.lut_report_k > 0 )
  {
    if ( params.lut_report_k >= 2 && params.lut_report_k <= 8 )
    {
      print_k_lut_stats_runtime( original_mig, params.lut_report_k, "Original MIG" );
    }
  }
  const auto original_inv_edges = count_inverted_edges( original_mig );

  // Adaptive window downsizing for large designs to keep runtime reasonable.
  const auto initial_gate_count = original_mig.num_gates();
  if ( initial_gate_count > 10000 )
  {
    params.rs_window_pre = std::min( params.rs_window_pre, 4 );
  }
  else if ( initial_gate_count > 4000 )
  {
    params.rs_window_pre = std::min( params.rs_window_pre, 5 );
  }
  else if ( initial_gate_count > 1500 )
  {
    params.rs_window_pre = std::min( params.rs_window_pre, 6 );
  }

  if ( auto pre_override = parse_int_option( argc, argv, "--window-pre=" ); pre_override )
  {
    params.rs_window_pre = std::clamp( *pre_override, 2, 32 );
  }

  params.rs_window_post = std::max( params.rs_window_pre, 4 ) + 2;
  if ( auto post_override = parse_int_option( argc, argv, "--window-post=" ); post_override )
  {
    params.rs_window_post = std::clamp( *post_override, 2, 40 );
  }

  std::cout << "[MIG] mode=" << ( params.mode == pipeline_mode::area ? "area" : "depth" )
            << " flow=" << to_string( params.flow )
            << " rs_cut=" << params.rs_max_pis
            << " rs_insert=" << params.rs_max_inserts
            << " rw_cut=" << params.cut_k
            << " rw_cut_limit=" << params.cut_limit
            << " zero_gain=" << ( params.flow_allow_zero_gain ? "on" : "off" )
            << " balance=" << ( params.flow_enable_balance ? "on" : "off" )
            << " exact=" << ( params.flow_enable_exact ? "on" : "off" ) << "\n";

  const auto majority_wall_begin = std::chrono::steady_clock::now();
  auto majority_result = run_pipeline( original_mig, params, "Majority", "Optimized MIG", npn_resyn );
  const double majority_wall_ms = to_milliseconds( majority_wall_begin, std::chrono::steady_clock::now() );

  print_simple_summary( majority_result, original_inv_edges );
  std::cout << "\nEnd-to-end mapping time (wall clock):\n";
  std::cout << "  Majority flow:             " << format_seconds_value( majority_wall_ms ) << " s\n";
  const std::string out_base = argv[2];
  if ( params.lut_report_k > 0 )
  {
    const std::string majority_lut_base = out_base + "_maj_opt";
    print_k_lut_stats_runtime(
        majority_result.network, params.lut_report_k, "Majority-only optimized MIG",
        params.write_lut_netlists ? &majority_lut_base : nullptr );
  }

  write_verilog( majority_result.network, out_base + "_maj_opt.v" );
  write_blif( majority_result.network, out_base + "_maj_opt.blif" );

  std::cout << "\nWrote " << out_base << "_maj_opt.{v,blif}\n";

  if ( enable_mmig )
  {
    mig_network mmig_ntk = majority_result.network;
    const auto min_before = count_minority_nodes( mmig_ntk );
    const bool mmig_use_advanced = mmig_advanced || mmig_flow != mmig_advanced_flow::round_robin;
    mmig_optimizer_params mmig_ps;
    mmig_ps.max_iterations = static_cast<uint32_t>( mmig_max_iters );
    mmig_ps.enable_minority_seeding = mmig_seed_enabled;
    mmig_ps.seeding_mode = mmig_seed_mode_ps;
    mmig_ps.seeding_max_candidates = static_cast<uint32_t>( mmig_seed_budget );
    mmig_ps.seeding_max_rounds = static_cast<uint32_t>( mmig_seed_rounds );
    mmig_ps.seeding_incremental_reranking = mmig_seed_incremental;
    mmig_ps.enable_sr5 = mmig_enable_sr5;
    mmig_ps.stage_mode = mmig_stage;
    mmig_ps.run_cec = mmig_cec_enabled;
    mmig_ps.enable_advanced = mmig_use_advanced;
    mmig_ps.advanced_rounds = static_cast<uint32_t>( mmig_advanced_rounds );
    mmig_ps.advanced_flow = mmig_flow;
    mmig_ps.advanced_allow_zero_gain = mmig_allow_zero_gain;
    mmig_ps.advanced_run_resub = !mmig_adv_no_resub;
    mmig_ps.advanced_run_cut = !mmig_adv_no_cut;
    mmig_ps.advanced_run_refactor = !mmig_adv_no_refac;
    mmig_ps.advanced_run_exact = !mmig_adv_no_exact;
    mmig_ps.advanced_run_balance = !mmig_adv_no_balance;
    mmig_ps.advanced_rollback_on_objective_regression = !mmig_adv_no_objective_guard;
    mmig_ps.advanced_stop_on_stagnation = !mmig_adv_no_stagnation_stop;
    mmig_ps.advanced_max_stagnation_rounds = static_cast<uint32_t>( mmig_adv_stagnation );
    mmig_ps.advanced_tuned_policy = !mmig_adv_no_tuned_policy;
    mmig_ps.enable_dont_cares = mmig_use_dont_cares;                   /* P2 */
    mmig_ps.dont_care_gate_threshold = static_cast<uint32_t>( mmig_dc_threshold );
    mmig_ps.enable_interleaved_seeding = mmig_interleaved_seeding;     /* P3 */
    mmig_ps.interleaved_seeding_interval = static_cast<uint32_t>( mmig_interleaved_seed_interval );
    mmig_ps.enable_dual_inv_propagation = !mmig_no_dual_inv;           /* P1 */
    mmig_ps.auto_sr5_with_cec = !mmig_no_auto_sr5;                    /* P7 */
    mmig_ps.rewrite_rank_candidates = !mmig_no_rewrite_ranking;
    mmig_ps.rewrite_normalize_complemented_inner = !mmig_no_normalize_inner;
    mmig_ps.resub_rank_candidates = !mmig_no_resub_ranking;
    mmig_ps.verbose = true;

    mmig_optimizer_stats mmig_st{};
    const auto mmig_wall_begin = std::chrono::steady_clock::now();
    if ( params.mode == pipeline_mode::area )
    {
      mmig_area_optimization( mmig_ntk, mmig_ps, &mmig_st );
    }
    else
    {
      mmig_depth_optimization( mmig_ntk, mmig_ps, &mmig_st );
    }
    const double mmig_wall_ms = to_milliseconds( mmig_wall_begin, std::chrono::steady_clock::now() );

    mmig_ntk = cleanup_dangling( mmig_ntk );
    const auto min_after = count_minority_nodes( mmig_ntk );

    std::cout << "\n[mMIG] enabled with stage=" << to_string( mmig_stage )
              << " max_iters=" << mmig_max_iters
              << " seed=" << ( mmig_seed_enabled ? to_string( mmig_seed_mode_ps ) : "off" )
              << " seed_budget=" << mmig_seed_budget
              << " seed_rounds=" << mmig_seed_rounds
              << " sr5=" << ( mmig_enable_sr5 ? "on" : "off" )
              << " advanced=" << ( mmig_use_advanced ? "on" : "off" )
              << " flow=" << to_string( mmig_flow )
              << " zero_gain=" << ( mmig_allow_zero_gain ? "on" : "off" )
              << " advanced_rounds=" << mmig_advanced_rounds
              << " adv_resub=" << ( mmig_adv_no_resub ? "off" : "on" )
              << " adv_cut=" << ( mmig_adv_no_cut ? "off" : "on" )
              << " adv_refac=" << ( mmig_adv_no_refac ? "off" : "on" )
              << " adv_exact=" << ( mmig_adv_no_exact ? "off" : "on" )
              << " adv_balance=" << ( mmig_adv_no_balance ? "off" : "on" )
              << " adv_obj_guard=" << ( mmig_adv_no_objective_guard ? "off" : "on" )
              << " adv_stag_stop=" << ( mmig_adv_no_stagnation_stop ? "off" : "on" )
              << " adv_stag_limit=" << mmig_adv_stagnation
              << " cec=" << ( mmig_cec_enabled ? "on" : "off" )
              << " dont_cares=" << ( mmig_use_dont_cares ? "on" : "off" )
              << " dc_threshold=" << mmig_dc_threshold
              << " dual_inv=" << ( mmig_no_dual_inv ? "off" : "on" )
              << " interleaved_seed=" << ( mmig_interleaved_seeding ? "on" : "off" )
              << " interleaved_seed_interval=" << mmig_interleaved_seed_interval
              << " auto_sr5=" << ( mmig_no_auto_sr5 ? "off" : "on" )
              << " rewrite_rank=" << ( mmig_no_rewrite_ranking ? "off" : "on" )
              << " normalize_inner=" << ( mmig_no_normalize_inner ? "off" : "on" )
              << " resub_rank=" << ( mmig_no_resub_ranking ? "off" : "on" )
              << " adv_tuned=" << ( mmig_adv_no_tuned_policy ? "off" : "on" ) << "\n";
    std::cout << "[mMIG] minority nodes: " << min_before << " -> " << min_after << "\n";
    print_mmig_stage_summary( mmig_st );
    print_stats( mmig_ntk, "Optimized mMIG" );

    if ( params.lut_report_k > 0 )
    {
      const std::string mmig_lut_base = out_base + "_mmig_opt";
      print_k_lut_stats_runtime(
          mmig_ntk, params.lut_report_k, "mMIG optimized network",
          params.write_lut_netlists ? &mmig_lut_base : nullptr );
    }

    write_verilog( mmig_ntk, out_base + "_mmig_opt.v" );
    write_blif( mmig_ntk, out_base + "_mmig_opt.blif" );

    std::cout << "  mMIG flow:                 " << format_seconds_value( mmig_wall_ms ) << " s\n";
    std::cout << "Wrote " << out_base << "_mmig_opt.{v,blif}\n";
  }

  return 0;
}
#endif
