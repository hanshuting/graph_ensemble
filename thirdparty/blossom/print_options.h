#pragma once
#include <sstream>
#include <string>
#include "PerfectMatching.h"

inline std::string print_options(const PerfectMatching::Options &opts) {
  std::stringstream ss;
  ss << "fractional_jumpstart = " << opts.fractional_jumpstart << "\n"
     << "dual_greedy_update_option = " << opts.dual_greedy_update_option << "\n"
     << "dual_LP_threshold = " << opts.dual_LP_threshold << "\n"
     << "update_duals_before = " << opts.update_duals_before << "\n"
     << "update_duals_after = " << opts.update_duals_after << "\n"
     << "single_tree_threshold = " << opts.single_tree_threshold << "\n"
     << "verbose = " << opts.verbose << "\n"
     ;
  return ss.str();
}

inline std::string interrogate(const PerfectMatching &pm) {
  std::stringstream ss;
  ss << "node_num" << pm.node_num << "\n"
     << "edge_num" << pm.edge_num << "\n"
     << "edge_num_max" << pm.edge_num_max << "\n"
     << "tree_num" << pm.tree_num << "\n"
     << "tree_num_max" << pm.tree_num_max << "\n"
     ;
  return ss.str();
}
