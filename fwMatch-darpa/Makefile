# Master Makefile

all:
	$(MAKE) -f compile.mk
	$(MAKE) -f cook.mk
	$(MAKE) -f report.mk


# EXPERIMENTS (Put here so we can run the experiments on demand, one at a time)

# Learns a graph of 10 nodes and plots the true and estimated graph comparison
expt/structure_learning_graph_plot/config1_10_node_graph: expt/structure_learning_graph_plot/config1.m expt/structure_learning_graph_plot/run.m expt/structure_learning_graph_plot/plot_results.m
	./run_and_plot.sh structure_learning_graph_plot 1

# Compares the effect in the AUC if we increase the number of samples from daily to hourly
expt/structure_learning_daily_to_hourly_data/config1: expt/structure_learning_daily_to_hourly_data/config1.m expt/structure_learning_daily_to_hourly_data/run.m expt/structure_learning_daily_to_hourly_data/plot_results.m
	./run_and_plot.sh structure_learning_daily_to_hourly_data 1

# Compares the two ways of building an adjencency matrix from the lasso logistic regression output (min and max)
expt/structure_learning_min_vs_max_graph/config1: expt/structure_learning_min_vs_max_graph/config1.m expt/structure_learning_min_vs_max_graph/run.m expt/structure_learning_min_vs_max_graph/plot_results.m
	./run_and_plot.sh structure_learning_min_vs_max_graph 1

# Varies the regulatization parameters and see how AUC varies for structure learning
expt/structure_learning_regularization/config1_absolute_lambda: expt/structure_learning_regularization/config1.m expt/structure_learning_regularization/run.m expt/structure_learning_regularization/plot_results.m
	./run_and_plot.sh structure_learning_regularization 1