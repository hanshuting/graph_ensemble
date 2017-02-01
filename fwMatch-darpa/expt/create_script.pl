#!/usr/bin/perl
use Cwd;
$EXPT_NAME = "m21_d2_vis";
@EE = ("all_high_0.4_28","all_high_0.4_29","all_high_0.4_30","all_high_0.6_1","all_high_0.6_2","all_high_0.6_3","all_high_0.6_4","all_high_0.6_5","all_high_0.6_6","all_high_0.6_7","all_high_0.6_8","all_high_0.6_9","all_high_0.6_10","all_high_0.6_11","all_high_0.6_12","all_high_0.6_13","all_high_0.6_14","all_high_0.6_15","all_high_0.6_16","all_high_0.6_17","all_high_0.6_18","all_high_0.6_19","all_high_0.6_20","all_high_0.6_21","all_high_0.6_22","all_high_0.6_23","all_high_0.6_24","all_high_0.6_25","all_high_0.6_26","all_high_0.6_27","all_high_0.6_28","all_high_0.6_29","all_high_0.6_30","all_high_0.8_1","all_high_0.8_2","all_high_0.8_3","all_high_0.8_4","all_high_0.8_5","all_high_0.8_6","all_high_0.8_7","all_high_0.8_8","all_high_0.8_9","all_high_0.8_10","all_high_0.8_11","all_high_0.8_12","all_high_0.8_13","all_high_0.8_14","all_high_0.8_15","all_high_0.8_16","all_high_0.8_17","all_high_0.8_18","all_high_0.8_19","all_high_0.8_20","all_high_0.8_21","all_high_0.8_22","all_high_0.8_23","all_high_0.8_24","all_high_0.8_25","all_high_0.8_26","all_high_0.8_27","all_high_0.8_28","all_high_0.8_29","all_high_0.8_30");
#@EE = ("all_high_add_neuron_100", "all_high_add_neuron_200", "all_high_add_neuron_300", "all_high_add_neuron_400", "all_high_add_neuron_500", "all_high_add_neuron_600", "all_high_add_neuron_700", "all_high_add_neuron_800", "all_high_add_neuron_900", "all_high_add_neuron_1000", "all_high_add_neuron_1100", "all_high_add_neuron_1200", "all_high_add_neuron_1300", "all_high_add_neuron_1400", "all_high_add_neuron_1500", "all_high_add_neuron_1600", "all_high_add_neuron_1700", "all_high_add_neuron_1800", "all_high_add_neuron_1900", "all_high_add_neuron_2000");
$MODEL_TYPE = "loopy";
$DATA_DIR = "/vega/brain/users/sh3276/data/luis";

for( $i = 0; $i <= $#EE; $i++){
    $DATA_FILE = sprintf("%s_%s", $EXPT_NAME, $EE[$i]);
    $EXPERIMENT = sprintf("%s_%s_%s", $EXPT_NAME, $EE[$i], $MODEL_TYPE);
    mkdir($EXPERIMENT);
    #$scommand = sprintf("cp %s_%s_template/* %s/", $EXPT_NAME, $EE[$i],$EXPERIMENT);
    $scommand = sprintf("cp %s_template/* %s/", $EXPT_NAME, $EXPERIMENT);
    print "Running: ".$scommand."\n";
    ($status, $result) = system($scommand);

    print "file: get_real_data.m\n";
    open(my $FID, ">", "$EXPERIMENT/get_real_data.m") 
	or die "cannot open < $!";
    print $FID "function [data,variable_names] = get_real_data()\n";
    print $FID "load(['".$DATA_DIR."' ...\n";
    print $FID "          '/".$DATA_FILE.".mat']);\n";
    print $FID "fprintf('Loaded: %s\\n', ['".$DATA_DIR."' ...\n".
	"'/".$DATA_FILE.".mat']);\n";
    print $FID "%data is time_frames by number_of_neurons\n";
    print $FID "data = full(data);\n";
    print $FID "N = size(data,2);\n";
    print $FID "variable_names = {};\n";
    print $FID "for i = 1:N\n";
    print $FID "\tvariable_names(end+1) = {int2str(i)};\n";
    print $FID "end\n";
    print $FID "fprintf('data is : %d, %d\\n', size(data,1), size(data,2));\n";
    print $FID "end\n";
    close($FID);
    
    print "done writing get_real_data.m\n";

    if($MODEL_TYPE eq "loopy"){
    open(my $FID, ">", "$EXPERIMENT/write_configs_for_loopy.m") 
	or die "cannot open < $!";    
    print $FID "create_config_files( ...\n";
    print $FID "    'experiment_name', '".$EXPERIMENT."', ...\n";
    print $FID "    'email_for_notifications', 'sh3276\@columbia.edu', ...\n";
    print $FID "    'yeti_user', 'sh3276', ...\n";
    print $FID "    'compute_true_logZ', false, ...\n";
    print $FID "    'reweight_denominator', 'mean_degree', ...\n";
    print $FID "    's_lambda_splits', 6, ...\n";
    print $FID "    's_lambdas_per_split', 1, ...\n";
    print $FID "    's_lambda_min', 2e-05, ...\n";
    print $FID "    's_lambda_max', 5e-03, ...\n";
    print $FID "    'density_splits', 1, ...\n";
    print $FID "    'densities_per_split', 6, ...\n";
    print $FID "    'density_min', 0.25, ...\n";
    print $FID "    'density_max', 0.3, ...\n";
    print $FID "    'p_lambda_splits', 5, ...\n";
    print $FID "    'p_lambdas_per_split', 1, ...\n";
    print $FID "    'p_lambda_min', 1e+01, ...\n";
    print $FID "    'p_lambda_max', 1e+04);\n";
    close($FID); 
    print "done writing write_configs_for_loopy.m\n";
    }else{
	#model is tree
    open(my $FID, ">", "$EXPERIMENT/write_configs_for_tree.m") 
	or die "cannot open < $!";    
    print $FID "create_config_files( ...\n";
    print $FID "    'experiment_name', '".$EXPERIMENT."', ...\n";
    print $FID "    'email_for_notifications', 'sh3276\@columbia.edu', ...\n";
    print $FID "    'yeti_user', 'sh3276', ...\n";
    print $FID "    'structure_type', 'tree', ...\n";
    print $FID "    'compute_true_logZ', true, ...\n";
    print $FID "    'p_lambda_splits', 20, ...\n";
    print $FID "    'p_lambdas_per_split', 1, ...\n";
    print $FID "    'p_lambda_min', 1e+01, ...\n";
    print $FID "    'p_lambda_max', 1e+07);\n";
    close($FID);

    print "done writing write_configs_for_tree.m\n";
    }

    $curr_dir = cwd();
    print "curr_dir = ".$curr_dir."\n";
    chdir($EXPERIMENT);
    print "changed into dir: ".cwd()."\n";
    $scommand = "matlab -nodesktop -nodisplay -r ".
	"\"try, write_configs_for_".$MODEL_TYPE.", catch, end, exit\"";
    print "About to run:\n".$scommand."\n";
    system($scommand);
    print "Done running system command\n";
    chdir($curr_dir);
    print "changed into dir: ".cwd()."\n";
}
