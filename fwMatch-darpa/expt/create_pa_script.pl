#!/usr/bin/perl
use Cwd;
$EXPT_NAME = "shuting_test";
@EE = ("0"); #, "45", "90", "135", "180", "225", "270", "315");
$DEPTH = 1;
$MODEL_TYPE = "tree";
$DATA_DIR = "~/data/steph_225/";

for( $i = 0; $i <= $#EE; $i++){
    $DATA_FILE = sprintf("%s_%d_%02d", $EXPT_NAME, $EE[$i], $DEPTH);
    $EXPERIMENT = sprintf("%s_%d_%02d_%s", $EXPT_NAME, $EE[$i], $DEPTH, $MODEL_TYPE);
    mkdir($EXPERIMENT);
    $scommand = sprintf("cp %s_template/* %s/", $EXPT_NAME, $EXPERIMENT);
    print "Running: ".$scommand."\n";
    ($status, $result) = system($scommand);

    print "file: get_real_data.m\n";
    open(my $FID, ">", "$EXPERIMENT/get_real_data.m") 
	or die "cannot open < $!";
    print $FID "function [data,variable_names] = get_real_data()\n";
    print $FID "load(['".$DATA_DIR."' ...\n";
    print $FID "          '/all_".$DATA_FILE.".mat']);\n";
    print $FID "fprintf('Loaded: %s\\n', ['".$DATA_DIR."' ...\n".
	"'/all_".$DATA_FILE.".mat']);\n";
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
    print $FID "    'email_for_notifications', '1live.life.queen.size1\@gmail.com', ...\n";
    print $FID "    'yeti_user', 'et2495', ...\n";
    print $FID "    'compute_true_logZ', false, ...\n";
    print $FID "    'reweight_denominator', 'mean_degree', ...\n";
    print $FID "    's_lambda_splits', 6, ...\n";
    print $FID "    's_lambdas_per_split', 1, ...\n";
    print $FID "    's_lambda_min', 2e-03, ...\n";
    print $FID "    's_lambda_max', 5e-01, ...\n";
    print $FID "    'density_splits', 1, ...\n";
    print $FID "    'densities_per_split', 4, ...\n";
    print $FID "    'density_min', 0.01, ...\n";
    print $FID "    'density_max', 0.04, ...\n";
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
    print $FID "    'email_for_notifications', '1live.life.queen.size1\@gmail.com', ...\n";
    print $FID "    'yeti_user', 'et2495', ...\n";
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
