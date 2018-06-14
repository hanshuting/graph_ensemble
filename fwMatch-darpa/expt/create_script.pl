#!/usr/bin/perl
if (not @ARGV){
    die "No experiment name(s) listed on command line\n";
}
use Cwd;
$EXPT_NAME = "temporal";
@EE = @ARGV;
$MODEL_TYPE = "loopy";
$DATA_DIR = "~/data/";
$USER = "jds2270";
$EMAIL = "jds2270\@columbia.edu";

for( $i = 0; $i <= $#EE; $i++){
    $DATA_FILE = sprintf("%s_%s", $EXPT_NAME, $EE[$i]);
    $EXPERIMENT = sprintf("%s_%s_%s", $EXPT_NAME, $EE[$i], $MODEL_TYPE);
    mkdir($EXPERIMENT);
    #$scommand = sprintf("cp %s_%s_template/* %s/", $EXPT_NAME, $EE[$i],$EXPERIMENT);
    $scommand = sprintf("cp %s_template/* %s", $EXPT_NAME, $EXPERIMENT);
    print "Running: ".$scommand."\n";
    ($status, $result) = system($scommand);

    if($MODEL_TYPE eq "loopy"){
    open(my $FID, ">", "$EXPERIMENT/write_configs_for_loopy.m")
        or die "cannot open < $!";
    print $FID "create_config_files( ...\n";
    print $FID "    'experiment_name', '".$EXPERIMENT."', ...\n";
    print $FID "    'email_for_notifications', '".$EMAIL."', ...\n";
    print $FID "    'yeti_user', '".$USER."', ...\n";
    print $FID "    'compute_true_logZ', false, ...\n";
    print $FID "    'reweight_denominator', 'mean_degree', ...\n";
    print $FID "    's_lambda_splits', 6, ...\n";
    print $FID "    's_lambdas_per_split', 1, ...\n";
    print $FID "    's_lambda_min', 2e-03, ...\n";
    print $FID "    's_lambda_max', 5e-01, ...\n";
    print $FID "    'density_splits', 1, ...\n";
    print $FID "    'densities_per_split', 6, ...\n";
    print $FID "    'density_min', 0.1, ...\n";
    print $FID "    'density_max', 0.3, ...\n";
    print $FID "    'p_lambda_splits', 5, ...\n";
    print $FID "    'p_lambdas_per_split', 1, ...\n";
    print $FID "    'p_lambda_min', 1e+01, ...\n";
    print $FID "    'p_lambda_max', 1e+04, ...\n";
    print $FID "    'time_span', 2);\n";
    close($FID);
    print "done writing write_configs_for_loopy.m\n";
    }else{
    #model is tree
    die "Only loopy model is supported.\n"
    # open(my $FID, ">", "$EXPERIMENT/write_configs_for_tree.m")
    #     or die "cannot open < $!";
    # print $FID "create_config_files( ...\n";
    # print $FID "    'experiment_name', '".$EXPERIMENT."', ...\n";
    # print $FID "    'email_for_notifications', '".$EMAIL."', ...\n";
    # print $FID "    'yeti_user', '".$USER."', ...\n";
    # print $FID "    'structure_type', 'tree', ...\n";
    # print $FID "    'compute_true_logZ', true, ...\n";
    # print $FID "    'p_lambda_splits', 20, ...\n";
    # print $FID "    'p_lambdas_per_split', 1, ...\n";
    # print $FID "    'p_lambda_min', 1e+01, ...\n";
    # print $FID "    'p_lambda_max', 1e+07);\n";
    # close($FID);

    # print "done writing write_configs_for_tree.m\n";
    }

    print "file: get_real_data.m\n";
    open(my $FID, ">", "$EXPERIMENT/get_real_data.m")
        or die "cannot open < $!";
    print $FID "function [data,variable_names, stimuli] = get_real_data()\n";
    print $FID "load(['".$DATA_DIR."' ...\n";
    print $FID "          '".$DATA_FILE.".mat']);\n";
    print $FID "fprintf('Loaded: %s\\n', ['".$DATA_DIR."' ...\n".
        "'".$DATA_FILE.".mat']);\n";
    print $FID "%data is time_frames by number_of_neurons\n";
    print $FID "data = full(data);\n";
    print $FID "N = size(data,2);\n";
    print $FID "fprintf('data is : %d, %d\\n', size(data,1), size(data,2));\n";
    print $FID "variable_names = {};\n";
    print $FID "for i = 1:N\n";
    print $FID "\tvariable_names(end+1) = {int2str(i)};\n";
    print $FID "end\n";
    print $FID "if exist('stimuli', 'var') == 1\n";
    print $FID "    stimuli = full(stimuli);\n";
    print $FID "    fprintf('stimuli is : %d, %d\\n', size(stimuli,1), size(stimuli,2));\n";
    print $FID "else\n";
    print $FID "    stimuli = [];\n";
    print $FID "end\n";
    print $FID "end\n";
    close($FID);

    print "done writing get_real_data.m\n";

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
