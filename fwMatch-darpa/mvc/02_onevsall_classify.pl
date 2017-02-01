#!/usr/bin/perl

use Env;

chdir("/u/4/e/et2495/src/crf/fwMatch/mnms/perf");

$this_hostname = `hostname`;
chomp($this_hostname);

$key = $ENV{'PBS_JOBID'};
#################################
#$key = "999999.warp.hpc1.cs.cmu.edu";

($job_id, $junk) = split(/\./, $key);
print "Job ID is: ".$job_id." from: ".$key." :) \n";
print "Running on ".$this_hostname."\n";

#$EXPT." ".$ORIENTS[$o]." ".$DEPTHS[$d]." ".$SUFFIX." ".$CPATH."\n";
#$filename = "/u/4/e/et2495/src/crf/fwMatch/mnms/perf/".
$filename = "s_".$ARGV[0]."_".$job_id;

print "File is: ".$filename."\n";

open(FID, ">".$filename.".m") or die("Error HERE $!\n");
print FID "EXPT = '".$ARGV[0]."';\n";
print FID "ORIENT = ".$ARGV[1].";\n";
print FID "DEPTH = ".$ARGV[2].";\n";
print FID "MODEL_TYPE = '".$ARGV[3]."';\n";
print FID "LPATH = '".$ARGV[4]."';\n";
print FID "MODEL_PREFIX = '".$ARGV[5]."';\n";
print FID "SCRIPT = 1;\n";
print FID "run('a_031_onevsall_classify_per_model.m');\n";
print FID "exit(0);\n";
close(FID);

print "about to system\n";
#print "/usr/local/bin/matlab -nodisplay -nojvm -nosplash -r $filename \n";
print "/usr/local/bin/matlab -nodisplay -nojvm -nosplash -r \"$filename; exit;\" \n";
system("/usr/local/bin/matlab -nodisplay -nojvm -nosplash -r \"$filename; exit;\" \n");

print "done system\n";

