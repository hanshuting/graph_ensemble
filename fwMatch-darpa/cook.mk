# Makefile to cook raw data into algorithm inputs (cooked data).
#
# Important: All targets should list the SOURCE FILE of any programs run in
# addition to the input data files. This will cause any update of the program
# source to regenerate any working files that depend on the program.
#
# Also add configuration files!
#
# A useful idiom is
#
# target: script input
# 	interpreter $^ $@
#
# Will invoke interpreter with the command line `interpreter script input target`
# which is often exactly what you want.
#
# TODO: Add a downloader/checker step to enable an empty start.

# Handy variable reference (to implement Don't Repeat Yourself)
# $@ -- name of the target
# $< -- name of the first prerequisite
# $^ -- name of all prerequisites, separated by spaces

# The runDataParsing script will generate the .mat dataset in the cp folder
# namely All_Data.mat (The whole data), Clean_Data (Due to the mismatch of rows 
# between NLP and Financial datasets, there will be NaN in the All_Data set. 
# Clean_Data.mat keeps all the rows without NaN), and Meta_Data.mat (contains all
# the meta_info of the whole datasets)

cp/All_Data.mat cp/Clean_Data.mat cp/Meta_Data.mat: src/DataParsing/runDataParsing.m 
	./runDataParsing.sh

all: data/fin_data/CookedDataCombined.xlsx

data/fin_data/CookedDataCombined.xlsx: src/FinMatlabFunctions/finBnry.m data/fin_data/RawDataCombined.xlsx
	matlab -nojvm -r "finBnry('data/fin_data/RawDataCombined.xlsx', 'daily', 'EFI', 'data/fin_data/CookedDataCombined.xlsx')"

# if somehow you need java, then switch to -nodesktop -nosplash
# 

