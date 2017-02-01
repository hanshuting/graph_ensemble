fwMatch -- Frank-Wolfe Bethe Matching
===
**Note to CU/DARPA/Bloomberg collaborators: You should be working in the `darpa` branch, NOT master. Do not push changes to this branch; we are working on other projects here!** I will occasionally move changes between the two branches in a way that does not break anything.

Please read about [git branches](Git Branches -- http://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging) if you are not yet familiar with this.

## Financial Data

*"rawData" has been cleaned (checked for matching of trading dates, missing data, etc.), and the values are adjusted closed prices. 

*For financial data, except for weekends, the market is also closed on public holidays, and occasionally on natural disasters (hurricanes, for example). 

*The binarization is based on the equal frequency bucketing (deciles, as of Feb 24 2015), applied on the daily percentage change: r_t = \frac{I_t - I_{t-1}}{I_{t-1}}, where I_t is the index value on t-th day. For employment-related statistics, first order differences (as opposed to percentage changes) are used.

 For volatility index and employment-related statistics, the original values are also maintained (and discretized).

*The binarized financial time series is inside  "EFI" folder, in .csv files

*Each .csv file corresponds to one asset class\sub-class.

*Each column represents one node. For example, for sticker GSPC (i.e. S&P500), it has 10 nodes (at this moment), GSPC-1, GSPC-2
....GSPC-10. For macroeconomic indicators, "NR" means "Not a Release date". For quarterly-updated indicators (like GDP), 2 binary states are used (above/below median). For monthly-updated indicators, quartiles are used.

*The pdf document is updated to reflect the current index list; an Excel file containing documentation for the Financial/Economic binary series is added in "EFI" folder.

=====The following is description of the consolidated program contained in "FinMatlabFunctions" folder======

*The function finBnry.m takes one Excel file containing the original time series data as input, and outputs the binary states time series.

*To use the function, call:   [allBnrySeries,nanInfo] = finBnry('RawDataCombined.xlsx','daily','EFI',true)

*The input raw data file is as specified in 'RawDataCombined.xlsx', with the first sheet being the raw time series, and the second
sheet specifying the type and updating frequency of each index. The order of the indices across the two sheets must be the same.


Getting Started
===
 - Compile all dependencies and processed data: `make`. You will need Matlab and the `mex` compiler set up with a recent clang or g++.
 - To run an experiment, type `./run.sh NAME ID`. NAME is a subdirectory of `expt/` and ID is a number denoting which experiment to run.
	 - An example demonstrating the features of the framework is in `expt/test`.
	 - An example showing how to call the graphical model inference code is in `expt/binaryDenoise`.
 - To run Matlab interactively, start Matlab in this directory and type `startup`.

Organization
===
 - cp/ -- checkpoint for running algorithms and storage for intermediate data. Avoid including these files in the git repository, although small files are OK. Eventually we should track this stuff with git-annex and synchronize the data with some shared server everyone has access to.
 - data/ -- raw data, NOT included in the git repository.
 - fig/ -- store pdf or eps figures here. Do NOT use raster file formats like jpg, png.
 - src/ -- library source code. Algorithms and functions go here.
 - run/ -- scripts. Test codes, exploratory programs, interactive programs, etc. go here. Do NOT include full experiments with any results we may need to reproduce -- use the experiment framework instead.
 - results/ -- store output .txt, .csv, and .mat files here. Do NOT include in the git repository.
 - tex/ -- writeups.

Experiment Framework
===
Please use the experiment framework to organize and run any experiment. This keeps everything organized, allows easy perusal of results, and allows us to quickly submit jobs to clusters (by having just a single command line).

Experiment configuration and driver code are stored in `expt/`. Each experiment lives in a single subdirectory and has a single driver script named `run.m`, which takes a structure of parameters. These parameters are stored in files named `config1.m` ... `configN.m`. Each config file defines a `params` struct, which is passed to the experiment driver (`run.m`). Additional fields of `saveFile`, `checkpointFile`, and `logFile` will be inserted into `params`, which by default will be named `result<n>.mat`, `config<n>_cp.mat`, and `config<n>_log.txt`.

Save both the final answers and any intermediate data we may want to analyze later in the provided `saveFile`.

You can also write a report template, which given configuration and results for each experimental run, generates a data table. See the example for more details.

Make sure to commit the code before running the experiment, so that a snapshot of what you ran is recorded in git! We want to be able to reproduce any result later.

I will eventually write some logic in the shell script to detect whether the experiment you are trying to run has any uncommitted files and give an error, but in the meantime, you should be sure to check that yourself.

Makefiles
===
 - _Every_ generated file should be tracked using Makefiles. This includes not only binary executables, but all intermediate data, figures, and other files. Anybody should be able to type a single make command to generate your figure or report output. There should be _no_ manual command link invocations to remember or communicate.
 - We may eventually move to more fully-featured solution like [Luigi](https://github.com/spotify/luigi), if anybody volunteers to set it up. For now, Makefiles are easier to get started with so we will use those.
 - The master `Makefile` simply includes three Makefiles, one for each type of file.
	 - `compile.mk` lists how to compile C/C++ code portions of the codebase.
	 - `cook.mk` lists how to convert raw data to every stage of intermediate data.
	 - `report.mk` lists how to figures, tables, and other reports from some intermediate data result, and how to compile the LaTeX files that include these autogenerated data.
 - To the extent practical, _avoid copying and pasting result data_. Write a script to read an intermediate data file and write out a csv or text file. This way, if we ever update any part of the pipeline, typing `make report` will automatically run the steps necessary to 
 - There is a script `thirdparty/matrix2latex.m` which creates a LaTeX table directly from a Matlab matrix, which you can then save and include in your writeup. Generating this table could then be another make rule.

Data Format for Wainwright's Structure Learning
===
To figure out.

Data Format for Parameter Learning
===
I recommend reading through `tex/mrf.pdf` to get a better sense of the notation which is a bit non-standard and which I copy in the code. Don't worry about following the derivations; just understand what each symbol means and what is its physical shape as a matrix..

The constructor in `src/Ising.m` looks like this:

        function th = Ising(YN, YE, Ut, Vt, Ns, edges, lambda, varargin)
            % th = Ising(y, U, V, edges, lambda, opts)
            %
            %   Initialize the objective. 
            %   YN, YE  : Overcomplete node and edge labels
            %   Ut, Vt  : TRANSPOSED vercomplete node and edge feature matrices
            %             Row index feature and columns index sample/node.
            %             This layout makes our products more efficient.
            %   Ns      : M vector listing number of nodes in each sample            
            %   edges   : M cell vector of 2 x nEdges(m) edge list
            %   lambda  : Regularization constant
            %   opts    : Tuning options blah blah blah

`YN` has 2 columns -- column 1 denotes whether the node is zero, 2 denotes whether it is 1, and `YE` has 4 columns,
denoting whether node (i,j) with i < j is 00, 01, 10, 11. That is, the first column corresponds to variables i and j both being 0, the second to variable i being zero but j being one, the third to variable i being one and j being zero, and the fourth to both being 1.

(Note: for libDAI, the order of the middle two potentials is flipped.)

`Ut` and `Vt` are the tranposed feature matrices. The columns of `Ut` should equal the rows of `YN`, and the columns of `Vt` should equal the rows of `YE`. These features are used to learn conditional random fields. However, for now we will not be using features, so we will set them to be long row vectors of ones. Please see Section 1 of `tex/mrf.pdf` if this confuses you.

We have stacked the data for all samples into long matrices. This is for speed and memory. The `.ixN(m)` and `.ixE(m)` return the row indices corresponding to sample `m`, which you can use for selection.


In more detail,
 - 


Data
===
...

Matching Solvers
===
 - We use Komogorov's code QPBO and BlossomV. See the files in `thirdparty` for details.

