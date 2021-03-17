
The performance data frame extracted from the ATC simulator logs
 is available in experiment_data.csv, and in img/dats_raw.RData. 

 1-extract.R contains the code that was used to extract this data from the
 (very large) raw XML logs. This script does not need to be run, as the extracted
 data is already provided in the repository. 

In the extracted data (experiment_data.csv and img/dats_raw.RData), the relevant columns are ppt (participant number),
 sess (session number), block (block number within a session),
 cond (experimental condition), stimulus (was the presented
 aircraft pair a conflict - c- or nonconflict - n),
 failtrial (did the automation make an incorrect recommendation, i.e., fail, on that trial,
 true or false), RT (response time in milliseconds), and R (did the participant
 respond conflict - C - or non-conflict - N). The 'response', 'score', and 'cumulative score'
 columns can be ignored, they are internal variables from the simulator that were only extracted for checking purposes.

All the paper's analyses are included in scripts numbered 2-5. Note that before the data was analysed,
non-responses and RTs <200ms were removed (as described in the manuscript). This 'cleaned' data is provided
in img/dats_clean.RData.

 Scripts 2 & 2.5 (standard analyses) stand alone just run from this data.
LBA analyses (scripts 3-5) were performed with a bespoke version of the 
dynamic models of choice R suite. The particular version
used is provided in the 'dmc' submodule of this repo. Standard dmc releases can be found at https://osf.io/pbwx8/, 
along with lessons and more information about the software. Script 3-model_specification.R
creates the LBA models reported in the paper. The most important model (reported in text) is the first one, which is stored in 
CA_top_samples. After models are created they must be sampled to obtain 
posterior distributions of parameter estimates. R/grid_dispatch.R dispatches sampling on a grid system
running pbs pro. Sampling could be performed on an individual computer by using h.RUN.dmc (note that run time would likely 
be overnight or longer). Once the model is fitted, the modeller should create a subdirectory samples_data/,
and the samples object should be stored inside. Other objects, such as posterior predictive data,
will also be stored in samples_data in later scripts (4-5). 


                 