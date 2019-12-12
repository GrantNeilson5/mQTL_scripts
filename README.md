# mQTL_scripts

These scripts will allow an mQTL analysis to be carried out,

1) start with filtering the genotypes for use in plink for use in mQTL analysis 
2) Then reformat the matrix files so they are in the correct formatt for analysis
3) then run the model, this should be ran twice once with and once without the principle components, notes on how to do this can be found
in the script,
It is sometimes best to run only on chromosome 22 first to make sure it works properly, this shouldnt take too long as chr 22 is the smallest

4) Then run Summarize mQTL data, again this needs to be run twice, once for the analysis with the 10PC and ocne without
5) finally visualise the data 
