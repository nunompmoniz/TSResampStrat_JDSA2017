Code used to load, model, run experiments and obtain results.<br/>
<br/>
<b>Exps.R</b>: After loading the respective data set (instructions here), this is the code to run the experiments. They use Monte Carlo estimates, with 50 repetitions using 50% of the data as training and the subsequent 25% as test data. As referenced in the article, for data sets 21 and 22 we used 10%/5% and for data sets 23 and 24, 20%/10%<br/><br/>
<b>GetResults.R</b>: After running the experiments, execute this code to obtain a table with the overall results.<br/><br/>
<b>PairedComparisons.R</b>: After obtaining the overall results, this code is executed to perform the paired comparisons using Wilcoxon signed rank tests in order to infer the statistical significance (with p-value < 0.05). This is used to prove the hypothesis set forth in the article submitted.

