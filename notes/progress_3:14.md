## Last Week

1. Redid MSE and BIC calculations
   1. Now showing to use many more basis functions (28x30 lon x lat)
      1. Actually I sped up the computations enough to make 30x30 a very manageble default. I need to remake the graphs though.
2. Wrote code to compare horizontally (i.e. for each ensemble member)
   1. Many ED values are > 0.05 so in the majority of cases the prior is **within **the 95% central region of the posterior. It seems that for a given ensemble member the posteriors are not that different from the prior.
   2. Also tested running ED on the entire spatial fields (i.e. not on the basis coefficients) to see how results differed from results on the basis coefficients
      1. Even for a small number of basis (say 50) the ED doesn't change all that much from the ED on the spatial field itself. More basis functions —> closer "basis" ED is to the "spatial field" ED. 30x30 basis was very close.
3. Wrote code to compare vertically (i.e. for each time point)
   1. I implemented the data depth based Wilcoxon sign rank test from the Romo and Pintado paper to compare ensemble to ensemble
      1. Without correcting for multiple comparisons there are only a few time points where the prior and posterior test as being significantly different. It seems that for a time point the posterior is that not different form the prior. When doing an FDR correction there are no differences. 
      2. I didnt see an analytical way to get the distribution of W (the test statistic) under H_0 so I simulated it based on the description.
         1. Ideally I would like an explicit form, but the p values are stable out to 3 or 4 decimal places when I resimulate the distribution over and over.
         2. In our case the distribution under H_0 is the sum of 100 integers sampled without replacement from 1 to 200.
   2. Instead of doing a wilcoxon test I think it also might be worthwhile to count the number of posteriors (or priors) outside the X% confidence intervals of the prior (posterior). 
   3. I'm less confident in the "vertical" results than the "horizontal" ones since I'm not sure that its valid to just adapt the Wilcoxon test to use ED.
4. Read more data depth
   1. The Jornsten, Vardi, Zhang paper on data depth clustering seems to just be a robust verion of k-means using data depth instead of euclidean distance. 

## This Week

1. Try to figure out if the Wilcoxon sign rank test is the way to go or not. 
   1. If so make sure that the implementation is correct.
      1. In particular make sure that ties between the ranks are being handled appropriately.
   2. Find a non Monte Carlo based way of getting the distribution under H_0?
   3. Meet with Naveen?
   4. Alternatively use or adapt a Kolmogorov-Smirnov like test?
      1. infact this may even be more appropriate since we're trying to compare distributions
2. Create MSE and BIC plots
3. Create plots for both vertical and horizontal comparisons
   1. Horizontal (first priority)
      1. How often each location is outside a given central region
         1. We're going to need a smaller number than 95% here since many fields have high ED. Maybe 50% will give us enough.
      2. Average distance outside the X% region
   2. Vertical 
      1. Number of prior spatial fields outside the posterior X% central region? Or vice versa
      2. more?
4. Start compiling results into a latex document.
5. start reading Riemannian papers
   1. CSDA_fdaTB_final —> CSDA_PCA_12 —> Rosenthal





