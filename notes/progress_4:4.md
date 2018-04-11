## Last Week

1. Worked on implementing the wilcoxon permutation test
   1. On hold until after projection depth and limiting p values papers are read 
2. Read the limiting p values paper
   1. Basic idea is that given x  ~ F and some functional g(x) we use the bootstrap distribution of g*(x) to see how likely g(x) was with data depth
3. Worked through (some of) the background material in the "Functional Shape Analysis" book
4. Started making presentation on the Manifold paper

## This Week

1. ~~Read the projection depth paper~~
2. Update/fix/rewrite the wilcoxon permutation test 
   1. ~~Verify projection depth is being calculated correctly~~
      1. How will using a sampling plan instead of random projections work?
   2. ~~limiting p values paper uses bootstrap distribution. Will permutation distribution also work?~~
3. Figure out answers to Nathans questions
   1. I'm wondering why the frequency of exceeding a particular central region is highest in locations without many proxies. If I understand things correctly, my intuition would be that you would see the opposite result.
   2. I'm also not sure what to make of the circular nature of all the features in the spatial plots. I assume that the features are due somehow to the basis knots. Patterns in climate data usually aren't so circular, so I'm wondering if they might be a by-product of the analysis?
4. Read the Xie paper on "A Geometric Approach to Visualization of Variability in Functional Data"
5. simulation study for naveen method 
   1. Same dist on prior and posterior vs different dist	
   2. Interesting results on which "way" to permute
      1. pointwise permutation breaks the "continuity" of the functions and so I probably shouldnt
      2. functionwise permutation seems to be the way to go. Also makes sense since our data is the functions themselves and not the observed points.
6. Think about how to combine ED with Manifolds

