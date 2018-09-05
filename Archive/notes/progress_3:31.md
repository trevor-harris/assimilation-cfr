##This Week

1. Got the new 'many to many' wilcoxon test working. This is how it would work for 16 regions, 1000 projections, and 100 permutations.
   1. Algorithm (part 1) - Wilcoxon
      1. Divide each prior ensemble field and posterior ensemble field into 16 disjoint evenly sized regions
      2. generate 1000 random projections (using rnorm()). These are fixed so each region and ensemble member gets projected with the SAME set of projection vectors.
      3. for region 1:
         1. project each prior and post ensemble member into R (results in 2 sets of 100 real numbers) using the first random projection.
         2. calculate the wilcoxon rank sum test between the two sets of numbers and store the statistic
         3. repeat step 3.1 and 3.2 for the next 999 random projections.
         4. This results in 1000 Wilcoxon rank sum test statistics. Take the mean.
      4. Repeat step 3 for the next 15 regions.
      5. After all 16 regions are completed we get a 4x4 field of (observed) mean Wilcoxon statistics
   2. Algorithm (part 2) - Permutation
      1. Generate a random permutation matrix.
      2. For each region and each ensemble member permute the prior and the posterior region using the same permutation matrix generated in step 1.
      3. Run Algorithm (part 1) on the newly permuted prior and posterior to get a new field of (permuted) mean wilcoxon statistics.
      4. Repeat step 1, 2, and 3 100 times to get a distribution of permuted wilcoxon fields.
   3. Use Extremal Depth to find the depth of the observed Wilcoxon field agaist the set of permuted Wilcoxon fields.
2. Some unverified assumptions I made:
   1. The same set of projections should be used for each region and each ensemble member. 
   2. the same permutation should be applied to each region and each ensemble member
   3. Projections are RANDOMLY generated using a N(0,1) distribution.
      1. Still looking into space filling designs.
3. The resulting wilcoxon fields are far smoother than the original basis coefficients. Also there is greater variability between the wilcoxon fields than the basis coefficients so the "slightly peaking out from the domain" problem looks to be potentially resolved.
4. Tested between the mean and the median for the wilcoxon tests. It appears the mean is more stable and has a tighter range. Even for low numbers of resamples. Weird?