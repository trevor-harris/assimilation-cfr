# Proxy Influence

- Introduction
  - Paleoclimate and Data Assimilation
  - Formulate the problem and our solution using words
  - Functional Data
  - Data Depth
- Method
  - formulate the problem into comparing two ensembles of spatial random fields, and further treat the spatial random fields as two dimensional functional data. Then will use the functional data depth to evaluate the difference between the distributions of two sets of functions. The advantage of using functional data depth is that it is fully nonparametric and no assumpution of any sepecific distribution of either the marginal or joint distribuiton of those random fields.
  - Notations, form the hypothesis  
    - Define the data depth we use
    - Describe the test statistics to use. 
  - Talk about p values
    - small sample (permutation)
    - large sample (asymptotic)

- Simulation Study
  - Size
  - Power
- Application to Data Assimilation
- Conclusion



Questions:

Why did we need a new depth thats only a slight modification of the old one? Because ED is NAB which is too sensitive for us and XD is much faster to compute.

Why did we use a KS test instead of a t-test? Because the depths dont change predictably in terms of mean or variance with the data so a comprehensive nonparametric test is needed to find distributional differences

How can you get a pvalue when the depths arent being generated independently? We can still do a permutation test to attain a pvalue. We also show that they're almost independent and that the KS distribution provides a close approximation as n,m go to infinity

Why do you do compare D(x | x) with D(y | x) insead of with D(x | y)? I dont know.

