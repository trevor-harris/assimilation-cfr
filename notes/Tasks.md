# Tasks

- Move jobs over to the campus cluster 

  - ~~Request access~~
  - Get access
  - Move jobs

- Implement basis function system

  1. ~~Scrap the Eigendecomp method~~
  2. ~~Create radial (bisquare) basis functions~~
     1. Check that these are accurate and everything lines up (semi done)

- Implement a (modified) version of Extremal data depth for handling noisy functional data

  1. ~~Implement standard ED~~
     1. Test these against easy functions (semi done)
  2. ~~Implement central regions~~
     1. Test these against easy functions (semi done)
  3. Modify ED to handle points being slightly outside 
     1. Think about what sort of tolerances should be allowed
     2. Ask Naveen for papers related to this
     3. Setup meeting with Naveen and Bo 

- Add Parallel processing for looping over time

  - Make sure results are not effected

- Compare ensemble with prior

  1. Establish method and output for single time point

  - What should results here look like? Graphs? a summary statistic? Check the Smerdon paper
  - Based on the above output, what does further investigation look like? How important are differences at a single site at a single time point? Are we more interested in aggregate effects?

- Compare ensemble with target ensemble

  - Need to get data for this
  - Implement or find FDR methods for these sort of tests
  - Rest same as above

