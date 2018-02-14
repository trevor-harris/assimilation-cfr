# Tasks

#### Implement basis function system

1. ~~Scrap the Eigendecomp method~~
2. ~~Create radial (bisquare) basis functions~~
   1. Check that these are accurate and everything lines up (semi done)
3. Add multi-resultion basis functions
   1. ~~Do you do one layer at a time? I.e. is res 2 fit on the residuals from res 1? Etc~~
   2. ~~Compare against lots of single resolution functions - Doesnt seem to help~~
4. Coefficients are very noisy - this probably needs to be handled in some way.
   1. Option 1: Modify ED (below) to allow tolerances
   2. Option 2: Modify depth function used to construct ED 
      1. Could allow tolerances or be based on local aggregation
   3. Option 3: Don't do anything - not ideal since prior always has ED = 0
   4. ~~Option 4: Shrink coefficients? terrible idea~~
5. Selecting an optimal number of basis functions.
   1. How does adding more basis functions minimize the MSE of smoothed field vs the original?
   2. min vs max on the radius

#### Implement a (modified) version of Extremal data depth for handling noisy functional data

1. ~~Implement standard ED~~
   1. Test these against easy functions (semi done)
2. ~~Implement central regions~~
   1. Test these against easy functions (semi done)
3. Modify ED to handle points being slightly outside 
   1. Think about what sort of tolerances should be allowed
   2. Ask Naveen for papers related to this
   3. Setup meeting with Naveen and Bo 

#### Compare ensemble with prior

1. Establish method and output for single time point
   - What should results here look like? Graphs? a summary statistic? Check the Smerdon paper
   - Based on the above output, what does further investigation look like? How important are differences at a single site at a single time point? Are we more interested in aggregate effects?
   - Check like 5 to 10 years evenly spaced over the 998 years 
   - turn coefficients back to map and color code

#### Compare ensemble with target ensemble

- Need to get data for this
- Implement or find FDR methods for these sort of tests
  - Benjamini-Yekuteli
- Rest same as above

#### Statistical Contribution to this Project

- What is it?

#### Misc

- ~~Add Parallel processing for looping over time~~
  - ~~make sure results not effected~~
- Move jobs over to the campus cluster
  - ~~Request access and get access~~
  - Move jobs
- setup meeting with collab next next week wednesday before noon or tuesday/thursday morning or friday afternoon

#### 

