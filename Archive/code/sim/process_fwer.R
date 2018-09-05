rm(list = ls())

matsplitter<-function(M, r, c) {
  # splits 1 matrix into c MxR matricies
  # I have no idea how this works
  
  rg <- (row(M)-1)%/%r+1
  cg <- (col(M)-1)%/%c+1
  rci <- (rg-1)*max(cg) + cg
  N <- prod(dim(M))/r/c
  cv <- unlist(lapply(1:N, function(x) M[rci==x]))
  dim(cv)<-c(r,c,N)
  cv
}

regions = 64
sims = 10
batches = 50

data.dir = paste0("/Users/trevh/research/assimilation-cfr/simdata/alpha6/")
list.files(data.dir)
files = list.files(data.dir)

dep = array(0, dim=c(regions, sims, batches))
bonf = array(0, dim=c(regions, sims, batches))
pw = array(0, dim=c(regions, sims, batches))
val = array(0, dim=c(regions, sims, batches))
pval = matrix(0, sims, batches)

d = 1; b = 1; p = 1; v = 1; e = 1
for(f in files) {
  if(grepl("Depth", f)) {
    dep[,,d] = readRDS(paste0(data.dir, f))
    d = d+1
  }
  if(grepl("Bonferroni", f)) {
    bonf[,,b] = readRDS(paste0(data.dir, f))
    b = b+1
  }
  if(grepl("Pointwise", f)) {
    pw[,,p] = readRDS(paste0(data.dir, f))
    p = p+1
  }
  if(grepl("Values", f)) {
    val[,,v] = readRDS(paste0(data.dir, f))
    v = v+1
  }
  if(grepl("pvals", f)) {
    pval[,e] = readRDS(paste0(data.dir, f))
    e = e+1
  }
}
# s = 10
# dep.sub = dep[,,s]
# val.sub = val[,,s]
# mean(colSums(val.sub > dep.sub) > 0)

# reformat 
flatten = function(mat) {
  matrix(mat, dim(mat)[1], prod(dim(mat)[2:3]))
}

d_alpha = flatten(dep)
b_alpha = flatten(bonf)
p_alpha = flatten(pw)
ks_vals = flatten(val)


# plot(d_alpha[,1], type = "l", ylim = c(0, 0.4))
# for(i in 2:2500) {
#   lines(d_alpha[,i])
# }


# compare familywise error rate
fwer = function(alphas, obs) {
  mean(sapply(1:ncol(obs), function(r) any(obs[,r] > alphas[,r])))
}

cat(" Depth: ", fwer(d_alpha, ks_vals), "\n",
    "Bonfe: ", fwer(b_alpha, ks_vals), "\n",
    "Point: ", fwer(p_alpha, ks_vals), "\n")




