library(extdepth)
library(ncdf4)

source('/home/trevorh2/assimilation-cfr/code/prep_functions.R')
source('/home/trevorh2/assimilation-cfr/code/ks_field_functions.R')

args = commandArgs(TRUE)

batch = as.double(args[1])
times = 1:10 + 10*(batch - 1)
if(batch == 100) times = 991:998


# read ensembles and prior ncdf4 objects
nc.post = nc_open('tas_ens_da_hydro_r.1000-2000_d.16-Feb-2018.nc')
nc.prior = nc_open('tas_prior_da_hydro_r.1000-2000_d.16-Feb-2018.nc')

cat("Data read \n")

# prep prior
prior.sub = read.csv("prior_ens.txt", header = F)
prior.sub = as.vector(prior.sub[,1])
prior.ens = prep_prior(nc.prior)

regions = 216

KS = matrix(0, length(times), regions)
KU = matrix(0, length(times), regions)

cat("setup complete \n")

prior_ens = prior.ens[,,prior.sub]
for (t in times) {

  cat(paste0("starting iter ", t, "\t"))

  post_ens = prep_post_ens(nc.post, t)

  lil_priors = vapply(1:100, function(x) matsplitter(prior_ens[,,x], 8, 8),
                      FUN.VALUE = array(0, dim = c(8, 8, regions)))

  lil_posts  = vapply(1:100, function(x) matsplitter(post_ens[,,x], 8, 8),
                      FUN.VALUE = array(0, dim = c(8, 8, regions)))

  # find the observed kst field
  kol.field = kst.field(lil_priors, lil_posts, 100)

  # find the permutation distribution
  perm.fields = kst.permute(lil_priors, lil_posts, 1000, 5)

  perm.ed = edepth_set(perm.fields, depth_function = "rank")
  perm.cr = central_region(perm.fields, perm.ed)

  KS[which(times==t), ] = kol.field
  KU[which(times==t), ] = perm.cr[[2]]
  
  # perm.upper = sapply(1:length(kol.field), function(r) quantile(perm.fields[r,], 0.95))
  
  # KS[which(times==t), ] = kol.field
  # KU[which(times==t), ] = perm.upper

  cat(paste0("finishing iter ", t, "\n"))
}

saveRDS(KS, file = paste0("values", batch, ".rds"))
saveRDS(KU, file = paste0("upper", batch, ".rds"))





