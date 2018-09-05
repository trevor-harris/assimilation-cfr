##### VERTICAL #####

# define data and import fuctions
source('code/setup.R')

# setup list of fields
z = prior.ens[,,-prior.sub]
x = prior.ens[,,prior.sub]
y = prep_post_ens(nc.post, 60)

# switch to basis functions
z.coef = sapply(1:dim(z)[3], function(e) proj %*% as.vector(z[,,e]))
x.coef = sapply(1:dim(x)[3], function(e) proj %*% as.vector(x[,,e]))
y.coef = sapply(1:dim(y)[3], function(e) proj %*% as.vector(y[,,e]))

# calculate ed for each wrt to z.coef
z.ed = edepth_set(z.coef)

# sort z.coef by zed
z.sort = z.coef[, order(z.ed)]

# fast sort x and y
x.ed = edepth_multi_fast(x.coef, z.sort, depth_set(z.sort))
y.ed = edepth_multi_fast(y.coef, z.sort, depth_set(z.sort))

x.ed
y.ed
# sum(x.ed == 0)
sum(y.ed == 0)
