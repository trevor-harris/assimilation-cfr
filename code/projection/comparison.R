
prior.gp = sim_gp(mu = 0, l = 1)
post.gp = sim_gp(mu = 0.3, l = 1)

# split em all
gp.prior.split = vapply(1:100, function(x) matsplitter(prior.gp[,,x], 10, 10),
                         FUN.VALUE = array(0, dim = c(10, 10, 9)))
gp.post.split = vapply(1:100, function(x) matsplitter(post.gp[,,x], 10, 10),
                        FUN.VALUE = array(0, dim = c(10, 10, 9)))

### NEW TEST
kst.val = kst.field.eig(gp.prior.split, gp.post.split, eigen = 20)
perm.fields = kst.permute.eig(gp.prior.split, gp.post.split, eigen = 20)

# find the central regions
perm.ed = edepth_set(perm.fields)
perm.cr = central_region(perm.fields, perm.ed)

plot(kst.val, type = "l", col = "red", ylim=c(0, 1))
for(i in 1:dim(perm.fields)[2]) {
  lines(perm.fields[,i])
}

sapply(1:length(kst.val), function(x) 1-isbetween(kst.val[x], perm.cr[[1]][x], perm.cr[[2]][x]))



### OLD TEST
# find the observed kst field
kol.field = kst.field(gp.prior.split, gp.post.split)

# find the permutation distribution
perm.fields = kst.permute(gp.prior.split, gp.post.split, 100, 100)

# find the central regions
perm.ed = edepth_set(perm.fields)
perm.cr = central_region(perm.fields, perm.ed)

plot(kol.field, type = "l", col = "red", ylim=c(0, 1))
for(i in 1:dim(perm.fields)[2]) {
  lines(perm.fields[,i])
}

sapply(1:length(kol.field), function(x) 1-isbetween(kol.field[x], perm.cr[[1]][x], perm.cr[[2]][x]))







### L2 (actuall l inf)
l2.val = l2.field(gp.prior.split, gp.post.split)
l2.perm = l2.permute(gp.prior.split, gp.post.split)

# find the central regions
perm.ed = edepth_set(l2.perm)
perm.cr = central_region(l2.perm, perm.ed)

plot(l2.val, type = "l", col = "red", ylim=c(0, 1))
for(i in 1:dim(l2.perm)[2]) {
  lines(l2.perm[,i])
}

sapply(1:length(l2.val), function(x) 1-isbetween(l2.val[x], perm.cr[[1]][x], perm.cr[[2]][x]))

