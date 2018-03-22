e = 1

# slice prior to member e
prior = prior.ens[,,e]

# prep posterior at member e
post = prep_post_time(nc.post, e)
# post = post[,,1:1001]

##### FIT BASIS AND FIND ED #####
prior.coef = proj %*% as.vector(prior)
post.coef = sapply(1:dim(post)[3], function(x) proj %*% as.vector(post[,,x]))

# calculate extremal depth
edepth(prior.coef, post.coef)

ed = edepth_set(post.coef)
central = central_region(post.coef, ed, 0.05)
median = post.coef[,ed = 1]


# create the prior/posterior comparison graph
coef.df = data.frame(ind = 1:896, post.coef)
coef.df = coef.df[coef.df$ind <= 15,]
coef.df = melt(coef.df, id.vars = "ind")

prior.df = data.frame(ind = 1:896, prior=prior.coef)
prior.df = prior.df[prior.df$ind <= 15,]

ggplot() +
  geom_line(data = coef.df, aes(x = ind, y = value, group = variable)) +
  geom_line(data = prior.df, aes(x = ind, y = prior, col = "Prior"),
            show.legend = FALSE) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Coefficient",
       y = "Value",
       title = "Prior vs Posterior")
ggsave(paste0("paper/figures/prior_post_coef", ".png"), width = 5, height = 3.2)

# create the prior/CR regions comparison graph
cr.df = data.frame(ind = 1:896, lower=central[[1]], upper=central[[2]])
cr.df = cr.df[cr.df$ind <= 15,]
cr.df = melt(cr.df, id.vars = "ind")

med.df = data.frame(ind = 1:896, median=median)
med.df = med.df[med.df$ind <= 15,]

ggplot() +
  geom_line(data = cr.df, aes(x = ind, y = value, group = variable)) +
  geom_line(data = med.df, aes(x = ind, y = median, color = "#F8766D"),
            show.legend = FALSE) +
  geom_line(data = prior.df, aes(x = ind, y = prior, color = "#00BFC4"),
            show.legend = FALSE) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Coefficient",
       y = "Value",
       title = "Prior vs Central Regions",
       fill = "sdfa")

ggsave(paste0("paper/figures/prior_cr_coef", ".png"), width = 5, height = 3.2)


# compare prior with ensemble of coeff
plot(post.coef[1:15, 1], type = "l", xlab = "Coefficient", ylab = "Value")
for (i in 2:1001) {
  lines(post.coef[1:15,i])
}
lines(prior.coef[1:15], col = "red")



# compare prior with median and central region
plot(prior.coef[1:15], type = "l", col = "red", xlab = "Coefficient", ylab = "Value")
lines(median[1:15], col = "blue")
lines(central[[2]][1:15])
lines(central[[1]][1:15])

