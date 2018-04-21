size = sum(diff_null) / (200*9)
power = 1:10
for (j in 1:10) {
  power[j] = sum(diff_alt[,,j])/(200*9)
}

alt_mu = c(0, seq(0.1:1, by = 0.1))
plot(alt_mu, c(size, power), type = "b"
     , xlab="Alternative Mean", ylim = c(0, 1), main = "Simulated Power Analysis for Shifts in the Mean")
points(0, size, col = "red")