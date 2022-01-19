# Difference between E(1/X) and 1/E(X) when X ~ Bin(n, p)

n = 2000
phi = 0.95
alpha = 8
p = phi^alpha / 2
x = 1:n
1 / sum(dbinom(x, n, p) * x)
sum(dbinom(x, n, p) * 1 / x)

n_grid = seq(500, 20500, 1000)
p_grid = seq(0.05, 0.5, 0.05)
n_n = length(n_grid)
n_p = length(p_grid)
diff_grid = matrix(NA, n_n, n_p)
for (i in 1:n_n) {
  for (j in 1:n_p) {
    n = n_grid[i]
    p = p_grid[j]
    x = 1:n
    diff_grid[i, j] = 
      abs(1 / sum(dbinom(x, n, p) * x) - sum(dbinom(x, n, p) * 1 / x))
  }
}

persp(
  n_grid, p_grid, diff_grid, 
  main = "Difference between E(1/X) and 1/E(X) when X ~ Bin(n, p)"
)
