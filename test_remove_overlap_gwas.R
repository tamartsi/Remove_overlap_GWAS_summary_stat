

## let's test a bit:
beta1 <- c(1, 1, 2, 2, 4)
beta2 <- c(1, 1, 2, 4, -4)

v1 <- c(0.5, 0.5, 0.5, 1, 1)
v2 <- c(0.5, 0.1, 0.5, 1, 2)

w1 <- 1/v1
w2 <- 1/v2

beta_all <- (w1*beta1 + w2*beta2)/(w1 + w2)
beta_all # looks right!

v_all <- 1/(w1 + w2)
v_all # looks good!

# (we can use a meta-analysis to check if needed, this is fairly easy though)

se_all <- sqrt(v_all)
se1 <- sqrt(v1)
res <- remove_overlap_gwas(beta_all, se_all, beta1, se1)

all(abs(res$beta_keep - beta2) < 1e-10) # TRUE
all(abs(res$se_keep - sqrt(v2)) < 1e-10) # TRUE
