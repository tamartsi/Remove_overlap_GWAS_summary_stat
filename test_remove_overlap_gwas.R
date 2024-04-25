

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


## more sophisticated example: meta-analyze summary statistics from multiple 
## studies, then assess the "removal" of one of them. 
require(metafor)
asn <- read.csv("asn.csv")
eur <- read.csv("eur.csv")
afr <- read.csv("afr.csv")
his <- read.csv("his.csv")

res <- data.frame(SNP = asn$SNP, Est_all = NA, SE_all = NA, Est_asn_eur_afr = NA, SE_asn_eur_afr = NA)

for (i in 1:nrow(res)){
  cur_res <- rma(yi = c(asn$Est[i], eur$Est[i], afr$Est[i], his$Est[i]), 
             vi = c(asn$SE[i]^2, eur$SE[i]^2, afr$SE[i]^2, his$SE[i]^2), 
             method="FE")
  
  res[i, "Est_all"] <- cur_res$beta
  res[i, "SE_all"] <- cur_res$se
  
  cur_res <- rma(yi = c(asn$Est[i], eur$Est[i], afr$Est[i]), 
                 vi = c(asn$SE[i]^2, eur$SE[i]^2, afr$SE[i]^2), 
                 method="FE")
  
  res[i, "Est_asn_eur_afr"] <- cur_res$beta
  res[i, "SE_asn_eur_afr"] <- cur_res$se
  
}


# now add an estimate from the remove overlap function:
res_rm_his <- remove_overlap_gwas(res$Est_all, res$SE_all, his$Est, his$SE)


all(abs(res_rm_his$beta_keep - res$Est_asn_eur_afr) < 1e-10) # TRUE
all(abs(res_rm_his$se_keep - res$SE_asn_eur_afr) < 1e-10) # TRUE

