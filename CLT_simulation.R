library(boot)
library(moments)
library(data.table)
set.seed(1234)

rzilnorm <- function(N, meanlog = 0, sdlog = 1, p = 0.5){
  rzilnorm <- rlnorm(N, meanlog, sdlog)
  rzilnorm[rbinom(N, 1, p)==1] = 0
  return(rzilnorm)
}

boot_mean <- function(data, indices) mean(data[indices])
boot_sd <- function(data, indices) sd(data[indices])
boot_skew <- function(data, indices) skewness(data[indices])

N <- 100000
N_rep <- 5000
zilnorm_data <- rzilnorm(N, meanlog=11.5, sdlog=1.2, p=0.95)
zilnorm_mean <- boot(zilnorm_data, statistic = boot_mean, R=N_rep)$t0
zilnorm_sd <- boot(zilnorm_data, statistic = boot_sd, R=N_rep)$t0
zilnorm_skew <- boot(zilnorm_data, statistic = boot_skew, R=N_rep)$t0

sprintf("mean: %.2f, standard deviation: %.2f, skewness: %.2f", zilnorm_mean, zilnorm_sd, zilnorm_skew)  

# distributions with the same variance of the zero inflated log-normal model

# skewness and kurtosis of normal, uniform, lognormal
# normal 0, 0
# uniform 0, -6/5
# lognormal (\exp(\sigma ^{2})+2){\sqrt {\exp(\sigma ^{2})-1}},  \exp(4\sigma ^{2})+2\exp(3\sigma ^{2})+3\exp(2\sigma ^{2})-6


sdlog <- sqrt(log(0.5*(1+sqrt(1+4*zilnorm_sd^2))))
lognormal_skew <- exp(sdlog^2+2)*(sqrt(exp(sdlog^2)-1))
# lognormal_kurtosis <- exp(4*sdlog^2)+2*exp(3*sdlog^2)+3*exp(2*sdlog^2)-6
  
sd(rnorm(N, mean = 0, sd = zilnorm_sd))
sd(runif(N, min = 0, max = sqrt(12)*zilnorm_sd))
sd(rlnorm(N, meanlog = 0, sdlog = sqrt(log(0.5*(1+sqrt(1+4*zilnorm_sd^2))))))

logzilnorm_data <- log1p(zilnorm_data)
logzilnorm_mean <- boot(logzilnorm_data, statistic = boot_mean, R=N_rep)$t0
logzilnorm_sd <- boot(logzilnorm_data, statistic = boot_sd, R=N_rep)$t0
logzilnorm_skew <- boot(logzilnorm_data, statistic = boot_skew, R=N_rep)$t0
sprintf("mean: %.2f, standard deviation: %.2f, skewness: %.2f", logzilnorm_mean, logzilnorm_sd, logzilnorm_skew)  


#############################################
shapiro_test_pval <- function(x){
  return(shapiro.test(x)$p.value)
}

n_trials <- 1000

samp_mean_list <- list()
sample_stats <- list()
true_means <- list(normal=0, uniform=0.5*sqrt(12)*zilnorm_sd, lognorm=exp(0+0.5*log(0.5*(1+sqrt(1+4*zilnorm_sd^2)))), zilnorm=zilnorm_mean, logzilnorm=logzilnorm_mean)

for(n_size in c(100, 1000, 10000)){
  for(i in 1:n_trials){
    samp_df <- data.frame(normal = rnorm(n_size, mean=0, sd=zilnorm_sd)
                          ,uniform = runif(n_size, min=0, max=sqrt(12)*zilnorm_sd)
                          ,lognorm = rlnorm(n_size, meanlog = 0, sdlog = sqrt(log(0.5*(1+sqrt(1+4*zilnorm_sd^2)))))
                          ,zilnorm = rzilnorm(n_size, meanlog=11.5, sdlog=1.2, p=0.95)
                          ,logzilnorm = log1p(rzilnorm(n_size, meanlog=11.5, sdlog=1.2, p=0.95)))
    
    samp_mean_list[[i]] <- apply(X = samp_df, MARGIN = 2, FUN = mean)
  }
  samp_mean_df <- data.frame(do.call('rbind', samp_mean_list))
  SE_df <- apply(samp_mean_df, MARGIN = 2, FUN = sd)
  pval_df <- apply(samp_mean_df, MARGIN=2, FUN = shapiro_test_pval)
  
  sample_stats[[paste0('N',n_size)]] <- list(means = samp_mean_df, SEs = SE_df, pvals = pval_df)
}


# plot histogram of sample means
for(distr in c('normal','uniform','lognorm','zilnorm','logzilnorm')){
  par(mfrow=c(1,3))
  for(n_size in c(100, 1000, 10000)){
    means <- sample_stats[[paste0('N',n_size)]][['means']][[distr]]
    hist(means, freq=FALSE, breaks=50, xlab='sample mean', ylab='frequency', main=paste(distr,'N =',n_size))
    abline(v = true_means[[distr]], col="red", lty=2, lwd=1)
  }
}


pval_df <- cbind(data.frame(sample_stats$N100$pvals), data.frame(sample_stats$N1000$pvals), data.frame(sample_stats$N10000$pvals))
colnames(pval_df)<-c('N100','N1000','N10000')


######################################
N = 1000000
p_t = 0.91
p_c = 0.9
meanlog_t = 9.05
meanlog_c = 9
sdlog_t = 1.45
sdlog_c = 1.4
df <- data.table(treatment=rep(c(1,0), each=N),value=c( rzilnorm(N, meanlog=meanlog_t, sdlog=sdlog_t, p=p_t), rzilnorm(N, meanlog=meanlog_c, sdlog=sdlog_c, p=p_c)) )

# true means
N_boot <- 5000
mean_t <- boot(df[treatment==1, value], statistic = boot_mean, R=N_boot)$t0
mean_c<- boot(df[treatment==0, value], statistic = boot_mean, R=N_boot)$t0


par(mfrow=c(1,2))
hist(df[treatment==1, value], breaks=100, xlab='value', ylab='frequency', main='treatment')
abline(v = mean_t, col="red", lty=2, lwd=1)
hist(df[treatment==1 & value>0, value],  breaks=100, xlab='value', ylab='frequency', main='treatment (value>0)')

par(mfrow=c(1,2))
hist(df[treatment==0, value],  breaks=100, xlab='value', ylab='frequency', main='control')
abline(v = mean_c, col="red", lty=2, lwd=1)
hist(df[treatment==0 & value>0, value], breaks=100, xlab='value', ylab='frequency', main='treatment (value>0)')


N_sim <- 1000
N_samp_size <- 80000
sim_pval_list <- vector("list", N_sim)
for(i in 1:N_sim){
  df_sub <- rbind(df[treatment==1,][sample(N, size=N_samp_size, replace=FALSE),], df[treatment==0,][sample(N, size=N_samp_size, replace=FALSE),])
  df_sub[, value_log := log1p(value)]
  pvals <- data.table(test=c('t_test','log_t_test','MWU_test')
                      , pval=c(t.test(value~treatment, data=df_sub)$p.value
                               , t.test(value_log~treatment, data=df_sub)$p.value
                               , wilcox.test(value~treatment, data=df_sub, paired=FALSE)$p.value))
  sim_pval_list[[i]] <- pvals
}
sim_pval_df <- rbindlist(sim_pval_list)

print(sim_pval_df[, .(sig005_prop=sum(pval<0.05)/.N, sig01_prop=sum(pval<0.1)/.N), by=test])


# A/A test
sim_pval_list <- vector("list", N_sim)
for(i in 1:N_sim){
  df_sub <- rbind(df[treatment==0,][sample(N, size=N_samp_size, replace=FALSE),][, treatment:=1], df[treatment==0,][sample(N, size=N_samp_size, replace=FALSE),])
  df_sub[, value_log := log1p(value)]
  pvals <- data.table(test=c('t_test','log_t_test','MWU_test')
                      , pval=c(t.test(value~treatment, data=df_sub)$p.value
                               , t.test(value_log~treatment, data=df_sub)$p.value
                               , wilcox.test(value~treatment, data=df_sub, paired=FALSE)$p.value))
  sim_pval_list[[i]] <- pvals
}
sim_pval_df <- rbindlist(sim_pval_list)


print(sim_pval_df[, .(sig005_prop=sum(pval<0.05)/.N, sig01_prop=sum(pval<0.1)/.N), by=test])






