
#Eva Schemmel
#Von B growth curve
#Reference: Smart and Grammer 2021 Modernizing fish and shark growth curves with 
#           Bayesian length at age models


library(rstan)
stan_version()
options(mc.cores = parallel::detectCores()) 
rstan_options(auto_write = TRUE) 

library(loo)
library(rethinking)
library(dplyr)
library(boot)


##relationship between L0 and t0/A0
#L0 = Linf(1-exp(A0*K))
#A0 = log(1-L0/Linf)/k

#model
#Von B model : length ~ Linf-(Linf-L0)exp(-k*age)

##prior checks
prior_k<-runif(1e4,0.0,0.6)
prior_Linf<-rnorm(1e4,50,5) 

prior_L0<-rnorm(1e4,0,1)
prior_L0<-rlnorm(1e4,0,2) #set to only positive #(1e4,0,2) too limiting? Changed to 0,1
sample_sigma<-runif(1e4, 0, 2)

prior_Linf<-rnorm(4e3, prior_Linf, sample_sigma)
prior_k<-rnorm(4e3, prior_k, sample_sigma)
prior_L0<-rnorm(4e3, prior_L0, sample_sigma)

dens(prior_Linf)
dens(prior_k) 
dens(prior_L0, xlim=c(-5,15))
dens(sample_sigma)

#alternative L0 prior
prior_L01<-rlnorm(1e4,0,1) 
prior_L01<-rnorm(4e3, prior_L01, sample_sigma)

#posterior check check after you get your posterior distributions
plot(prior_L0, post$L0, ylim=c(0,20), xlim=c(-20,20), col="red",main="black=lnorm(0,1), red=lnorm(0,2)")
points(prior_L0, post$L0, col="black")
#checked L0 priors for lnorm(0,1) and lnorm(0,2) produce different growth curves lnorm(0,2) is less regulating and model fit is better

plot(prior_Linf, post$Linf, ylim=c(30,60), xlim=c(30,60), main="black=norm(45,5)")
#changing Linf did not change my model and results once I got a good L0 prior

write("// Stan model for basic vonB with L0 instead of t0

data {
 int < lower = 1 > N; // sample size
 vector[N] length; // response
 vector[N] age; // predictor
}

parameters {
  real alpha; // intercept
  real < lower = 0 > Linf; // Linf
  real < lower = 0 > k; // K
  real < lower = 0 > L0; // L0
  real < lower = 0, upper = 5 > sigma; // error
}

transformed parameters {
  vector[N] mu;
  mu = Linf-((Linf-L0)*exp(-k*age));
}

model {
  length ~ normal(mu, sigma);
  Linf ~ normal(50,5); 
  k ~ uniform(0.0,0.6);
  L0 ~ normal(0,2);//(0,2) works best here 
  sigma ~ uniform(0,5);
}

generated quantities {
  vector[N] log_lik;
  vector[N] yrep;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(length[i] | mu[i], sigma);
    yrep[i] = normal_rng(mu[i],sigma); //posterior for each prediction for each data point
  }
} 

",

"stan_vonBL0.stan")

stanc("stan_vonBL0.stan")



#data
valo_age<-read.csv("Age_VALO_March_2023.csv")
valo_age<-valo_age[complete.cases(valo_age$reage),]
valo_age<-valo_age[complete.cases(valo_age$Length.cm.),]
nrow(valo_age)

##glm von B for comparision

modvb<- nls(Length.cm. ~ Linf-(Linf-L0)*exp(-k*reage),
            data= valo_age, 
            start = list (Linf = 50,
                          k = 0.5,
                          L0=0))

summary(modvb)

#bootstrap glm model for mean and CIs
bmodvb<-bootCase(modvb, R=1000)
bmodvb<-as.data.frame(bmodvb)
nrow(bmodvb)
head(bmodvb)
new_glm<-as.data.frame(age_p)
pred_glm<-matrix(NA, nrow=nrow(bmodvb), ncol=nrow(new_glm)) #rows equal number of posterior samples , columns = a value for every data point in new data - prediction for the new data - predict posterior
for(i in 1:nrow(new_glm)){
  pred_glm[,i]<-bmodvb$Linf-(bmodvb$Linf - bmodvb$L0)*exp(-bmodvb$k*(new_glm$age_p[i]))
}


new_glm$mean<-NA
new_glm$up80<-NA
new_glm$down80<-NA
for(i in 1:nrow(new_glm)){
  new_glm$mean[i]<-mean(pred_glm[,i]) #mean of every column gives a new row in the new data frame
  new_glm$up80[i]<-quantile(pred_glm[,i], 0.10) #mean of every column gives a new row in the new data frame
  new_glm$down80[i]<-quantile(pred_glm[,i], 0.90)
}



#select data for stan model
standata <- list(N=nrow(valo_age),age=valo_age$reage,length=valo_age$Length.cm.)
#run stan vonBL0 model
vonbL0<-stan(file="stan_vonBL0.stan", data=standata)
print(vonbL0, pars = c('Linf','k','L0'))

post<-rstan::extract(vonbL0, pars=c('Linf','k','L0')) 

pairs(vonbL0, pars=c('Linf','k','L0'))

hist(post$Linf)
hist(post$k)
hist(post$L0) 

#compare prior and posterior
dens(post$L0, xlim=c(0,20), ylim=c(0,.6), col="blue")
dens(prior_L0, add=TRUE)
plot(post$L0, prior_L0, ylab="Simulated Data", xlab="Observed L0", pch=1, ylim=c(0,100))
points(post$L0, prior_L01, col="red")

plot(post$Linf, prior_Linf, ylab="Simulated Data", xlab="Observed Linf", pch=1)

plot(post$k, prior_k, ylab="Simulated Data", xlab="Observed k", pch=1, ylim=c(0,6))


log_lik_1 <- extract_log_lik(vonbL0, merge_chains = FALSE) #log lik needs to be in generated quantities
r_eff_1 <- relative_eff(exp(log_lik_1)) # A vector of relative effective sample sizes from the MCMC
loo_1 <- loo(log_lik_1, r_eff=r_eff_1, save_psis = T)
loo_1

# check for influential observations
plot(loo_1) #looks good

# save posterior predictions
yrep <- rstan::extract(vonbL0)["yrep"]

# posterior predictive check
#reference: Gabry et al Visualization in Bayesian workflow
#leave one out cross validation loo-pit - influential observations check - points of leverage and outliers
png("bayesplots_pcc_loo_pit_overlay_valo.png", width=1100, height=800, res=350)
par(mfrow=c(1,3),mar=c(2,2.5,1.5
                       ,1),mgp=c(2,.7,0))

bayesplot::ppc_loo_pit_overlay(valo_age$Length.cm., as.array(yrep)$yrep, lw = weights(loo_1$psis_object))  #this looks good
dev.off()

png("bayesplots_pcc_dens_overlay_valo.png", width=1100, height=800, res=350)
par(mfrow=c(1,3),mar=c(2,2.5,1.5
                       ,1),mgp=c(2,.7,0))

bayesplot::ppc_dens_overlay(valo_age$Length.cm.,as.array(yrep)$yrep) #this looks good
dev.off()


#model statistic to data statistic
png("bayesplots_pcc_stat_valo.png", width=1100, height=800, res=350)
par(mfrow=c(1,3),mar=c(2,2.5,1.5
                       ,1),mgp=c(2,.7,0))
bayesplot::ppc_stat(valo_age$Length.cm., yrep$yrep, "mean") #this looks good
dev.off()

png("bayesplots_pcc_stat_valo.png", width=1100, height=800, res=350)
par(mfrow=c(1,3),mar=c(2,2.5,1.5
                       ,1),mgp=c(2,.7,0))
bayesplot::ppc_stat(valo_age$Length.cm., yrep$yrep, "mean") #this looks good
dev.off()


#posterior predictive checks P-value fit and R2
png("traceplot_Linf.png", width=1600, height=700, res=250)
par(mfrow=c(1,3),mar=c(2,2.5,1.5
                       ,1),mgp=c(2,.7,0))
traceplot(vonbL0, pars=c('Linf'))
dev.off()
png("traceplot_k.png", width=1600, height=700, res=250)
par(mfrow=c(1,3),mar=c(2,2.5,1.5
                       ,1),mgp=c(2,.7,0))
traceplot(vonbL0, pars=c('k'))
dev.off()

png("traceplot_L0.png", width=1600, height=700, res=250)
par(mfrow=c(1,3),mar=c(2,2.5,1.5
                       ,1),mgp=c(2,.7,0))
traceplot(vonbL0, pars=c('L0'))
dev.off()

pairs(vonbL0, pars=c('Linf','k')) #red points are divergent transitions - dropped MCMC samples
#approximation of the posterior is not good enough if many divergent transitions


#Graph it vonB L0
age_p<-seq(0,15, by =1)

new_data<-as.data.frame(age_p)
pred_new<-matrix(NA, nrow=nrow(post$Linf), ncol=nrow(new_data)) #rows equal number of posterior samples , columns = a value for every data point in new data - prediction for the new data - predict posterior
for(i in 1:nrow(new_data)){
  pred_new[,i]<-post$Linf-(post$Linf - post$L0)*exp(-post$k*(new_data$age_p[i]))
}

new_data$mean<-NA
new_data$up80<-NA
new_data$down80<-NA
for(i in 1:nrow(new_data)){
  new_data$mean[i]<-mean(pred_new[,i]) #mean of every column gives a new row in the new data frame
  new_data$up80[i]<-quantile(pred_new[,i], 0.25) #mean of every column gives a new row in the new data frame
  new_data$down80[i]<-quantile(pred_new[,i], 0.975)
}

#plot sexes and stages on figure
valo_ar<-read.csv("~/Documents/Courses_trainings/ASU_Bayesian/Bayesian_ASU/valo_bayes/valo_ar.csv")


valo_ar$H_Sex<-as.factor(valo_ar$H_Sex)
valo_ar_f<-valo_ar %>%
  subset(H_Sex=="F")
valo_ar_m<-valo_ar %>%
  subset(H_Sex=="M" )
valo_ar_po<-valo_ar %>%
  subset(H_Sex=="M_PO" )
valo_ar_t<-valo_ar %>%
  subset(H_Sex=="T" )

#plot it
png(file="vbgf_valo_bayesL0_jan2022.png",width=2000,height=2000,res=300)
par(mfrow=c(1,1))
plot(valo_age$FinalAge, valo_age$Length.cm., pch = 19,bty="l", main="", col="grey", yaxt="n",yaxs="i",xaxs="i", xlim=c(0,15), ylim=c(0,55), xlab = 'Age (years)', ylab = 'Fork Length (cm)', las=1, tcl  = -0.3 )
plot(valo_ar_f$FinalAge, valo_ar_f$Length.cm., pch = 19,bty="l", main="", col="grey", yaxt="n",yaxs="i",xaxs="i", xlim=c(0,15), ylim=c(0,55), xlab = 'Age (years)', ylab = 'Fork Length (cm)', las=1, tcl  = -0.3 )
points(valo_ar_m$FinalAge, valo_ar_m$Length.cm., col="black",pch = 19)
points(valo_ar_po$FinalAge, valo_ar_po$Length.cm., col="black",pch = 19)
points(valo_ar_t$FinalAge, valo_ar_t$Length.cm., col="black",pch = 2)
#bayesian
points(new_data$age_p, new_data$mean, type="l", col='black', lty=1, lwd=2)
points(new_data$age_p, new_data$up80, type="l", col='black', lty=2)
points(new_data$age_p, new_data$down80, type="l", col='black', lty=2)
#glm
#points(new_glm$age_p, new_glm$mean, type="l", col='green', lty=1, lwd=2)
#points(new_glm$age_p, new_glm$up80, type="l", col='green', lty=2)
#points(new_glm$age_p, new_glm$down80, type="l", col='green', lty=2)
#lines(pred.data.b, col="black", lwd=2)
#lines(pred.data, col="green", lwd=2)
axis(2, at=c(seq(0,55,by=5)),labels=c(seq(0,55,by=5)), tcl  = -0.3 )
dev.off()




###estimated t0 from L0
Linf_mean<-mean(post$Linf)
k_mean<-mean(post$k)
L0_mean<-mean(post$L0)
log(1-L0_mean/Linf_mean)/k_mean #~-1 so pretty good, better than frequentist glm



