
##### Interaction  simulations#####

### 5x2 cross-validation t-test four covariates 
library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(rms)
library(Runuran)
library(cmprsk)
library(survsim)
library(prodlim)
library(randomForestSRC)
library("timereg")
library(riskRegression)
library(pec)
library(reshape2)

simulate_data <- function(N, seed) {
  # create covariates
  x1 <-rnorm(N, mean=0, sd=1) 
  x2 <-rnorm(N, mean=0, sd=1)
  x3 <-rnorm(N, mean=0, sd=1)
  x4 <-rnorm(N, mean=0, sd=1)
  x5 <-rnorm(N, mean=0, sd=1)
  x6 <-rnorm(N, mean=0, sd=1)
  x7 <- rbinom(N,1,0.5)
  x8 <- rbinom(N,1,0.5)
  x9 <- rbinom(N,1,0.5)
  x10 <- rbinom(N,1,0.5)
  x11 <- rbinom(N,1,0.5)
  x12 <- rbinom(N,1,0.5)
  rate = 0.01
  b1 = log(2)
  betasI1<-c(-b1,b1,0,0,-b1,b1)
  betasI2<-c(0,0,-b1,b1,-b1,b1)
  # create censoring times (unrelated to survival time)
  cens <- rexp(N,0.01)
  # create the two hazards for the events 
  lambda1<-rate*exp(x7*betasI1[1]*I(x1>0)+x8* betasI1[2]*I(x2>0)+ 
                      x9*betasI1[3]*I(x3>0)+x10*betasI1[4]*I(x4>0)+
                      x11*betasI1[5]*I(x5>0)+ 
                      x12*betasI1[6]*I(x6>0))
  lambda2<-rate*exp(x7*betasI2[1]*I(x1>0)+ x8*betasI2[2]*I(x2>0)+ 
                      x9*betasI2[3]*I(x3>0)+
                      x10*betasI2[4]*I(x4>0)+ x11*betasI2[5]*I(x5>0)+
                      x12*betasI2[6]*I(x6>0)
          )
  # create the overall hazard 
  h <- lambda1+lambda2
  # Create survival time based on  Bender formula for times T ~ Exponential(baseline_hazard)
  dt <- -log(runif(N)) / h 
  # Create event indicator 
  fcause<-rbinom(N, size=1, prob=lambda1/h)
  cause<-ifelse(fcause==0,2,1)
  e <- c(dt<=cens)*cause
  dt <- pmin(dt, cens)
  # Create a dataframe of the simulated dataset
  dat1<-data.frame(dt,e,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12)
  k<-which(dat1$e==1)
  # Model predictive performance
  ##Forest
  Results = list()
  folds <- 2
  for (i  in 1:5) {
    idxsplit <- split(sample(nrow(dat1), replace=FALSE), 1:folds)
    IBSF1<- sapply(idxsplit, function(x) {
      ##Train and test data
      test_data <- dat1[x, ]
      train_data <- dat1[-x, ]
      timepoints <- sort(unique(train_data$dt))
      ###grow a forest
      FG.fit <- FGR(Hist(dt, e)~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12,data=train_data,cause=1)
      forestC       <-rfsrc(Surv(dt, e)~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12, data=train_data,
                            ntree = 500,nsplit = 2,mtry =4,cause=c(1,1))
      # F1<- as.formula(Hist(dt, e)~1)
      ## Compute IBS
      predictionerror <-crps(pec(list(FG.fit,forestC),
                                 formula = as.formula(Hist(dt, e)~1),data=test_data,cause=1, times = timepoints, 
                                 cens.model = "marginal", reference = FALSE))[]
      diff = predictionerror[1]-predictionerror[2]
      rbind(predictionerror,diff )
    })
    Results[[i]] = IBSF1
    rownames(Results[[i]])<-c("FG","RSF","Diff")
  }
  Results<-rbind(Results[[1]],Results[[2]],Results[[3]],Results[[4]],
                 Results[[5]])
  Newresults<-data.frame(Results)
  colnames(Newresults)<-c("fold1","fold2")
  score1 <- Newresults[c(seq(1,15,3)),]
  score2 <- Newresults[c(seq(2,15,3)),]
  FG_scores<- rbind(score1[,1],score1[,2])
  RSF_scores<-rbind(score2[,1],score2[,2])
  diff   <- Newresults[c(seq(0,15,3)),]
  d_1_1 = diff[1,1]
  # 5 x2 cv t-test
  diff_i_bar = (diff[,1] +diff[,2]) / 2
  diff_i_sqr = (diff[,1] - diff_i_bar) ** 2 + (diff[,2] - diff_i_bar) ** 2 
  diff_sqr =sum(diff_i_sqr)
  t_bar = d_1_1 / ((diff_sqr / 5) ** .5) 
  p_value_t = 2*pt(-abs(t_bar),df=5)
  # Combined 5 x2 cv F-test
  M = diff_sqr
  N1 = sum((diff[,1]**2 +diff[,2]**2) / 2)
  F_bar = N1/(2*M)
  p_value_f=pf(F_bar, 10, 5, lower.tail = FALSE)
  c(t_stat=t_bar,
    p_value_tstat= p_value_t,
    F_stat = F_bar, 
    p_value_fstat= p_value_f,
    mean_FG_sc=mean(FG_scores),
    mean_RSF_sc=mean(RSF_scores),
    max_FG_sc= max(FG_scores),
    max_RSF_sc= max(RSF_scores),
    SD_FG_Sc = sd(FG_scores),
    SD_RSF_Sc = sd(RSF_scores),
    samples = N,
    seed = seed)
}
simulate_trial <- function(n_sims, N, seed) {
  set.seed(seed)
  results <- replicate(n_sims, simulate_data(N, seed))
  data.frame(t(results))
}
build_design_matrix <- function(initial_seed, num_seeds, sample_sizes) {
  set.seed(initial_seed)
  seeds <- sample.int(100000, num_seeds)
  design <- expand.grid(
    N = sample_sizes,
    seed = seeds
  )
}
n_sims <- 100 # At each sample size we repeat n_sims times
sample_sizes <- c(200, 300,400, 500,2000,3000) # We are considering four sample sizes
# These are the coefficients for the covariates
setup_seed <- 42 # 
n_seeds <- 10# trying the experiment at several number of seeds
design <- build_design_matrix(setup_seed, n_seeds, sample_sizes)
results <- design %>%
  rowwise() %>%do(simulate_trial(n_sims, .$N,.$seed))

write.csv(results,file="Interaction_results_new.csv")


