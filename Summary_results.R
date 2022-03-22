
library(foreign)
library(xtable)
library(stargazer)
library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(rms)
library(survival)
library(Runuran)
library(cmprsk)
library(survsim)
library(prodlim)
library(randomForestSRC)
library("timereg")
library(riskRegression)
library(pec)
library(dplyr)
library(reshape2)


## Summary

### Linear Simulations#####
results_summary <-read.csv("first_set_of_results_summary .csv")
New<-data.frame(results_summary[,c(6,7,12)])
df <- melt(New, id="samples")
df$samples<- as.factor(df$samples)
p <- ggplot(df, aes(x = samples, y =value,fill = variable ))+geom_boxplot(outlier.colour = "red", # Outliers color
                                                                          alpha = 0.9)+stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = c("#E3256B", "#0641A5"),name = "Models", labels = c("FG", "RSF"))
p+
  labs(x = "Sample sizes")+ 
  labs(title = "Cross-validated integrated brier scores for the linear simulations")+
  labs(y = " Integrated brier scores")+
  theme(legend.position="top")+theme(text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),
                                     axis.line = element_line(colour = "black", # Theme customization
                                                              size = 1))
#Stats
results_summary$samples<- as.factor(results_summary$samples)
summary_statistics <-results_summary %>%
  group_by(samples) %>%
  summarise(mean_p_value_tstat = sum(p_value_tstat<= .05),
            mean_p_value_fstat = sum(p_value_fstat<= .05),
            mean_CV_IBS_FG = mean(mean_FG_sc), 
            mean_CV_IBS_RSF = mean(mean_RSF_sc),
            n = n())
View(summary_statistics)



xtable(summary_statistics, type = "latex", file = "filename2.tex")

## Type I error 
results_summary <-read.csv("Type_one_error_simulations.csv")
data1<-data.frame(results_summary[,c(3,5,12)])
data1 <- melt(data1, id="samples")
data1$samples<- as.numeric(as.vector(data1$samples))

ggplot(data1, aes(x=value, color=variable)) +
  geom_histogram( bins=25)+facet_wrap(samples ~ .)+
  scale_x_continuous(limits=c(0, 1), breaks=seq(0,1, 0.2))+
  labs(y = " Frequency",title = "Observed Type I error rates with 1000 simulations",x = "Observed Type I error rates")+
  labs(color='p-values') +scale_colour_manual(labels = c("5 x 2-fold cv  paired t-test", "Combined 5 x 2-fold cv F-test"), 
                                              values = c("#a51606", "#06a5a2"))


ggplot(data=data1, aes(x=samples, y=value, group=variable)) +
  stat_summary(aes(color=variable), geom="line")


results_summary$samples<- as.factor(results_summary$samples)
summary_statistics <-results_summary %>%
  group_by(samples) %>%
  summarise(FG = sum(mean_FG_sc)/n(),
            RSF = sum(mean_RSF_sc)/n())
View(summary_statistics)



data2 <- melt(summary_statistics, id="samples")
xtable(summary_statistics, type = "latex", file = "filename2.tex")

## Line graphs 

ggplot(data2, aes(x=samples, y=value,colour = variable,group=variable)) +geom_point(aes(color= variable)) +labs(x = "Sample size")+
  labs(title = " Mean integrated Brier scores for 1000 simulations for the FG and RSF model")+
  labs(y = "Crossvalidated integrated brier scores")+geom_path()

## The distribution F and T statistics

# T stat
results_summary <-read.csv("Type_one_error_simulations.csv")
data1<-data.frame(results_summary[,c(2,12)])
data1 <- melt(data1, id="samples")
data1$samples<- as.numeric(as.vector(data1$samples))

ggplot(data1, aes(x=value, color=variable)) +
  geom_histogram( bins=25)+facet_wrap(samples ~ .)+labs(x = "t- Statistics values")+
  labs(title = "Observed t- Statistics under the null hypothesis for 1000 simulations")+
  labs(y = " Frequency")+
  scale_colour_manual(values="red", name="t - Statistics under the null", labels ="t - Distribution")

# F stat

results_summary <-read.csv("Type_one_error_simulations.csv")
data1<-data.frame(results_summary[,c(4,12)])
data1 <- melt(data1, id="samples")
data1$samples<- as.numeric(as.vector(data1$samples))

ggplot(data1, aes(x=value, color=variable)) +
  geom_histogram( bins=25)+facet_wrap(samples ~ .)+labs(x = "F-Statistics values")+
  labs(title = "Observed F-statistics under the null hypothesis for 1000 simulations")+labs(y = " Frequency")+
  scale_colour_manual(values="green", name="F - Statistics under the null", labels ="F - Distribution")


# Bar graph
library(devtools)
library(easyGgplot2)

ggplot2.barplot(data=data2, xName='samples', yName="value",
                groupName='variable', position=position_dodge(),
                groupColors=c("#E3256B", "#0641A5"),xtitle="Sample sizes", 
                ytitle="Integrated Brier scores", 
                mainTitle="Crossvalidated integrated Brier scores for 1000 simulations for the FG and RSF model under the null hypothesis",legendTitle="Models")





