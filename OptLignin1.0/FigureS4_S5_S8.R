# R Script Copyright (c) 2023 Arjun Chakrawal
# This script is intended for creating # Figure S4,S5, S8
# and Tables S1, S2, S3.

library("readxl")
library("ggplot2")
library("GGally")
library(effects)
library(modelsummary)
library(lmerTest)
library(lme4)
library(sjPlot)
library(sjmisc)
library(MuMIn)
library(performance)
library(dplyr)
library(latex2exp)
library(corrplot)
library(boot)
library(glmmTMB)

# get data and apply transformations--------

my_data = read_excel("data_forRStudio.xlsx")

my_data$studyID = as.factor(my_data$studyID)
my_data$LitterType = as.factor(my_data$LitterType)
my_data$Ntreat = as.factor(my_data$Ntreat)

c= quantile(my_data$tau, .5)^2/quantile(my_data$tau, .75)
c= min(my_data$tau[my_data$tau>0])/2

my_data$log_tau_C = log(my_data$tau+c)
my_data$logrO = log(my_data$ro)
my_data$logvh= log(my_data$vhmax)
my_data$logavgvo = log(my_data$avg_vo)
my_data$logmaxvo = log(my_data$max_vo)

my_data$MATC=my_data$MATC+273
my_data$sTemp =  datawizard::standardize(my_data$MATC)
my_data$cMAP =  datawizard::standardize(my_data$MAPMm)
my_data$sCN0 =  datawizard::standardize((my_data$CN0))
my_data$sARC0 =  datawizard::standardize(my_data$ARC0)
my_data$cAN0 =  datawizard::standardize(my_data$AN0)
str(my_data)
summary(my_data)

# correlation and histogram -----------
col=c('sTemp','cMAP','sCN0','sARC0','cAN0','log_tau_C','logrO','logvh','logavgvo','logmaxvo')
M = cor(my_data[col])
colnames(M)= c('$MAT[0][S]','$MAP[0][S]','$CN[0][S]','$ARC[0][S]','$AN[0][S]','$log(tau+K)','$log(r[O])',
               '$log(v[H])','$log(bar(v)[O])','$log(max((v[O])))')
rownames(M)=colnames(M)
png("results/FigureS4.png", width = 12, height = 8, units = "in", res = 300)  # Adjust width and height as needed
corr_plot=corrplot(M, method = 'circle', type = 'lower',diag = FALSE,
                   addCoef.col = 'black',
                   tl.srt = 45,tl.col = 'black'
)
dev.off()
graphics.off()


par(mfrow = c(2,4))
hist(my_data$vhmax)
hist(my_data$avg_vo)
hist(my_data$max_vo)
hist(my_data$ro)
hist(log(my_data$vhmax))
hist(log(my_data$avg_vo))
hist(log(my_data$max_vo))
hist(log(my_data$ro))



# LMM fitting ------------
## lmer ------------
models <- list(
  "log(tau+K)" = lmer(log(tau+c)~ (sTemp +sCN0+sARC0)^2 + (1|studyID), data = my_data,REML = FALSE),
  "log(r_O)" = lmer(log(ro)~ (sTemp+sCN0+sARC0)^2 + (1|studyID),data = my_data,REML = FALSE),
  "log(v_H)"     = lmer(log(vhmax)~ (sTemp+sCN0+sARC0)^2 + (1|studyID),data = my_data,REML = FALSE),
  "log(avg. v_O)" = lmer(log(avg_vo)~ (sTemp+sCN0+sARC0)^2 + (1|studyID),data = my_data,REML = FALSE),
  "log(max v_O)" = lmer(log(max_vo)~ (sTemp+sCN0+sARC0)^2 + (1|studyID),data = my_data,REML = FALSE)
  # "log(voAtTau)" = lmer(log(voAtTau)~ (sTemp+sCN0+sARC0)^2 + (1|studyID),data = my_data,REML = FALSE)
)
## glmmTMB -------

glmm_mod<- glmmTMB(tau ~ (sTemp +sCN0+sARC0)^2 + (1|studyID),
             family = ziGamma(link = "log"),
             ziformula = ~ (sTemp +sCN0+sARC0)^2,
             dispformula = ~1,
             data = my_data)
modelsummary(glmm_mod,fmt = 2,
             estimate  = "{estimate} ({std.error}){stars}",
             statistic = NULL,
             notes = list('Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1'),
             gof_omit = 'ICC|RMSE'
             # output = "table.docx"
             )
png("results/FigureS8.png", width = 6, height = 7, units = "in", res = 300)  # Adjust width and height as needed
plot_model(glmm_mod, grid = TRUE, show.p=TRUE,show.values=TRUE)
dev.off()

plot_models(models, grid = TRUE, show.p=TRUE,show.values=TRUE)
plot_models(models$`log(r_O)`, grid = TRUE, show.p=TRUE,show.values=TRUE)
## modelsummary ------
modelsummary(models,fmt = 2,
             estimate  = "{estimate} ({std.error}){stars}",
             statistic = NULL,
             notes = list('Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1'),
             gof_omit = 'ICC|RMSE',
             output = "table.docx"
             )

bbmle::AICtab(models,base=TRUE,logLik=TRUE)
# check_model -------
check_model(models$`log(v_H)`)
check_model(models$`log(avg. v_O)`)
check_model(models$`log(r_O)`)
check_model(models$`log(tau+K)`)


# Bootstrapping for Temperature sensitivity Q10 -------------
## Non Parametric bootstrap ----------- 

getEstimate <- function(model){
  beta_T=fixef(model)[2]/sd(my_data$MATC)
  exp(10*beta_T)
}

getBoot_vh = function(data, k){
  fitNew = lmer(log(vhmax)~ (sTemp+sCN0+sARC0)^2 + (1|studyID),data = data[k,],REML = FALSE)
  return(getEstimate(fitNew))
}
b_vh= boot(data = my_data, statistic = getBoot_vh, R = 1000)
cat("Mean:", mean(b_vh$t),"Standard Deviation:", sd(b_vh$t))
# png("b_vh.png", width = 10*1.66, height = 8, units = "cm", res = 300)
plot(b_vh)
# dev.off()
boot.ci(b_vh,type=c("norm","basic","perc"))



getBoot_avg_vo = function(data, k){
  fitNew = lmer(log(avg_vo)~ (sTemp+sCN0+sARC0)^2 + (1|studyID),data = data[k,],REML = FALSE)
  return(getEstimate(fitNew))
}
b_avg_vo= boot(data = my_data, statistic = getBoot_avg_vo, R = 1000)
cat("Mean:", mean(b_avg_vo$t),"Standard Deviation:", sd(b_avg_vo$t))
# png("b_avg_vo.png", width = 10*1.66, height = 8, units = "cm", res = 300)
plot(b_avg_vo)
# dev.off()
boot.ci(b_avg_vo,type=c("norm","basic","perc"))


getBoot_max_vo = function(data, k){
  fitNew = lmer(log(max_vo)~ (sTemp+sCN0+sARC0)^2 + (1|studyID),data = data[k,],REML = FALSE)
  return(getEstimate(fitNew))
}
b_max_vo= boot(data = my_data, statistic = getBoot_max_vo, R = 1000)
cat("Mean:", mean(b_max_vo$t),"Standard Deviation:", sd(b_max_vo$t))

# png("b_max_vo.png", width = 10*1.66, height = 8, units = "cm", res = 300)
plot(b_max_vo)
# dev.off()
boot.ci(b_max_vo,type=c("norm","basic","perc"))


Q10 = data.frame('$log(v[H])' = b_vh$t, '$log(bar(v)[O])' = b_avg_vo$t, '$log(max((v[O])))' = b_max_vo$t)
png("results/FigureS5.png", width = 10*1.66, height = 10, units = "cm", res = 300)
par(mar = c(2, 4, 2, 2)) # c(bottom, left, top, right)
# par(oma = c(0, 0, 0, 0))
label <- c(expression(paste("log ", italic((v[H])))), expression(paste("log ", italic( (bar(v)[O]) ))),
           expression(paste("log ", italic(max((v[O]))))))
boxplot(Q10, names = label, ylab = expression(Q[10]), ylim = c(0, 2.5))
grid()
dev.off()


# Table S3 
# t.test(b_vh$t, b_avg_vo$t,mu= 0, paired = FALSE, var.equal = FALSE)
t.test(b_vh$t, b_avg_vo$t,alternative = "greater", paired = FALSE, var.equal = FALSE)

# t.test(b_vh$t, b_max_vo$t,mu= 0, paired = FALSE, var.equal = FALSE)
t.test( b_max_vo$t,b_vh$t,"greater", paired = FALSE, var.equal = FALSE)

# t.test(b_avg_vo$t, b_max_vo$t,mu= 0, paired = FALSE, var.equal = FALSE)
t.test( b_max_vo$t,b_avg_vo$t,"greater", paired = FALSE, var.equal = FALSE)




