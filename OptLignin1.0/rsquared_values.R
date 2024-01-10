# R Script Copyright (c) 2023 Arjun Chakrawal

library("readxl")
library(modelsummary)
library(lmerTest)
library(lme4)
library(performance)

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

## modelsummary ------
modelsummary(models,fmt = 2,
             estimate  = "{estimate} ({std.error}){stars}",
             statistic = NULL,
             notes = list('Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1'),
             gof_omit = 'ICC|RMSE',
             output = "table.docx"
)
rsquared <-r2(models$`log(v_H)`)

rsquared$R2_conditional

# Function to calculate and extract R-squared values
get_rsquared <- function(model) {
  rsquared <- r2(model)
  return(c(marginal = rsquared$R2_marginal, conditional = rsquared$R2_conditional))
}

get_rsquared(models$`log(v_H)`)
# Apply the function to all models
rsquared_values <- sapply(models, get_rsquared)

# Convert the result to a data frame
rsquared_df <- as.data.frame(t(rsquared_values))

# Add model names as a column
rsquared_df$model <- rownames(rsquared_df)

# Write the result to a CSV file
write.csv(rsquared_df, "rsquared_values.csv", row.names = FALSE)

