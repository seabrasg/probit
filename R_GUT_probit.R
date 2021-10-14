###### PROBIT ANALYSIS
# pipeline adapted from https://stats.idre.ucla.edu/r/dae/probit-regression/

# install these packages first 
require(dplyr)
require(ggplot2)
require(aod)
require(lmtest)

setwd("/Users/seabrasg/Dropbox/0_Coenosia_GUT/Submission_Github")
m=read.table("./matrix_GUT_detection_time.csv",header=TRUE,sep=",")
summary(m)


#####################################
# Species with 2 primer sets (Bim, Lhu, Dmrc, Dis)
######################################

for(i in c("Bim","Lhu","Dmrc","Dis")) {

prey<-m[which(m$Prey==i),]

###### PROPORTION OF POSITIVES ON EACH TIME - FOR EACH PREY

# FOR PCR1 
prey_PCR1<-prey[which(prey$PCR=="PCR1"),]
table_prey<-table(prey[which(prey$PCR=="PCR1"),c(4,5)])
table_proportion<-prop.table(table_prey,2)
df_proportion<-as.data.frame(table_proportion)
table_positive_PCR1<-df_proportion[which(df_proportion$Detection=="1"),]
time_variable_PCR1<-as.numeric(as.numeric(levels(table_positive_PCR1$Time)))
freq_variable_PCR1<-table_positive_PCR1$Freq
samplesize_PCR1<-colSums(table_prey)

total_positives_PCR1<-as.vector(table_prey[2,])

# FOR PCR2 
prey_PCR2<-prey[which(prey$PCR=="PCR2"),]
table_prey<-table(prey[which(prey$PCR=="PCR2"),c(4,5)])
table_proportion<-prop.table(table_prey,2)
df_proportion<-as.data.frame(table_proportion)
table_positive_PCR2<-df_proportion[which(df_proportion$Detection=="1"),]
time_variable_PCR2<-as.numeric(as.numeric(levels(table_positive_PCR2$Time)))
freq_variable_PCR2<-table_positive_PCR2$Freq
samplesize_PCR2<-colSums(table_prey)

total_positives_PCR2<-as.vector(table_prey[2,])


#### PROBIT ANALYSIS
sink(file="./output_probit_R.txt", append=TRUE, type="output")

print("PROBIT ANALYSIS")
print(i)

table_prey_total_positives<-data.frame(Time=time_variable_PCR1, 
                            N1=samplesize_PCR1,PCR1=total_positives_PCR1,
                            N2=samplesize_PCR2,PCR2=total_positives_PCR2)
print(table_prey_total_positives)

table_prey_prop<-data.frame(Time=time_variable_PCR1, 
                       N1=samplesize_PCR1,PCR1=freq_variable_PCR1,
                       N2=samplesize_PCR2,PCR2=freq_variable_PCR2)
print(table_prey_prop)


# FOR PCR1 ########################################################
print("PCR1")
myprobit_PCR1 <- glm(Detection ~ Time, family = binomial(link = "probit"), 
                data = prey_PCR1)
print(summary(myprobit_PCR1))
print(confint(myprobit_PCR1)) 

# Test whether the model with predictors fits significantly better than a model with just an intercept (i.e. a null model). 
# The test statistic is the difference between the residual deviance for the model with predictors and the null model.

het=myprobit_PCR1$deviance/myprobit_PCR1$df.residual #heterogeneity factor (chisquare divided by degrees of freedom)
changedev<-with(myprobit_PCR1, null.deviance - deviance)
changedf<-with(myprobit_PCR1, df.null - df.residual)
pquisquare<-with(myprobit_PCR1, pchisq(null.deviance - deviance, df.null - df.residual,lower.tail = FALSE))
LL<-logLik(myprobit_PCR1)

stat_results<-data.frame(HeterogeneityFactor=het, ChangeDeviance=changedev,ChangeDF=changedf,
                         pChiSquare=pquisquare,LogLikelihood=LL)
print(stat_results)

## PREDICT
newdata_PCR1<- data.frame(Time=seq(0,48,0.25))
newdata_PCR1[, c("p", "se")] <- predict(myprobit_PCR1, newdata_PCR1, type = "response", se.fit = TRUE)[-3]


# FOR PCR2 ########################################################
print("PCR2")
myprobit_PCR2 <- glm(Detection ~ Time, family = binomial(link = "probit"), 
                     data = prey_PCR2)
print(summary(myprobit_PCR2))
print(confint(myprobit_PCR2)) 

with(myprobit_PCR2, pchisq(deviance, df.residual))

# Test whether the model with predictors fits significantly better than a model with just an intercept (i.e. a null model). 
# The test statistic is the difference between the residual deviance for the model with predictors and the null model.
het=myprobit_PCR2$deviance/myprobit_PCR2$df.residual #heterogeneity factor (chisquare divided by degrees of freedom)
changedev<-with(myprobit_PCR2, null.deviance - deviance)
changedf<-with(myprobit_PCR2, df.null - df.residual)
pquisquare<-with(myprobit_PCR2, pchisq(null.deviance - deviance, df.null - df.residual,lower.tail = FALSE))
LL<-logLik(myprobit_PCR2)

stat_results<-data.frame(HeterogeneityFactor=het, ChangeDeviance=changedev,ChangeDF=changedf,
                         pChiSquare=pquisquare,LogLikelihood=LL)
print(stat_results)

## PREDICT
newdata_PCR2<- data.frame(Time=seq(0,48,0.25))
newdata_PCR2[, c("p", "se")] <- predict(myprobit_PCR2, newdata_PCR2, type = "response", se.fit = TRUE)[-3]



print("Fit separate lines") ############################

mod1=glm(Detection ~ PCR + PCR:Time, family = binomial(link = "probit"), 
         data = prey)
print(summary(mod1))

het=mod1$deviance/mod1$df.residual #heterogeneity factor

print("Heterogeneity factor")
print(het)

print("Test hypothesis of parallel lines") ############################

mod2=glm(Detection ~ PCR + Time, family = binomial(link = "probit"), 
         data = prey)
print(summary(mod2))

df1=mod1$df.residual
df2=mod2$df.residual
tstat=(mod2$deviance-mod1$deviance)/(het*(df2-df1))
testp=1-pf(tstat,df2-df1,df1)
print(paste("test deviance=",tstat,", ","Pvalue=",testp))

## PLOT VALUES AND MODEL PREDICTION #################################
name_prey<-i
namefile<-paste0("./plot_",i,".jpeg")
jpeg(file = namefile, bg="white", antialias = "default",width = 6, height = 5,  
     units = "in", res = 600)

plot(c(0,48), c(0,1), type="n", xlab="Time since feeding (hours)", ylab="Probability of detection by PCR", main=name_prey)
ld<-seq(0,48,0.25)
lines(ld, newdata_PCR1$p, lty=1)
points(time_variable_PCR1,freq_variable_PCR1, pch=20)

lines(ld, newdata_PCR2$p, lty=2)
points(time_variable_PCR2,freq_variable_PCR2, pch=4)

dev.off()

sink()

}

#################################
# Species with 3 primer sets (Tva)
################################

for(i in c("Tva")) {
  
  prey<-m[which(m$Prey==i),]
  
  ###### PROPORTION OF POSITIVES ON EACH TIME - FOR EACH PREY
  
  # FOR PCR1 
  prey_PCR1<-prey[which(prey$PCR=="PCR1"),]
  table_prey<-table(prey[which(prey$PCR=="PCR1"),c(4,5)])
  table_proportion<-prop.table(table_prey,2)
  df_proportion<-as.data.frame(table_proportion)
  table_positive_PCR1<-df_proportion[which(df_proportion$Detection=="1"),]
  time_variable_PCR1<-as.numeric(as.numeric(levels(table_positive_PCR1$Time)))
  freq_variable_PCR1<-table_positive_PCR1$Freq
  samplesize_PCR1<-colSums(table_prey)
  
  total_positives_PCR1<-as.vector(table_prey[2,])
  
  # FOR PCR2 
  prey_PCR2<-prey[which(prey$PCR=="PCR2"),]
  table_prey<-table(prey[which(prey$PCR=="PCR2"),c(4,5)])
  table_proportion<-prop.table(table_prey,2)
  df_proportion<-as.data.frame(table_proportion)
  table_positive_PCR2<-df_proportion[which(df_proportion$Detection=="1"),]
  time_variable_PCR2<-as.numeric(as.numeric(levels(table_positive_PCR2$Time)))
  freq_variable_PCR2<-table_positive_PCR2$Freq
  samplesize_PCR2<-colSums(table_prey)
  
  total_positives_PCR2<-as.vector(table_prey[2,])
  
  # FOR PCR3 
  prey_PCR3<-prey[which(prey$PCR=="PCR3"),]
  table_prey<-table(prey[which(prey$PCR=="PCR3"),c(4,5)])
  table_proportion<-prop.table(table_prey,2)
  df_proportion<-as.data.frame(table_proportion)
  table_positive_PCR3<-df_proportion[which(df_proportion$Detection=="1"),]
  time_variable_PCR3<-as.numeric(as.numeric(levels(table_positive_PCR3$Time)))
  freq_variable_PCR3<-table_positive_PCR3$Freq
  samplesize_PCR3<-colSums(table_prey)
  
  total_positives_PCR3<-as.vector(table_prey[2,])
  
  
  #### PROBIT ANALYSIS
  sink(file="./output_probit_R_Tva.txt", append=TRUE, type="output")
  
  print("")
  print("PROBIT ANALYSIS")
  print(i)
  
  table_prey_total_positives<-data.frame(Time=time_variable_PCR1, 
                                         N1=samplesize_PCR1,PCR1=total_positives_PCR1,
                                         N2=samplesize_PCR2,PCR2=total_positives_PCR2,
                                         N3=samplesize_PCR3,PCR3=total_positives_PCR3)
  print(table_prey_total_positives)
  
  table_prey_prop<-data.frame(Time=time_variable_PCR1, 
                              N1=samplesize_PCR1,PCR1=freq_variable_PCR1,
                              N2=samplesize_PCR2,PCR2=freq_variable_PCR2,
                              N3=samplesize_PCR3,PCR3=freq_variable_PCR3)
  print(table_prey_prop)
  
  
  # FOR PCR1 ########################################################
  print("PCR1")
  myprobit_PCR1 <- glm(Detection ~ Time, family = binomial(link = "probit"), 
                       data = prey_PCR1)
  print(summary(myprobit_PCR1))
  print(confint(myprobit_PCR1)) 
  
  
  # Test whether the model with predictors fits significantly better than a model with just an intercept (i.e. a null model). 
  # The test statistic is the difference between the residual deviance for the model with predictors and the null model.
  het=myprobit_PCR1$deviance/myprobit_PCR1$df.residual #heterogeneity factor (chisquare divided by degrees of freedom)
  changedev<-with(myprobit_PCR1, null.deviance - deviance)
  changedf<-with(myprobit_PCR1, df.null - df.residual)
  pquisquare<-with(myprobit_PCR1, pchisq(null.deviance - deviance, df.null - df.residual,lower.tail = FALSE))
  LL<-logLik(myprobit_PCR1)
  
  stat_results<-data.frame(HeterogeneityFactor=het, ChangeDeviance=changedev,ChangeDF=changedf,
                           pChiSquare=pquisquare,LogLikelihood=LL)
  print(stat_results)
  
  ## PREDICT
  newdata_PCR1<- data.frame(Time=seq(0,48,0.25))
  newdata_PCR1[, c("p", "se")] <- predict(myprobit_PCR1, newdata_PCR1, type = "response", se.fit = TRUE)[-3]
  
  
  # FOR PCR2 ########################################################
  print("PCR2")
  myprobit_PCR2 <- glm(Detection ~ Time, family = binomial(link = "probit"), 
                       data = prey_PCR2)
  print(summary(myprobit_PCR2))
  print(confint(myprobit_PCR2)) 
  
  with(myprobit_PCR2, pchisq(deviance, df.residual))
  
  # Test whether the model with predictors fits significantly better than a model with just an intercept (i.e. a null model). 
  # The test statistic is the difference between the residual deviance for the model with predictors and the null model.
  het=myprobit_PCR2$deviance/myprobit_PCR2$df.residual #heterogeneity factor (chisquare divided by degrees of freedom)
  changedev<-with(myprobit_PCR2, null.deviance - deviance)
  changedf<-with(myprobit_PCR2, df.null - df.residual)
  pquisquare<-with(myprobit_PCR2, pchisq(null.deviance - deviance, df.null - df.residual,lower.tail = FALSE))
  LL<-logLik(myprobit_PCR2)
  
  stat_results<-data.frame(HeterogeneityFactor=het, ChangeDeviance=changedev,ChangeDF=changedf,
                           pChiSquare=pquisquare,LogLikelihood=LL)
  print(stat_results)
  
  ## PREDICT
  newdata_PCR2<- data.frame(Time=seq(0,48,0.25))
  newdata_PCR2[, c("p", "se")] <- predict(myprobit_PCR2, newdata_PCR2, type = "response", se.fit = TRUE)[-3]
  
  # FOR PCR3 ########################################################
  print("PCR3")
  myprobit_PCR3 <- glm(Detection ~ Time, family = binomial(link = "probit"), 
                       data = prey_PCR3)
  print(summary(myprobit_PCR3))
  print(confint(myprobit_PCR3)) 
  
  # Test whether the model with predictors fits significantly better than a model with just an intercept (i.e. a null model). 
  # The test statistic is the difference between the residual deviance for the model with predictors and the null model.
  het=myprobit_PCR3$deviance/myprobit_PCR3$df.residual #heterogeneity factor (chisquare divided by degrees of freedom)
  changedev<-with(myprobit_PCR3, null.deviance - deviance)
  changedf<-with(myprobit_PCR3, df.null - df.residual)
  pquisquare<-with(myprobit_PCR3, pchisq(null.deviance - deviance, df.null - df.residual,lower.tail = FALSE))
  LL<-logLik(myprobit_PCR3)
  
  stat_results<-data.frame(HeterogeneityFactor=het, ChangeDeviance=changedev,ChangeDF=changedf,
                           pChiSquare=pquisquare,LogLikelihood=LL)
  print(stat_results)
  
  ## PREDICT
  newdata_PCR3<- data.frame(Time=seq(0,48,0.25))
  newdata_PCR3[, c("p", "se")] <- predict(myprobit_PCR3, newdata_PCR3, type = "response", se.fit = TRUE)[-3]
  
  
  
  
  print("Fit separate lines") ############################
  
  mod1=glm(Detection ~ PCR + PCR:Time, family = binomial(link = "probit"), 
           data = prey)
  print(summary(mod1))
  
  het=mod1$deviance/mod1$df.residual #heterogeneity factor
  
  print("Heterogeneity factor")
  print(het)
  
  print("Test hypothesis of parallel lines") ############################
  
  mod2=glm(Detection ~ PCR + Time, family = binomial(link = "probit"), 
           data = prey)
  print(summary(mod2))
  
  df1=mod1$df.residual
  df2=mod2$df.residual
  tstat=(mod2$deviance-mod1$deviance)/(het*(df2-df1))
  testp=1-pf(tstat,df2-df1,df1)
  print(paste("test deviance=",tstat,", ","Pvalue=",testp))
  
  ## PLOT VALUES AND MODEL PREDICTION #################################
  name_prey<-i
  namefile<-paste0("./plot_",i,".jpeg")
  jpeg(file = namefile, bg="white", antialias = "default",width = 6, height = 5,  
       units = "in", res = 600)
  
  plot(c(0,48), c(0,1), type="n", xlab="Time after ingestion (hours)", ylab="Probability of detection by PCR", main=name_prey)
  ld<-seq(0,48,0.25)
  lines(ld, newdata_PCR1$p, lty=1)
  points(time_variable_PCR1,freq_variable_PCR1, pch=20)
  #text(time_variable_PCR1,freq_variable_PCR1, labels=samplesize_PCR1, cex= 0.7, pos=2)
  
  lines(ld, newdata_PCR2$p, lty=2)
  points(time_variable_PCR2,freq_variable_PCR2, pch=4)
  #text(time_variable_PCR2,freq_variable_PCR2, labels=samplesize_PCR2, cex= 0.7, pos=2)
  
  lines(ld, newdata_PCR3$p, lty=3)
  points(time_variable_PCR3,freq_variable_PCR3, pch=5)
  #text(time_variable_PCR2,freq_variable_PCR2, labels=samplesize_PCR2, cex= 0.7, pos=2)
  
  dev.off()
  
  sink()
  
}


##################################################################
# PLOT ALL IN ONE FIGURE
##################################################################
# FUNCTION for 2 PCR
f.plot_2PCR<-function(i){
  prey<-m[which(m$Prey==i),] 
  
  ###### PROPORTION OF POSITIVES ON EACH TIME - FOR EACH PREY
  
  # FOR PCR1 
  prey_PCR1<-prey[which(prey$PCR=="PCR1"),]
  table_prey<-table(prey[which(prey$PCR=="PCR1"),c(4,5)])
  table_proportion<-prop.table(table_prey,2)
  df_proportion<-as.data.frame(table_proportion)
  table_positive_PCR1<-df_proportion[which(df_proportion$Detection=="1"),]
  time_variable_PCR1<-as.numeric(as.numeric(levels(table_positive_PCR1$Time)))
  freq_variable_PCR1<-table_positive_PCR1$Freq
  samplesize_PCR1<-colSums(table_prey)
  
  total_positives_PCR1<-as.vector(table_prey[2,])
  
  # FOR PCR2 
  prey_PCR2<-prey[which(prey$PCR=="PCR2"),]
  table_prey<-table(prey[which(prey$PCR=="PCR2"),c(4,5)])
  table_proportion<-prop.table(table_prey,2)
  df_proportion<-as.data.frame(table_proportion)
  table_positive_PCR2<-df_proportion[which(df_proportion$Detection=="1"),]
  time_variable_PCR2<-as.numeric(as.numeric(levels(table_positive_PCR2$Time)))
  freq_variable_PCR2<-table_positive_PCR2$Freq
  samplesize_PCR2<-colSums(table_prey)
  
  total_positives_PCR2<-as.vector(table_prey[2,])
  
  #### PROBIT ANALYSIS
  
  table_prey_total_positives<-data.frame(Time=time_variable_PCR1, 
                                         N1=samplesize_PCR1,PCR1=total_positives_PCR1,
                                         N2=samplesize_PCR2,PCR2=total_positives_PCR2)
  
  table_prey_prop<-data.frame(Time=time_variable_PCR1, 
                              N1=samplesize_PCR1,PCR1=freq_variable_PCR1,
                              N2=samplesize_PCR2,PCR2=freq_variable_PCR2)
  
  
  # FOR PCR1 ########################################################
  myprobit_PCR1 <- glm(Detection ~ Time, family = binomial(link = "probit"), 
                       data = prey_PCR1)
  
  # Test whether the model with predictors fits significantly better than a model with just an intercept (i.e. a null model). 
  # The test statistic is the difference between the residual deviance for the model with predictors and the null model.
  changedev<-with(myprobit_PCR1, null.deviance - deviance)
  df<-with(myprobit_PCR1, df.null - df.residual)
  quisquare<-with(myprobit_PCR1, pchisq(null.deviance - deviance, df.null - df.residual,lower.tail = FALSE))
  LL<-logLik(myprobit_PCR1)
  
  stat_results<-data.frame(ChangeDeviance=changedev,ChangeDF=df,
                           pChiSquare=quisquare,LogLikelihood=LL)
  
  ## PREDICT
  newdata_PCR1<- data.frame(Time=seq(0,48,0.25))
  newdata_PCR1[, c("p", "se")] <- predict(myprobit_PCR1, newdata_PCR1, type = "response", se.fit = TRUE)[-3]
  
  # Confidence intervals 
  ##########################
  ## grad the inverse link function
  ilink <- family(myprobit_PCR1)$linkinv
  ## add fit and se.fit on the **link** scale
  newdata_PCR1 <- bind_cols(newdata_PCR1, setNames(as_tibble(predict(myprobit_PCR1, newdata_PCR1, se.fit = TRUE)[1:2]),
                                                   c('fit_link','se_link')))
  ## create the interval and backtransform
  newdata_PCR1 <- mutate(newdata_PCR1,
                         fit_resp  = ilink(fit_link),
                         right_upr = ilink(fit_link + (2 * se_link)),
                         right_lwr = ilink(fit_link - (2 * se_link)))
  
  # FOR PCR2 ########################################################
  myprobit_PCR2 <- glm(Detection ~ Time, family = binomial(link = "probit"), 
                       data = prey_PCR2)
  with(myprobit_PCR2, pchisq(deviance, df.residual))
  
  # Test whether the model with predictors fits significantly better than a model with just an intercept (i.e. a null model). 
  # The test statistic is the difference between the residual deviance for the model with predictors and the null model.
  changedev<-with(myprobit_PCR2, null.deviance - deviance)
  df<-with(myprobit_PCR2, df.null - df.residual)
  quisquare<-with(myprobit_PCR2, pchisq(null.deviance - deviance, df.null - df.residual,lower.tail = FALSE))
  LL<-logLik(myprobit_PCR2)
  
  stat_results<-data.frame(ChangeDeviance=changedev,ChangeDF=df,
                           pChiSquare=quisquare,LogLikelihood=LL)
  ## PREDICT
  newdata_PCR2<- data.frame(Time=seq(0,48,0.25))
  newdata_PCR2[, c("p", "se")] <- predict(myprobit_PCR2, newdata_PCR2, type = "response", se.fit = TRUE)[-3]
  
  # Confidence intervals 
  ##########################33
  ## grad the inverse link function
  ilink <- family(myprobit_PCR2)$linkinv
  ## add fit and se.fit on the **link** scale
  newdata_PCR2 <- bind_cols(newdata_PCR2, setNames(as_tibble(predict(myprobit_PCR2, newdata_PCR2, se.fit = TRUE)[1:2]),
                                                   c('fit_link','se_link')))
  ## create the interval and backtransform
  newdata_PCR2 <- mutate(newdata_PCR2,
                         fit_resp  = ilink(fit_link),
                         right_upr = ilink(fit_link + (2 * se_link)),
                         right_lwr = ilink(fit_link - (2 * se_link)))
  
  
  # MODELs
  mod1=glm(Detection ~ PCR + PCR:Time, family = binomial(link = "probit"), 
           data = prey)
  
  het=mod1$deviance/mod1$df.residual #heterogeneity factor
  
  mod2=glm(Detection ~ PCR + Time, family = binomial(link = "probit"), 
           data = prey)
  
  df1=mod1$df.residual
  df2=mod2$df.residual
  tstat=(mod2$deviance-mod1$deviance)/(het*(df2-df1))
  testp=1-pf(tstat,df2-df1,df1)
  
  ## PLOT VALUES AND MODEL PREDICTION #################################
  name_prey=i
  plot(c(0,48), c(0,1), type="n", xlab="Time after ingestion (hours)", ylab="Probability of detection", main=name_prey)
  ld<-seq(0,48,0.25)
  lines(ld, newdata_PCR1$p, lty=1,col= "#E69F00",lwd=2)
  points(time_variable_PCR1,freq_variable_PCR1, pch=20,col= "#E69F00")
  #text(time_variable_PCR1,freq_variable_PCR1, labels=samplesize_PCR1, cex= 0.7, pos=2)
  lines(ld, newdata_PCR1$right_upr, lty=2,col= "#E69F00",lwd=2)
  lines(ld, newdata_PCR1$right_lwr, lty=2,col= "#E69F00",lwd=2)
  
  lines(ld, newdata_PCR2$p, lty=1,col= "#0072B2",lwd=2)
  points(time_variable_PCR2,freq_variable_PCR2, pch=4,col= "#0072B2")
  #text(time_variable_PCR2,freq_variable_PCR2, labels=samplesize_PCR2, cex= 0.7, pos=2)
  lines(ld, newdata_PCR2$right_upr, lty=2,col= "#0072B2",lwd=2)
  lines(ld, newdata_PCR2$right_lwr, lty=2,col= "#0072B2",lwd=2)
  
}

# FUNCTION for 3 PCR
f.plot_3PCR<-function(i){
  prey<-m[which(m$Prey==i),] 
  
  ###### PROPORTION OF POSITIVES ON EACH TIME - FOR EACH PREY
  
  # FOR PCR1 
  prey_PCR1<-prey[which(prey$PCR=="PCR1"),]
  table_prey<-table(prey[which(prey$PCR=="PCR1"),c(4,5)])
  table_proportion<-prop.table(table_prey,2)
  df_proportion<-as.data.frame(table_proportion)
  table_positive_PCR1<-df_proportion[which(df_proportion$Detection=="1"),]
  time_variable_PCR1<-as.numeric(as.numeric(levels(table_positive_PCR1$Time)))
  freq_variable_PCR1<-table_positive_PCR1$Freq
  samplesize_PCR1<-colSums(table_prey)
  
  total_positives_PCR1<-as.vector(table_prey[2,])
  
  # FOR PCR2 
  prey_PCR2<-prey[which(prey$PCR=="PCR2"),]
  table_prey<-table(prey[which(prey$PCR=="PCR2"),c(4,5)])
  table_proportion<-prop.table(table_prey,2)
  df_proportion<-as.data.frame(table_proportion)
  table_positive_PCR2<-df_proportion[which(df_proportion$Detection=="1"),]
  time_variable_PCR2<-as.numeric(as.numeric(levels(table_positive_PCR2$Time)))
  freq_variable_PCR2<-table_positive_PCR2$Freq
  samplesize_PCR2<-colSums(table_prey)
  
  total_positives_PCR2<-as.vector(table_prey[2,])
  
  # FOR PCR3 
  prey_PCR3<-prey[which(prey$PCR=="PCR3"),]
  table_prey<-table(prey[which(prey$PCR=="PCR3"),c(4,5)])
  table_proportion<-prop.table(table_prey,2)
  df_proportion<-as.data.frame(table_proportion)
  table_positive_PCR3<-df_proportion[which(df_proportion$Detection=="1"),]
  time_variable_PCR3<-as.numeric(as.numeric(levels(table_positive_PCR3$Time)))
  freq_variable_PCR3<-table_positive_PCR3$Freq
  samplesize_PCR3<-colSums(table_prey)
  
  total_positives_PCR3<-as.vector(table_prey[2,])
  
  
  #### PROBIT ANALYSIS
  
  table_prey_total_positives<-data.frame(Time=time_variable_PCR1, 
                                         N1=samplesize_PCR1,PCR1=total_positives_PCR1,
                                         N2=samplesize_PCR2,PCR2=total_positives_PCR2)
  
  table_prey_prop<-data.frame(Time=time_variable_PCR1, 
                              N1=samplesize_PCR1,PCR1=freq_variable_PCR1,
                              N2=samplesize_PCR2,PCR2=freq_variable_PCR2)
  
  
  # FOR PCR1 ########################################################
  myprobit_PCR1 <- glm(Detection ~ Time, family = binomial(link = "probit"), 
                       data = prey_PCR1)
  
  # Test whether the model with predictors fits significantly better than a model with just an intercept (i.e. a null model). 
  # The test statistic is the difference between the residual deviance for the model with predictors and the null model.
  changedev<-with(myprobit_PCR1, null.deviance - deviance)
  df<-with(myprobit_PCR1, df.null - df.residual)
  quisquare<-with(myprobit_PCR1, pchisq(null.deviance - deviance, df.null - df.residual,lower.tail = FALSE))
  LL<-logLik(myprobit_PCR1)
  
  stat_results<-data.frame(ChangeDeviance=changedev,ChangeDF=df,
                           pChiSquare=quisquare,LogLikelihood=LL)
  
  ## PREDICT
  newdata_PCR1<- data.frame(Time=seq(0,48,0.25))
  newdata_PCR1[, c("p", "se")] <- predict(myprobit_PCR1, newdata_PCR1, type = "response", se.fit = TRUE)[-3]
  
  # Confidence intervals 
  ##########################
  ## grad the inverse link function
  ilink <- family(myprobit_PCR1)$linkinv
  ## add fit and se.fit on the **link** scale
  newdata_PCR1 <- bind_cols(newdata_PCR1, setNames(as_tibble(predict(myprobit_PCR1, newdata_PCR1, se.fit = TRUE)[1:2]),
                                                   c('fit_link','se_link')))
  ## create the interval and backtransform
  newdata_PCR1 <- mutate(newdata_PCR1,
                         fit_resp  = ilink(fit_link),
                         right_upr = ilink(fit_link + (2 * se_link)),
                         right_lwr = ilink(fit_link - (2 * se_link)))
  
  # FOR PCR2 ########################################################
  myprobit_PCR2 <- glm(Detection ~ Time, family = binomial(link = "probit"), 
                       data = prey_PCR2)
  with(myprobit_PCR2, pchisq(deviance, df.residual))
  
  # Test whether the model with predictors fits significantly better than a model with just an intercept (i.e. a null model). 
  # The test statistic is the difference between the residual deviance for the model with predictors and the null model.
  changedev<-with(myprobit_PCR2, null.deviance - deviance)
  df<-with(myprobit_PCR2, df.null - df.residual)
  quisquare<-with(myprobit_PCR2, pchisq(null.deviance - deviance, df.null - df.residual,lower.tail = FALSE))
  LL<-logLik(myprobit_PCR2)
  
  stat_results<-data.frame(ChangeDeviance=changedev,ChangeDF=df,
                           pChiSquare=quisquare,LogLikelihood=LL)
  ## PREDICT
  newdata_PCR2<- data.frame(Time=seq(0,48,0.25))
  newdata_PCR2[, c("p", "se")] <- predict(myprobit_PCR2, newdata_PCR2, type = "response", se.fit = TRUE)[-3]
  
  # Confidence intervals 
  ##########################
  ## grad the inverse link function
  ilink <- family(myprobit_PCR2)$linkinv
  ## add fit and se.fit on the **link** scale
  newdata_PCR2 <- bind_cols(newdata_PCR2, setNames(as_tibble(predict(myprobit_PCR2, newdata_PCR2, se.fit = TRUE)[1:2]),
                                                   c('fit_link','se_link')))
  ## create the interval and backtransform
  newdata_PCR2 <- mutate(newdata_PCR2,
                         fit_resp  = ilink(fit_link),
                         right_upr = ilink(fit_link + (2 * se_link)),
                         right_lwr = ilink(fit_link - (2 * se_link)))
  
  
  # FOR PCR3 ########################################################
  myprobit_PCR3 <- glm(Detection ~ Time, family = binomial(link = "probit"), 
                       data = prey_PCR3)
  # Test whether the model with predictors fits significantly better than a model with just an intercept (i.e. a null model). 
  # The test statistic is the difference between the residual deviance for the model with predictors and the null model.
  changedev<-with(myprobit_PCR3, null.deviance - deviance)
  df<-with(myprobit_PCR3, df.null - df.residual)
  quisquare<-with(myprobit_PCR3, pchisq(null.deviance - deviance, df.null - df.residual,lower.tail = FALSE))
  LL<-logLik(myprobit_PCR3)
  
  stat_results<-data.frame(ChangeDeviance=changedev,ChangeDF=df,
                           pChiSquare=quisquare,LogLikelihood=LL)
  
  ## PREDICT
  newdata_PCR3<- data.frame(Time=seq(0,48,0.25))
  newdata_PCR3[, c("p", "se")] <- predict(myprobit_PCR3, newdata_PCR3, type = "response", se.fit = TRUE)[-3]
  
  # Confidence intervals 
  ##########################
  ## grad the inverse link function
  ilink <- family(myprobit_PCR3)$linkinv
  ## add fit and se.fit on the **link** scale
  newdata_PCR3 <- bind_cols(newdata_PCR3, setNames(as_tibble(predict(myprobit_PCR3, newdata_PCR3, se.fit = TRUE)[1:2]),
                                                   c('fit_link','se_link')))
  ## create the interval and backtransform
  newdata_PCR3 <- mutate(newdata_PCR3,
                         fit_resp  = ilink(fit_link),
                         right_upr = ilink(fit_link + (2 * se_link)),
                         right_lwr = ilink(fit_link - (2 * se_link)))
  
  
  # MODELs
  mod1=glm(Detection ~ PCR + PCR:Time, family = binomial(link = "probit"), 
           data = prey)
  
  het=mod1$deviance/mod1$df.residual #heterogeneity factor
  
  mod2=glm(Detection ~ PCR + Time, family = binomial(link = "probit"), 
           data = prey)
  
  df1=mod1$df.residual
  df2=mod2$df.residual
  tstat=(mod2$deviance-mod1$deviance)/(het*(df2-df1))
  testp=1-pf(tstat,df2-df1,df1)
  
  ## PLOT VALUES AND MODEL PREDICTION #################################
  name_prey=i
  plot(c(0,48), c(0,1), type="n", xlab="Time after ingestion (hours)", ylab="Probability of detection", main=name_prey)
  ld<-seq(0,48,0.25)
  lines(ld, newdata_PCR1$p, lty=1,col= "#E69F00",lwd=2)
  points(time_variable_PCR1,freq_variable_PCR1, pch=20,col= "#E69F00")
  #text(time_variable_PCR1,freq_variable_PCR1, labels=samplesize_PCR1, cex= 0.7, pos=2)
  lines(ld, newdata_PCR1$right_upr, lty=2,col= "#E69F00",lwd=2)
  lines(ld, newdata_PCR1$right_lwr, lty=2,col= "#E69F00",lwd=2)
  
  lines(ld, newdata_PCR2$p, lty=1,col= "#0072B2",lwd=2)
  points(time_variable_PCR2,freq_variable_PCR2, pch=4,col= "#0072B2")
  #text(time_variable_PCR2,freq_variable_PCR2, labels=samplesize_PCR2, cex= 0.7, pos=2)
  lines(ld, newdata_PCR2$right_upr, lty=2,col= "#0072B2",lwd=2)
  lines(ld, newdata_PCR2$right_lwr, lty=2,col= "#0072B2",lwd=2)
  
  lines(ld, newdata_PCR3$p, lty=1,col= "#CC79A7",lwd=2)
  points(time_variable_PCR3,freq_variable_PCR3, pch=5,col= "#CC79A7")
  #text(time_variable_PCR2,freq_variable_PCR2, labels=samplesize_PCR2, cex= 0.7, pos=2)
  lines(ld, newdata_PCR3$right_upr, lty=2,col= "#CC79A7",lwd=2)
  lines(ld, newdata_PCR3$right_lwr, lty=2,col= "#CC79A7",lwd=2)
  
}


jpeg(file = "./plot_all.jpeg", bg="white", antialias = "default",width = 8, height = 10,  
     units = "in", res = 600)

par(fig = c(0,0.5,0.7,1),oma=c(1,1,0.8,1),mar=c(4,4,2,3),xpd=TRUE,mgp=c(3,1,0))
f.plot_3PCR("Tva")

par(fig = c(0.5,1,0.7,1),oma=c(1,1,0.8,1),mar=c(4,4,2,3),xpd=TRUE,mgp=c(3,1,0), new = T)
f.plot_2PCR("Lhu")

par(fig = c(0,0.5,0.35,0.65),oma=c(1,1,0.8,1),mar=c(4,4,2,3),xpd=TRUE,mgp=c(3,1,0), new = T)
f.plot_2PCR("Dis")

par(fig = c(0.5,1,0.35,0.65),oma=c(1,1,0.8,1),mar=c(4,4,2,3),xpd=TRUE,mgp=c(3,1,0), new = T)
f.plot_2PCR("Bim")

par(fig = c(0,0.5,0,0.3),oma=c(1,1,0.8,1),mar=c(4,4,2,3),xpd=TRUE,mgp=c(3,1,0), new = T)
f.plot_2PCR("Dmrc")

dev.off()


