###### PROBIT ANALYSIS
# pipeline adapted from https://stats.idre.ucla.edu/r/dae/probit-regression/
  
require(ggplot2)
require(aod)
require(lmtest)

setwd("/User_directory")
m=read.table("./matrix_GUT_detection_time.csv",header=TRUE,sep=",")
summary(m)

table(m)

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

changedev<-with(myprobit_PCR1, null.deviance - deviance)
df<-with(myprobit_PCR1, df.null - df.residual)
quisquare<-with(myprobit_PCR1, pchisq(null.deviance - deviance, df.null - df.residual,lower.tail = FALSE))
LL<-logLik(myprobit_PCR1)

stat_results<-data.frame(ChangeDeviance=changedev,ChangeDF=df,
                         pChiSquare=quisquare,LogLikelihood=LL)
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

changedev<-with(myprobit_PCR2, null.deviance - deviance)
df<-with(myprobit_PCR2, df.null - df.residual)
quisquare<-with(myprobit_PCR2, pchisq(null.deviance - deviance, df.null - df.residual,lower.tail = FALSE))
LL<-logLik(myprobit_PCR2)

stat_results<-data.frame(ChangeDeviance=changedev,ChangeDF=df,
                         pChiSquare=quisquare,LogLikelihood=LL)
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


