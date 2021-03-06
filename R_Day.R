library(survival)
library(rms)
#Load car package for outlier analysis
require(car)
pharynx <- read.csv("C://Users//lenovo//Desktop//R Day//pharynx.csv"
                    , header=TRUE)
head(pharynx)
summary(pharynx)

#Cox PH with COND, T_STAGE, and N_STAGE as our predictors and TX stratified
cox.ph_1 <- coxph(Surv(x$TIME, x$STATUS)~x$COND + x$T_STAGE + x$N_STAGE + strata(x$TX), 
                  data = pharynx)
summary(cox.ph_1)

#Cox PH with COND, T_STAGE, and N_STAGE as our predictors
cox.ph_2 <- coxph(Surv(x$TIME, x$STATUS)~x$COND+ x$T_STAGE + x$N_STAGE, 
                  data = pharynx)
summary(cox.ph_2)

#Create initial regression model
REG <- lm(TIME ~ COND + T_STAGE, data = pharynx)
plot(REG)

#Outlier Tests
outlierTest(REG)
qqPlot(REG, main = "QQ Plot")
leveragePlots(REG)
avPlots(REG)

fit <- survreg(Surv(TIME, STATUS)~1, pharynx)
#Cook's Distance
cutoff <- 4/((nrow(REG)-length(fit$coefficients)-2))
plot(REG, which=4, cook.levels=cutoff)

#Influence plot (Click bubbles to identify numbers)
influencePlot(REG, id.method="identify", 
              main="Figure 1: Influence Plot", 
              sub="Circle size is proportial to Cook's Distance" )

##Removes outlier
x <- pharynx[which(pharynx$CASE != 159),]
nrow(x)


#Create initial regression model without outlier
REG_1 <- lm(TIME ~ COND + T_STAGE, data = x)
plot(REG_1)

#Outlier Tests
outlierTest(REG_1)
qqPlot(REG_1, main = "QQ Plot")
leveragePlots(REG_1)
avPlots(REG_1)

fit1 <- survreg(Surv(TIME, STATUS)~1, pharynx)
#Cook's Distance
cutoff1 <- 4/((nrow(REG_1)-length(fit1$coefficients)-2))
plot(REG_1, which=4, cook.levels=cutoff1)

#Influence plot (Click bubbles to identify numbers)
influencePlot(REG_1, id.method="identify", 
              main="Figure 2: Influence Plot After Removal Of Outlier", 
              sub="Circle size is proportial to Cook's Distance")

#Add survival object. status == 1 is dead
x$SurvObj <- with(x, Surv(x$TIME, x$STATUS == 1))

#Kaplan-Meier estimator. The "log-log" CI is preferred.
km.as.one <- survfit(SurvObj ~ 1, data = x, conf.type = "log-log")

#Plot
plot(km.as.one, conf = F, mark.time = F, xlab= "TIME", ylab= "Survival probabilty",
     main="Figure 3: Survival Estimate")

#Cox PH with COND and T_STAGE as predictors and TX stratified
cox.ph_3 <- coxph(Surv(x$TIME, x$STATUS)~x$COND + x$T_STAGE + strata(x$TX), data = x)
summary(cox.ph_3)

#Cox PH with COND and T_STAGE as predictors
cox.ph_4 <- coxph(Surv(x$TIME, x$STATUS)~x$COND+ x$T_STAGE, data = x)
summary(cox.ph_4)


plot(survfit(formula = Surv(x$TIME, x$STATUS)~ x$COND + x$T_STAGE, 
             data = x, conf.type="none"), 
     lty=0, xlab="Time", ylab="Survival Probability" )

cox.ph_4.zph <- cox.zph(time.dep, transform = 'log')
time.dep.zph
plot(time.dep.zph[1])
abline (h=0, lty=2)
