library(AER)
library(tidyverse)
library(caret)
library(pscl)
library(DescTools)
source('linktest.R')
source('marginaleffects.R')
library(ResourceSelection)
library(stargazer)
library(LogisticDx)
library(mfx)
library(ggplot2)
data <- read.csv('mushroom.csv')
data <- subset(data, select = c(class, gill.color, cap.color, habitat, season, cap.diameter, stem.height, stem.width))
###Numerical var
hist(data$cap.diameter) #highly left-skewed, lets log it
hist(log(data$cap.diameter))
data$logcap <- log(data$cap.diameter)
hist(data$stem.height)
hist(log(data$stem.height))
data$logsth <- log(data$stem.height)
hist(data$logsth)
hist(data$stem.width)
hist(log(data$stem.width))
data$logstw <- log(data$stem.width)
data$logsth[is.infinite(data$logsth)] <- 0
data$logstw[is.infinite(data$logstw)] <- 0
hist(data$logcap)
hist(data$logsth)
hist(data$logstw)
# normality test
shapiro.test(sample(data$cap.diameter, 5000)) #H0 rejected - not normal
shapiro.test(sample(data$logcap, 5000)) #H0 rejected - still not normal but looks more like it.
shapiro.test(sample(data$logsth, 5000))
### Transformation of categorical vars into factors
data$class <- ifelse(data$class == "p", 1, 0) #this ones a binary
data$gill.color <- as.factor(data$gill.color)
data$cap.color <- as.factor(data$cap.color)
data$habitat <- as.factor(data$habitat)
data$season <- as.factor(data$season)
# Frequency table for each categorical variable
table(data$gill.color)
table(data$cap.color)
table(data$habitat)
table(data$season)
ggplot(data, aes(x = gill.color, fill = as.factor(class))) + 
  geom_bar(position = "fill") + 
  labs(title = "Gill Color vs Class", y = "Proportion")
ggplot(data, aes(x = cap.color, fill = as.factor(class))) + 
  geom_bar(position = "fill") + 
  labs(title = "Cap Color vs Class", y = "Proportion")
ggplot(data, aes(x = habitat, fill = as.factor(class))) + 
  geom_bar(position = "fill") + 
  labs(title = "Habitat vs Class", y = "Proportion")
ggplot(data, aes(x = season, fill = as.factor(class))) + 
  geom_bar(position = "fill") + 
  labs(title = "Season vs Class", y = "Proportion")
# chisq test on association
chisq.test(table(data$class, data$gill.color))
chisq.test(table(data$class, data$cap.color))
chisq.test(table(data$class, data$habitat))
chisq.test(table(data$class, data$season))
### Creation of dummies
data_transformed <- model.matrix(~ . - 1, data = data)
data_transformed <- as.data.frame(data_transformed)
data_transformed <- dplyr::select(data_transformed, -gill.colorb) #we remove this one manually, as it is not removed correctly by the code
# View the structure of the transformed data
str(data_transformed)
###Multicolinearity - we construct a linear model and do some pseudo GETS with vif
linmod <- lm(class ~ ., data = data_transformed)
vif(linmod)
#remove highly colinear variables from the dataset
data_transformed <- dplyr::select(data_transformed, -c(gill.colorw, stem.height, cap.colorn, stem.width, logcap)) #we remove cap.diameter instead of logcap
linmod <- lm(class ~ ., data = data_transformed)
vif(linmod)
#### Constructing a model

#General model for GETS
#general logit
genmod <- glm(class ~ ., data = data_transformed, 
              family = binomial(link = 'logit'))
#general probit
probmod <- glm(class ~ ., data = data_transformed, family = binomial(link = 'probit'))
summary(genmod)
summary(probmod) #stargazer would be optimal but this is nicer for screenshots
#Blank model to check if theres any sense
blankmod <- glm(class ~ 1, data = data_transformed, family = binomial(link = 'logit'))
summary(blankmod)
lrtest(genmod, blankmod) #the model does matter - the vars are jointly significant
#iter0 - check if all the insignificant variables are jointly insignificant
insmod <- glm(class ~ gill.colorg + cap.colorl + habitath + habitatp + habitatu + habitatw + seasonu, data = data_transformed, family = binomial(link = 'logit'))
summary(insmod)
anova(blankmod, insmod, test = 'LRT') #all the insignificant vars at once are jointly singificant
### iter1
summary(genmod)
linktest(genmod) #genmod passes the linktest
##iter2 - we start with removing habitatu at pval 90
protomod1 <- glm(class ~ . - habitatu, data = data_transformed, family = binomial(link = 'logit'))
linktest(protomod1)
anova(genmod, protomod1, test = 'LRT') #we reject H0, thus indicating that the variable was important after all
##iter3 - instead we remove habitatw at pval 83
protomod2 <- glm(class ~ . - habitatw , data = data_transformed, family = binomial(link = 'logit'))
anova(genmod, protomod2, test = 'LRT') #we reject H0, thus indicating that the model is different than the general model
##iter4 - instead we remove habitatp with pval 83
protomod3 <- glm(class ~ . - habitatp , data = data_transformed, family = binomial(link = 'logit'))
anova(genmod, protomod3, test = 'LRT') #we reject H0, thus indicating that the model is different than the general model
##iter5 - instead we remove habitath with pval 72
protomod4 <- glm(class ~ . - habitath , data = data_transformed, family = binomial(link = 'logit'))
anova(genmod, protomod4, test = 'LRT') #H0 stands, indicating this variable was in fact insignificant
##iter6 - seasonu is removed with pval 53
protomod5 <- glm(class ~ . - habitath - seasonu , data = data_transformed, family = binomial(link = 'logit'))
anova(genmod, protomod5, test = 'LRT') #H0 stands, indicating this variable was in fact insignificant
summary(protomod5)
##iter7 - gill.colorg is removed with pval 50
protomod6 <- glm(class ~ . - habitath - seasonu - gill.colorg , data = data_transformed, family = binomial(link = 'logit'))
anova(genmod, protomod6, test = 'LRT') #H0 stands, indicating this variable was in fact insignificant
summary(protomod6)
##iter8 - cap.colorl is removed with pval 16
protomod7 <- glm(class ~ . - habitath - seasonu - gill.colorg - cap.colorl, data = data_transformed, family = binomial(link = 'logit'))
anova(genmod, protomod7, test = 'LRT') #H0 stands, indicating this variable was in fact insignificant
summary(protomod7) #no insignificant vars
#gof tests
finmod <- glm(class ~ gill.colore + gill.colorf + gill.colork + gill.colorn + gill.coloro + cap.colorl + gill.colorp + gill.colorr + gill.coloru + gill.colory + cap.colore + cap.colorg + cap.colork + cap.coloro + cap.colorp + cap.colorr + cap.coloru + cap.colorw + cap.colory + habitatg + habitatl + habitatm + habitatp + habitatu + habitatw + seasons + seasonw + cap.diameter + logsth + logstw
, data = data_transformed, family = binomial(link = 'logit')) #for logisticDx it has to be this way
loggen_gof <- glm(class~ gill.colore + gill.colorf + gill.colorg + gill.colork + gill.colorn + gill.coloro + gill.colorp + gill.colorr + gill.coloru + gill.colory + cap.colore + cap.colorg + cap.colork + cap.colorl + cap.coloro + cap.colorp + cap.colorr + cap.coloru + cap.colorw + cap.colory + habitatg + habitath + habitatl + habitatm + habitatp + habitatu + habitatw + seasons + seasonu + seasonw + cap.diameter + logsth + logstw
                  , data = data_transformed, family = binomial(link = 'probit')) #here i change between gen logit and probit for logisticDx
predicted_probs <- predict(protomod7, type = "response")
hoslem.test(data_transformed$class, predicted_probs) #i change between models here
a <- LogisticDx::gof(loggen_gof, g = 9) #i change between models here
PseudoR2(probmod, c("McKelveyZavoina"))#i change bet... you get it
###count r2 from stack overflow
countR2<-function(m) mean(m$y==round(m$fitted.values))
countR2(genmod)
countR2(probmod)
adj.countR2<-function(m) {
  n<-length(m$y)
  k<-max(table(m$y))
  correct<-table(m$y==round(m$fitted.values))[["TRUE"]]
  (correct-k)/(n-k)
}
adj.countR2(probmod)
stargazer(protomod7, probmod, type = 'text')
##marginal effects
probitmfx(formula = class ~ . - habitath - seasonu - gill.colorg - cap.colorl, data = data_transformed, atmean = T)
##wald test
waldtest(protomod7, genmod, test = 'Chisq') #the test gives us no reason to drop the null hypothesis of omitted coefs being different from 0





