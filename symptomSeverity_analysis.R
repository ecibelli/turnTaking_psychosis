#..............................................................
#
# Preliminary analysis of turn-taking data
# Laura Sichlinger, Emily Cibelli, Matt Goldrick, Vijay Mittal
#
# Last updated: 4.11.18
#
#..............................................................

# Specify the dir where the data lives, and change to that directory
#wd = "/Users/laurasichlinger1/Desktop/turntaking_stats"
wd = "C:/Users/esc642/Box/schizophrenia/turnTaking/data"
setwd(wd)

# Load packages
library(ggplot2)      # plotting
library(lattice)      # plotting
library(gridExtra)    # allows multiple ggplots on one figure
library(lme4)         # mixed models (lmer)
library(RePsychLing)  # assess random effects in mixed models
# (install: https://github.com/dmbates/RePsychLing)

### Read in and set up data ------------------------------------

# Read in data that was set up in setupData.R
btpDF1 = read.csv(sprintf("%s/btpDF.csv", wd), 
                  head = T, row.names = 1)

length(unique(btpDF1$speaker))
# 72

# Read in data from BDI data frame
bdiDF = read.csv(sprintf("%s/BDI.APs.csv", wd), 
                 head = T, row.names = 1)

# Rename "ID" to "speaker" for merge
names(bdiDF)[names(bdiDF)=="ID"] <- "speaker"

# Fix speaker column in bdiDF to only be first 4 digits
bdiDF$speaker = substr(bdiDF$speaker, 1, 4)

# How many speakers in the current dataset are on antipsychotics?
table(bdiDF[bdiDF$speaker %in% btpDF1$speaker,]$antipsychotic_current)
# no yes 
# 70   2 

# Which speakers are on antipsychotics?
bdiDF[bdiDF$speaker %in% btpDF1$speaker &
        bdiDF$antipsychotic_current == "yes",]$speaker
# 1] "1005" "1006"

# Merge dataframes
btpDF = merge(btpDF1, bdiDF)

# Make sure we still have 72 speakers
length(unique(btpDF$speaker))
# 72

# Restrict data frame to rows with between-turn pauses only and UHR only
btpDF = btpDF[btpDF$intervalType == "BTP" & btpDF$group == "UHR",]

# ..............................

# Log-transform dependent variable (duration of BTPs)

btpDF$logDur = log(btpDF$dur)

## Visualize difference between linear and log duration 
## of the dependent variable

# Set up two panels for plotting
par(mfrow=c(1,2))
# linear BTP plot, with main plot title and x-axis title
plot(density(btpDF$dur), main = "BTP duration (linear scale)", 
     xlab = "Dur (ms)")
# log BTP plot, with main plot title and x-axis title
plot(density(btpDF$logDur), main = "BTP duration (log scale)",
     xlab = "Dur (log-transformed)")
# Re-set plot panels to 1
par(mfrow = c(1,1))

# .....................................

# Summarize reponse complexity

# Try a three-way version:
# low: yes/no, yes/no+, single-word
# medium: question, short sentence, yes/no++
# high: complex
btpDF$respComplexity1 = ifelse(btpDF$response %in% c("yn", "yn+", "sw"),
                               "low", ifelse(btpDF$response %in% 
                                               c("qu", "ss", "yn++"), "medium", "high"))

# Plot this grouping to evaluate if all three levels are needed
# Set up the parameter(s) to be plotted - logDur as the sole variable
ggplot(btpDF, aes(x = logDur)) + 
  # set up the type of plot (density) and groups to get their own shape/color
  # alpha - transparency of each density
  geom_density(aes(group=respComplexity1, 
                   fill = respComplexity1), alpha = 0.3) + 
  # split the plot into panels by group, with two columns
  facet_wrap(~group, ncol = 2) +
  # aesthetics
  theme_bw()
# Looks like high and medium group together

# Reconstruct another summarization which treats both high and medium = "high"
btpDF$respComplexity = as.factor(ifelse(
  btpDF$respComplexity1 == "low", "low", "high"))

# Re-plot
ggplot(btpDF, aes(x = logDur)) + 
  geom_density(aes(group=respComplexity, 
                   fill = respComplexity), alpha = 0.3) + 
  facet_wrap(~group, ncol = 2) + 
  theme_bw()

# Final coding of response complexity:
# low: yes/no, yes/no+, single-word
# high: question ,short sentence, yes/no++, complex

# Summarize symptom severity

# Symptom scale in general:
# 0-2: absent or mildly present
# 3-4: moderate 
# 5-6: severe

# Positive symptoms
# Five question blocks P1-P5 - each is rated according to symtpom scale above (0-6)
# Positive symtpom scale in total:
# mild: 0-10
# moderate: 11-20
# severe: 21-30
btpDF$pos = ifelse(btpDF$posT1 %in% 0:10, "mild",
                   ifelse(btpDF$posT1 %in% 11:20, "moderate", "severe"))

# Plot this grouping 
# Set up the parameter(s) to be plotted - logDur as the sole variable
ggplot(btpDF, aes(x = logDur)) + 
  # set up the type of plot (density) and groups to get their own shape/color
  # alpha - transparency of each density
  geom_density(aes(group=pos, 
                   fill = pos), alpha = 0.3) + 
  # split the plot into panels by group, with two columns
  # aesthetics
  theme_bw()

# bar chart
ggplot(data=btpDF, aes(x=pos, y=logDur, fill=pos)) +
  geom_col()

# Negative symptoms
# Six question blocks N1-N6 - each is rated according to symtpom scale above (0-6)
# Negative symptom scale in total:
# mild: 0-12
# moderate: 13-24
# severe: 25-36
btpDF$neg = ifelse(btpDF$negT1 %in% 0:12, "mild",
                   ifelse(btpDF$negT1 %in% 13:24, "moderate", "severe"))

# Plot
ggplot(btpDF, aes(x = logDur)) + 
  geom_density(aes(group=negT1, 
                   fill = negT1), alpha = 0.3) + 
  facet_wrap(~group, ncol = 2) +
  theme_bw()

# bar chart
ggplot(data=btpDF, aes(x=neg, y=logDur, fill = neg)) +
  geom_col()

### Dependent variables - description ----------

# Fixed effects structure:

# Column:         Description/explanation:

# age             Age of the participant
# sex             Male or female
# antipsychotics  Antipsychotic treatment, yes/no
# respDur         Duration of the response after the pause
# questDur        Duration of the question preceding the pause
# pauseType       What portion of the interview did the question come from?
#                 (background, SIPS, or "after" questions)
# respComplexity  Complexity of the response (high/low)
# pos             Severity of positive symptoms (mild/moderate)
# neg             Severity of negative symtoms (mild/moderate)

# Two-way interactions of severity of symptoms and linguistic variables:
# Positive:
# pos * questDur
# pos * respDur
# pos * pauseType
# pos * respComplexity

# Negative:
# neg * questDur
# neg * respDur
# neg * pauseType
# neg * respComplexity

# Two way interaction:
# questDur * respDur

# Three way interaction:
# pos * questDur * respDur
# neg * questDur * respDur

# .........................

# Random effects structure: 
# intercept for speaker
# slopes for all variables that the data will support

### Center variables and numerically code ---------------------------------

# Centering variables will allow us to interpret each individual effect at the mean value of the others. For numeric predictors, this means that the mean value will be set to 0. 

# Categorical predictors will be numerically coded as well. for these, groups will be assigned to a binary split. For a perfectly balanced binary factor (e.g. 50% male, 50% female in the data set), one category would be set to -0.5 (the one that comes first alphabetically), and one to 0.5 (the one that is second alphabetically). This ensures that they sum to zero. If there are uneven numbers in the groups, the myCenter() function below will balance the weights so that zero is still the center - in this case, the values may not be precisely 0.5 and -0.5.  

# We will indicate these centered, numerically-coded variables with a "C" at the end of the column name.

# You can check to ensure that each column is properly centered by running:
# summary(btpDF) 
# The mean for each centered column should be 0. 

# .......

# To make sure we are centering uneven groups properly, we'll use the function described here:
# https://hlplab.wordpress.com/2009/04/27/centering-several-variables/

myCenter= function(x) {
  if (is.numeric(x)) { return(x - mean(x, na.rm=T)) }
  if (is.factor(x)) {
    x= as.numeric(x)
    return(x - mean(x, na.rm=T))
  }
  if (is.data.frame(x) || is.matrix(x)) {
    m= matrix(nrow=nrow(x), ncol=ncol(x))
    colnames(m)= paste("c", colnames(x), sep="")
    for (i in 1:ncol(x)) {
      m[,i]= myCenter(x[,i])
    }
    return(as.data.frame(m))
  }
}

# ........

btpDF$ageC = myCenter(btpDF$age)
btpDF$sexC = myCenter(btpDF$sex) # female =~ -0.50, male =~ +0.50
btpDF$negC = myCenter(btpDF$negT1) 
btpDF$posC = myCenter(btpDF$posT1) 
btpDF$antipsychoticC = myCenter(btpDF$antipsychotic_current) # no = ~ -0.16, yes = ~ +0.84

# also scale these numeric predictors so that they are not on a very different scale from other predictors
btpDF$respDurC = scale(btpDF$respDur, center = T, scale = T)
btpDF$questDurC = scale(btpDF$questDur, center = T, scale = T)

# Reorder response complexity so that low = -0.5
btpDF$respComplexity = factor(btpDF$respComplexity,
                              levels(btpDF$respComplexity)[c(2:1)])
btpDF$respComplexC = myCenter(btpDF$respComplexity)

# .........

# pauseType is a three-way predictor
# We'll set up one predictor to compare background to the others.
# We'll set up a 2nd predictor to compare SIPS to after (with background = 0).

# Predictor 1: background/others
btpDF$bkgdC = ifelse(btpDF$pauseType == "background", -0.5, 0.5)
btpDF$bkgdC = myCenter(btpDF$bkgdC)

# Predictor 2: sips/after
# We need to find the proportion of each in the dataset, excluding background, so as to know what weights to assign.
sipsProp = nrow(btpDF[btpDF$pauseType == "SIPS",])/nrow(btpDF[btpDF$pauseType %in% c("SIPS", "after"),])
afterProp = nrow(btpDF[btpDF$pauseType == "after",])/nrow(btpDF[btpDF$pauseType %in% c("SIPS", "after"),])

# assign the weight of one category to the other, to balance them
# set the first to negative as well
# together, this will ensure they average to 0
btpDF$sipsAfterC = ifelse(btpDF$pauseType == "SIPS", (afterProp*-1), 
                          ifelse(btpDF$pauseType == "after", sipsProp, 
                                 0)) # background = 0
# ...............

# Fix some NAs by setting to 0
btpDF[is.na(btpDF$bkgdC),]$bkgdC = 0
btpDF[is.na(btpDF$sipsAfterC),]$sipsAfterC = 0
btpDF[is.na(btpDF$questDurC),]$questDurC = 0

### Model selection procedure description ----------------------------------

# To select the random effects structure of the model, we will use the procedure recommended in Bates et al. (2015), "Parsimonious Mixed Models".

# https://arxiv.org/abs/1506.04967

# This procedure allows us to balance the principle of maximal models - including all fixed effects as random slopes - with the fact that these models often fail to converge, and even when they do, there may not be sufficient data to support them. 

# Assuming all models do converge, they will be fit as follows:

# First model: maximal random effects structure - all fixed effects included as slopes - including correlations of intercept(s) and slopes. If this converges, check to see if there are parameters in the random effects structure which have zero variance; the existence of these suggests that there are too many parameters in the model, and that some are not contributing any explanatory power.

# Second model: if there are zero-variance parameters in model 1, or if it does not converge, fit the same model without correlations in the random effects structure. This will allow us to determine which components do not explain any variance in the model. (The zero-correlation model has fewer parameters, so it is also more likely to converge than the first model.)

# Third model, if the second model does *not* converge: remove the most complex components in the random slopes. So first, remove any three-way interactions if they are present, and re-fit. If this does not converge, then remove all two-way interactions, re-fit, and so on. Then proceed with the procedure below.

# Next model, if the second model *does* converge, or once the third model above does: remove all the zero-variance components identified in that model. Run an anova() to make sure that this model is not a worse fit than the model that includes those zero-variance parameters.

# Final model: once there is a stable modle that converges and has no zero-variance components, add the correlations back in to the model. If it converges, run an anova to compare the model with and without the correlations. If there is a significant difference in fit, select the model with correlations. If not, select the model without (prefer a simpler model when there is no evidence for improvement with a more complex model.)

# Trim residuals: once we have a final model structure, we will trim the points in the model with high residuals (> 2.5) and re-fit, to ensure that we are not over-fitting to extreme data points.

# ..................

# A final note: only predictors which vary within a speaker can be used as random slopes for the speaker intercept. Factors such as sex, group, and age do not vary for an individual, and so will not be included in the by-speaker random slopes.

### Model selection --------------------------------------------------

# Model 1: maximal random effects structure
group.lmer1 = lmer(logDur ~ 
                     # -- Fixed Effects -- #
                     # Simple effects: speaker variables
                     ageC + sexC + antipsychoticC + negC + posC + 
                     # Simple effects: linguistic variables
                     respDurC + questDurC + respComplexC +
                     bkgdC + sipsAfterC + 
                     # Two-way interactions: neg * linguistic vars
                     negC:questDurC + negC:respDurC +
                     negC:respComplexC + negC:bkgdC + negC:sipsAfterC +
                     # Two-way interactions: pos * linguistic vars
                     posC:questDurC + posC:respDurC +
                     posC:respComplexC + posC:bkgdC + posC:sipsAfterC +
                     # Two-way interaction: quest * resp durations
                     questDurC:respDurC +
                     # Three-way: neg * resp * group 
                     negC:questDurC:respDurC +
                     # Three-way: pos * resp * group 
                     posC:questDurC:respDurC +
                     ## -- Random Effects -- ##
                     # Only variables that vary within-speaker
                     (1 + respDurC + questDurC + respComplexC +
                        bkgdC + sipsAfterC +
                        questDurC:respDurC
                      | speaker),
                   data = btpDF,
                   REML = FALSE, # don't use REML when comparing fixed effects
                   verbose = 2) # Track progress in fitting model

# Model did converge

# .............

# Model 2: zero correlation model
group.lmer2 = lmer(logDur ~ 
                     ageC + sexC + antipsychoticC + negC + posC +
                     respDurC + questDurC + respComplexC +
                     bkgdC + sipsAfterC + 
                     negC:questDurC + negC:respDurC +
                     negC:respComplexC + negC:bkgdC + negC:sipsAfterC +
                     posC:questDurC + posC:respDurC +
                     posC:respComplexC + posC:bkgdC + posC:sipsAfterC +
                     questDurC:respDurC +
                     negC:questDurC:respDurC +
                     posC:questDurC:respDurC +
                     (1 + respDurC + questDurC + respComplexC +
                        bkgdC + sipsAfterC +
                        questDurC:respDurC || speaker), # || = no correlation
                   data = btpDF,
                   REML = FALSE, verbose = 2) 

# Check for zero correlation parameters in the model
summary(rePCA(group.lmer2))
# Importance of components:
#   [,1]   [,2]   [,3]    [,4] [,5] [,6] [,7]
# Standard deviation     0.4966 0.4055 0.2607 0.12208    0    0    0
# Proportion of Variance 0.4993 0.3329 0.1376 0.03018    0    0    0
# Cumulative Proportion  0.4993 0.8322 0.9698 1.00000    1    1    1

# 3 components explain no variance (5, 6, 7)

# Check to see which components 
summary(group.lmer2, corr = FALSE)
# Random effects:
#   Groups    Name               Variance Std.Dev.
# speaker   (Intercept)        0.104882 0.32385 
# speaker.1 respDurC           0.000000 0.00000 
# speaker.2 questDurC          0.009508 0.09751 
# speaker.3 respComplexC       0.043366 0.20824 
# speaker.4 bkgdC              0.000000 0.00000 
# speaker.5 sipsAfterC         0.157322 0.39664 
# speaker.6 respDurC:questDurC 0.000000 0.00000 
# Residual                     0.637976 0.79873 
# Number of obs: 1636, groups:  speaker, 36

# ..................

# Model 3: remove the zero variance random slope
group.lmer3 = lmer(logDur ~ 
                     ageC + sexC + antipsychoticC + negC + posC +
                     respDurC + questDurC + respComplexC +
                     bkgdC + sipsAfterC + 
                     negC:questDurC + negC:respDurC +
                     negC:respComplexC + negC:bkgdC + negC:sipsAfterC +
                     posC:questDurC + posC:respDurC +
                     posC:respComplexC + posC:bkgdC + posC:sipsAfterC +
                     questDurC:respDurC +
                     negC:questDurC:respDurC +
                     posC:questDurC:respDurC +
                     (1 + questDurC + respComplexC + sipsAfterC || speaker), # || = no correlation
                   data = btpDF, REML = FALSE, verbose = 2) 

# Check to see if removing this component negatively impacts the model 
anova(group.lmer2, group.lmer3)
# Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# group.lmer3 29 4082.7 4239.3 -2012.3   4024.7                        
# group.lmer2 32 4088.7 4261.5 -2012.3   4024.7     0      3          1

# p = 1, no evidence for a reduction in model explanation as a function of removing respDurC, bkgdC,  and respDurC:questDurC from the random slopes
# prefer the simpler model (lmer3)

# Check to ensure there are no other zero-variance components
summary(rePCA(group.lmer3))
# $speaker
# Importance of components:
#   [,1]   [,2]   [,3]    [,4]
# Standard deviation     0.4966 0.4055 0.2607 0.12208
# Proportion of Variance 0.4993 0.3329 0.1376 0.03018
# Cumulative Proportion  0.4993 0.8322 0.9698 1.00000

# ................

# Model 4: add correlations back in
group.lmer4 = lmer(logDur ~ 
                     ageC + sexC + antipsychoticC + negC + posC +
                     respDurC + questDurC + respComplexC +
                     bkgdC + sipsAfterC + 
                     negC:questDurC + negC:respDurC +
                     negC:respComplexC + negC:bkgdC + negC:sipsAfterC +
                     posC:questDurC + posC:respDurC +
                     posC:respComplexC + posC:bkgdC + posC:sipsAfterC +
                     questDurC:respDurC +
                     negC:questDurC:respDurC +
                     posC:questDurC:respDurC +
                     (1 + questDurC + respComplexC + sipsAfterC | speaker), 
                   data = btpDF,
                   REML = FALSE, verbose = 2,
                   control=lmerControl(optCtrl=list(maxfun=1e06)))

# Compare zero correlation and correlation model
anova(group.lmer3, group.lmer4)
# p = 0.05078
# marginal but close, prefer lmer4

# Final model: lmer4

### Refit final model and assess fixed effects ----------------------

# Remove data points with extreme residuals (> 2.5) in the final model
btpDF.trim = btpDF[abs(scale(resid(group.lmer4)))<2.5,]

# How much data was retained?
nrow(btpDF.trim)/nrow(btpDF)
#  0.9871717

# Re-fit model to this trimmed data set
group.lmer.trim = lmer(logDur ~ 
                         ageC + sexC + antipsychoticC + negC + posC +
                         respDurC + questDurC + respComplexC +
                         bkgdC + sipsAfterC + 
                         negC:questDurC + negC:respDurC +
                         negC:respComplexC + negC:bkgdC + negC:sipsAfterC +
                         posC:questDurC + posC:respDurC +
                         posC:respComplexC + posC:bkgdC + posC:sipsAfterC +
                         questDurC:respDurC +
                         negC:questDurC:respDurC +
                         posC:questDurC:respDurC +
                         (1 + questDurC + respComplexC + sipsAfterC | speaker), 
                   data = btpDF.trim,
                   REML = FALSE,
                   control=lmerControl(optCtrl=list(maxfun=1e06))) 

summary(group.lmer.trim)
# AIC      BIC   logLik deviance df.resid 
# 4024.2   4212.7  -1977.1   3954.2     1580 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -3.7837 -0.6244  0.0766  0.6686  2.9734 
# 
# Random effects:
#   Groups   Name         Variance Std.Dev. Corr             
# speaker  (Intercept)  0.10454  0.3233                    
# questDurC    0.01302  0.1141   -0.63            
# respComplexC 0.03157  0.1777    0.63 -0.86      
# sipsAfterC   0.81778  0.9043   -0.96  0.39 -0.40
# Residual              0.63376  0.7961                    
# Number of obs: 1615, groups:  speaker, 36
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)              5.9763862  0.0663198   90.11
# ageC                    -0.0889454  0.0320015   -2.78
# sexC                     0.3648249  0.0945784    3.86
# antipsychoticC           0.2220802  0.1543969    1.44
# negC                     0.0140397  0.0108142    1.30
# posC                     0.0243913  0.0203702    1.20
# respDurC                 0.0377350  0.0245446    1.54
# questDurC                0.0563734  0.0305055    1.85
# respComplexC             0.3474535  0.0566898    6.13
# bkgdC                    0.1158841  0.1234300    0.94
# sipsAfterC               0.4188833  0.2588457    1.62
# negC:questDurC           0.0036874  0.0048573    0.76
# negC:respDurC           -0.0061591  0.0039705   -1.55
# negC:respComplexC       -0.0076840  0.0090371   -0.85
# negC:bkgdC               0.0196259  0.0145513    1.35
# negC:sipsAfterC         -0.0982901  0.0359200   -2.74
# posC:questDurC          -0.0030418  0.0085914   -0.35
# posC:respDurC            0.0118071  0.0076102    1.55
# posC:respComplexC       -0.0008390  0.0159212   -0.05
# posC:bkgdC              -0.0596637  0.0425452   -1.40
# posC:sipsAfterC          0.1683878  0.0663721    2.54
# respDurC:questDurC      -0.0059455  0.0210880   -0.28
# negC:respDurC:questDurC -0.0009374  0.0036675   -0.26
# posC:respDurC:questDurC  0.0003409  0.0061575    0.06

# ..........................

## Assess significance of fixed effects

# P-values are not reliable for mixed models and not provided by lmer. In lieu of them, we will use chi-square tests to compare versions of the model with and without a particular predictor of interest. If the model fit significantly improves with the inclusion of the predictor, we consider it to be significant. 

# This function runs the test and provides a summary of the chi-square statistic and significance.

#lme convenience function
chiReport.func <- function(a){
  ifelse (a$"Pr(>Chisq)"[2] > .0001,
          return(paste("chisq(",a$"Chi Df"[2],")=",round(a$Chisq[2],2)," p = ",round(a  $"Pr(>Chisq)"[2],4),sep="")),
          return(paste("chisq(",a$"Chi Df"[2],")=",round(a$Chisq[2],2),", p < .0001")))
}

chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-ageC)))
# "chisq(1)=2.28 p = 0.1311"
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-sexC)))
# "chisq(1)=5.58 p = 0.0182" ***
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-antipsychoticC)))
# "chisq(1)=0.82 p = 0.3665"
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-negC)))
# "chisq(1)=1.89 p = 0.1697"
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-posC)))
# "chisq(1)=1.41 p = 0.235"
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-respDurC)))
# "chisq(1)=3.55 p = 0.0594" ***
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-questDurC)))
# "chisq(1)=2.38 p = 0.1232"
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-respComplexC)))
# "chisq( 1 )= 20.32 , p < .0001" ***
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-bkgdC)))
# "chisq(1)=0.6 p = 0.4374"
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-sipsAfterC)))
# "chisq(1)=2.42 p = 0.1196"
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-negC:questDurC)))
# "chisq(1)=0.28 p = 0.5943"
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-negC:respDurC)))
# "chisq(1)=4.93 p = 0.0264" ***
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-negC:respComplexC)))
# "chisq(1)=0.84 p = 0.3606"
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-negC:bkgdC)))
# "chisq(1)=2.48 p = 0.1153"
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-negC:sipsAfterC)))
# "chisq(1)=1.77 p = 0.183"
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-posC:questDurC)))
# "chisq(1)=0.1 p = 0.7533"
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-posC:respDurC)))
# "chisq(1)=2.29 p = 0.1303"
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-posC:respComplexC)))
# "chisq(1)=0.27 p = 0.6052"
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-posC:bkgdC)))
# "chisq(1)=3.47 p = 0.0624" ***
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-posC:sipsAfterC)))
# "chisq(1)=3.14 p = 0.0765" ***
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-respDurC:questDurC)))
# "chisq(1)=0.13 p = 0.7201"
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-negC:respDurC:questDurC)))
# "chisq(1)=0 p = 0.9678"
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-posC:respDurC:questDurC)))
# "chisq(1)=0.09 p = 0.7641"
