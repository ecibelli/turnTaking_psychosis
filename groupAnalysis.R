#..............................................................
#
# Analysis of turn-taking data: group analysis
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
bdiDF[bdiDF$speaker %in% btpDF1$speaker & bdiDF$antipsychotic_current == "yes",]$speaker
# 1] "1005" "1006"

# Merge dataframes
btpDF = merge(btpDF1, bdiDF)

# Make sure we still have 72 speakers
length(unique(btpDF$speaker))
# 72

# Restrict data frame to rows with between-turn pauses only
btpDF = btpDF[btpDF$intervalType == "BTP",]

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

# Summarize response complexity

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
# high: question ,short sentence, yes/no++

### Dependent variables - description ----------

# Fixed effects structure:

# Column:         Description/explanation:

# age             Age of the participant
# sex             Male or female
# bdi             Depression symptoms score
# antipsychotics  Antipsychotic treatment, yes/no
# respDur         Duration of the response after the pause
# questDur        Duration of the question preceding the pause
# pauseType       What portion of the interview did the question come from?
#                 (background, SIPS, or "after" questions)
# respComplexity  Complexity of the response (high/low)
# group           UHR (ultra-high risk) or NCP (non-clinical psychosis)

# Two-way interactions of group and linguistic variables:
# group * questDur
# group * respDur
# group * pauseType
# group * respComplexity

# Two way interaction:
# questDur * respDur

# Three way interaction:
# group * questDur * respDur

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
btpDF$bdiC = myCenter(btpDF$base_BDI_total)
btpDF$sexC = myCenter(btpDF$sex) # female =~ -0.50, male =~ +0.50
btpDF$groupC = myCenter(btpDF$group) # NCP =~ -0.50, UHR =~ +0.50
btpDF$antipsychoticC = myCenter(btpDF$antipsychotic_current) # no = ~ -0.075, yes = ~ 0.925

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

# Final model: once there is a stable model that converges and has no zero-variance components, add the correlations back in to the model. If it converges, run an anova to compare the model with and without the correlations. If there is a significant difference in fit, select the model with correlations. If not, select the model without (prefer a simpler model when there is no evidence for improvement with a more complex model.)

# Trim residuals: once we have a final model structure, we will trim the points in the model with high residuals (> 2.5) and re-fit, to ensure that we are not over-fitting to extreme data points.

# ..................

# A final note: only predictors which vary within a speaker can be used as random slopes for the speaker intercept. Factors such as sex, group, and age do not vary for an individual, and so will not be included in the by-speaker random slopes.

### Model selection --------------------------------------------------

# Model 1: maximal random effects structure
group.lmer1 = lmer(logDur ~ 
                     # -- Fixed Effects -- #
                     # Simple effects: speaker variables
                     ageC + sexC + groupC + bdiC + antipsychoticC +
                     # Simple effects: linguistic variables
                     respDurC + questDurC + respComplexC +
                     bkgdC + sipsAfterC + 
                     # Simple effects: other variables
                     # interviewerC +
                     # Two-way interactions: group * linguistic vars
                     groupC:questDurC + groupC:respDurC +
                     groupC:respComplexC + groupC:bkgdC + groupC:sipsAfterC +
                     # Two-way interaction: quest * resp durations
                     questDurC:respDurC +
                     # Three-way: quest * resp * group 
                     groupC:questDurC:respDurC +
                     ## -- Random Effects -- ##
                     # Only variables that vary within-speaker
                     (1 + respDurC + questDurC + respComplexC + 
                       bkgdC + sipsAfterC +
                        questDurC:respDurC
                      | speaker),
                   data = btpDF,
                   REML = FALSE, # don't use REML when comparing fixed effects
                   verbose = 2) # Track progress in fitting model

# Model did not converge
# Warning messages:
# 1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#   Model failed to converge with max|grad| = 1.0149 (tol = 0.002, component 1)
# 2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#   Model is nearly unidentifiable: large eigenvalue ratio
# - Rescale variables?
# .............

# Model 2: zero correlation model
group.lmer2 = lmer(logDur ~ 
                     ageC + sexC + groupC + bdiC + antipsychoticC +
                     respDurC + questDurC + respComplexC +
                     bkgdC + sipsAfterC + 
                     groupC:questDurC + groupC:respDurC +
                     groupC:respComplexC + groupC:bkgdC + groupC:sipsAfterC +
                     questDurC:respDurC +
                     groupC:questDurC:respDurC +
                     (1 + respDurC + questDurC + respComplexC +
                        bkgdC + sipsAfterC +
                        questDurC:respDurC || speaker), # || = no correlation
                   data = btpDF, REML = FALSE, verbose = 2) 

# Check for zero correlation parameters in the model
summary(rePCA(group.lmer2))
# $speaker
# Importance of components:
#   [,1]   [,2]    [,3]    [,4]    [,5]    [,6] [,7]
# Standard deviation     0.5449 0.3616 0.10361 0.09932 0.09597 0.02143    0
# Proportion of Variance 0.6484 0.2854 0.02344 0.02154 0.02011 0.00100    0
# Cumulative Proportion  0.6484 0.9339 0.95735 0.97889 0.99900 1.00000    1

# Component 7 explains no variance

# Check to see which component this is
summary(group.lmer2, corr = FALSE)
# Random effects:
#   Groups    Name               Variance  Std.Dev.
# speaker   (Intercept)        0.0841423 0.29007 
# speaker.1 respDurC           0.0063490 0.07968 
# speaker.2 questDurC          0.0059279 0.07699 
# speaker.3 respComplexC       0.0069095 0.08312 
# speaker.4 bkgdC              0.0000000 0.00000 
# speaker.5 sipsAfterC         0.1911417 0.43720 
# speaker.6 respDurC:questDurC 0.0002957 0.01720 
# Residual                     0.6436629 0.80229 
# Number of obs: 3481, groups:  speaker, 72


# exclude bkgdC

# ..................

# Model 3: remove the zero variance random slope
group.lmer3 = lmer(logDur ~ 
                     ageC + sexC + groupC + bdiC + antipsychoticC +
                     respDurC + questDurC + respComplexC +
                     bkgdC + sipsAfterC + 
                     groupC:questDurC + groupC:respDurC +
                     groupC:respComplexC + groupC:bkgdC + groupC:sipsAfterC +
                     questDurC:respDurC +
                     groupC:questDurC:respDurC +
                     (1 + respDurC + questDurC + respComplexC + sipsAfterC + 
                        questDurC:respDurC || speaker), # || = no correlation
                   data = btpDF, REML = FALSE, verbose = 2) 

# Check to see if removing this component negatively impacts the model 
anova(group.lmer2, group.lmer3)
# Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# group.lmer3 25 8617.8 8771.7 -4283.9   8567.8                        
# group.lmer2 26 8619.8 8779.8 -4283.9   8567.8     0      1          1
# p = 1, no evidence for a reduction in model explanation as a function of removing bkgdC from the random slopes
# prefer simpler model (lmer3)

# Check to ensure there are no other zero-variance components
summary(rePCA(group.lmer3))

#$speaker
# Importance of components:
#   [,1]   [,2]    [,3]    [,4]    [,5]    [,6]
# Standard deviation     0.5449 0.3616 0.10361 0.09932 0.09597 0.02143
# Proportion of Variance 0.6484 0.2854 0.02344 0.02154 0.02011 0.00100
# Cumulative Proportion  0.6484 0.9339 0.95735 0.97889 0.99900 1.00000


# ................

# Model 4: add correlations back in
group.lmer4 = lmer(logDur ~ 
                     ageC + sexC + groupC + bdiC + antipsychoticC +
                     respDurC + questDurC + respComplexC +
                     bkgdC + sipsAfterC + 
                     groupC:questDurC + groupC:respDurC +
                     groupC:respComplexC + groupC:bkgdC + groupC:sipsAfterC +
                     questDurC:respDurC +
                     groupC:questDurC:respDurC +
                     (1 + respDurC + questDurC + respComplexC + sipsAfterC +
                        questDurC:respDurC | speaker), # | = correlation  
                   data = btpDF,
                   REML = FALSE, verbose = 2,
                   control=lmerControl(optCtrl=list(maxfun=1e06)))


# Compare zero correlation and correlation model
anova(group.lmer3, group.lmer4)
# p = 0.3366
# model4 has no significantly better explanatory power; prefer lmer3

# Final model: lmer3

### Refit final model and assess fixed effects ----------------------

# Remove data points with extreme residuals (> 2.5) in the final model
btpDF.trim = btpDF[abs(scale(resid(group.lmer3)))<2.5,]

# How much data was retained?
nrow(btpDF.trim)/nrow(btpDF)
# 0.987027

# Re-fit model to this trimmed data set
group.lmer.trim =  lmer(logDur ~ 
                          ageC + sexC + groupC + bdiC + antipsychoticC +
                          respDurC + questDurC + respComplexC +
                          bkgdC + sipsAfterC + 
                          groupC:questDurC + groupC:respDurC +
                          groupC:respComplexC + groupC:bkgdC + groupC:sipsAfterC +
                          questDurC:respDurC +
                          groupC:questDurC:respDurC +
                          (1 + respDurC + questDurC + respComplexC + sipsAfterC + 
                             questDurC:respDurC || speaker), # || = no correlation  
                        data = btpDF.trim,
                        REML = FALSE, #verbose = 2,
                        control=lmerControl(optCtrl=list(maxfun=1e06)))

summary(group.lmer.trim)

# AIC      BIC   logLik deviance df.resid 
# 8169.5   8323.0  -4059.7   8119.5     3411 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -2.6601 -0.6606  0.0411  0.6977  2.5234 
# 
# Random effects:
#   Groups    Name               Variance Std.Dev.
# speaker   (Intercept)        0.088736 0.29789 
# speaker.1 respDurC           0.007725 0.08789 
# speaker.2 questDurC          0.008668 0.09310 
# speaker.3 respComplexC       0.013210 0.11493 
# speaker.4 sipsAfterC         0.222706 0.47192 
# speaker.5 respDurC:questDurC 0.000000 0.00000 
# Residual                     0.577043 0.75963 
# Number of obs: 3436, groups:  speaker, 72
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)                6.0003799  0.0411789  145.72
# ageC                      -0.0042584  0.0263343   -0.16
# sexC                       0.1874443  0.0847402    2.21
# groupC                     0.1356252  0.1029983    1.32
# bdiC                      -0.0076660  0.0049842   -1.54
# antipsychoticC             0.5809789  0.2423000    2.40
# respDurC                   0.0538923  0.0213311    2.53
# questDurC                  0.0591083  0.0192698    3.07
# respComplexC              -0.3461357  0.0343227  -10.08
# bkgdC                      0.3220585  0.0517920    6.22
# sipsAfterC                -0.0335824  0.1090884   -0.31
# groupC:questDurC          -0.0292264  0.0383615   -0.76
# groupC:respDurC           -0.0006061  0.0424108   -0.01
# groupC:respComplexC       -0.0188353  0.0687179   -0.27
# groupC:bkgdC               0.1103375  0.1044071    1.06
# groupC:sipsAfterC          0.2949733  0.2215059    1.33
# respDurC:questDurC         0.0161013  0.0161711    1.00
# groupC:respDurC:questDurC -0.0333539  0.0317334   -1.05
# ..........................

## Assess significance of fixed effects

# P-values are not reliable for mixed models and not provided by lmer. In lieu of them, we will use chi-square tests to compare versions of the model with and without a particular predictor of interest. If the model fit significantly improves with the inclusion of the predictor, we consider it to be significant. 

# This function runs the test and provides a summary of the chi-square statistic and significance.

#lme convenience function
chiReport.func <- function(a){
  ifelse (a$"Pr(>Chisq)"[2] > .0001,
          return(paste("chisq(",a$"Chi Df"[2],")=",round(a$Chisq[2],2)," p = ",round(a	$"Pr(>Chisq)"[2],4),sep="")),
          return(paste("chisq(",a$"Chi Df"[2],")=",round(a$Chisq[2],2),", p < .0001")))
}

chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-ageC)))
# "chisq(1)=0.03 p = 0.8716"
 
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-sexC)))
# "chisq(1)=4.72 p = 0.0298"
 
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-groupC)))
# "chisq(1)=1.7 p = 0.1924"
 
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-bdiC)))
# "chisq(1)=2.28 p = 0.1308"

chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-antipsychoticC)))
# "chisq(1)=5.44 p = 0.0197" ***
 
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-respDurC)))
# "chisq(1)=6.18 p = 0.0129" ***
 
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-questDurC)))
# "chisq(1)=8.86 p = 0.0029" ***
 
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-respComplexC)))
# "chisq( 1 )= 58.09 , p < .0001" ***
 
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-bkgdC)))
# "chisq( 1 )= 29.98 , p < .0001" ***
 
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-sipsAfterC)))
# "chisq(1)=0.09 p = 0.7602"
 
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-groupC:questDurC)))
# "chisq(1)=0.58 p = 0.4477"
 
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-groupC:respDurC)))
# "chisq(1)=0 p = 0.9889"
 
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-groupC:respComplexC)))
# "chisq(1)=0.07 p = 0.7856"
 
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-groupC:bkgdC)))
# "chisq(1)=1.07 p = 0.3009"
 
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-groupC:sipsAfterC)))
# "chisq(1)=1.69 p = 0.1933"
 
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-respDurC:questDurC)))
# "chisq(1)=0.99 p = 0.3203"
 
chiReport.func(anova(group.lmer.trim,update(group.lmer.trim,.~.-groupC:respDurC:questDurC)))
# "chisq(1)=1.1 p = 0.2945"

## Visualization: do those on antipsychotics have distinct BTP durations? --------

# Density plot of BTP durations by group, for those not on antipsychotics (N = 70)
# Lines show individual BTP durations for those on antipsychotics (N = 2, both UHR)

ggplot(btpDF[btpDF$antipsychotic_current == "no",],
                   aes(x = logDur)) +
  geom_vline(xintercept=btpDF[btpDF$antipsychotic_current == "yes",]$logDur,
             linetype="dotted", color = "navy")+
  geom_density(aes(group = group, fill = group, color = group), alpha = 0.3) +
  theme_bw() 
