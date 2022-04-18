SEMfinalCSV$NumCrabs <- scale(SEMfinalCSV$logCrabs)
SEMfinalCSV$PneumatophoreCount <- scale(SEMfinalCSV$lnPneumatophore)
SEMfinalCSV$Growing2020 <- scale(SEMfinalCSV$Growing2020)
SEMfinalCSV$KDecomp1 <- scale(SEMfinalCSV$logKDecomp1)




install.packages(c("lavaan", "semPlot", "MPsychoR", "corrplot"))
library(MPsychoR)
library(lavaan)
library(semPlot)
library(corrplot)
data("SEMfinalCSV")
View(CSVSEM)
attach(SEMfinalCSV)
#sample size
nrow(SEMfinalCSV)
#create mean scores
#Model 1 
#Model specification

-------------------------ABGBIO------------------------
  model1<- '
        KDecomp1 ~ NumCrabs
        KDecomp1 ~ PneumatophoreCount
        Growing2020 ~ KDecomp1
        
'
#model estimation using maximum likelihood which is default, turn off means
model1.fit <-sem(model1,
                 data=SEMfinalCSV, 
                 meanstructure = FALSE,
                 estimator = "MLR")
#in summary, regression is just telling us the measure of correlation, want standardized coefficients, fit measures if the model fits
summary(model1.fit, rsquare=TRUE, standardized = TRUE, fit.measures = TRUE)

#visualize path model
semPaths(model1.fit, rotation = 2, layout = "tree2", what = "std", posCol = "black", edge.width = 0.5, style = "Lisrel", fade = T, edge.label.position = 0.55)


------------------Correlations-------------------------------

