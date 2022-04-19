SEMfinalCSV$NumCrabs <- scale(lnCrabs)
SEMfinalCSV$PneumatophoreCount <- scale(lnpneu)
SEMfinalCSV$Growing2020 <- scale(SEMfinalCSV$Growing2020)
SEMfinalCSV$KDecomp1 <- scale(lnKDecomp)




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
?sem
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


------------------Correlations------------------------------
  rma<-function(x,y) {
    
    # Compile the bivariate data into a workable format.
    data <- as.data.frame(cbind(x,y))
    colnames(data) <- c("x","y")
    
    # Calculate sample size, mean of x, mean of y, SD of x, SD of y, correltation between x and y.
    n <- length(data$x)
    mx <- mean(data$x)
    my <- mean(data$y)
    sx <- sqrt(var(data$x))
    sy <- sqrt(var(data$y))
    r <- cor(data$x,data$y)
    
    # Calculate the RMA slope, then the RMA intercept, combine into object.
    slp <- sy/sx
    int <- my-(slp*mx)
    rma <- c(int,slp)
    names(rma) <- c("int","slp")
    
    # Calculate slope and intercept SEs and add to RMA object.
    # Not needed for current application. Not done.
    
    # Calculated predicted values for x and y based on rma line.
    x_hat <- rep(0,n)
    y_hat <- rep(0,n)
    for (k in 1:n)
    {x_hat[k] <- (data$y[k]-int)/slp
    y_hat[k] <- int+(slp*data$x[k])}
    
    # Calculate the legs of the triangle formed by the RMA line and each point - x and y distances.
    x_dist <- rep(0,n)
    y_dist <- rep(0,n)
    for (k in 1:n)
    {x_dist[k]<-abs(data$x[k]-x_hat[k])
    y_dist[k]<-abs(data$y[k]-y_hat[k])}
    
    # Calculate the hypothenuse of the triangle, which is the length of the RMA line that is part of the triangle.
    h_dist <- rep(0,n)
    for (k in 1:n)
    {h_dist[k]<-sqrt((x_dist[k]^2)+(y_dist[k]^2))}
    
    # Make a vector that tells you whether each point is above or below the RMA line.
    sign <- rep(0,n)
    for (k in 1:n)
    {ifelse(data$y[k] < y_hat[k],sign[k] <- -1,sign[k] <- 1)}
    
    # Calculate classic RMA residuals, as the area of the triangle formed by the point and the line and apply the correct sign.
    # Also calculate the square roots of the RMA residuals so that the units are the same as other residuals.
    rma_resid <- rep(0,n)
    rma_resid_sqrt <- rep(0,n)
    for (k in 1:n)
    {rma_resid[k] <- 0.5*x_dist[k]*y_dist[k]*sign[k]
    rma_resid_sqrt[k] <- sqrt(0.5*x_dist[k]*y_dist[k])*sign[k]}
    
    # Calculate residuals as for ma regression, as the perpendicular distance from the RMA line to each point.
    ma_resid <- rep(0,n)
    for (k in 1:n)
    {ma_resid[k] <- y_dist[k]*x_dist[k]/h_dist[k]*sign[k]}
    
    # Test for normality of each set of residuals using a Shapiro-Wilks test and a Kolmogorov-Smirnov test.
    sw_rma_res <- shapiro.test(rma_resid)
    sw_rma_res_sqrt <- shapiro.test(rma_resid_sqrt)
    sw_ma_res <- shapiro.test(ma_resid)
    ks_rma_res <- ks.test(rma_resid,"pnorm",mean=mean(rma_resid),sd=sqrt(var(rma_resid)))
    ks_rma_res_sqrt <- ks.test(rma_resid_sqrt,"pnorm",mean=mean(rma_resid_sqrt),sd=sqrt(var(rma_resid_sqrt)))
    ks_ma_res <- ks.test(ma_resid,"pnorm",mean=mean(ma_resid),sd=sqrt(var(ma_resid)))
    
    # Compile results into a list that contains the x y data and residuals, the rma coefficients, and the residual test results.
    # First compile the normality test results
    sw_rma <- c(sw_rma_res$stat,sw_rma_res$p)
    sw_rma_sqrt <- c(sw_rma_res_sqrt$stat,sw_rma_res_sqrt$p)
    sw_ma <- c(sw_ma_res$stat,sw_ma_res$p)
    ks_rma <- c(ks_rma_res$stat,ks_rma_res$p)
    ks_rma_sqrt <- c(ks_rma_res_sqrt$stat,ks_rma_res_sqrt$p)
    ks_ma <- c(ks_ma_res$stat,ks_ma_res$p)
    res_norm_tests <- rbind(sw_rma,sw_rma_sqrt,sw_ma,ks_rma,ks_rma_sqrt,ks_ma)
    colnames(res_norm_tests) <- c("stat","p")
    
    # Second compile that data and residuals.
    data_res <- cbind(data,rma_resid,rma_resid_sqrt,ma_resid)
    
    # Finally put it all together.
    rma_results <- list(data_and_resids=data_res,correlation=r,rma_coefs=rma,resid_norm_tests=res_norm_tests)
    rma_results
  }




logRedoxx<-log10(abs(SEMfinalCSV$Redox1))
shapiro.test(logRedoxx)
lnRedoxx<-log(SEMfinalCSV$Redox1+1)
shapiro.test(lnRedoxx)
SqRedoxx<-sqrt(abs(SEMfinalCSV$Redox1))
shapiro.test(SqRedoxx)
ArcRedoxx<-asin(abs(SEMfinalCSV$Redox1))

cor.test(lnCrabs,SEMfinalCSV$Redox1, method="spearman", exact=FALSE)

shapiro.test(SEMfinalCSV$KDecomp1)
logKDecomp<-log10(SEMfinalCSV$KDecomp1)
shapiro.test(logKDecomp)
cor.test(logNumCrabs, logKDecomp, method="pearson")


logNAb<-log10(SEMfinalCSV$ActualAbsorb+1)
shapiro.test(logNAb)
lnNab<-log(SEMfinalCSV$ActualAbsorb)
shapiro.test(lnNab)


lnCrabs<-log(SEMfinalCSV$NumCrabs)
shapiro.test(lnCrabs)


cor.test(lnCrabs, lnKDecomp, method="pearson")
lnKDecomp<-log(SEMfinalCSV$KDecomp1)
shapiro.test(lnKDecomp)

cor.test(lnCrabs, lnNab, method="pearson")

shapiro.test(SEMfinalCSV$ProductivityWETFEETWay)
lnProd<-log(SEMfinalCSV$ProductivityWETFEETWay)
shapiro.test(lnProd)
cor.test(lnCrabs, lnProd, method="pearson")

cor.test(lnCrabs, SEMfinalCSV$Growing2020, method="pearson")

cor.test(SEMfinalCSV$Growing2020, lnKDecomp, method="pearson")


#pneumatophore and fiddler crab correlatiom
shapiro.test(SEMfinalCSV$PneumatophoreCount)
lnpneu<-log(SEMfinalCSV$PneumatophoreCount)
shapiro.test(lnpneu)
shapiro.test(SEMfinalCSV$lnPneumatophore)
logpneu<-sqrt(SEMfinalCSV$PneumatophoreCount)
shapiro.test(logpneu)


cor.test(lnCrabs,SEMfinalCSV$PneumatophoreCount, method="spearman", exact=FALSE)
cor.test(SEMfinalCSV$NumCrabs,SEMfinalCSV$PneumatophoreCount, method="spearman", exact=FALSE)

cor.test(lnKDecomp,SEMfinalCSV$PneumatophoreCount, method="spearman", exact=FALSE)



cor.test(SEMfinalCSV$Growing2020,SEMfinalCSV$Redox1, method="spearman", exact=FALSE)
