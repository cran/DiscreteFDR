## ----format-am----------------------------------------------------------------
library(DiscreteFDR)
library(DiscreteDatasets)
data(amnesia)
amnesia.formatted <- generate.pvalues(amnesia, "fisher", list(alternative = "greater"), reconstruct_two)

## -----------------------------------------------------------------------------
library(DiscreteDatasets)
library(DiscreteTests)
amnesia.formatted <- amnesia |>
  reconstruct_two() |>
  fisher.test.pv(alternative = "greater")

## ----DFDR-Pharmacovigilance---------------------------------------------------
DBH.su  <-  DBH(amnesia.formatted, ret.crit.consts = TRUE)
DBH.sd  <-  DBH(amnesia.formatted, direction = "sd", ret.crit.consts = TRUE)
ADBH.su <- ADBH(amnesia.formatted, ret.crit.consts = TRUE)
ADBH.sd <- ADBH(amnesia.formatted, direction = "sd", ret.crit.consts = TRUE)
DBR     <-  DBR(amnesia.formatted, ret.crit.consts = TRUE)

## ----Hist-Pharmaco------------------------------------------------------------
hist(DBH.sd)

## ----Plot-Pharmaco------------------------------------------------------------
raw.pvalues <- amnesia.formatted$get_pvalues()
m <- length(raw.pvalues)
crit.values.BH <- 0.05 * seq_len(m) / m
scale.points <- 0.7

plot(DBH.su, col = c("black", "black", "orange"), pch = NA, type.crit = 'p', xlim = c(1, 100),
     ylim = c(0, DBH.su$Critical.values[100]), ylab = "critical values", cex = scale.points, main = "")

points(crit.values.BH[1:105],          col = "green",  pch = 19, cex = scale.points)
points(DBH.sd$Critical.values[1:105],  col = "red",    pch = 19, cex = scale.points)
points(ADBH.su$Critical.values[1:105], col = "blue",   pch = 19, cex = scale.points)
points(ADBH.sd$Critical.values[1:105], col = "purple", pch = 19, cex = scale.points)
points(DBR$Critical.values[1:105],     col = "yellow", pch = 19, cex = scale.points)
points(sort(raw.pvalues),                              pch = 4,  cex = scale.points)
mtext("Figure 2", 1, outer = TRUE, line = -1)

## ----Reject-Pharmaco----------------------------------------------------------
rej.BH      <- length(which(p.adjust(raw.pvalues, method = "BH") <= 0.05))
rej.DBH.su  <- length(DBH.su$Indices)
rej.DBH.sd  <- length(DBH.sd$Indices)
rej.ADBH.su <- length(ADBH.su$Indices)
rej.ADBH.sd <- length(ADBH.sd$Indices)
rej.DBR     <- length(DBR$Indices)
c(rej.BH, rej.DBH.su, rej.DBH.sd, rej.ADBH.su, rej.ADBH.sd, rej.DBR)

## ----Poisson-Setup------------------------------------------------------------
lambda.vector <- c(0.6, 1.2, 0.7, 1.3, 1.0, 0.2, 0.8, 1.3, 0.9)
observations <- c(3, 3, 1, 2, 3, 3, 1, 2, 4)
configuration <- cbind(observations, lambda.vector)
alpha <- 0.05
m <- length(observations)
print(configuration)

## ----Poisson-RawPValues-------------------------------------------------------
raw.pvalues <- ppois(observations - 1, lambda.vector, lower.tail = FALSE)
poisson.p <- poisson.test.pv(observations, lambda.vector, "greater")
raw.pvalues.2 <- poisson.p$get_pvalues()
print(raw.pvalues.2)

## ----Poisson-tauMin-----------------------------------------------------------
y.min <- alpha/m * (1 + alpha/m)^(-1)
n.max <- qpois(y.min,     lambda.vector, lower.tail = FALSE) + 1
t.min <- ppois(n.max - 1, lambda.vector, lower.tail = FALSE)
s.min <- min(t.min)
print(s.min)

## -----------------------------------------------------------------------------
sapply(poisson.p$get_pvalue_supports(), min)

## ----Poisson-Supports---------------------------------------------------------
supports <- lapply(1:m, function(w){sort(ppois(0:n.max[w] - 1, lambda.vector[w], lower.tail = FALSE))})
DBH.sd <- DBH(raw.pvalues, supports, direction = "sd", ret.crit.consts = TRUE)
print(DBH.sd)

## ----Poisson-Supports-with-object---------------------------------------------
DBH.sd.2 <- DBH(poisson.p, direction = "sd", ret.crit.consts = TRUE)
print(DBH.sd.2)

## ----Poisson-BH---------------------------------------------------------------
p.adjust(raw.pvalues, method = "BH")

## ----Poisson-DBH--------------------------------------------------------------
DBH.sd$Adjusted

## ----Poisson-Print------------------------------------------------------------
print(DBH.sd)
summary(DBH.sd)

## ----Poisson-Rejections-------------------------------------------------------
stepf <- lapply(supports, function(x) stepfun(x, c(0, x)))
par(mfcol = c(1, 3), mai = c(1, 0.5, 0.3, 0.1))

plot(stepf[[1]], xlim = c(0,1), ylim = c(0,1), do.points = FALSE, lwd = 1, lty = 1, ylab = "F(x)", 
     main = "(a)")
for(i in (2:9)){
  plot(stepf[[i]], add = TRUE, do.points = FALSE, lwd = 1, col = i)
}
segments(0, 0, 1, 1, col = "grey", lty = 2)

#   Plot xi
support <- sort(unique(unlist(supports)))
components <- lapply(stepf, function(s){s(support) / (1 - s(support))}) 
xi.values <- 1/9 * Reduce('+', components)
xi <- stepfun(support, c(0, xi.values))
plot(xi, xlim = c(0, 0.10), ylim = c(0, 0.10), do.points = FALSE, ylab = expression(xi), main = "(b)")
segments(0, 0, 0.1, 0.1, col = "grey", lty = 2)

#   Plot discrete critical values as well a BH constants
DBH.sd <- DBH(raw.pvalues, supports, direction = "sd", ret.crit.consts = TRUE)
plot(DBH.sd, col = c("black", "black", "red"), pch = c(4, 4, 19), type.crit = 'p', ylim = c(0, 0.15),
     cex = 1.3, main = "(c)", ylab = "Critical Values")
points(1:9, 0.05 * (1:9) / 9, col = "green", pch = 19, cex = 1.3)

mtext("Figure 3", 1, outer = TRUE, line = -2)

## -----------------------------------------------------------------------------
library(DiscreteDatasets)
library(DiscreteFDR)

lambda.0 <- disorderdetection[,2]
lambda.vector <- lambda.0
observations <- disorderdetection[,1]
configuration <- cbind(observations, lambda.0)
alpha <- 0.05
m <- length(observations)

raw.pvalues <- sapply(1:m,function(i){ppois(observations[i] , lambda.vector[i], lower.tail = TRUE)})

pvalues_a3=p.adjust(raw.pvalues, method="BH",n=m)
BH_NumRejected <- length(pvalues_a3[pvalues_a3<0.05])
BH_NumRejected

#pvalues_a=SB_SU(raw.pvalues,hu,gg,alt="del")#AHSU correction
#pvalues_a$Num.rejected


## -----------------------------------------------------------------------------

n.max <- qpois(2^-1074, lambda.vector, FALSE)

supports <- lapply(1:m, function(i){sort(ppois(0:n.max[i], lambda.vector[i]))})

## ----other-procedures---------------------------------------------------------
# rc = real case

library(DiscreteFDR)
raw.pvalues.rc <- raw.pvalues
pCDFlist.rc <- supports

DBH.su.rc  <-  DBH(raw.pvalues.rc, pCDFlist.rc, ret.crit.consts = TRUE)
DBH.sd.rc  <-  DBH(raw.pvalues.rc, pCDFlist.rc, direction = "sd", ret.crit.consts = TRUE)
ADBH.su.rc <- ADBH(raw.pvalues.rc, pCDFlist.rc, ret.crit.consts = TRUE)
ADBH.sd.rc <- ADBH(raw.pvalues.rc, pCDFlist.rc, direction = "sd", ret.crit.consts = TRUE)
DBR.su.rc  <-  DBR(raw.pvalues.rc, pCDFlist.rc, ret.crit.consts = TRUE)

## ----plot---------------------------------------------------------------------


num.plot<-60
index.plot<-1:num.plot
m <- length(raw.pvalues.rc)
crit.values.BH <- 0.05 * (1:m) / m

plot(index.plot,sort(raw.pvalues.rc)[index.plot],pch=3,lwd=1.5,ylim=c(0,0.02),ylab="critical values / p-values",col="black")
points(index.plot, ADBH.su.rc$Critical.values[index.plot],pch = 19, col="blue")
points(crit.values.BH[1:105],             col = "green",  pch = 19)
points(DBH.sd.rc$Critical.values[1:105],  col = "red",   pch = 19)
points(ADBH.sd.rc$Critical.values[1:105], col = "purple", pch = 19)
points(DBR.su.rc$Critical.values[1:105],  col = "cyan", pch = 19)

legend("topleft",inset=0.05, legend=c("raw p-values (sorted)",   "[A-DBH-SU]","[A-DBH-SD]","[DBH-SU]","[DBH-SD]", "[DBR-SU]","BH"),
       col=c("black","blue","purple","orange", "red", "cyan", "green"), pch=c(3, 19, 19, 19, 19, 19, 19), cex=0.75)
title(expression("NGS Disorder Detection: critical values for " ~  "FDR" <=0.05),line=1)


## ----table--------------------------------------------------------------------
df.rc <- c("FDR<=0.05")

table_rc <- cbind("FDR <= 0.05",ADBH.su.rc$Num.rejected, ADBH.sd.rc$Num.rejected, DBH.su.rc$Num.rejected, DBH.sd.rc$Num.rejected, DBR.su.rc$Num.rejected, BH_NumRejected )


## ----rc-table, echo=FALSE-----------------------------------------------------
library(knitr)
library(kableExtra)
kable(table_rc, col.names = c("procedure controls", "[A-DBH-SU]", "[A-DBH-SD]", "[DBH-SU]", "[DBH-SD]", "[DBR-SU]", "[BH]" ),caption = "*Number of rejections for the disorder detection data.*")%>%
  kable_styling(latex_options = "striped")


