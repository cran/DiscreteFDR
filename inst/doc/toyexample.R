## ----toy-example-data, results='asis'-----------------------------------------
library(knitr)
X1 <- c(4, 2, 2, 14, 6, 9, 4, 0, 1)
X2 <- c(0, 0, 1, 3, 2, 1, 2, 2, 2)
N1 <- rep(148, 9)
N2 <- rep(132, 9)
Y1 <- N1 - X1
Y2 <- N2 - X2
df <- data.frame(X1, Y1, X2, Y2)
kable(df, caption = "Toy Example")

## ----toy-example-direct-------------------------------------------------------
library(DiscreteFDR)
DBH.sd.fast <- direct.discrete.BH(df, "fisher", direction = "sd")
print(DBH.sd.fast)
summary(DBH.sd.fast)

## ----toy-example-summary, results='asis'--------------------------------------
DBH.sd.fast.summary <- summary(DBH.sd.fast)
summary.table <- DBH.sd.fast.summary$Table
kable(summary.table, caption = "Summary table")

## ----toy-example-generate-fisher----------------------------------------------
# compute results of Fisher's exact test for each row of 'df'
fisher.p <- generate.pvalues(df, "fisher", list(alternative = "two.sided"))
# extract raw observed p-values
raw.pvalues <- fisher.p$get_pvalues()

# or
library(DiscreteTests)
fisher.p.2 <- fisher.test.pv(df, "two.sided")
raw.pvalues.2 <- fisher.p.2$get_pvalues()

# or:
raw.pvalues.3 <- DBH.sd.fast$Data$raw.pvalues

all(raw.pvalues == raw.pvalues.2) && all(raw.pvalues == raw.pvalues.3)
p.adjust(raw.pvalues, method = "BH")

## ----toy-example-crit---------------------------------------------------------
# extract raw observed p-values
raw.pvalues <- fisher.p$get_pvalues()
# extract p-value support sets
pCDFlist <- fisher.p$get_pvalue_supports()

# use raw pvalues and list of supports:
DBH.sd.crit <- DBH(raw.pvalues, pCDFlist, 0.05, "sd", TRUE)
crit.vals.BH.disc <- DBH.sd.crit$Critical.values

# use test results object directly
DBH.sd.crit.2 <- DBH(fisher.p, 0.05, "sd", TRUE)
crit.vals.BH.disc.2 <- DBH.sd.crit.2$Critical.values

# compare
all(crit.vals.BH.disc == crit.vals.BH.disc.2)

## ----toy-example-pipe---------------------------------------------------------
df |>
  fisher.test.pv(alternative = "two.sided") |>
  DBH(0.05, "sd", TRUE) |>
  with(Critical.values)

## ----toy-example-crit-table, results='asis'-----------------------------------
crit.vals.BH.cont <- 1:9 * 0.05/9
tab <- data.frame(raw.pvalues = sort(raw.pvalues), crit.vals.BH.disc, crit.vals.BH.cont)
kable(tab, caption = "Observed p-values and critical values")

## ----toy-example-plot-1-------------------------------------------------------
plot(DBH.sd.crit, col = c("red", "blue", "green"), pch = c(4, 2, 19), lwd = 2, type.crit = 'o',
     legend = "topleft", cex = 1.3)

## ----toy-example-plot-2-------------------------------------------------------
plot(DBH.sd.crit, col = c("red", "blue", "green"), pch = c(4, 2, 19), lwd = 2, type.crit = 'o',
     cex = 1.3, ylim = c(0, 0.25), main = "Comparison of discrete and continuous BH procedures")
points(crit.vals.BH.cont, pch = 19, cex = 1.3, lwd = 2)
legend("topright", c("Rejected", "Accepted", "Critical Values (disc.)", "Critical Values (cont.)"),
       col = c("red", "blue", "green", "black"), pch = c(4, 2, 19, 19), lwd = 2, lty = 0)

## ----toy-example-3------------------------------------------------------------
# extract p-value support sets
pCDFlist <- fisher.p$get_pvalue_supports()
# extract 1st and 5th support set
pCDFlist[c(1,5)]

## ----toy-example-4------------------------------------------------------------
stepf <- lapply(pCDFlist, function(x) stepfun(x, c(0, x)))
par(mfcol = c(1, 3), mai = c(1, 0.5, 0.3, 0.1))
plot(stepf[[1]], xlim = c(0, 1), ylim = c(0, 1), do.points = FALSE, lwd = 1, lty = 1, ylab = "F(x)", 
     main = "(a)")
for(i in (2:9)){
  plot(stepf[[i]], add = TRUE, do.points = FALSE, lwd = 1, col = i)
}
segments(0, 0, 1, 1, col = "grey", lty = 2)

#   Plot xi
support <- sort(unique(unlist(pCDFlist)))
components <- lapply(stepf, function(s){s(support) / (1 - s(support))}) 
xi.values <- 1/9 * Reduce('+', components)
xi <- stepfun(support, c(0, xi.values))
plot(xi, xlim = c(0, 0.10), ylim = c(0, 0.10), do.points = FALSE, ylab = expression(xi), main = "(b)")
segments(0, 0, 0.1, 0.1, col = "grey", lty = 2)

#   Plot discrete critical values as well a BH constants and raw p-values
DBH.sd <- DBH(fisher.p, direction = "sd", ret.crit.consts = TRUE)
plot(DBH.sd, col = c("black", "black", "red"), pch = c(4, 4, 19), type.crit = 'p', ylim = c(0, 0.15),
     cex = 1.3, main = "(c)", ylab = "Critical Values")
points(1:9, 0.05 * (1:9) / 9, col = "green", pch = 19, cex = 1.3)

mtext("Figure 1", 1, outer = TRUE, line = -2)

