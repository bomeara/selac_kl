setwd("/Users/bomeara/Documents/MyDocuments/GitClones/selac_kl/scripts")
load("KLrun2.rda")
library(MuMIn)
results$neg2lnL <- results$AIC - 2
results$logNtax <- log(results$ntax)
results$logNchar <- log(results$nchar)
results$NcharTimesNtax <- results$ntax * results$nchar
results$LogNcharTimesNtax <- results$ntax * log(results$nchar)
results$NcharTimesLogNtax <- log(results$ntax) * results$nchar

K=1
#results$AICc_ntax <- results$neg2lnL + 2 * K * (results$ntax / ( results$ntax - K - 1))
#results$AICc_nchar <- results$neg2lnL + 2 * K * (results$nchar / ( results$nchar - K - 1))
#results$AICc_ncharTimesntax <- results$neg2lnL + 2 * K * (results$nchar*results$ntax / ( results$nchar*results$ntax - K - 1))
results$AICc_ntax_penalty <-  2 * K * (results$ntax / ( results$ntax - K - 1))
results$AICc_nchar_penalty <- 2 * K * (results$nchar / ( results$nchar - K - 1))
results$AICc_nchar_penalty[which(!is.finite(results$AICc_nchar_penalty ))] <- 1e6
results$AICc_ncharTimesntax_penalty <-  2 * K * (results$nchar*results$ntax / ( results$nchar*results$ntax - K - 1))
results$KLminusNeg2lnL <- results$KL - results$neg2lnL


results.for.model <- results
results.for.model$rep <- NULL
results.for.model$AIC <- NULL
results.for.model$K <- rep(1, nrow(results.for.model))
results.for.model$KL <- NULL
results.for.model$neg2lnL <- NULL


options(na.action = "na.fail")
global.model <- lm(KLminusNeg2lnL ~ ., data=results.for.model)
dredged.model <- MuMIn::dredge(global.model, beta="sd", m.lim=c(1,2))
pdf(file="BestModels.pdf")
#plot(results$AIC, results$KLminusNeg2lnL, type="p", bty="n", pch=16, col=rgb(1,0,0,.5), log="xy",asp=1, xlab="estimator", ylab="KL")

best.model <- get.models(dredged.model, 1)[[1]]
best.model.predictions <- predict(best.model)
plot(best.model.predictions, best.model$model$KLminusNeg2lnL, pch=16, col=rgb(0,0,0,.5))
points(results$AICc_nchar_penalty, best.model$model$KLminusNeg2lnL, pch=16, col=rgb(1,0,0,.5))
abline(a=0, b=1)
abline(h=2, col="red")
dev.off()
write.csv(dredged.model, file="AllModels.csv")
# all.tests <- list(
# cor.test.AIC = cor.test(results$KL, results$AIC),
# cor.test.AICc_ntax = cor.test(results$KL, results$AICc_ntax),
# cor.test.AICc_nchar =  cor.test(results$KL, results$AICc_nchar),
# cor.test.best.model = cor.test(best.model$model$KL, predict(best.model))
# )

#print(lapply(all.tests, '[[', 'p.value'))
print(best.model)
