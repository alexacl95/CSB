library('fda.usc')
setwd("~/MATLAB-Drive/PaperCSB/Repositorio/FunctionalIntervals")
Estimations <- read.csv(file = 'DengueEstimations.csv')
fdaEst = fdata(Estimations, argvals = NULL, rangeval = NULL, names = NULL, fdata2d = FALSE)

out.boot=fdata.bootstrap(fdaEst,statistic=func.med.FM,nb=1000,alpha=0.05,draw=TRUE)



################################################################
## Paquetes requeridos
my_packages = c("expm","gridExtra","tidyverse","knitr","mrfDepth","e1071","pracma","MASS","mixtools","rrcov","Rfast","depthTools", "mvoutlier","fds")
not_installed = my_packages[!(my_packages %in% installed.packages()[ , "Package"])]
if (length(not_installed)) install.packages(not_installed, dependencies = TRUE)
for (q in 1:length(my_packages)) {
  library(my_packages[q], character.only = TRUE)
}
## Carga de datos
#data(Octanespectrum)
setwd("~/MATLAB-Drive/PaperCSB/Repositorio/FunctionalIntervals")
Estimations <- read.csv(file = 'DengueEstimations.csv', header = FALSE)
## Computo de observación más profunda e IQR funcional
wdbc.b.z =matrix(unlist(Estimations),ncol =13, nrow =1000)
#wdbc.b.z = scale(wdbc.b.z)
fda.wdbc = depthTools::MBD(wdbc.b.z, band = T, plotting = F)
p.fda.wdbc = fda.wdbc$ordering[1:ceiling(0.5*nrow(wdbc.b.z))]
p2.fda.wdbc = fda.wdbc$MBD[p.fda.wdbc]
#d50.wdbc.b.z = wdbc.b.z[p.fda.wdbc,]
d50.wdbc.b.z =matrix(unlist(out.boot$resample$data),ncol =13, nrow =1000)
max.d50.wdbc.b.z = apply(d50.wdbc.b.z, 2, max)
min.d50.wdbc.b.z = apply(d50.wdbc.b.z, 2, min)
# Método minimum diagonal product estimator
W1 = as.numeric(rmdp(wdbc.b.z, itertime = 500)$wei)
W1[which(W1 == 1)] = 2
W1[which(W1 == 0)] = 1
W1[which(W1 == 2)] = 0
Index2 = which(W1 == 1)
# Método minimum regularized covariance estimator
W2 = as.numeric(getFlag(CovMrcd(wdbc.b.z)))
W2[which(W2 == 1)] = 2
W2[which(W2 == 0)] = 1
W2[which(W2 == 2)] = 0
Index3 = which(W2 == 1)
# Gráfico 1
Index = Index2
wdbc.outl = wdbc.b.z[Index,]
wdbc.no.outl = wdbc.b.z[-Index,]
par(mar = c(5.1, 4.1, 4.1, 10), xpd = TRUE)
plot(wdbc.outl[1,], col = "gray", lwd = 2, ylim = c(min(wdbc.outl),max(wdbc.outl)),
     type = "l", xlim = c(0, ncol(wdbc.b.z)),
     main = "MDP Multivariate Outlier Detection Parallel Coordinates: OctaneSpectrum Data",
     xlab = "Dimension (p)",
     ylab = "Observations (n)")
for (i in 2:nrow(wdbc.outl)) {
  lines(wdbc.outl[i,], col = "gray", lwd = 2)
}
for (i in 1:nrow(wdbc.no.outl)) {
  lines(wdbc.no.outl[i,], col = "black", lwd = 2)
}
#Curva mas profunda
lines(wdbc.b.z[p.fda.wdbc[1],], col = "cyan", lwd = 4)
lines(max.d50.wdbc.b.z, col = "deeppink", lwd = 4)
lines(min.d50.wdbc.b.z, col = "deeppink", lwd = 4)
legend("topright", inset = c(-0.22, 0),
       legend = c("Outliers","Non-Outliers", "Deepest\nObservation","f-IQR Borders"),
       col = c("gray","black","cyan","deeppink"),
       lty = 1, lwd = 4, box.lty = 1, horiz = FALSE, text.font = 1)
## Gráfico 2
Index = Index3
wdbc.outl = wdbc.b.z[Index,]
wdbc.no.outl = wdbc.b.z[-Index,]
par(mar = c(5.1, 4.1, 4.1, 10), xpd = TRUE)
plot(wdbc.outl[1,], col = "gray", lwd = 2, ylim = c(min(wdbc.outl),max(wdbc.outl)),
     type = "l", xlim = c(0, ncol(wdbc.b.z)),
     main = "MRCD Multivariate Outlier Detection Parallel Coordinates: OctaneSpectrum Data",
     xlab = "Dimension (p)",
     ylab = "Observations (n)")
for (i in 2:nrow(wdbc.outl)) {
  lines(wdbc.outl[i,], col = "gray", lwd = 2)
}
for (i in 1:nrow(wdbc.no.outl)) {
  lines(wdbc.no.outl[i,], col = "black", lwd = 2)
}
lines(wdbc.b.z[p.fda.wdbc[1],], col = "cyan", lwd = 4)
lines(max.d50.wdbc.b.z, col = "deeppink", lwd = 4)
lines(min.d50.wdbc.b.z, col = "deeppink", lwd = 4)
legend("topright", inset = c(-0.22, 0),
       legend = c("Outliers","Non-Outliers", "Deepest\nObservation","f-IQR Borders"),
       col = c("gray","black","cyan","deeppink"),
       lty = 1, lwd = 4, box.lty = 1, horiz = FALSE, text.font = 1)


DataSave = cbind(wdbc.b.z[p.fda.wdbc[1],],min.d50.wdbc.b.z,max.d50.wdbc.b.z)
write.csv(DataSave, "DengueIntervals.csv", row.names = FALSE)