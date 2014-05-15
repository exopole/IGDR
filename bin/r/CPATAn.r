#!/local/R/3.0.1/bin/Rscript

#
# program: CPATAn.r
# 
# aim: Analysis of the CPAT result from several test.
#
# input:
#	several file of table that contain result
#
# output:
#	a plot that let us to compare result
#
###############################################################


args <- commandArgs(trailingOnly = TRUE)
len <- length(args)


png("result.png")
par(mfrow=c(3,2),mar=c(5,4,2,2),cex.axis=1.2, cex.lab=1.2)



A <- read.table(args[1], header=T, sep="\t")
B <- read.table(args[2], header=T, sep="\t")
C <- read.table(args[3], header=T, sep="\t")
D <- read.table(args[4], header=T, sep="\t")
# E <- read.table(args[5], header=T, sep="\t")
# F <- read.table(args[6], header=T, sep="\t")

x= 1:length(A$Sensibility)
print(A$Sensibility)

plot(x, sort(A$Sensibility), type ='l', col="red",ylab="Sensibility-specificity", xlab="shuffle", ylim = c(0.90,1))
# lines(x, sort(B$Sensibility), type ='l',lty=2, col="red")
# lines(x, sort(C$Sensibility), type ='l',lty=1, col="green")
# lines(x, sort(D$Sensibility), type ='l',lty=2, col="green")
# lines(x, sort(E$Sensibility), type ='l',lty=1, col="blue")
# lines(x, sort(F$Sensibility), type ='l',lty=2, col="blue")


plot(x, sort(A$Specificity), type ='l', col="red",ylab="Specificity",xlab="shuffle", ylim = c(0,1))
# lines(x, sort(A$Specificity), type ='l', lty=2,col="blue")
# lines(x, sort(B$Specificity), type ='l',lty=2, col="red")
# lines(x, sort(C$Specificity), type ='l',lty=1, col="green")
# lines(x, sort(D$Specificity), type ='l',lty=2, col="green")
# lines(x, sort(E$Specificity), type ='l',lty=1, col="blue")
# lines(x, sort(F$Specificity), type ='l',lty=2, col="blue")k
# legend(0.5,0.7,col=c("red","blue"),lwd=2,legend=c("sensibility","specificity"))


plot(x, sort(A$Precision), type ='l', col="red", ylab="Precision",xlab="shuffle", ylim = c(0.90,1))
# lines(x, sort(B$Precision), type ='l',lty=2, col="red")
# lines(x, sort(C$Precision), type ='l',lty=1, col="green")
# lines(x, sort(D$Precision), type ='l',lty=2, col="green")
# lines(x, sort(E$Precision), type ='l',lty=1, col="blue")
# lines(x, sort(F$Precision), type ='l',lty=2, col="blue")



plot(x, sort(A$Accuracy), type ='l', col="red", ylab="Accuracy",xlab="shuffle", ylim = c(0.90,1))
# lines(x, sort(B$Accuracy), type ='l',lty=2, col="red")
# lines(x, sort(C$Accuracy), type ='l',lty=1, col="green")
# lines(x, sort(D$Accuracy), type ='l',lty=2, col="green")
# lines(x, sort(E$Accuracy), type ='l',lty=1, col="blue")
# lines(x, sort(F$Accuracy), type ='l',lty=2, col="blue")

boxplot(A$cuttof_performance, A$cuttof_accuracy, col=c("green", "red"),ylab="cutoff value distribution", xlab="performance & accuracy 10x")


dev.off()
