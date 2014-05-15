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

# take 4 file in line commander
args <- commandArgs(trailingOnly = TRUE)


# create of the image to stocke the result
# name: result.png
# size: 1100of width and 1000 of height
png("result.png", width=1200, height=1300)

# allow to fix the disposition of image, size of marge, size of label and title
par(mfrow=c(3,4),mar=c(5,5,4,1),cex.axis=2.2, cex.lab=2.2, cex.main=2.2)

#use of a library to have nice color
library("RColorBrewer")

mypalette1=brewer.pal(8,"Set1")

mypalette2=brewer.pal(8,"Set2")

# stocke 4 file in 4 table
A <- read.table(args[1], header=T, sep="\t")
B <- read.table(args[2], header=T, sep="\t")
C <- read.table(args[3], header=T, sep="\t")
D <- read.table(args[4], header=T, sep="\t")
E <- read.table(args[5], header=T, sep="\t")
F <- read.table(args[6], header=T, sep="\t")
G <- read.table(args[7], header=T, sep="\t")
H <- read.table(args[8], header=T, sep="\t")


# create of 4 boxplot for the sensibility, specificity, the precision, and cutoff of accuracy and performance
# notch do a student test to known is the difference is significative

boxplot(A$Sensibility, B$Sensibility,C$Sensibility,D$Sensibility, col=mypalette1,ylab="cutoff value distribution", xlab="Sensibility, PredictV1", ylim=c(0.95,0.965), notch=T, main = "evolve of the sensibility")

boxplot(E$Sensibility, F$Sensibility,G$Sensibility,H$Sensibility, col=mypalette1,ylab="cutoff value distribution", xlab="Sensibility, PredictV2", ylim=c(0.95,0.965), notch=T, main = "evolve of the sensibility")

boxplot(A$Specificity, B$Specificity,C$Specificity,D$Specificity, col=mypalette1,ylab="cutoff value distribution", xlab="Specificity, PredictV1",ylim=c(0.955,0.965),main = "evolve of the specificity", notch=T)

boxplot(E$Specificity, F$Specificity,G$Specificity,H$Specificity, col=mypalette1,ylab="cutoff value distribution", xlab="Specificity, PredictV2",ylim=c(0.955,0.965),main = "evolve of the specificity", notch=T)

boxplot(A$Precision, B$Precision,C$Precision,D$Precision, col=mypalette1,ylab="cutoff value distribution", xlab="Precision, PredictV1",ylim=c(0.955,0.965), main = "evolve of the precision",notch= T)

boxplot(E$Precision, F$Precision,G$Precision,H$Precision, col=mypalette1,ylab="cutoff value distribution", xlab="Precision, PredictV2",ylim=c(0.955,0.965), main = "evolve of the precision",notch= T)

boxplot(A$Accuracy, B$Accuracy,C$Accuracy,D$Accuracy, col=mypalette1,ylab="cutoff value distribution", xlab="Accuracy, PredictV1", ylim=c(0.955,0.965), notch= T,main = "evolve of the cutoff calculate")

boxplot(E$Accuracy, F$Accuracy,G$Accuracy,H$Accuracy, col=mypalette1,ylab="cutoff value distribution", xlab="Accuracy, PredictV2", ylim=c(0.955,0.965), notch= T,main = "evolve of the accuracy")

boxplot(A$cuttof_performance, B$cuttof_performance,C$cuttof_performance,D$cuttof_performance, col=mypalette1,ylab="cutoff value distribution", xlab="performance, PredictV1", ylim=c(0.40,0.60), main = "evolve of the accuracy", notch=T)

boxplot(E$cuttof_performance, F$cuttof_performance, G$cuttof_performance, 	H$cuttof_performance, col=mypalette1,ylab="cutoff value distribution", xlab="performance, PredictV2", ylim=c(0.40,0.60), main = "evolve of the cutoff of performance", notch=T)

# create an empty plot to put extra legend to explain what we do and the correspondance of each colors
# box.lwd=0 allow to eliminate all border
# plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10), axes = F)
# legend("top",fill=mypalette1,lwd=2,legend=c("size of 50","size of 100","size of 200","ORF size"), cex=2)
# legend("center", legend= "In these graphic we can see the evolve\n of several statisctic test in term of:\n\n- if we use the fonction of predict size with:\n   + 50 like minimal length\n   + 100 like minimal length\n   + 200 like minimal length\n\n- or if we use ORF sequences to predict\n the size of pseudo-lncRNA", cex=2, box.lwd =0)

dev.off()
