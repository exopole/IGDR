#!/usr/bin/Rscript


args <- commandArgs(trailingOnly = TRUE)
monfichier <-args[1]

dat = read.table(monfichier, header=T) 
summary(dat$coding_prob)
image=paste(monfichier,"png", sep=".")
png(image)
boxplot(dat$coding_prob,outline =FALSE)
dev.off()
