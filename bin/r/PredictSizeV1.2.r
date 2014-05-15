#!/usr/bin/Rscript

# arguments extraction
args <- commandArgs(trailingOnly = TRUE)
mBp <-as.numeric(args[1])
sdBp <-as.numeric(args[2])
n <- as.numeric(args[3])


x = n*2

normBp <-rnorm(n, mean=mBp, sd=sdBp)
normBp[normBp<100]=100
print(normBp)
