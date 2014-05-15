#!/local/R/3.0.1/bin/Rscript

# arguments extraction
args <- commandArgs(trailingOnly = TRUE)
mBp <-as.numeric(args[1])
sdBp <-as.numeric(args[2])
n <- as.numeric(args[3])
min <- as.numeric(args[5])
max <- as.numeric(args[6])

normBp <-rnorm(n, mean=mBp, sd=sdBp)
normBp[normBp<min]=min
normBp[normBp>max]=max
print(normBp)