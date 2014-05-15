#!/usr/bin/Rscript



A <- read.table("size_ORF.1", header=T,sep="\t")
Aprime <- read.table("size_ORF.2", header=T, sep="\t")
B <- read.table("size_NONCODE.csv", header=T, sep="\t")

# A <- read.table("size_NONCODE.csv", header=T,sep="\t")
# Aprime <- read.table("size_NONCODE.csv", header=T, sep="\t")
# B <- read.table("size_ORF.csv", header=T, sep="\t")

png("size.png")
par(mfrow=c(3,4),mar=c(5,4,2,2),cex.axis=1.2, cex.lab=1.2)

plot(density(A[,2]),ylim=c(0,0.002),xlim=c(0,2000))
lines(density(B[,2]))
lines(density(Aprime[,2]),col=2)


### ks test est super sensible aux grandes tailles d'echantillon
ks.test(A[,2],Aprime[,2])

### Vision graphique de la normalite: Si t'as une droite c'est que ca suit à peu près une loi normale...

logA <- log(A[,2])
logAp <- log(Aprime[,2])

qqnorm((logA-mean(logA))/sd(logA))
abline(c(0,1))
## Plutot bon pour log A

qqnorm((logAp-mean(logAp))/sd(logAp))
abline(c(0,1))
## Plutot bon pour log Aprime


logB <- log(B[,2])
qqnorm((logB-mean(logB))/sd(logB))
abline(c(0,1))
## Plutot moyen pour log B
 
## Pour B prime voila ce que je ferai

mLogBp <- mean(logB)*mean(logAp)/mean(logA)
# vLogBp <- var(logB)*var(logAp)/var(logA)


########## moi ##################
sdLogBp <- sd(logB)*sd(logAp)/sd(logA)

mBp = exp(mLogBp)
sdLogBp <- sd(logB)*sd(logAp)/sd(logA)
sdBp <- exp(sdLogBp)
sdBp2 <- sd(exp(logB))*sd(exp(logAp))/sd(exp(logA))


sdA =exp(sd(logA))
moyA = mean(A[,2])

maxA <- max(A[,2])
denA <- density(B[,2])

toto[toto<50]=mBp

# Je te joins une commande R qui permet d'extraire des valeurs aléatoires d'une distribution étant donné sa moyenne (moyB) et son écart-type (sdB).
# Ce sera utile pour extraire des tailles de lncRNAs (Bprime).
# Tu peux utiliser :
set.seed(1) # pour intiliaser la graine du random et reproduire la même valeur
normBp <-rnorm(100, mean=mBp, sd=sdBp) # avec nbvaleur ==1 si tu veux seulement 1 valeur
plot(density(normBp))

normBp2 <-rnorm(1000, mean=mBp, sd=sdBp2)
plot(density(normBp2),xlim=c(0,4000))


toto <- normBp2
# toto[toto<50]=mBp
qqnorm((log(toto)-mean(log(toto)))/sd(log(toto)))
abline(c(0,1))

plot(density(toto))

dev.off()


