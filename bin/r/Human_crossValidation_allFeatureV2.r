#!/local/R/3.0.1/bin/Rscript

library(ROCR)
args <- commandArgs(trailingOnly = TRUE)
monfichier <-args[1]
name_png <-args[2]
data=read.table(file=monfichier,header=T,sep="\t")
#data=read.table(file= "cpat2.out",header=T,sep="\t")
attach(data)
names(data)



number_row = nrow(data)
d1 =seq(1,as.integer(number_row/10))
d2 =seq(as.integer(number_row/10)+1,as.integer(2*number_row/10))
d3 = seq(as.integer(2*number_row/10)+1,as.integer(3*number_row/10))
d4 = seq(as.integer(3*number_row/10)+1,as.integer(4*number_row/10))
d5 = seq(as.integer(4*number_row/10)+1,as.integer(5*number_row/10))
d6 = seq(as.integer(5*number_row/10)+1,as.integer(6*number_row/10))
d7 = seq(as.integer(6*number_row/10)+1,as.integer(7*number_row/10))
d8 = seq(as.integer(7*number_row/10)+1,as.integer(8*number_row/10))
d9 = seq(as.integer(8*number_row/10)+1,as.integer(9*number_row/10))
d10 = seq(as.integer(9*number_row/10)+1,number_row)
 

#d1
vlabel = Label[-d1]
vmrna = mRNA[-d1]
vorf = ORF[-d1]
vfickett = Fickett[-d1]
vhexamer = Hexamer[-d1]


mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d1], vorf = ORF[d1],vfickett = Fickett[d1], vhexamer = Hexamer[d1], vlabel=Label[d1])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind(	"mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file="test1.xls",quote=F,sep="\t",row.names=ID[d1])


#d2
vlabel = Label[-d2]
vmrna = mRNA[-d2]
vorf = ORF[-d2]
vfickett = Fickett[-d2]
vhexamer = Hexamer[-d2]

mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d2], vorf = ORF[d2],vfickett = Fickett[d2], vhexamer = Hexamer[d2], vlabel=Label[d2])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file="test2.xls",quote=F,sep="\t",row.names=ID[d2])



#d3
vlabel = Label[-d3]
vmrna = mRNA[-d3]
vorf = ORF[-d3]
vfickett = Fickett[-d3]
vhexamer = Hexamer[-d3]

mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d3], vorf = ORF[d3],vfickett = Fickett[d3], vhexamer = Hexamer[d3], vlabel=Label[d3])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file="test3.xls",quote=F,sep="\t",row.names=ID[d3])



#d4
vlabel = Label[-d4]
vmrna = mRNA[-d4]
vorf = ORF[-d4]
vfickett = Fickett[-d4]
vhexamer = Hexamer[-d4]

mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d4], vorf = ORF[d4],vfickett = Fickett[d4], vhexamer = Hexamer[d4], vlabel=Label[d4])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file="test4.xls",quote=F,sep="\t",row.names=ID[d4])



#d5
vlabel = Label[-d5]
vmrna = mRNA[-d5]
vorf = ORF[-d5]
vfickett = Fickett[-d5]
vhexamer = Hexamer[-d5]

mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d5], vorf = ORF[d5],vfickett = Fickett[d5], vhexamer = Hexamer[d5], vlabel=Label[d5])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file="test5.xls",quote=F,sep="\t",row.names=ID[d5])



#d6
vlabel = Label[-d6]
vmrna = mRNA[-d6]
vorf = ORF[-d6]
vfickett = Fickett[-d6]
vhexamer = Hexamer[-d6]

mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d6], vorf = ORF[d6],vfickett = Fickett[d6], vhexamer = Hexamer[d6], vlabel=Label[d6])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file="test6.xls",quote=F,sep="\t",row.names=ID[d6])



#d7
vlabel = Label[-d7]
vmrna = mRNA[-d7]
vorf = ORF[-d7]
vfickett = Fickett[-d7]
vhexamer = Hexamer[-d7]

mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d7], vorf = ORF[d7],vfickett = Fickett[d7], vhexamer = Hexamer[d7], vlabel=Label[d7])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file="test7.xls",quote=F,sep="\t",row.names=ID[d7])



#d8
vlabel = Label[-d8]
vmrna = mRNA[-d8]
vorf = ORF[-d8]
vfickett = Fickett[-d8]
vhexamer = Hexamer[-d8]

mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d8], vorf = ORF[d8],vfickett = Fickett[d8], vhexamer = Hexamer[d8], vlabel=Label[d8])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file="test8.xls",quote=F,sep="\t",row.names=ID[d8])



#d9
vlabel = Label[-d9]
vmrna = mRNA[-d9]
vorf = ORF[-d9]
vfickett = Fickett[-d9]
vhexamer = Hexamer[-d9]

mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d9], vorf = ORF[d9],vfickett = Fickett[d9], vhexamer = Hexamer[d9], vlabel=Label[d9])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file="test9.xls",quote=F,sep="\t",row.names=ID[d9])



#d10
vlabel = Label[-d10]
vmrna = mRNA[-d10]
vorf = ORF[-d10]
vfickett = Fickett[-d10]
vhexamer = Hexamer[-d10]

mylogit <- glm(vlabel ~ vmrna + vorf + vfickett + vhexamer, family=binomial(link="logit"), na.action=na.pass)
test <- data.frame(vmrna = mRNA[d10], vorf = ORF[d10],vfickett = Fickett[d10], vhexamer = Hexamer[d10], vlabel=Label[d10])
test$prob <- predict(mylogit,newdata=test,type="response")
output = cbind("mRNA"=test$vmrna, "ORF"=test$vorf, "Fickett"=test$vfickett, "Hexamer" = test$vhexamer, "Label"=test$vlabel,"Prob"=test$prob)
write.table(output,file="test10.xls",quote=F,sep="\t",row.names=ID[d10])



#ROC
test1=read.table(file="test1.xls",header=T,sep="\t")
test2=read.table(file="test2.xls",header=T,sep="\t")
test3=read.table(file="test3.xls",header=T,sep="\t")
test4=read.table(file="test4.xls",header=T,sep="\t")
test5=read.table(file="test5.xls",header=T,sep="\t")
test6=read.table(file="test6.xls",header=T,sep="\t")
test7=read.table(file="test7.xls",header=T,sep="\t")
test8=read.table(file="test8.xls",header=T,sep="\t")
test9=read.table(file="test9.xls",header=T,sep="\t")
test10=read.table(file="test10.xls",header=T,sep="\t")


Response = list(test1$Prob,test2$Prob,test3$Prob,test4$Prob,test5$Prob,test6$Prob,test7$Prob,test8$Prob,test9$Prob,test10$Prob)
Labls = list(test1$Label,test2$Label,test3$Label,test4$Label,test5$Label,test6$Label,test7$Label,test8$Label,test9$Label,test10$Label)
ROCR_data = list(predictions=Response,Labels=Labls)
pred <- prediction(ROCR_data$predictions, ROCR_data$Labels)
#perf <- performance(pred,"auc")
#avergae AUC = 0.9927


png(name_png)
par(mfrow=c(2,3),mar=c(5,4,2,2),cex.axis=1.2, cex.lab=1.2)
#ROC curve
#pdf("Human_10fold.ROC.pdf")
perf <- performance(pred,"tpr","fpr")
plot(perf,col="blue",lty=3,xlab="1-Specificity",ylab="Sensitivity",ylim=c(0.7,1),xlim=c(0,0.3),main="",cex.axis=1.5,cex.label=1.5)	#AUC = 0.9927 
plot(perf,lwd=2,avg="vertical",add=TRUE,col="red",xlab="1-specificity",ylab="sensitivity",main="",cex.axis=1.2,cex.label=1.2) 
abline(v=0,lty="dashed",lwd=0.5)
abline(h=1.0,lty="dashed",lwd=0.5)
abline(v=0.05,lty="dashed",lwd=0.5)
abline(h=0.95,lty="dashed",lwd=0.5)
#dev.off()

#precision
#pdf("Human_10fold.precision_vs_recall.pdf")
d=performance(pred,measure="prec", x.measure="rec")
plot(d,col="blue",lty=3,xlab="Recall (TPR)",ylab="Precision (PPV)",xlim=c(0.7,1),ylim=c(0.7,1),cex.axis=1.2,cex.label=1.2)
plot(d,lwd=2,avg="vertical",col="red",xlab="Recall (TPR)",ylab="Precision (PPV)",add=T,cex.axis=1.2,cex.label=1.2)
abline(v=1.0,lty="dashed",lwd=0.5)
abline(h=1.0,lty="dashed",lwd=0.5)
abline(v=0.95,lty="dashed",lwd=0.5)
abline(h=0.95,lty="dashed",lwd=0.5)
#dev.off()


#Accuracy
#pdf("Human_10fold.Accuracy.pdf")
perf <- performance(pred,"acc")
plot(perf,col="blue",lty=3,xlab="Coding probability cutoff",ylab="Accuracy",ylim=c(0.7,1),cex.axis=1.2,cex.label=1.2) 
plot(perf,lwd=2,avg="vertical",add=TRUE,col="red",cex.axis=1.2,cex.label=1.2)

############ ajout de AnCurv.r ############

# compute the list size 
len = lapply(slot(perf,"x.values"),length)
# print(len)

# take the size of the shortest list
minlen = min(unlist(len))
# print(minlen)

# number of list
nblist = length(slot(perf,"x.values"))
# print(nblist)
# attributes(perf)

# Initialization of a matrice of 10 row and 1622 column
mat <- matrix(data = NA, nrow = nblist, ncol = minlen , byrow = FALSE,dimnames = NULL)

# Assigning values to different lines (using of a for loop) 
for(elt in 1:nblist){mat[elt,] <- slot(perf, "y.values")[[elt]][1:minlen]}

# Compute the mean by column
moyC = apply(mat, 2,mean)

# maximum of moy
Maxi = max(moyC)

# find the matching indexs vallues indicate
coord = seq_along(moyC)[sapply(moyC,FUN=function(X) Maxi %in% X)]

# list length of x.values 
len2 = lapply(slot(perf,"x.values"),length)

# Spring size of the shortest list
minlen2 = min(unlist(len2))

# number of list
nblist2 = length(slot(perf,"x.values"))

# Initialization of a matrice of 10 row and 1622 column
mat2 <- matrix(data = NA, nrow = nblist2, ncol = minlen2, byrow = FALSE,dimnames = NULL)

# Assigning values to different lines (using of a for loop) 
for(elt in 1:nblist){mat2[elt,] <- slot(perf, "x.values")[[elt]][1:minlen2]}

# Compute the mean by column
moyC2 = apply(mat2, 2,mean)

# extract value of moyc2
xval = moyC2[coord]

# mean of coordonne find
xval_moy = mean(xval)

print (xval_moy)
############################################

abline(h=1,lty="dashed",lwd=0.5)
abline(h=0.95,lty="dashed",lwd=0.5)
#dev.off()


#sensitivity vs specificity
pred <- prediction(ROCR_data$predictions, ROCR_data$Labels)
S <- performance(pred,measure="sens")
P <- performance(pred,measure="spec")


#pdf("Human_10fold_sens_vs_spec.pdf")
plot(S,col="blue",lty=3,ylab="Zoom Performance",xlab="Coding Probability Cutoff",ylim=c(0.95,1),cex.axis=1.2,cex.label=1.2, xlim = c(0.3,0.8)) 
plot(S,lwd=2,avg="vertical",add=TRUE,col="blue") 
plot(P,col="red",lty=3, add=TRUE,) 
plot(P,lwd=2,avg="vertical",add=TRUE,col="red")

############ ajout de Ancurv.r ##############
### sensibility

# Calcul de la longueur des listes
lenS = lapply(slot(S,"y.values"),length)
# Ressort la taille de la liste la plus courte
minlenS = min(unlist(lenS))
# nombre de list
nblistS = length(slot(S,"y.values"))
# Initialisation d'une matrice de 10 lignes et de 1622 colonnes.
matS <- matrix(data = NA, nrow = nblist, ncol = minlenS , byrow = FALSE,dimnames = NULL)
# utilisation d'une boucle pour l'incorporation des données
for(elt in 1:nblistS){matS[elt,] <- slot(S, "y.values")[[elt]][1:minlenS]}
moyS = apply(matS, 2,mean)
#### Specificite #####

# Calcul de la longueur des listes
lenP = lapply(slot(P,"y.values"),length)
# Ressort la taille de la liste la plus courte
minlenP = min(unlist(lenP))
# nombre de list
nblistP = length(slot(P,"y.values"))
# Initialisation d'une matrice de nblist lignes et de minlenP colonnes.
matP <- matrix(data = NA, nrow = nblistP, ncol = minlenP , byrow = FALSE,dimnames = NULL)
# utilisation d'une boucle pour l'incorporation des données
for(elt in 1:nblistP){matP[elt,] <- slot(P, "y.values")[[elt]][1:minlenP]}
# Calcul de la moyenne par colonne
moyP = apply(matP, 2,mean)


### Calcule du point le plus proche ####

# calcule de la difference entre la matrice de moyenne de sensibilite et celle de specificite
moyDif = moyP-moyS
# calcul de la valeur absolu
moyDif2 = moyDif * moyDif
# prend la valeur la plus petite
minDiff = min(moyDif2)
# coordonne de minDiff
coordDiff = match(min(moyDif2),moyDif2)
# matrice des x.values de S ou P (ce sont les memes)
matPx <- matrix(data = NA, nrow = nblistP, ncol = minlenP , byrow = FALSE,dimnames = NULL)
for(elt in 1:nblistP){matPx[elt,] <- slot(P, "x.values")[[elt]][1:minlenS]}
# Moyenne par colonne
moyPSx = apply(matPx, 2, mean)
# Coordonne de x de l'intersection
coordPSx = moyPSx[coordDiff]
print (coordPSx)

#############################################
abline(h=moyP[coordDiff],lty="dashed",lwd=0.5)
abline(v=coordPSx,lty="dashed",lwd=0.5)
legend(0.4,0.85,col=c("blue","red"),lwd=2,legend=c("Sensitivity","Specificity"))

## without zoom
plot(S,col="blue",lty=3,ylab="Performance",xlab="Coding Probability Cutoff",ylim=c(0.8,1),cex.axis=1.2,cex.label=1.2) 
plot(S,lwd=2,avg="vertical",add=TRUE,col="blue") 
plot(P,col="red",lty=3, add=TRUE,) 
plot(P,lwd=2,avg="vertical",add=TRUE,col="red")
abline(h=moyP[coordDiff],lty="dashed",lwd=0.5)
abline(v=coordPSx,lty="dashed",lwd=0.5)
legend(0.4,0.85,col=c("blue","red"),lwd=2,legend=c("Sensitivity","Specificity"))

dev.off()
