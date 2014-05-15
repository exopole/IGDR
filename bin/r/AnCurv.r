

#########################################################################################
#
#
#
#	Calcule des coordonnees du maximum de la courbe accuracy



# Calcul de la longueur des listes
len = lapply(slot(perf,"x.values"),length)

# Ressort la taille de la liste la plus courte
minlen = min(unlist(len))

# nombre de list
nblist2 = length(slot(perf,"x.values"))

# Initialisation d'une matrice de 10 lignes et de 1622 colonnes.
mat <- matrix(data = NA, nrow = nblist, ncol = minlen2 , byrow = FALSE,dimnames = NULL)

# Assignation des valeurs aux differentes lignes
#mat[1,] <- slot(perf, "y.values")[[1]][1:minlen2]
#mat[2,] <- slot(perf, "y.values")[[2]][1:minlen2]
#mat[3,] <- slot(perf, "y.values")[[3]][1:minlen2]
#mat[4,] <- slot(perf, "y.values")[[4]][1:minlen2]
#mat[5,] <- slot(perf, "y.values")[[5]][1:minlen2]
#mat[6,] <- slot(perf, "y.values")[[6]][1:minlen2]
#mat[7,] <- slot(perf, "y.values")[[7]][1:minlen2]
#mat[8,] <- slot(perf, "y.values")[[8]][1:minlen2]
#mat[9,] <- slot(perf, "y.values")[[9]][1:minlen2]
#mat[10,] <- slot(perf, "y.values")[[10]][1:minlen2]

# utilisation d'une boucle pour l'incorporation des données
for(elt in 1:nblist){mat[elt,] <- slot(perf, "y.values")[[elt]][1:minlen2]}

# Calcul de la moyenne par colonne
moyC = apply(mat, 2,mean)

# Maximun de moyC
Maxi = max(moyC)

# Trouver les indexs correspondant au valeurs indiquer
coord = seq_along(moyC)[sapply(moyC,FUN=function(X) Maxi %in% X)]

# longueur des liste de x.values
len2 = lapply(slot(perf,"x.values"),length)

# Ressort la taille de la liste la plus courte
minlen2 = min(unlist(len2))

# nombre de list
nblist2 = length(slot(perf,"x.values"))

# Initialisation d'une matrice de 10 lignes et de 1622 colonnes.
mat2 <- matrix(data = NA, nrow = nblist2, ncol = minlen2, byrow = FALSE,dimnames = NULL)

# Assignation des valeurs aux differentes lignes
#mat2[1,] <- slot(perf, "x.values")[[1]][1:minlen2]
#mat2[2,] <- slot(perf, "x.values")[[2]][1:minlen2]
#mat2[3,] <- slot(perf, "x.values")[[3]][1:minlen2]
#mat2[4,] <- slot(perf, "x.values")[[4]][1:minlen2]
#mat2[5,] <- slot(perf, "x.values")[[5]][1:minlen2]
#mat2[6,] <- slot(perf, "x.values")[[6]][1:minlen2]
#mat2[7,] <- slot(perf, "x.values")[[7]][1:minlen2]
#mat2[8,] <- slot(perf, "x.values")[[8]][1:minlen2]
#mat2[9,] <- slot(perf, "x.values")[[9]][1:minlen2]
#mat2[10,] <- slot(perf, "x.values")[[10]][1:minlen2]

# loop for
for(elt in 1:nblist){mat2[elt,] <- slot(perf, "x.values")[[elt]][1:minlen2]}

# Calcul de la moyenne par colonne
moyC2 = apply(mat2, 2,mean)

# extractions des valeurs dans moyC2
xval = moyC2[coord]

# moyennes des coordonnees trouver
xval_moy = mean(xval)



##########################################################################################
#
#
#
#
#		Calcule du point d'intersection entre la courbe de specificite et de sensibilite


#### Sensibilite #####

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

# Calcul de la moyenne par colonne
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
moyDif2 = moyDif %*% moyDif

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


#########################################################################################
#
#
#
#
#		Afficher les resultats



######## Accurancy  curve #######
# Affiche la meilleur accurancy (axe des Y)
Maxi 

# Affiche le cuttof du potentiel codant associe (axe des X)
xval_moy 

######## Performance ########
# Affiche la valeur de cutoff correspondant à la distance la plus courte entre P.y et S.y 
coordPSx

