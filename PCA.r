library(ade4)
library(rgl)

#--------- On importe les données ---------
X <- read.csv2("data_PDE20.csv") # Matrice des individus
X <- X[,-1]    # On supprime les colonnes qui ne nous intéressent pas
X <- X[,-9]
X <- as.matrix(X)

n <- dim(X)[1]
p <- dim(X)[2]

#----------------------------
#--------- Partie I ---------
#----------------------------

#--------- Fonction qui centre notre nuage de points ---------
centrage <- function(X) # ACP centrée
{
  X0 <- matrix(0, dim(X)[1], dim(X)[2])
  for (j in 1:length(X))
  {
    moy <- mean(X[,j])
    for (i in 1:length(X[,1]))
    {
      X0[i,j] <- X[i,j]-moy
    }
  }
  X0
}
#--------- Fonction qui centre et réduit notre nuage de points ---------
centrageRed <- function(X) # ACP centrée et réduite
{
  X0 <- matrix(0, dim(X)[1], dim(X)[2])
  for (j in 1:dim(X)[2])
  {
    moy <- mean(X[,j])
    var <- var(X[,j])
    for (i in 1:dim(X)[1])
    {
      X0[i,j] <- (X[i,j]-moy)/sqrt(var)
    }
  }
  X0
}
#--------- On recherche les hyperplans de projections ---------
XCR <- centrageRed(X)
A <- cov(XCR) # Matrice de variance/covariance
diag <- eigen(A) # On diagonalise la matrice en base orthonormée
values <- diag$values # Valeurs propres (ordonnées par ordre décroissant)
vector <- diag$vectors # Vecteurs propres
#--------- Inertie selon chaque composante ---------
barplot(values,
        main = "Valeurs propres de la matrice de covariance",
        xlab = "Lambda",
        ylab = "Inertie",
        names.arg = c("1", "2", "3", "4", "5", "6", "7", "8"),
        col = "darkred")
#--------- Inertie Cumulée ---------
IC <- matrix(0, nrow = 1, ncol = p)
for ( k in 1:p)
{
  IC[k] <- sum(values[1:k])
}
barplot(IC,
        main = "Inertie cumulée",
        xlab = "Nombre de composantes comptabilisées",
        ylab = "Inertie",
        names.arg = c("1", "2", "3", "4", "5", "6", "7", "8"),
        col = "darkred")
#--------- Taux d'inertie expliquée par composantes ---------

TI <- matrix(values/sum(values), nrow=1, ncol=p)
x <- seq(1,p)

plot(x, TI, type = "b",
     pch = 19, col = 4,
     ylim=c(0,1),
     main = "Taux d'inertie expliquée par chaque composante",
     xlab = "Composante",
     ylab = "Taux d'inertie expliqué")
points(x, TI, pch = 19)

#--------- Taux d'inertie expliquée par k composantes ---------
TIk <- matrix(0, nrow=1, ncol=p)
for (k in 1:p)
{
  TIk[k] <- sum(TI[1:k])
}

plot(x, TIk, type = "b",
     pch = 19, col = 4,
     ylim=c(0,1),
     main = "Taux d'inertie expliquée par k composante",
     xlab = "Nb composantes",
     ylab = "Taux d'inertie expliqué")
points(x, TIk, pch = 19)

#--------- On choisi de garder les 3 premières composantes principales ---------
p2 <- 3
HP <- vector[,1:p2] # Hyperplan retenu

#--------- Fonction qui calcule le produit scalaire de deux vecteurs ---------

scal <- function(a,b)
{
  c <- 0
  for (k in 1:length(a))
  {
    c <- c +a[k]*b[k]
  }
  as.numeric(c)
}

#--------- On calcule les nouvelles coordonnées des individus ---------

Xproj <- XCR%*%vector # On projette les individus dans la bon des vecteurs propres
X2 <- Xproj[,1:p2] # On ne garde que les trois premières composantes principales

#--------- Fonction qui calcul la qualité de la projection ---------

Qk <- function(C,k=3)
{
  sum(C[1:k]**2)/sum(C**2)
}

#--------- On calcul la qualité de la projection, et on la représente ---------

Qtot <- vector(length = n)
for (i in 1:n)
{
  Qtot[i] <- Qk(Xproj[i,])
}

boxplot(Qtot,
        main = "Qualité de la projection",
        ylim=c(0,1), xaxt="n" , col="red" , frame=F)


#--------- On calcul la contribution des individus par rapport aux axes ---------

gamma <- matrix(0, nrow = n, ncol = p2)
for (k in 1:n)
{
  for (i in 1:p2)
  {
    gamma[k,i] <- (X2[k,i]**2)/(n*values[i])
  }
}

barplot(t(gamma),
        main = "Contribution des individus par rapport aux trois axes",
        col=c("blue","red", "green"),
        space=c(0,1),
        xlab="Individus",ylim = c(0,0.2),legend = c("Axe 1", "Axe 2", "Axe 3"), beside=TRUE)

#--------- Comparaison avec ade4 ---------

X3 <- dudi.pca(df = X, scannf = FALSE, nf = 3)

print(X3$c1-HP)


#--------- On représente les individus dans le nouveau sous espace formée par les deux premiers axes retenus ---------
CP1 <- X2[,1]
CP2 <- X2[,2]
CP3 <- X2[,3]
plot(CP1,CP2,
     main = "Individus dans le nouveau sous espace")

plot3d(CP1,CP2,CP3, main="Titre") 

#-----------------------------
#--------- Partie II ---------
#-----------------------------

#--------- Nuage isotrope ---------
#--------- On génére les données ---------

données.gen <- function(n)
{
  V <- matrix(c(rnorm(n),rnorm(n),rnorm(n)), nrow = n, ncol = 3)
  for (i in 1:n)
  {
    V[i,] <- V[i,]/sqrt(scal(V[i,],V[i,]))
  }
  V
}

#--------- Fonction qui réalise l'ACP ---------
ACP <- function(data, plot = TRUE, plotbar = TRUE, plotq = TRUE, plotvp = FALSE)
{
  n <- dim(data)[1]
  p <- dim(data)[2]
  if (plot)
  {
    plot3d(data[,1], data[,2], data[,3],xlab="x", ylab="y", zlab="z", col = rainbow(n))
  }
  DCR <- centrageRed(data)
  G <- cov(DCR) # Matrice de variance/covariance
  G.diag <- eigen(G) # On diagonalise la matrice en base orthonormée
  G.values <- G.diag$values # Valeurs propres (ordonnées par ordre décroissant)
  G.vector <- G.diag$vectors # Vecteurs propres
  if (plotbar)
  {
    barplot(G.values,main = "Valeurs propres de la matrice de covariance",xlab = "Lambda",ylab = "Inertie",names.arg = c("1", "2", "3"),col = "darkred") 
  }
  if (plotvp)
  {
    x <- seq(1,p)
    
    plot(x, G.values, type = "b", col = 4, ylim = c(0, max(G.values)), 
         main = "Valeurs propres de la matice de covariance",)
  }
  Dproj <- DCR%*%G.vector # On projette les individus dans la bon des vecteurs propres
  D2 <- Dproj[,1:2] # On ne garde que les trois premières composantes principales
  Qtot2 <- vector(length = n)
  for (i in 1:n)
  {
    Qtot2[i] <- Qk(Dproj[i,],2)
  }
  if (plotq)
  {
    boxplot(Qtot2, ylim = c(0,1),
            main = "Qualité de la projection", xaxt="n" , col="red" , frame=F)
  }
  plot(D2[,1], D2[,2], col = rainbow(n))
}

#--------- On fait l'ACP ---------
for (n in c(10,100,1000))
{
  ACP(données.gen(n))
}


#--------- Nuage non isotrope ---------
#--------- On génére les nouvelles données ---------
données.gen2 <- function(n, normé =TRUE)
{
  x <- sort(rnorm(n))
  y <- rnorm(n)
  z <- rnorm(n) + atan2(x, y)
  V <- matrix(c(x,y,z), nrow = n, ncol = 3)
  if (normé)
  {
    for (i in 1:n)
    {
      V[i,] <- V[i,]/sqrt(scal(V[i,],V[i,]))
    }
  }
  V
}

#--------- On fait l'ACP ---------
for (i in c(10,100,1000))
{
  ACP(données.gen2(n,FALSE))
}

for (i in c(10,100,1000))
{
  ACP(données.gen2(n))
}

#--------- On ajoute des points extrémaux, en choissisant de normer nos données ---------
n <- 1000
Dex <- données.gen2(n)

extr <- function(data, ratio)
{
  n <- dim(Dex)[1]
  m <- max(abs(data))
  nex <- ratio*n
  ntot <- n + nex
  data2 <- matrix(0, nrow = ntot, ncol = 3)
  for (i in 1:n)
  {
    for (k in 1:3)
    {
      data2[i,k] <- data[i,k]
    }
  }
  
  for (i in (n+1):ntot)
  {
    for (k in 1:3)
    {
      data2[i,k] <- sample(x = c(-1,1), size = 1)*runif(1,min=0, max=2*m)
    }
  }
  data2
}
for (i in c(0.1,0.2,0.5))
{
  ACP(extr(Dex, i))
}

#--------- Etape 2 ---------
#--------- 1er nuage de points ---------
DRp <- données.gen(1000)
DRn <- t(DRp)
ACP(DRp, FALSE, FALSE, FALSE, TRUE)
ACP(DRn, FALSE, FALSE, FALSE, TRUE)

#--------- 2eme nuage de points ---------
DRp2 <- données.gen2(1000)
DRn2 <- t(DRp2)
ACP(DRp2, FALSE, FALSE, FALSE, TRUE)
ACP(DRn2, FALSE, FALSE, FALSE, TRUE)

#--------- Vérification des relations de passage ---------

X <- données.gen2(1000)
lbda1 <- eigen(X%*%t(X))$values
lbda1 <- lbda1[which(lbda1 > 10**(-5))]
lbda2 <- eigen(t(X)%*%X)$values
print(lbda1-lbda2)

u <- eigen(t(X)%*%X)$vector
v <- eigen(X%*%t(X))$vector
v <- v[,which(lbda1 > 10**(-5))]

result <- matrix(0, nrow = dim(v)[1], ncol = dim(v)[2])
for (i in 1:length(lbda1))
{
  result[,i] <- v[,i] - X%*%u[,i]/(lbda1[i])
}
boxplot(result,
        main = "Vérification relations de passages" , col="red" , frame=F)

#------------------------------
#--------- Partie III ---------
#------------------------------
#install.packages("corrplot")
#install.packages("FactoMineR")
library(corrplot)
library(FactoMineR)
#--------- 1) ---------
X <- read.csv2("TP4_covC1234_DS19_20.csv") # Matrice des individus
X <- X[,-1]
X <- X[-141,] 
#--------- Pré-traitement des données ---------
#--------- On supprime chaque ligne qui comporte un 0 ou un ? ---------
ind0 <- c()
for (k in 1:dim(X)[1])
{
  for (i in 1:dim(X)[2])
  {
    if ((X[k,i]==0 | X[k,i]=="?") & !(k %in% ind0))
    {
      ind0 <- c(ind0, k)
    }
  }
}
X <- X[-ind0,]
#--------- Analyse des données ---------
mat.Xquanti <- as.matrix(X[,1:14])

boxplot(mat.Xquanti[,13:14],
        main = "Résumé des données" , col="red" , frame=F)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(cor(mat.Xquanti), method="color", col=col(200), 
         type="upper",
         addCoef.col = "black")

compte <- function(mat.cov)
{
  inf <- 0
  sup <- 0
  for (i in 1:dim(mat.cov)[1])
  {
    for (j in 1:dim(mat.cov)[2])
    {
      if (i!=j & abs(mat.cov[i,j])>0.9)
      {
        sup <- sup+1
      }
      if (i!=j & abs(mat.cov[i,j])<0.1)
      {
        inf <- inf +1
      }
    }
  }
  c(inf/2,sup/2)
}

compte.X <- compte(cor(mat.Xquanti))
compte.BF2 <- compte(cor(X[which((X$Campagne=="BF2")),1:14]))
compte.BF3 <- compte(cor(X[which((X$Campagne=="BF3")),1:14]))
compte.CA1 <- compte(cor(X[which((X$Campagne=="CA1")),1:14]))
compte.CA2 <- compte(cor(X[which((X$Campagne=="CA2")),1:14]))
compte.CA3 <- compte(cor(X[which((X$Campagne=="CA3")),1:14]))
compte.CA4 <- compte(cor(X[which((X$Campagne=="CA4")),1:14]))
compte.ete <- compte(cor(X[which((X$SAISON=="ete")),1:14]))
compte.hiver <- compte(cor(X[which((X$SAISON=="hiver")),1:14]))
compte.BF <- compte(cor(rbind(mat.BF2, mat.BF3)))
compte.CA <- compte(cor(rbind(mat.CA1, mat.CA2, mat.CA3, mat.CA4)))

#--------- ACP ---------

mat.covdiag <- eigen(cov(centrageRed(Xquanti)))
values <- mat.covdiag$values
p <- dim(Xquanti)[2]
TI <- matrix(values/sum(values), nrow=1, ncol=p)
x <- seq(1,p)

plot(x, TI, type = "b",
     pch = 19, col = 4,
     ylim=c(0,1),
     main = "Taux d'inertie expliquée par chaque composante",
     xlab = "Composante",
     ylab = "Taux d'inertie expliqué")
points(x, TI, pch = 19)

resultat<-PCA(Xquanti)

pie(resultat$var$contrib[,1], main = "Composante 1") 
pie(resultat$var$contrib[,2], main = "Composante 2") 
