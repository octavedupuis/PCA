#----------------------------------------------------------------------------------
#----------------- TP 2 : Analyse factorielle des correspondances -----------------
#----------------------------------------------------------------------------------
library(FactoMineR)


X <- read.csv2("TP_AFC_data.csv")
X <- X[,-1]
X <- X[-251,]
X <- X[-250,]
X <- X[-249,]


#----------------- TP 2 : Analyse factorielle des correspondances -----------------
tabl_corres <- function(X) #Prend en argument un tableau de données
{
  ndata <- dim(X)[1]
  lign <- unique(X[,1])
  col <- unique(X[,2])
  n <- length(lign)
  p <- length(col)
  T <- matrix(0, nrow = n, ncol = p)
  rownames(T) <- lign
  colnames(T) <- col
  for (i in 1:ndata)
  {
    T[as.character(X[i,1]),as.character(X[i,2])] <- T[as.character(X[i,1]),as.character(X[i,2])]+1
  }
  T
}

freq_relat <- function(T) # Prend en argument un tableau de correspondance
{
  F <- T/sum(T)
  F
}

profil_ligne <- function(F) # Prend en arguement un tableau de fréquences
{
  n <- dim(F)[1]
  sum_i <- vector(length=n)
  for (i in 1:n)
  {
    sum_i[i] <- sum(F[i,])
  }
  PL <- F
  for (i in 1:n)
  {
    PL[i,] <- PL[i,]/sum_i[i]
  }
  PL
}

profil_colonne <- function(F) # Prend en arguement un tableau de fréquences
{
  p <- dim(F)[2]
  sum_j <- vector(length=p)
  for (j in 1:p)
  {
    sum_j[j] <- sum(F[,j])
  }
  PC <- F
  for (j in 1:p)
  {
    PC[,j] <- PC[,j]/sum_j[j]
  }
  PC
}

AFC <- function(F, plot_hist=FALSE) # Prend en argument un tableau de fréquences
{
  n <- dim(F)[1]
  p <- dim(F)[2]
  Dn <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n)
  {
    Dn[i,i] <- sum(F[i,])
  }
  Dp <- matrix(0, nrow = p, ncol = p)
  for (j in 1:p)
  {
    Dp[j,j] <- sum(F[,j])
  }
  A <- sqrt(solve(Dp))%*%t(F)%*%solve(Dn)%*%F%*%sqrt(solve(Dp))
  diagA <- eigen(A)
  
  coord_fac <- solve(Dn)%*%F%*%solve(Dp)%*%sqrt(Dp)%*%diagA$vectors
  
  diagA$values <- diagA$values[-1]
  diagA$vectors <- diagA$vectors[, -1]
  coord_fac <- coord_fac[, -1]
  
  if (plot_hist)
  {
    v_non_nulles <- diagA$values[which(diagA$values>10**(-10))]
    
    TI <- matrix(v_non_nulles/sum(v_non_nulles), nrow=1, ncol=length(v_non_nulles))
    
    barplot(TI,
            main = "Taux d'inertie expliquée par chaque composante",
            xlab = "Composante principale",
            ylab = "Taux d'inertie expliquée",
            col = "darkred")
  }
  coord_fac                
}

F <- freq_relat(tabl_corres(X))
AFC1 <- AFC(F, TRUE)
AFC2 <- AFC(t(F))

#--------- On projette sur deux axes ---------

nproj <- 2
C1 <- AFC1[,1:nproj]
C2 <- AFC2[,1:nproj]
mat=rbind(C1,C2)
plot(c(),
     xlim=range(mat[,1]),
     ylim=range(mat[,2]),
     main="AFC",
     bty="l", tcl= -.25,
     xlab="Dim 1", ylab="Dim 2")

grid()
points(C1, col="blue", pch=23, bg="red")
points(C2, col="red", pch=25, bg="yellow") 

AFC_facto <- CA(F)

#--------- On arrange les signes de nos coordonnées ---------

C1[,1] <- -C1[,1]
C1[,2] <- -C1[,2]
C2[,2] <- -C2[,2]
mat=rbind(C1,C2)
plot(c(),
     xlim=range(mat[,1]),
     ylim=range(mat[,2]),
     main="AFC",
     bty="l", tcl= -.25,
     xlab="Dim 1", ylab="Dim 2"
)
grid()
points(C1, col="blue", pch=23, bg="red")
points(C2, col="red", pch=25, bg="yellow") 


#--------- Fonction qui calcul la qualité de la projection ---------

Qk <- function(C,k)
{
  sum(C[1:k]**2)/sum(C**2)
}

#--------- On calcul la qualité de la projection, et on la représente ---------

Qualité <- function(coord, nproj = 2)
{
  n <- dim(coord)[1]
  Qtot <- vector(length = n)
  for (i in 1:n)
  {
    Qtot[i] <- Qk(coord[i,], nproj)
  }

  boxplot(Qtot,
          main = "Qualité de la projection",
          ylim=c(0,1), xaxt="n" , col="red" )
}

