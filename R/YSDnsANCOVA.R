#' Neutrosophic Analysis of Covariance for Youden Square Design
#'
#' Performs Neutrosophic Analysis of Covariance (NANCOVA) for Youden
#' square designs using interval-valued observations and covariates.
#' The function computes neutrosophic sums of squares, mean squares,
#' interval-valued F-statistics, significance tests, and Least
#' Significant Difference (LSD)-based treatment comparisons.
#' For crisp data, enter identical lower and upper values to obtain
#' the corresponding classical ANCOVA results.
#'
#' @usage
#' YSDnsANCOVA(Lower_y, Upper_y, Lower_z, Upper_z, design, alpha = 0.05, verbose = FALSE)
#'
#' @param Lower_y Matrix containing lower bounds of response observations.
#' @param Upper_y Matrix containing upper bounds of response observations.
#' @param Lower_z Matrix containing lower bounds of covariate observations.
#' @param Upper_z Matrix containing upper bounds of covariate observations.
#' @param design Matrix representing Youden square treatment allocation.
#' @param alpha Significance level for the F-test and LSD test. Default is 0.05.
#' @param verbose Logical. If TRUE, prints the ANCOVA table,
#' LSD interval, treatment comparisons and significance codes.
#' Default is FALSE.
#'
#' @details
#' Input matrix structure:
#' \itemize{
#'   \item Rows represent blocks of the design.
#'   \item Columns represent treatment positions within each row.
#'   \item `Lower_y` and `Upper_y` must have the same dimensions as
#'   the design matrix.
#'   \item `Lower_z` and `Upper_z` must have the same dimensions as
#'   the design matrix.
#' }
#'
#' @return
#' A list containing:
#' \itemize{
#' \item \code{nancova_table}: Neutrosophic ANCOVA table.
#' \item \code{comparison}: LSD treatment comparisons, if performed.
#' \item \code{LSD}: Lower and upper limits of the LSD interval.
#' }
#'
#' @examples
#' Lower_y <- matrix(c(
#' 2.06,2.34,0.13,15.63,15.72,8.18,292.48,
#' 15.39,218.67,10.02,13.67,8.90,25.23,24.86,
#' 35.78,24.90,308.36,19.96,22.01,20.17,28.99,
#' 45.05,42.44,33.19,280.58,31.40,31.82,28.41
#' ), nrow = 4, byrow = TRUE)
#'
#' Upper_y <- matrix(c(
#' 5.94,8.26,2.33,18.17,18.08,12.42,295.52,
#' 19.61,221.33,14.38,17.33,13.10,27.77,29.54,
#' 38.22,27.10,311.64,25.44,26.39,22.63,33.60,
#' 48.55,45.96,35.41,283.42,36.00,35.58,32.59
#' ), nrow = 4, byrow = TRUE)
#'
#' Lower_z <- matrix(c(
#' 10.57,269.93,224.08,260.05,257.81,257.58,13.76,
#' 246.72,7.42,216.00,257.15,232.38,237.91,215.32,
#' 246.50,215.17,2.63,257.26,267.90,212.10,254.15,
#' 262.67,250.51,250.45,6.09,245.01,228.31,258.90
#' ), nrow = 4, byrow = TRUE)
#'
#' Upper_z <- matrix(c(
#' 13.43,272.07,227.92,265.95,260.19,260.42,18.24,
#' 249.28,12.58,220.00,262.85,235.62,240.09,218.68,
#' 251.50,220.83,7.37,260.74,272.10,215.90,257.85,
#' 265.33,253.49,253.55,11.91,250.99,231.69,263.10
#' ), nrow = 4, byrow = TRUE)
#'
#' design <- matrix(c(
#' 2,3,4,5,6,7,1,
#' 7,1,2,3,4,5,6,
#' 6,7,1,2,3,4,5,
#' 5,6,7,1,2,3,4
#' ), nrow = 4, byrow = TRUE)
#'
#' YSDnsANCOVA(Lower_y, Upper_y, Lower_z, Upper_z, design, alpha = 0.05,  verbose = TRUE)
#'
#' @importFrom MASS ginv
#' @importFrom stats qf qt
#' @export

YSDnsANCOVA <- function(Lower_y, Upper_y, Lower_z, Upper_z, design, alpha = 0.05,  verbose = FALSE){

  if(!is.matrix(Lower_y) || !is.matrix(Upper_y) || !is.matrix(Lower_z) || !is.matrix(Upper_z) || !is.matrix(design)){
    stop("Inputs must be matrices.")
  }

  if(any(dim(Lower_y) != dim(Upper_y))){
    stop("Lower_y and Upper_y must have same dimensions.")
  }

  if(any(dim(Lower_z) != dim(Upper_z))){
    stop("Lower_z and Upper_z must have same dimensions.")
  }

  if(any(dim(Lower_y) != dim(design))){
    stop("Design and response matrices must have same dimensions.")
  }

  if(any(dim(Lower_z) != dim(design))){
    stop("Design and covariate matrices must have same dimensions.")
  }

  v <- max(design)
  b <- nrow(design)
  k <- ncol(design)

  rep_counts <- table(design)

  r <- as.numeric(unique(rep_counts))

  interval <- function(L,U){
    list(
      d = (L+U)/2,
      i = abs(U-L)/2
    )
  }

  n <- b*k

  #### X matrix

trt_vec <- as.vector(t(design))

Tmat <- matrix(0,n,v)

for(i in 1:n){
  Tmat[i,trt_vec[i]] <- 1
}

mu <- matrix(1,n,1)

blk <- rep(1:k,each=b)
Bmat <- diag(v)[rep(1:k, times = b), ]
rep_vec <- rep(1:r,each=b/r*k)

Rmat <- matrix(0,n,r)

for(i in 1:n){
  Rmat[i,rep_vec[i]] <- 1
}

X <- cbind(Tmat,mu,Bmat,Rmat)
d_y<-(Lower_y+Upper_y)/2
i_y<-abs(Upper_y-Lower_y)/2
d_z<-(Lower_z+Upper_z)/2
i_z<-abs(Upper_z-Lower_z)/2
ybar_l<-mean(Lower_y)
ybar_u<-mean(Upper_y)
ybar_d<-mean(d_y)
ybar_i<-abs(ybar_u-ybar_l)/2
zbar_l<-mean(Lower_z)
zbar_u<-mean(Upper_z)
zbar_d<-(zbar_l+zbar_u)/2
zbar_i<-abs(zbar_u-zbar_l)/2
yminusbar_l<-d_y-ybar_d
yminusbar_u<-yminusbar_l+i_y-ybar_i
yminusbar_d<-(yminusbar_l+yminusbar_u)/2
yminusbar_i<-abs(yminusbar_u-yminusbar_l)/2
zminusbar_l<-d_z-zbar_d
zminusbar_u<-zminusbar_l+i_z-zbar_i
zminusbar_d<-(zminusbar_l+zminusbar_u)/2
zminusbar_i<-abs(zminusbar_u-zminusbar_l)/2
dy<- as.vector(t(yminusbar_d))
iy<-as.vector(t(yminusbar_i))
dz<- as.vector(t(zminusbar_d))
iz<- as.vector(t(zminusbar_i))
zpz_l<-sum(dz*dz)
zpz_u<-zpz_l+t(dz)%*%iz+t(iz)%*%dz+t(iz)%*%iz
zpz_d<-(zpz_l+zpz_u)/2
zpz_i<-abs(zpz_u-zpz_l)/2
zpy_l<-t(dz)%*%dy
zpy_u<-zpy_l+t(dz)%*%iy+t(iz)%*%dy+t(iz)%*%iy
zpy_d<-(zpy_l+zpy_u)/2
zpy_i<-abs(zpy_u-zpy_l)/2
ypy_l<-t(dy)%*%dy
ypy_u<-ypy_l+t(dy)%*%iy+t(iy)%*%dy+t(iy)%*%iy
ypy_d<-(ypy_l+ypy_u)/2
ypy_i<-abs(ypy_u-ypy_l)/2

# Response#treatment
# Response#treatment
##Trt unadj
Ct_unadj <- t(Tmat) %*% Tmat
bcdy_t_un_d<-t(cbind(Tmat,Bmat,Rmat))%*%mu%*%ginv(t(mu)
                                                  %*%mu)%*%t(mu)%*%dy
bcdy_t_un_i<-t(cbind(Tmat,Bmat,Rmat))%*%mu%*%ginv(t(mu)
                                                  %*%mu)%*%t(mu)%*%iy
Qy_un_l<- (t(cbind(Tmat,Bmat,Rmat))%*%dy-bcdy_t_un_d)
Qy_un_u<- Qy_un_l+t(cbind(Tmat,Bmat,Rmat))%*%iy-bcdy_t_un_i
Qy_t_un_d<-(Qy_un_l[1:v,]+Qy_un_u[1:v,])/2
Qy_t_un_i<-abs(Qy_un_l[1:v,]-Qy_un_u[1:v,])/2
##trt adj
Z <- cbind(mu, Bmat,Rmat)
M <- diag(nrow(X)) - Z %*% ginv(t(Z) %*% Z) %*% t(Z)
Ct_adj <- t(Tmat) %*% M %*% Tmat
bcdy_t_ad_d<-t(Tmat)%*%cbind(mu,Bmat,Rmat)%*%ginv(t(cbind(mu,Bmat,Rmat))
                                                  %*%cbind(mu,Bmat,Rmat))%*%t(cbind(mu,Bmat,Rmat))%*%dy
bcdy_t_ad_i<-t(Tmat)%*%cbind(mu,Bmat,Rmat)%*%ginv(t(cbind(mu,Bmat,Rmat))
                                                  %*%cbind(mu,Bmat,Rmat))%*%t(cbind(mu,Bmat,Rmat))%*%iy
Qy_ad_l<- t(Tmat)%*%dy-bcdy_t_ad_d
Qy_ad_u<- Qy_ad_l+t(Tmat)%*%iy-bcdy_t_ad_i
Qy_t_ad_d<-(Qy_ad_l+Qy_ad_u)/2
Qy_t_ad_i<-abs(Qy_ad_l-Qy_ad_u)/2
ay_t_un_l <- ginv(Ct_unadj) %*% Qy_t_un_d
ay_t_un_u <- ay_t_un_l+ginv(Ct_unadj)%*%Qy_t_un_i
ay_t_un_d<-(ay_t_un_l +ay_t_un_u)/2
ay_t_un_i<-abs(ay_t_un_l -ay_t_un_u)/2
ay_t_ad_l <- ginv(Ct_adj) %*% Qy_t_ad_d
ay_t_ad_u <- ay_t_ad_l+ginv(Ct_adj)%*%Qy_t_ad_i
ay_t_ad_d<-(ay_t_ad_l +ay_t_ad_u)/2
ay_t_ad_i<-abs(ay_t_ad_l -ay_t_ad_u)/2
##b unadj y
Cb_unadj <- t(Bmat) %*% (Bmat)

bcdy_b_un_d<-t(cbind(Bmat,Tmat,Rmat))%*%mu%*%ginv(t(mu)
                                                  %*%mu)%*%t(mu)%*%dy
bcdy_b_un_i<-t(cbind(Bmat,Tmat,Rmat))%*%mu%*%ginv(t(mu)
                                                  %*%mu)%*%t(mu)%*%iy
Qy_un_l<- (t(cbind(Bmat,Tmat,Rmat))%*%dy-bcdy_b_un_d)
Qy_un_u<- Qy_un_l+t(cbind(Bmat,Tmat,Rmat))%*%iy-bcdy_b_un_i
Qy_b_un_d<-(Qy_un_l[1:k,]+Qy_un_u[1:v,])/2
Qy_b_un_i<-abs(Qy_un_l[1:k,]-Qy_un_u[1:v,])/2
##B adj y
Z <- cbind(mu, Tmat,Rmat)
M <- diag(nrow(X)) - Z %*% ginv(t(Z) %*% Z) %*% t(Z)
Cb_adj <- t(Bmat) %*% M %*% Bmat
bcdy_b_ad_d<-t(Bmat)%*%cbind(mu,Tmat,Rmat)%*%ginv(t(cbind(mu,Tmat,Rmat))
                                                  %*%cbind(mu,Tmat,Rmat))%*%t(cbind(mu,Tmat,Rmat))%*%dy
bcdy_b_ad_i<-t(Bmat)%*%cbind(mu,Tmat,Rmat)%*%ginv(t(cbind(mu,Tmat,Rmat))
                                                  %*%cbind(mu,Tmat,Rmat))%*%t(cbind(mu,Tmat,Rmat))%*%iy
Qy_ad_l<- t(Bmat)%*%dy-bcdy_b_ad_d
Qy_ad_u<- Qy_ad_l+t(Bmat)%*%iy-bcdy_b_ad_i
Qy_b_ad_d<-(Qy_ad_l+Qy_ad_u)/2
Qy_b_ad_i<-abs(Qy_ad_l-Qy_ad_u)/2
ay_b_un_l <- ginv(Cb_unadj) %*% Qy_b_un_d
ay_b_un_u <- ay_b_un_l+ginv(Cb_unadj)%*%Qy_b_un_i
ay_b_un_d<-(ay_b_un_l +ay_b_un_u)/2
ay_b_un_i<-abs(ay_b_un_l-ay_b_un_u)/2
ay_b_ad_l <- ginv(Cb_adj) %*% Qy_b_ad_d
ay_b_ad_u <- ay_b_ad_l+ginv(Cb_adj)%*%Qy_b_ad_i
ay_b_ad_d<-(ay_b_ad_l +ay_b_ad_u)/2
ay_b_ad_i<-abs(ay_b_ad_l-ay_b_ad_u)/2
###rep
Cr_unadj <- t(Rmat) %*% Rmat
bcdy_r_un_d<-t(cbind(Rmat,Bmat,Tmat))%*%mu%*%ginv(t(mu)
                                                  %*%mu)%*%t(mu)%*%dy
bcdy_r_un_i<-t(cbind(Rmat,Bmat,Tmat))%*%mu%*%ginv(t(mu)
                                                  %*%mu)%*%t(mu)%*%iy
Qy_un_l<- (t(cbind(Rmat,Bmat,Tmat))%*%dy-bcdy_r_un_d)
Qy_un_u<- Qy_un_l+t(cbind(Rmat,Bmat,Tmat))%*%iy-bcdy_r_un_i
Qy_r_un_d<-(Qy_un_l[1:r,]+Qy_un_u[1:r,])/2
Qy_r_un_i<-abs(Qy_un_l[1:r,]-Qy_un_u[1:r,])/2

ay_r_l <- ginv(Cr_unadj) %*% Qy_r_un_d
ay_r_u <- ay_r_l+ginv(Cr_unadj)%*%Qy_r_un_i
ay_r_d<-(ay_r_l +ay_r_u)/2
ay_r_i<-abs(ay_r_l-ay_r_u)/2
##mue
Qy_mu_l<-sum(round(dy,6))
Qy_mu_u<-Qy_mu_l+sum(round(iy,7))
Qy_mu_d<-(Qy_mu_l+Qy_mu_u)/2
Qy_mu_i<-abs(Qy_mu_l-Qy_mu_u)/2
ay_mu_l<-Qy_mu_d*1/(b*k)
ay_mu_u<-ay_mu_l+Qy_mu_d/(b*k)
ay_mu_d<-(ay_mu_l+ay_mu_u)/2
ay_mu_i<-abs(ay_mu_l-ay_mu_u)/2
Qy_all_ad_d<-c(Qy_t_ad_d,Qy_mu_d,Qy_b_un_d,Qy_r_un_d)
Qy_all_ad_i<-c(Qy_t_ad_i,Qy_mu_i,Qy_b_un_i,Qy_r_un_i)
Qy_all_un_d<-c(Qy_t_un_d,Qy_mu_d,Qy_b_ad_d,Qy_r_un_d)
Qy_all_un_i<-c(Qy_t_un_i,Qy_mu_i,Qy_b_ad_i,Qy_r_un_i)
Qy_tmueb_ad_d<-c(Qy_t_ad_d,Qy_mu_d,Qy_b_un_d)
Qy_tmueb_ad_i<-c(Qy_t_ad_i,Qy_mu_i,Qy_b_un_i)
Qy_tmueb_un_d<-c(Qy_t_un_d,Qy_mu_d,Qy_b_ad_d)
Qy_tmueb_un_i<-c(Qy_t_un_i,Qy_mu_i,Qy_b_ad_i)
Qy_tmuer_ad_d<-c(Qy_t_ad_d,Qy_mu_d,Qy_r_un_d)
Qy_tmuer_ad_i<-c(Qy_t_ad_i,Qy_mu_i,Qy_r_un_i)
Qy_tmuer_un_d<-c(Qy_t_un_d,Qy_mu_d,Qy_r_un_d)
Qy_tmuer_un_i<-c(Qy_t_un_i,Qy_mu_i,Qy_r_un_i)
Qy_muebr_ad_d<-c(Qy_mu_d,Qy_b_ad_d,Qy_r_un_d)
Qy_muebr_ad_i<-c(Qy_mu_i,Qy_b_ad_i,Qy_r_un_i)
Qy_muebr_un_d<-c(Qy_mu_d,Qy_b_un_d,Qy_r_un_d)
Qy_muebr_un_i<-c(Qy_mu_i,Qy_b_un_i,Qy_r_un_i)
ay_all_ad_d<-c(ay_t_ad_d,ay_mu_d,ay_b_un_d,ay_r_d)
ay_all_ad_i<-c(ay_t_ad_i,ay_mu_i,ay_b_un_i,ay_r_i)
ay_all_un_d<-c(ay_t_un_d,ay_mu_d,ay_b_ad_d,ay_r_d)
ay_all_un_i<-c(ay_t_un_i,ay_mu_i,ay_b_ad_i,ay_r_i)
ay_tmueb_ad_d<-c(ay_t_ad_d,ay_mu_d,ay_b_un_d)
ay_tmueb_ad_i<-c(ay_t_ad_i,ay_mu_i,ay_b_un_i)
ay_tmueb_un_d<-c(ay_t_un_d,ay_mu_d,ay_b_ad_d)
ay_tmueb_un_i<-c(ay_t_un_i,ay_mu_i,ay_b_ad_i)
ay_tmuer_ad_d<-c(ay_t_ad_d,ay_mu_d,ay_r_d)
ay_tmuer_ad_i<-c(ay_t_ad_i,ay_mu_i,ay_r_i)
ay_tmuer_un_d<-c(ay_t_un_d,ay_mu_d,ay_r_d)
ay_tmuer_un_i<-c(ay_t_un_i,ay_mu_i,ay_r_i)
ay_muebr_ad_d<-c(ay_mu_d,ay_b_ad_d,ay_r_d)
ay_muebr_ad_i<-c(ay_mu_i,ay_b_ad_i,ay_r_i)
ay_muebr_un_d<-c(ay_mu_d,ay_b_un_d,ay_r_d)
ay_muebr_un_i<-c(ay_mu_i,ay_b_un_i,ay_r_i)
# Covariate#
##Trt unadj
Ct_unadj <- t(Tmat) %*% Tmat
bcdz_t_un_d<-t(cbind(Tmat,Bmat,Rmat))%*%mu%*%ginv(t(mu)
                                                  %*%mu)%*%t(mu)%*%dz
bcdz_t_un_i<-t(cbind(Tmat,Bmat,Rmat))%*%mu%*%ginv(t(mu)
                                                  %*%mu)%*%t(mu)%*%iz
Qz_un_l<- (t(cbind(Tmat,Bmat,Rmat))%*%dz-bcdz_t_un_d)
Qz_un_u<- Qz_un_l+t(cbind(Tmat,Bmat,Rmat))%*%iz-bcdz_t_un_i
Qz_t_un_d<-(Qz_un_l[1:v,]+Qz_un_u[1:v,])/2
Qz_t_un_i<-abs(Qz_un_l[1:v,]-Qz_un_u[1:v,])/2
##trt adj
Z <- cbind(mu, Bmat,Rmat)
M <- diag(nrow(X)) - Z %*% ginv(t(Z) %*% Z) %*% t(Z)
Ct_adj <- t(Tmat) %*% M %*% Tmat
bcdz_t_ad_d<-t(Tmat)%*%cbind(mu,Bmat,Rmat)%*%ginv(t(cbind(mu,Bmat,Rmat))
                                                  %*%cbind(mu,Bmat,Rmat))%*%t(cbind(mu,Bmat,Rmat))%*%dz
bcdz_t_ad_i<-t(Tmat)%*%cbind(mu,Bmat,Rmat)%*%ginv(t(cbind(mu,Bmat,Rmat))
                                                  %*%cbind(mu,Bmat,Rmat))%*%t(cbind(mu,Bmat,Rmat))%*%iz
Qz_ad_l<- t(Tmat)%*%dz-bcdz_t_ad_d
Qz_ad_u<- Qz_ad_l+t(Tmat)%*%iz-bcdz_t_ad_i
Qz_t_ad_d<-(Qz_ad_l+Qz_ad_u)/2
Qz_t_ad_i<-abs(Qz_ad_l-Qz_ad_u)/2
az_t_un_l <- ginv(Ct_unadj) %*% Qz_t_un_d
az_t_un_u <- az_t_un_l+ginv(Ct_unadj)%*%Qz_t_un_i
az_t_un_d<-(az_t_un_l +az_t_un_u)/2
az_t_un_i<-abs(az_t_un_l -az_t_un_u)/2
az_t_ad_l <- ginv(Ct_adj) %*% Qz_t_ad_d
az_t_ad_u <- az_t_ad_l+ginv(Ct_adj)%*%Qz_t_ad_i
az_t_ad_d<-(az_t_ad_l +az_t_ad_u)/2
az_t_ad_i<-abs(az_t_ad_l -az_t_ad_u)/2

##b unadj y
Cb_unadj <- t(Bmat) %*% Bmat
X1<-cbind(Bmat,Tmat)
X2 <- cbind(mu, Rmat)

bcdz_b_un_d <- t(Bmat) %*% X2 %*%ginv(t(X2)%*%X2) %*%t(X2) %*% dz
bcdz_b_un_i<-t(Bmat) %*% X2 %*%ginv(t(X2)%*%X2) %*%t(X2) %*% iz
Qz_b_un_l<-t(Bmat)%*%dz-bcdz_b_un_d
Qz_b_un_u<-Qz_b_un_l+ t(Bmat)%*%iz-bcdz_b_un_i
Qz_b_un_d<-(Qz_b_un_l+Qz_b_un_u)/2
Qz_b_un_i<-abs(Qz_b_un_u-Qz_b_un_l)/2
##B adj z
Z <- cbind(mu, Tmat,Rmat)
M <- diag(nrow(X)) - Z %*% ginv(t(Z) %*% Z) %*% t(Z)
Cb_adj <- t(Bmat) %*% M %*% Bmat
bcdz_b_ad_d<-t(Bmat)%*%cbind(mu,Tmat,Rmat)%*%ginv(t(cbind(mu,Tmat,Rmat))
                                                  %*%cbind(mu,Tmat,Rmat))%*%t(cbind(mu,Tmat,Rmat))%*%dz
bcdz_b_ad_i<-t(Bmat)%*%cbind(mu,Tmat,Rmat)%*%ginv(t(cbind(mu,Tmat,Rmat))
                                                  %*%cbind(mu,Tmat,Rmat))%*%t(cbind(mu,Tmat,Rmat))%*%iz
Qz_ad_l<- t(Bmat)%*%dz-bcdz_b_ad_d
Qz_ad_u<- Qz_ad_l+t(Bmat)%*%iz-bcdz_b_ad_i
Qz_b_ad_d<-(Qz_ad_l+Qz_ad_u)/2
Qz_b_ad_i<-abs(Qz_ad_l-Qz_ad_u)/2
az_b_ad_l <- ginv(Cb_adj) %*% Qz_b_ad_d
az_b_ad_u <- az_b_ad_l+ginv(Cb_adj)%*%Qz_b_ad_i
az_b_ad_d<-(az_b_ad_l +az_b_ad_u)/2
az_b_ad_i<-abs(az_b_ad_l-az_b_ad_u)/2
az_b_un_l <- ginv(Cb_unadj) %*% Qz_b_un_d
az_b_un_u <- az_b_un_l+ginv(Cb_unadj)%*%Qz_b_un_i
az_b_un_d<-(az_b_un_l +az_b_un_u)/2
az_b_un_i<-abs(az_b_un_l-az_b_un_u)/2
###rep
Cr_unadj <- t(Rmat) %*% Rmat
bcdz_r_un_d<-t(cbind(Rmat,Bmat,Tmat))%*%mu%*%ginv(t(mu)
                                                  %*%mu)%*%t(mu)%*%dz
bcdz_r_un_i<-t(cbind(Rmat,Bmat,Tmat))%*%mu%*%ginv(t(mu)
                                                  %*%mu)%*%t(mu)%*%iz
Qz_un_l<- (t(cbind(Rmat,Bmat,Tmat))%*%dz-bcdz_r_un_d)
Qz_un_u<- Qz_un_l+t(cbind(Rmat,Bmat,Tmat))%*%iz-bcdz_r_un_i
Qz_r_un_d<-(Qz_un_l[1:r,]+Qz_un_u[1:r,])/2
Qz_r_un_i<-abs(Qz_un_l[1:r,]-Qz_un_u[1:r,])/2

az_r_l <- ginv(Cr_unadj) %*% Qz_r_un_d
az_r_u <- az_r_l+ginv(Cr_unadj)%*%Qz_r_un_i
az_r_d<-(az_r_l +az_r_u)/2
az_r_i<-abs(az_r_l-az_r_u)/2
##mue
Qz_mu_l<-sum(round(dz,6))
Qz_mu_u<-Qz_mu_l+sum(round(iz,7))
Qz_mu_d<-(Qz_mu_l+Qz_mu_u)/2
Qz_mu_i<-abs(Qz_mu_l-Qz_mu_u)/2
az_mu_l<-Qz_mu_d*1/(b*k)
az_mu_u<-az_mu_l+Qz_mu_d/(b*k)
az_mu_d<-(az_mu_l+az_mu_u)/2
az_mu_i<-abs(az_mu_l-az_mu_u)/2
####qz,az
Qz_all_ad_d<-c(Qz_t_ad_d,Qz_mu_d,Qz_b_un_d,Qz_r_un_d)
Qz_all_ad_i<-c(Qz_t_ad_i,Qz_mu_i,Qz_b_un_i,Qz_r_un_i)
Qz_all_un_d<-c(Qz_t_un_d,Qz_mu_d,Qz_b_ad_d,Qz_r_un_d)
Qz_all_un_i<-c(Qz_t_un_i,Qz_mu_i,Qz_b_ad_i,Qz_r_un_i)
Qz_tmueb_ad_d<-c(Qz_t_ad_d,Qz_mu_d,Qz_b_un_d)
Qz_tmueb_ad_i<-c(Qz_t_ad_i,Qz_mu_i,Qz_b_un_i)
Qz_tmueb_un_d<-c(Qz_t_un_d,Qz_mu_d,Qz_b_ad_d)
Qz_tmueb_un_i<-c(Qz_t_un_i,Qz_mu_i,Qz_b_ad_i)
Qz_tmuer_ad_d<-c(Qz_t_ad_d,Qz_mu_d,Qz_r_un_d)
Qz_tmuer_ad_i<-c(Qz_t_ad_i,Qz_mu_i,Qz_r_un_i)
Qz_tmuer_un_d<-c(Qz_t_un_d,Qz_mu_d,Qz_r_un_d)
Qz_tmuer_un_i<-c(Qz_t_un_i,Qz_mu_i,Qz_r_un_i)
Qz_muebr_ad_d<-c(Qz_mu_d,Qz_b_ad_d,Qz_r_un_d)
Qz_muebr_ad_i<-c(Qz_mu_i,Qz_b_ad_i,Qz_r_un_i)
Qz_muebr_un_d<-c(Qz_mu_d,Qz_b_un_d,Qz_r_un_d)
Qz_muebr_un_i<-c(Qz_mu_i,Qz_b_un_i,Qz_r_un_i)
az_all_ad_d<-c(az_t_ad_d,az_mu_d,az_b_un_d,az_r_d)
az_all_ad_i<-c(az_t_ad_i,az_mu_i,az_b_un_i,az_r_i)
az_all_un_d<-c(az_t_un_d,az_mu_d,az_b_ad_d,az_r_d)
az_all_un_i<-c(az_t_un_i,az_mu_i,az_b_ad_i,az_r_i)
az_tmueb_ad_d<-c(az_t_ad_d,az_mu_d,az_b_un_d)
az_tmueb_ad_i<-c(az_t_ad_i,az_mu_i,az_b_un_i)
az_tmueb_un_d<-c(az_t_un_d,az_mu_d,az_b_ad_d)
az_tmueb_un_i<-c(az_t_un_i,az_mu_i,az_b_ad_i)
az_tmuer_ad_d<-c(az_t_ad_d,az_mu_d,az_r_d)
az_tmuer_ad_i<-c(az_t_ad_i,az_mu_i,az_r_i)
az_tmuer_un_d<-c(az_t_un_d,az_mu_d,az_r_d)
az_tmuer_un_i<-c(az_t_un_i,az_mu_i,az_r_i)
az_muebr_ad_d<-c(az_mu_d,az_b_ad_d,az_r_d)
az_muebr_ad_i<-c(az_mu_i,az_b_ad_i,az_r_i)
az_muebr_un_d<-c(az_mu_d,az_b_un_d,az_r_d)
az_muebr_un_i<-c(az_mu_i,az_b_un_i,az_r_i)

##adj trt unadj blk   qz*az
QZAZ_all_ad_l<-c(Qz_t_ad_d*az_t_ad_d, Qz_mu_d*az_mu_d,Qz_b_un_d*az_b_un_d,Qz_r_un_d*az_r_d)
QZAZ_all_ad_u<-QZAZ_all_ad_l+c(Qz_t_ad_d*az_t_ad_i,Qz_mu_d*az_mu_i,Qz_b_un_d*az_b_un_i,Qz_r_un_d*az_r_i
)+c(Qz_t_ad_i*az_t_ad_d, Qz_mu_i*az_mu_d,Qz_b_un_i*az_b_un_d,Qz_r_un_i*az_r_d
)+c(Qz_t_ad_i*az_t_ad_i, Qz_mu_i*az_mu_i,Qz_b_un_i*az_b_un_i,Qz_r_un_i*az_r_i
)
QZAZ_all_ad_d<-sum((QZAZ_all_ad_l+QZAZ_all_ad_u)/2)
QZAZ_all_ad_i<-sum(abs(QZAZ_all_ad_u-QZAZ_all_ad_l)/2)
QZAZ_tmueb_ad_l<-c(Qz_t_ad_d*az_t_ad_d, Qz_mu_d*az_mu_d,Qz_b_un_d*az_b_un_d)
QZAZ_tmueb_ad_u<-QZAZ_tmueb_ad_l+c(Qz_t_ad_d*az_t_ad_i,Qz_mu_d*az_mu_i,Qz_b_un_d*az_b_un_i)+
  c(Qz_t_ad_i*az_t_ad_d, Qz_mu_i*az_mu_d,Qz_b_un_i*az_b_un_d)+
  c(Qz_t_ad_i*az_t_ad_i, Qz_mu_i*az_mu_i,Qz_b_un_i*az_b_un_i)
QZAZ_tmueb_ad_d<-sum((QZAZ_tmueb_ad_l+QZAZ_tmueb_ad_u)/2)
QZAZ_tmueb_ad_i<-sum(abs((QZAZ_tmueb_ad_u-QZAZ_tmueb_ad_l)/2))
QZAZ_tmuer_ad_l<-c(Qz_t_ad_d*az_t_ad_d, Qz_mu_d*az_mu_d,Qz_r_un_d*az_r_d)
QZAZ_tmuer_ad_u<-QZAZ_tmuer_ad_l+c(Qz_t_ad_d*az_t_ad_i,Qz_mu_d*az_mu_i,Qz_r_un_d*az_r_i)+
  c(Qz_t_ad_i*az_t_ad_d, Qz_mu_i*az_mu_d,Qz_r_un_i*az_r_d)+
  c(Qz_t_ad_i*az_t_ad_i, Qz_mu_i*az_mu_i,Qz_r_un_i*az_r_i)
QZAZ_tmuer_ad_d<-sum((QZAZ_tmuer_ad_l+QZAZ_tmuer_ad_u)/2)
QZAZ_tmuer_ad_i<-sum(abs((QZAZ_tmuer_ad_u-QZAZ_tmuer_ad_l)/2))
QZAZ_muebr_un_l<-c(Qz_mu_d*az_mu_d,Qz_b_un_d*az_b_un_d,Qz_r_un_d*az_r_d)
QZAZ_muebr_un_u<-QZAZ_muebr_un_l+c(Qz_mu_d*az_mu_i,Qz_b_un_d*az_b_un_i,Qz_r_un_d*az_r_i
)+c( Qz_mu_i*az_mu_d,Qz_b_un_i*az_b_un_d,Qz_r_un_i*az_r_d
)+c(Qz_mu_i*az_mu_i,Qz_b_un_i*az_b_un_i,Qz_r_un_i*az_r_i
)
QZAZ_muebr_un_d<-sum(QZAZ_muebr_un_l+QZAZ_muebr_un_u)/2
QZAZ_muebr_un_i<-sum(abs(QZAZ_muebr_un_u-QZAZ_muebr_un_l)/2)

##adj trt unadj blk  ######qz*ay
QZAy_all_ad_l<-c(Qz_t_ad_d*ay_t_ad_d, Qz_mu_d*ay_mu_d,Qz_b_un_d*ay_b_un_d,Qz_r_un_d*ay_r_d)
QZAy_all_ad_u<-QZAy_all_ad_l+c(Qz_t_ad_d*ay_t_ad_i,Qz_mu_d*ay_mu_i,Qz_b_un_d*ay_b_un_i,Qz_r_un_d*ay_r_i
)+c(Qz_t_ad_i*ay_t_ad_d, Qz_mu_i*ay_mu_d,Qz_b_un_i*ay_b_un_d,Qz_r_un_i*ay_r_d
)+c(Qz_t_ad_i*ay_t_ad_i, Qz_mu_i*ay_mu_i,Qz_b_un_i*ay_b_un_i,Qz_r_un_i*ay_r_i
)
QZAy_all_ad_d<-sum((QZAy_all_ad_l+QZAy_all_ad_u)/2)
QZAy_all_ad_i<-sum(abs(QZAy_all_ad_u-QZAy_all_ad_l)/2)
QZAy_tmueb_ad_l<-c(Qz_t_ad_d*ay_t_ad_d, Qz_mu_d*ay_mu_d,Qz_b_un_d*ay_b_un_d)
QZAy_tmueb_ad_u<-QZAy_tmueb_ad_l+c(Qz_t_ad_d*ay_t_ad_i,Qz_mu_d*ay_mu_i,Qz_b_un_d*ay_b_un_i)+
  c(Qz_t_ad_i*ay_t_ad_d, Qz_mu_i*ay_mu_d,Qz_b_un_i*ay_b_un_d)+
  c(Qz_t_ad_i*ay_t_ad_i, Qz_mu_i*ay_mu_i,Qz_b_un_i*ay_b_un_i)
QZAy_tmueb_ad_d<-sum((QZAy_tmueb_ad_l+QZAy_tmueb_ad_u)/2)
QZAy_tmueb_ad_i<-sum(abs((QZAy_tmueb_ad_u-QZAy_tmueb_ad_l)/2))

QZAy_tmuer_ad_l<-c(Qz_t_ad_d*ay_t_ad_d, Qz_mu_d*ay_mu_d,Qz_r_un_d*ay_r_d)
QZAy_tmuer_ad_u<-QZAy_tmuer_ad_l+c(Qz_t_ad_d*ay_t_ad_i,Qz_mu_d*ay_mu_i,Qz_r_un_d*ay_r_i)+
  c(Qz_t_ad_i*ay_t_ad_d, Qz_mu_i*ay_mu_d,Qz_r_un_i*ay_r_d)+
  c(Qz_t_ad_i*ay_t_ad_i, Qz_mu_i*ay_mu_i,Qz_r_un_i*ay_r_i)
QZAy_tmuer_ad_d<-sum((QZAy_tmuer_ad_l+QZAy_tmuer_ad_u)/2)
QZAy_tmuer_ad_i<-sum(abs((QZAy_tmuer_ad_u-QZAy_tmuer_ad_l)/2))
QZAy_muebr_un_l<-c(Qz_mu_d*ay_mu_d,Qz_b_un_d*ay_b_un_d,Qz_r_un_d*ay_r_d)
QZAy_muebr_un_u<-QZAy_muebr_un_l+c(Qz_mu_d*ay_mu_i,Qz_b_un_d*ay_b_un_i,Qz_r_un_d*ay_r_i
)+c( Qz_mu_i*ay_mu_d,Qz_b_un_i*ay_b_un_d,Qz_r_un_i*ay_r_d
)+c(Qz_mu_i*ay_mu_i,Qz_b_un_i*ay_b_un_i,Qz_r_un_i*ay_r_i
)
QZAy_muebr_un_d<-sum(QZAy_muebr_un_l+QZAy_muebr_un_u)/2
QZAy_muebr_un_i<-sum(abs(QZAy_muebr_un_u-QZAy_muebr_un_l)/2)

##unadj trt adj blk  ######qz*az
QZAZ_all_un_l<-c(Qz_t_un_d*az_t_un_d, Qz_mu_d*az_mu_d,Qz_b_ad_d*az_b_ad_d,Qz_r_un_d*az_r_d)
QZAZ_all_un_u<-QZAZ_all_un_l+c(Qz_t_un_d*az_t_un_i,Qz_mu_d*az_mu_i,Qz_b_ad_d*az_b_ad_i,Qz_r_un_d*az_r_i
)+c(Qz_t_un_i*az_t_un_d, Qz_mu_i*az_mu_d,Qz_b_ad_i*az_b_ad_d,Qz_r_un_i*az_r_d
)+c(Qz_t_un_i*az_t_un_i, Qz_mu_i*az_mu_i,Qz_b_ad_i*az_b_ad_i,Qz_r_un_i*az_r_i
)
QZAZ_all_un_d<-sum((QZAZ_all_un_l+QZAZ_all_un_u)/2)
QZAZ_all_un_i<-sum(abs(QZAZ_all_un_u-QZAZ_all_un_l)/2)
QZAZ_tmueb_un_l<-c(Qz_t_un_d*az_t_un_d, Qz_mu_d*az_mu_d,Qz_b_ad_d*az_b_ad_d)
QZAZ_tmueb_un_u<-QZAZ_tmueb_un_l+c(Qz_t_un_d*az_t_un_i,Qz_mu_d*az_mu_i,Qz_b_ad_d*az_b_ad_i)+
  c(Qz_t_un_i*az_t_un_d, Qz_mu_i*az_mu_d,Qz_b_ad_i*az_b_ad_d)+
  c(Qz_t_un_i*az_t_un_i, Qz_mu_i*az_mu_i,Qz_b_ad_i*az_b_ad_i)
QZAZ_tmueb_un_d<-sum((QZAZ_tmueb_un_l+QZAZ_tmueb_un_u)/2)
QZAZ_tmueb_un_i<-sum(abs((QZAZ_tmueb_un_u-QZAZ_tmueb_un_l)/2))

QZAZ_tmuer_un_l<-c(Qz_t_un_d*az_t_un_d, Qz_mu_d*az_mu_d,Qz_r_un_d*az_r_d)
QZAZ_tmuer_un_u<-QZAZ_tmuer_un_l+c(Qz_t_un_d*az_t_un_i,Qz_mu_d*az_mu_i,Qz_r_un_d*az_r_i)+
  c(Qz_t_un_i*az_t_un_d, Qz_mu_i*az_mu_d,Qz_r_un_i*az_r_d)+
  c(Qz_t_un_i*az_t_un_i, Qz_mu_i*az_mu_i,Qz_r_un_i*az_r_i)
QZAZ_tmuer_un_d<-sum((QZAZ_tmuer_un_l+QZAZ_tmuer_un_u)/2)
QZAZ_tmuer_un_i<-sum(abs((QZAZ_tmuer_un_u-QZAZ_tmuer_un_l)/2))
QZAZ_muebr_ad_l<-c(Qz_mu_d*az_mu_d,Qz_b_ad_d*az_b_ad_d,Qz_r_un_d*az_r_d)
QZAZ_muebr_ad_u<-QZAZ_muebr_ad_l+c(Qz_mu_d*az_mu_i,Qz_b_ad_d*az_b_ad_i,Qz_r_un_d*az_r_i
)+c( Qz_mu_i*az_mu_d,Qz_b_ad_i*az_b_ad_d,Qz_r_un_i*az_r_d
)+c(Qz_mu_i*az_mu_i,Qz_b_ad_i*az_b_ad_i,Qz_r_un_i*az_r_i
)
QZAZ_muebr_ad_d<-sum(QZAZ_muebr_ad_l+QZAZ_muebr_ad_u)/2
QZAZ_muebr_ad_i<-sum(abs(QZAZ_muebr_ad_u-QZAZ_muebr_ad_l)/2)
##unadj trt adj blk  ######qz*ay
QZAy_all_un_l<-c(Qz_t_un_d*ay_t_un_d, Qz_mu_d*ay_mu_d,Qz_b_ad_d*ay_b_ad_d,Qz_r_un_d*ay_r_d)
QZAy_all_un_u<-QZAy_all_un_l+c(Qz_t_un_d*ay_t_un_i,Qz_mu_d*ay_mu_i,Qz_b_ad_d*ay_b_ad_i,Qz_r_un_d*ay_r_i
)+c(Qz_t_un_i*ay_t_un_d, Qz_mu_i*ay_mu_d,Qz_b_ad_i*ay_b_ad_d,Qz_r_un_i*ay_r_d
)+c(Qz_t_un_i*ay_t_un_i, Qz_mu_i*ay_mu_i,Qz_b_ad_i*ay_b_ad_i,Qz_r_un_i*ay_r_i
)
QZAy_all_un_d<-sum((QZAy_all_un_l+QZAy_all_un_u)/2)
QZAy_all_un_i<-sum(abs(QZAy_all_un_u-QZAy_all_un_l)/2)

QZAy_tmueb_un_l<-c(Qz_t_un_d*ay_t_un_d, Qz_mu_d*ay_mu_d,Qz_b_ad_d*ay_b_ad_d)
QZAy_tmueb_un_u<-QZAy_tmueb_un_l+c(Qz_t_un_d*ay_t_un_i,Qz_mu_d*ay_mu_i,Qz_b_ad_d*ay_b_ad_i)+
  c(Qz_t_un_i*ay_t_un_d, Qz_mu_i*ay_mu_d,Qz_b_ad_i*ay_b_ad_d)+
  c(Qz_t_un_i*ay_t_un_i, Qz_mu_i*ay_mu_i,Qz_b_ad_i*ay_b_ad_i)
QZAy_tmueb_un_d<-sum((QZAy_tmueb_un_l+QZAy_tmueb_un_u)/2)
QZAy_tmueb_un_i<-sum(abs((QZAy_tmueb_un_u-QZAy_tmueb_un_l)/2))

QZAy_tmuer_un_l<-c(Qz_t_un_d*ay_t_un_d, Qz_mu_d*ay_mu_d,Qz_r_un_d*ay_r_d)
QZAy_tmuer_un_u<-QZAy_tmuer_un_l+c(Qz_t_un_d*ay_t_un_i,Qz_mu_d*ay_mu_i,Qz_r_un_d*ay_r_i)+
  c(Qz_t_un_i*ay_t_un_d, Qz_mu_i*ay_mu_d,Qz_r_un_i*ay_r_d)+
  c(Qz_t_un_i*ay_t_un_i, Qz_mu_i*ay_mu_i,Qz_r_un_i*ay_r_i)
QZAy_tmuer_un_d<-sum((QZAy_tmuer_un_l+QZAy_tmuer_un_u)/2)
QZAy_tmuer_un_i<-sum(abs((QZAy_tmuer_un_u-QZAy_tmuer_un_l)/2))
QZAy_muebr_ad_l<-c(Qz_mu_d*ay_mu_d,Qz_b_ad_d*ay_b_ad_d,Qz_r_un_d*ay_r_d)
QZAy_muebr_ad_u<-QZAy_muebr_ad_l+c(Qz_mu_d*ay_mu_i,Qz_b_ad_d*ay_b_ad_i,Qz_r_un_d*ay_r_i
)+c( Qz_mu_i*ay_mu_d,Qz_b_ad_i*ay_b_ad_d,Qz_r_un_i*ay_r_d
)+c(Qz_mu_i*ay_mu_i,Qz_b_ad_i*ay_b_ad_i,Qz_r_un_i*ay_r_i
)
QZAy_muebr_ad_d<-sum(QZAy_muebr_ad_l+QZAy_muebr_ad_u)/2
QZAy_muebr_ad_i<-sum(abs(QZAy_muebr_ad_u-QZAy_muebr_ad_l)/2)
######GAMMA
EXX_all_ad_l<-zpz_d-QZAZ_all_ad_d
EXX_all_ad_u<-EXX_all_ad_l+zpz_i-QZAZ_all_ad_i
EXX_all_ad_d<-(EXX_all_ad_l+EXX_all_ad_u)/2
EXX_all_ad_i<-abs(EXX_all_ad_u-EXX_all_ad_l)/2
EXX_tmueb_ad_l<-zpz_d-QZAZ_tmueb_ad_d
EXX_tmueb_ad_u<-EXX_tmueb_ad_l+zpz_i-QZAZ_tmueb_ad_i
EXX_tmueb_ad_d<-(EXX_tmueb_ad_l+EXX_tmueb_ad_u)/2
EXX_tmueb_ad_i<-abs(EXX_tmueb_ad_u-EXX_tmueb_ad_l)/2
EXX_tmuer_ad_l<-zpz_d-QZAZ_tmuer_ad_d
EXX_tmuer_ad_u<-EXX_tmuer_ad_l+zpz_i-QZAZ_tmuer_ad_i
EXX_tmuer_ad_d<-(EXX_tmuer_ad_l+EXX_tmuer_ad_u)/2
EXX_tmuer_ad_i<-abs(EXX_tmuer_ad_u-EXX_tmuer_ad_l)/2
EXX_muebr_un_l<-zpz_d-QZAZ_muebr_un_d
EXX_muebr_un_u<-EXX_muebr_un_l+zpz_i-QZAZ_muebr_un_i
EXX_muebr_un_d<-(EXX_muebr_un_l+EXX_muebr_un_u)/2
EXX_muebr_un_i<-abs(EXX_muebr_un_u-EXX_muebr_un_l)/2
EXy_all_ad_l<-zpy_d-QZAy_all_ad_d
EXy_all_ad_u<-EXy_all_ad_l+zpy_i-QZAy_all_ad_i
EXy_all_ad_d<-(EXy_all_ad_l+EXy_all_ad_u)/2
EXy_all_ad_i<-abs(EXy_all_ad_u-EXy_all_ad_l)/2
EXy_tmueb_ad_l<-zpy_d-QZAy_tmueb_ad_d
EXy_tmueb_ad_u<-EXy_tmueb_ad_l+zpy_i-QZAy_tmueb_ad_i
EXy_tmueb_ad_d<-(EXy_tmueb_ad_l+EXy_tmueb_ad_u)/2
EXy_tmueb_ad_i<-abs(EXy_tmueb_ad_u-EXy_tmueb_ad_l)/2
EXy_tmuer_ad_l<-zpy_d-QZAy_tmuer_ad_d
EXy_tmuer_ad_u<-EXy_tmuer_ad_l+zpy_i-QZAy_tmuer_ad_i
EXy_tmuer_ad_d<-(EXy_tmuer_ad_l+EXy_tmuer_ad_u)/2
EXy_tmuer_ad_i<-abs(EXy_tmuer_ad_u-EXy_tmuer_ad_l)/2
EXy_muebr_un_l<-zpy_d-QZAy_muebr_un_d
EXy_muebr_un_u<-EXy_muebr_un_l+zpy_i-QZAy_muebr_un_i
EXy_muebr_un_d<-(EXy_muebr_un_l+EXy_muebr_un_u)/2
EXy_muebr_un_i<-abs(EXy_muebr_un_u-EXy_muebr_un_l)/2
EXX_all_un_l<-zpz_d-QZAZ_all_un_d
EXX_all_un_u<-EXX_all_un_l+zpz_i-QZAZ_all_un_i
EXX_all_un_d<-(EXX_all_un_l+EXX_all_un_u)/2
EXX_all_un_i<-abs(EXX_all_un_u-EXX_all_un_l)/2
EXX_tmueb_un_l<-zpz_d-QZAZ_tmueb_un_d
EXX_tmueb_un_u<-EXX_tmueb_un_l+zpz_i-QZAZ_tmueb_un_i
EXX_tmueb_un_d<-(EXX_tmueb_un_l+EXX_tmueb_un_u)/2
EXX_tmueb_un_i<-abs(EXX_tmueb_un_u-EXX_tmueb_un_l)/2
EXX_tmuer_un_l<-zpz_d-QZAZ_tmuer_un_d
EXX_tmuer_un_u<-EXX_tmuer_un_l+zpz_i-QZAZ_tmuer_un_i
EXX_tmuer_un_d<-(EXX_tmuer_un_l+EXX_tmuer_un_u)/2
EXX_tmuer_un_i<-abs(EXX_tmuer_un_u-EXX_tmuer_un_l)/2
EXX_muebr_ad_l<-zpz_d-QZAZ_muebr_ad_d
EXX_muebr_ad_u<-EXX_muebr_ad_l+zpz_i-QZAZ_muebr_ad_i
EXX_muebr_ad_d<-(EXX_muebr_ad_l+EXX_muebr_ad_u)/2
EXX_muebr_ad_i<-abs(EXX_muebr_ad_u-EXX_muebr_ad_l)/2
EXy_all_un_l<-zpy_d-QZAy_all_un_d
EXy_all_un_u<-EXy_all_un_l+zpy_i-QZAy_all_un_i
EXy_all_un_d<-(EXy_all_un_l+EXy_all_un_u)/2
EXy_all_un_i<-abs(EXy_all_un_u-EXy_all_un_l)/2
EXy_tmueb_un_l<-zpy_d-QZAy_tmueb_un_d
EXy_tmueb_un_u<-EXy_tmueb_un_l+zpy_i-QZAy_tmueb_un_i
EXy_tmueb_un_d<-(EXy_tmueb_un_l+EXy_tmueb_un_u)/2
EXy_tmueb_un_i<-abs(EXy_tmueb_un_u-EXy_tmueb_un_l)/2
EXy_tmuer_un_l<-zpy_d-QZAy_tmuer_un_d
EXy_tmuer_un_u<-EXy_tmuer_un_l+zpy_i-QZAy_tmuer_un_i
EXy_tmuer_un_d<-(EXy_tmuer_un_l+EXy_tmuer_un_u)/2
EXy_tmuer_un_i<-abs(EXy_tmuer_un_u-EXy_tmuer_un_l)/2
EXy_muebr_ad_l<-zpy_d-QZAy_muebr_ad_d
EXy_muebr_ad_u<-EXy_muebr_ad_l+zpy_i-QZAy_muebr_ad_i
EXy_muebr_ad_d<-(EXy_muebr_ad_l+EXy_muebr_ad_u)/2
EXy_muebr_ad_i<-abs(EXy_muebr_ad_u-EXy_muebr_ad_l)/2

gamma_1_ad_l<-EXy_all_ad_d/EXX_all_ad_d
gamma_1_ad_u<-gamma_1_ad_l+(EXy_all_ad_d*EXX_all_ad_i+EXy_all_ad_i*EXX_all_ad_d+
                              EXy_all_ad_i*EXX_all_ad_i)/(EXX_all_ad_d*EXX_all_ad_d)
gamma_1_ad_d<-(gamma_1_ad_l+gamma_1_ad_u)/2
gamma_1_ad_i<-abs(gamma_1_ad_u-gamma_1_ad_l)/2
gamma_2_ad_l<-EXy_tmueb_ad_d/EXX_tmueb_ad_d
gamma_2_ad_u<-gamma_2_ad_l+(EXy_tmueb_ad_d*EXX_tmueb_ad_i+EXy_tmueb_ad_i*EXX_tmueb_ad_d+
                              EXy_tmueb_ad_i*EXX_tmueb_ad_i)/(EXX_tmueb_ad_d*EXX_tmueb_ad_d)
gamma_2_ad_d<-(gamma_2_ad_l+gamma_2_ad_u)/2
gamma_2_ad_i<-abs(gamma_2_ad_u-gamma_2_ad_l)/2
gamma_3_ad_l<-EXy_tmuer_ad_d/EXX_tmuer_ad_d
gamma_3_ad_u<-gamma_3_ad_l+(EXy_tmuer_ad_d*EXX_tmuer_ad_i+EXy_tmuer_ad_i*EXX_tmuer_ad_d+
                              EXy_tmuer_ad_i*EXX_tmuer_ad_i)/(EXX_tmuer_ad_d*EXX_tmuer_ad_d)
gamma_3_ad_d<-(gamma_3_ad_l+gamma_3_ad_u)/2
gamma_3_ad_i<-abs(gamma_3_ad_u-gamma_3_ad_l)/2
gamma_4_ad_l<-EXy_muebr_un_d/EXX_muebr_un_d
gamma_4_ad_u<-gamma_4_ad_l+(EXy_muebr_un_d*EXX_muebr_un_i+EXy_muebr_un_i*EXX_muebr_un_d+
                              EXy_muebr_un_i*EXX_muebr_un_i)/(EXX_muebr_un_d*EXX_muebr_un_d)
gamma_4_ad_d<-(gamma_4_ad_l+gamma_4_ad_u)/2
gamma_4_ad_i<-abs(gamma_4_ad_u-gamma_4_ad_l)/2

gamma_1_un_l<-EXy_all_un_d/EXX_all_un_d
gamma_1_un_u<-gamma_1_un_l+(EXy_all_un_d*EXX_all_un_i+EXy_all_un_i*EXX_all_un_d+
                              EXy_all_un_i*EXX_all_un_i)/(EXX_all_un_d*EXX_all_un_d)
gamma_1_un_d<-(gamma_1_un_l+gamma_1_un_u)/2
gamma_1_un_i<-abs(gamma_1_un_u-gamma_1_un_l)/2
gamma_2_un_l<-EXy_tmueb_un_d/EXX_tmueb_un_d
gamma_2_un_u<-gamma_2_un_l+(EXy_tmueb_un_d*EXX_tmueb_un_i+EXy_tmueb_un_i*EXX_tmueb_un_d+
                              EXy_tmueb_un_i*EXX_tmueb_un_i)/(EXX_tmueb_un_d*EXX_tmueb_un_d)
gamma_2_un_d<-(gamma_2_un_l+gamma_2_un_u)/2
gamma_2_un_i<-abs(gamma_2_un_u-gamma_2_un_l)/2
gamma_4_un_l<-EXy_tmuer_un_d/EXX_tmuer_un_d
gamma_4_un_u<-gamma_4_un_l+(EXy_tmuer_un_d*EXX_tmuer_un_i+EXy_tmuer_un_i*EXX_tmuer_un_d+
                              EXy_tmuer_un_i*EXX_tmuer_un_i)/(EXX_tmuer_un_d*EXX_tmuer_un_d)
gamma_4_un_d<-(gamma_4_un_l+gamma_4_un_u)/2
gamma_4_un_i<-abs(gamma_4_un_u-gamma_4_un_l)/2
gamma_3_un_l<-EXy_muebr_ad_d/EXX_muebr_ad_d
gamma_3_un_u<-gamma_3_un_l+(EXy_muebr_ad_d*EXX_muebr_ad_i+EXy_muebr_ad_i*EXX_muebr_ad_d+
                              EXy_muebr_ad_i*EXX_muebr_ad_i)/(EXX_muebr_ad_d*EXX_muebr_ad_d)
gamma_3_un_d<-(gamma_3_un_l+gamma_3_un_u)/2
gamma_3_un_i<-abs(gamma_3_un_u-gamma_3_un_l)/2
############gamma az (trt adj, b un)
gaz_all_ad_l<-az_all_ad_d%*%gamma_1_ad_d
gaz_all_ad_u<-gaz_all_ad_l+az_all_ad_d%*%gamma_1_ad_i+az_all_ad_i%*%gamma_1_ad_d+
  az_all_ad_i%*%gamma_1_ad_i
gaz_all_ad_d<-(gaz_all_ad_l+gaz_all_ad_u)/2
gaz_all_ad_i<-abs(gaz_all_ad_u-gaz_all_ad_l)/2
gaz_all_un_l<-az_all_un_d%*%gamma_1_un_d
gaz_all_un_u<-gaz_all_un_l+az_all_un_d%*%gamma_1_un_i+az_all_un_i%*%gamma_1_un_d+
  az_all_un_i%*%gamma_1_un_i
gaz_all_un_d<-(gaz_all_un_l+gaz_all_un_u)/2
gaz_all_un_i<-abs(gaz_all_un_u-gaz_all_un_l)/2
acap_all_ad_l<-ay_all_ad_d-gaz_all_ad_d
acap_all_ad_u<-acap_all_ad_l+ay_all_ad_i-gaz_all_ad_i
acap_all_ad_d<-(acap_all_ad_l+acap_all_ad_u)/2
acap_all_ad_i<-abs(acap_all_ad_u-acap_all_ad_l)/2
acap_all_un_l<-ay_all_un_d-gaz_all_un_d
acap_all_un_u<-acap_all_un_l+ay_all_un_i-gaz_all_un_i
acap_all_un_d<-(acap_all_un_l+acap_all_un_u)/2
acap_all_un_i<-abs(acap_all_un_u-acap_all_un_l)/2

gaz_tmueb_ad_l<-az_tmueb_ad_d%*%gamma_2_ad_d
gaz_tmueb_ad_u<-gaz_tmueb_ad_l+az_tmueb_ad_d%*%gamma_2_ad_i+az_tmueb_ad_i%*%gamma_2_ad_d+
  az_tmueb_ad_i%*%gamma_2_ad_i
gaz_tmueb_ad_d<-(gaz_tmueb_ad_l+gaz_tmueb_ad_u)/2
gaz_tmueb_ad_i<-abs(gaz_tmueb_ad_u-gaz_tmueb_ad_l)/2
gaz_tmueb_un_l<-az_tmueb_un_d%*%gamma_2_un_d
gaz_tmueb_un_u<-gaz_tmueb_un_l+az_tmueb_un_d%*%gamma_2_un_i+az_tmueb_un_i%*%gamma_2_un_d+
  az_tmueb_un_i%*%gamma_2_un_i
gaz_tmueb_un_d<-(gaz_tmueb_un_l+gaz_tmueb_un_u)/2
gaz_tmueb_un_i<-abs(gaz_tmueb_un_u-gaz_tmueb_un_l)/2
acap_tmueb_ad_l<-ay_tmueb_ad_d-gaz_tmueb_ad_d
acap_tmueb_ad_u<-acap_tmueb_ad_l+ay_tmueb_ad_i-gaz_tmueb_ad_i
acap_tmueb_ad_d<-(acap_tmueb_ad_l+acap_tmueb_ad_u)/2
acap_tmueb_ad_i<-abs(acap_tmueb_ad_u-acap_tmueb_ad_l)/2
acap_tmueb_un_l<-ay_tmueb_un_d-gaz_tmueb_un_d
acap_tmueb_un_u<-acap_tmueb_un_l+ay_tmueb_un_i-gaz_tmueb_un_i
acap_tmueb_un_d<-(acap_tmueb_un_l+acap_tmueb_un_u)/2
acap_tmueb_un_i<-abs(acap_tmueb_un_u-acap_tmueb_un_l)/2

gaz_tmuer_ad_l<-az_tmuer_ad_d%*%gamma_3_ad_d
gaz_tmuer_ad_u<-gaz_tmuer_ad_l+az_tmuer_ad_d%*%gamma_3_ad_i+az_tmuer_ad_i%*%gamma_3_ad_d+
  az_tmuer_ad_i%*%gamma_3_ad_i
gaz_tmuer_ad_d<-(gaz_tmuer_ad_l+gaz_tmuer_ad_u)/2
gaz_tmuer_ad_i<-abs(gaz_tmuer_ad_u-gaz_tmuer_ad_l)/2
gaz_tmuer_un_l<-az_tmuer_un_d%*%gamma_4_un_d
gaz_tmuer_un_u<-gaz_tmuer_un_l+az_tmuer_un_d%*%gamma_4_un_i+az_tmuer_un_i%*%gamma_4_un_d+
  az_tmuer_un_i%*%gamma_4_un_i
gaz_tmuer_un_d<-(gaz_tmuer_un_l+gaz_tmuer_un_u)/2
gaz_tmuer_un_i<-abs(gaz_tmuer_un_u-gaz_tmuer_un_l)/2
acap_tmuer_ad_l<-ay_tmuer_ad_d-gaz_tmuer_ad_d
acap_tmuer_ad_u<-acap_tmuer_ad_l+ay_tmuer_ad_i-gaz_tmuer_ad_i
acap_tmuer_ad_d<-(acap_tmuer_ad_l+acap_tmuer_ad_u)/2
acap_tmuer_ad_i<-abs(acap_tmuer_ad_u-acap_tmuer_ad_l)/2
acap_tmuer_un_l<-ay_tmuer_un_d-gaz_tmuer_un_d
acap_tmuer_un_u<-acap_tmuer_un_l+ay_tmuer_un_i-gaz_tmuer_un_i
acap_tmuer_un_d<-(acap_tmuer_un_l+acap_tmuer_un_u)/2
acap_tmuer_un_i<-abs(acap_tmuer_un_u-acap_tmuer_un_l)/2

gaz_muebr_ad_l<-az_muebr_ad_d%*%gamma_3_un_d
gaz_muebr_ad_u<-gaz_muebr_ad_l+az_muebr_un_d%*%gamma_3_un_i+az_muebr_un_i%*%gamma_3_un_d+
  az_muebr_un_i%*%gamma_3_un_i
gaz_muebr_ad_d<-(gaz_muebr_ad_l+gaz_muebr_ad_u)/2
gaz_muebr_ad_i<-abs(gaz_muebr_ad_u-gaz_muebr_ad_l)/2
gaz_muebr_un_l<-az_muebr_un_d%*%gamma_4_ad_d
gaz_muebr_un_u<-gaz_muebr_un_l+az_muebr_un_d%*%gamma_4_ad_i+az_muebr_un_i%*%gamma_4_ad_d+
  az_muebr_un_i%*%gamma_4_ad_i
gaz_muebr_un_d<-(gaz_muebr_un_l+gaz_muebr_un_u)/2
gaz_muebr_un_i<-abs(gaz_muebr_un_u-gaz_muebr_un_l)/2
acap_muebr_ad_l<-ay_muebr_ad_d-gaz_muebr_ad_d
acap_muebr_ad_u<-acap_muebr_ad_l+ay_muebr_ad_i-gaz_muebr_ad_i
acap_muebr_ad_d<-(acap_muebr_ad_l+acap_muebr_ad_u)/2
acap_muebr_ad_i<-abs(acap_muebr_ad_u-acap_muebr_ad_l)/2
acap_muebr_un_l<-ay_muebr_un_d-gaz_muebr_un_d
acap_muebr_un_u<-acap_muebr_un_l+ay_muebr_un_i-gaz_muebr_un_i
acap_muebr_un_d<-(acap_muebr_un_l+acap_muebr_un_u)/2
acap_muebr_un_i<-abs(acap_muebr_un_u-acap_muebr_un_l)/2

acapss_all_ad_l<-Qy_all_ad_d*acap_all_ad_d
acapss_all_ad_u<-acapss_all_ad_l+Qy_all_ad_d*acap_all_ad_i+
  Qy_all_ad_i*acap_all_ad_d+Qy_all_ad_i*acap_all_ad_i
acapss_all_ad_d<-(acapss_all_ad_u+acapss_all_ad_l)/2
acapss_all_ad_i<-abs(acapss_all_ad_u-acapss_all_ad_l)/2
acapss_all_un_l<-Qy_all_un_d*acap_all_un_d
acapss_all_un_u<-acapss_all_un_l+Qy_all_un_d*acap_all_un_i+
  Qy_all_un_i*acap_all_un_d+Qy_all_un_i*acap_all_un_i
acapss_all_un_d<-(acapss_all_un_u+acapss_all_un_l)/2
acapss_all_un_i<-abs(acapss_all_un_u-acapss_all_un_l)/2

acapss_tmueb_ad_l<-Qy_tmueb_ad_d*acap_tmueb_ad_d
acapss_tmueb_ad_u<-acapss_tmueb_ad_l+Qy_tmueb_ad_d*acap_tmueb_ad_i+
  Qy_tmueb_ad_i*acap_tmueb_ad_d+Qy_tmueb_ad_i*acap_tmueb_ad_i
acapss_tmueb_ad_d<-(acapss_tmueb_ad_u+acapss_tmueb_ad_l)/2
acapss_tmueb_ad_i<-abs(acapss_tmueb_ad_u-acapss_tmueb_ad_l)/2
acapss_tmueb_un_l<-Qy_tmueb_un_d*acap_tmueb_un_d
acapss_tmueb_un_u<-acapss_tmueb_un_l+Qy_tmueb_un_d*acap_tmueb_un_i+
  Qy_tmueb_un_i*acap_tmueb_un_d+Qy_tmueb_un_i*acap_tmueb_un_i
acapss_tmueb_un_d<-(acapss_tmueb_un_u+acapss_tmueb_un_l)/2
acapss_tmueb_un_i<-abs(acapss_tmueb_un_u-acapss_tmueb_un_l)/2

acapss_tmuer_ad_l<-Qy_tmuer_ad_d*acap_tmuer_ad_d
acapss_tmuer_ad_u<-acapss_tmuer_ad_l+Qy_tmuer_ad_d*acap_tmuer_ad_i+
  Qy_tmuer_ad_i*acap_tmuer_ad_d+Qy_tmuer_ad_i*acap_tmuer_ad_i
acapss_tmuer_ad_d<-(acapss_tmuer_ad_u+acapss_tmuer_ad_l)/2
acapss_tmuer_ad_i<-abs(acapss_tmuer_ad_u-acapss_tmuer_ad_l)/2
acapss_tmuer_un_l<-Qy_tmuer_un_d*acap_tmuer_un_d
acapss_tmuer_un_u<-acapss_tmuer_un_l+Qy_tmuer_un_d*acap_tmuer_un_i+
  Qy_tmuer_un_i*acap_tmuer_un_d+Qy_tmuer_un_i*acap_tmuer_un_i
acapss_tmuer_un_d<-(acapss_tmuer_un_u+acapss_tmuer_un_l)/2
acapss_tmuer_un_i<-abs(acapss_tmuer_un_u-acapss_tmuer_un_l)/2

acapss_muebr_ad_l<-Qy_muebr_ad_d*acap_muebr_ad_d
acapss_muebr_ad_u<-acapss_muebr_ad_l+Qy_muebr_ad_d*acap_muebr_ad_i+
  Qy_muebr_ad_i*acap_muebr_ad_d+Qy_muebr_ad_i*acap_muebr_ad_i
acapss_muebr_ad_d<-(acapss_muebr_ad_u+acapss_muebr_ad_l)/2
acapss_muebr_ad_i<-abs(acapss_muebr_ad_u-acapss_muebr_ad_l)/2
acapss_muebr_un_l<-Qy_muebr_un_d*acap_muebr_un_d
acapss_muebr_un_u<-acapss_muebr_un_l+Qy_muebr_un_d*acap_muebr_un_i+
  Qy_muebr_un_i*acap_muebr_un_d+Qy_muebr_un_i*acap_muebr_un_i
acapss_muebr_un_d<-(acapss_muebr_un_u+acapss_muebr_un_l)/2
acapss_muebr_un_i<-abs(acapss_muebr_un_u-acapss_muebr_un_l)/2

ss_all_ad_l<-sum(acapss_all_ad_l)
ss_all_ad_u<-sum(acapss_all_ad_u)
ss_all_ad_d<-(ss_all_ad_l+ss_all_ad_u)/2
ss_all_ad_i<-abs((ss_all_ad_u-ss_all_ad_l)/2)
ss_all_un_l<-sum(acapss_all_un_l)
ss_all_un_u<-sum(acapss_all_un_u)
ss_all_un_d<-sum(acapss_all_un_d)
ss_all_un_i<-abs(ss_all_un_l-ss_all_un_u)/2
ss_tmueb_ad_l<-sum(acapss_tmueb_ad_l)
ss_tmueb_ad_u<-sum(acapss_tmueb_ad_u)
ss_tmueb_ad_d<-(ss_tmueb_ad_l+ss_tmueb_ad_u)/2
ss_tmueb_ad_i<-abs(ss_tmueb_ad_l-ss_tmueb_ad_u)/2
ss_tmueb_un_l<-sum(acapss_tmueb_un_l)
ss_tmueb_un_u<-sum(acapss_tmueb_un_u)
ss_tmueb_un_d<-sum(acapss_tmueb_un_d)
ss_tmueb_un_i<-abs(ss_tmueb_un_l-ss_tmueb_un_u)/2

ss_tmuer_ad_l<-sum(acapss_tmuer_ad_l)
ss_tmuer_ad_u<-sum(acapss_tmuer_ad_u)
ss_tmuer_ad_d<-(ss_tmuer_ad_l+ss_tmuer_ad_u)/2
ss_tmuer_ad_i<-abs(ss_tmuer_ad_l-ss_tmuer_ad_u)/2
ss_tmuer_un_l<-sum(acapss_tmuer_un_l)
ss_tmuer_un_u<-sum(acapss_tmuer_un_u)
ss_tmuer_un_d<-sum(acapss_tmuer_un_d)
ss_tmuer_un_i<-abs(ss_tmuer_un_l-ss_tmuer_un_u)/2

ss_muebr_ad_l<-sum(acapss_muebr_ad_l)
ss_muebr_ad_u<-sum(acapss_muebr_ad_u)
ss_muebr_ad_d<-sum(acapss_muebr_ad_d)
ss_muebr_ad_i<-abs(ss_muebr_ad_l-ss_muebr_ad_u)/2
ss_muebr_un_l<-sum(acapss_muebr_un_l)
ss_muebr_un_u<-sum(acapss_muebr_un_u)
ss_muebr_un_d<-sum(acapss_muebr_un_d)
ss_muebr_un_i<-abs(ss_muebr_un_l-ss_muebr_un_u)/2
gzy_1_ad_l<-gamma_1_ad_d*zpy_d
gzy_1_ad_u<-gzy_1_ad_l+gamma_1_ad_d*zpy_i+gamma_1_ad_i*zpy_d+gamma_1_ad_i*zpy_i
gzy_1_ad_d<-(gzy_1_ad_l+gzy_1_ad_u)/2
gzy_1_ad_i<-abs(gzy_1_ad_u-gzy_1_ad_l)/2
gzy_1_un_l<-gamma_1_un_d*zpy_d
gzy_1_un_u<-gzy_1_un_l+gamma_1_un_d*zpy_i+gamma_1_un_i*zpy_d+gamma_1_un_i*zpy_i
gzy_1_un_d<-(gzy_1_un_l+gzy_1_un_u)/2
gzy_1_un_i<-abs(gzy_1_un_u-gzy_1_un_l)/2
gzy_2_ad_l<-gamma_2_ad_d*zpy_d
gzy_2_ad_u<-gzy_2_ad_l+gamma_2_ad_d*zpy_i+gamma_2_ad_i*zpy_d+gamma_2_ad_i*zpy_i
gzy_2_ad_d<-(gzy_2_ad_l+gzy_2_ad_u)/2
gzy_2_ad_i<-abs(gzy_2_ad_u-gzy_2_ad_l)/2
gzy_2_un_l<-gamma_2_un_d*zpy_d
gzy_2_un_u<-gzy_2_un_l+gamma_2_un_d*zpy_i+gamma_2_un_i*zpy_d+gamma_2_un_i*zpy_i
gzy_2_un_d<-(gzy_2_un_l+gzy_2_un_u)/2
gzy_2_un_i<-abs(gzy_2_un_u-gzy_2_un_l)/2
gzy_3_ad_l<-gamma_3_ad_d*zpy_d
gzy_3_ad_u<-gzy_3_ad_l+gamma_3_ad_d*zpy_i+gamma_3_ad_i*zpy_d+gamma_3_ad_i*zpy_i
gzy_3_ad_d<-(gzy_3_ad_l+gzy_3_ad_u)/2
gzy_3_ad_i<-abs(gzy_3_ad_u-gzy_3_ad_l)/2
gzy_3_un_l<-gamma_3_un_d*zpy_d
gzy_3_un_u<-gzy_3_un_l+gamma_3_un_d*zpy_i+gamma_3_un_i*zpy_d+gamma_3_un_i*zpy_i
gzy_3_un_d<-(gzy_3_un_l+gzy_3_un_u)/2
gzy_3_un_i<-abs(gzy_3_un_u-gzy_3_un_l)/2
gzy_4_ad_l<-gamma_4_ad_d*zpy_d
gzy_4_ad_u<-gzy_4_ad_l+gamma_4_ad_d*zpy_i+gamma_4_ad_i*zpy_d+gamma_4_ad_i*zpy_i
gzy_4_ad_d<-(gzy_4_ad_l+gzy_4_ad_u)/2
gzy_4_ad_i<-abs(gzy_4_ad_u-gzy_4_ad_l)/2
gzy_4_un_l<-gamma_4_un_d*zpy_d
gzy_4_un_u<-gzy_4_un_l+gamma_4_un_d*zpy_i+gamma_4_un_i*zpy_d+gamma_4_un_i*zpy_i
gzy_4_un_d<-(gzy_4_un_l+gzy_4_un_u)/2
gzy_4_un_i<-abs(gzy_4_un_u-gzy_4_un_l)/2


mdss_all_ad_l<-ss_all_ad_d+gzy_1_ad_d
mdss_all_ad_u<-mdss_all_ad_l+ss_all_ad_i+gzy_1_ad_i
mdss_all_ad_d<-(mdss_all_ad_u+mdss_all_ad_l)/2
mdss_all_ad_i<-abs(mdss_all_ad_u-mdss_all_ad_l)/2
mdss_all_un_l<-ss_all_un_d+gzy_1_un_d
mdss_all_un_u<-mdss_all_un_l+ss_all_un_i+gzy_1_un_i
mdss_all_un_d<-(mdss_all_un_u+mdss_all_un_l)/2
mdss_all_un_i<-abs(mdss_all_un_u-mdss_all_un_l)/2

mdss_tmueb_ad_l<-ss_tmueb_ad_d+gzy_2_ad_d
mdss_tmueb_ad_u<-mdss_tmueb_ad_l+ss_tmueb_ad_i+gzy_2_ad_i
mdss_tmueb_ad_d<-(mdss_tmueb_ad_u+mdss_tmueb_ad_l)/2
mdss_tmueb_ad_i<-abs(mdss_tmueb_ad_u-mdss_tmueb_ad_l)/2
mdss_tmueb_un_l<-ss_tmueb_un_d+gzy_2_un_d
mdss_tmueb_un_u<-mdss_tmueb_un_l+ss_tmueb_un_i+gzy_2_un_i
mdss_tmueb_un_d<-(mdss_tmueb_un_u+mdss_tmueb_un_l)/2
mdss_tmueb_un_i<-abs(mdss_tmueb_un_u-mdss_tmueb_un_l)/2

mdss_tmuer_ad_l<-ss_tmuer_ad_d+gzy_3_ad_d
mdss_tmuer_ad_u<-mdss_tmuer_ad_l+ss_tmuer_ad_i+gzy_3_ad_i
mdss_tmuer_ad_d<-(mdss_tmuer_ad_u+mdss_tmuer_ad_l)/2
mdss_tmuer_ad_i<-abs(mdss_tmuer_ad_u-mdss_tmuer_ad_l)/2
mdss_tmuer_un_l<-ss_tmuer_un_d+gzy_4_un_d
mdss_tmuer_un_u<-mdss_tmuer_un_l+ss_tmuer_un_i+gzy_4_un_i
mdss_tmuer_un_d<-(mdss_tmuer_un_u+mdss_tmuer_un_l)/2
mdss_tmuer_un_i<-abs(mdss_tmuer_un_u-mdss_tmuer_un_l)/2

mdss_muebr_ad_l<-ss_muebr_ad_d+gzy_3_un_d
mdss_muebr_ad_u<-mdss_muebr_ad_l+ss_muebr_ad_i+gzy_3_un_i
mdss_muebr_ad_d<-(mdss_muebr_ad_u+mdss_muebr_ad_l)/2
mdss_muebr_ad_i<-abs(mdss_muebr_ad_u-mdss_muebr_ad_l)/2
mdss_muebr_un_l<-ss_muebr_un_d+gzy_4_ad_d
mdss_muebr_un_u<-mdss_muebr_un_l+ss_muebr_un_i+gzy_4_ad_i
mdss_muebr_un_d<-(mdss_muebr_un_u+mdss_muebr_un_l)/2
mdss_muebr_un_i<-abs(mdss_muebr_un_u-mdss_muebr_un_l)/2
####error sum of sq
ess_all_ad_l<-ypy_d-mdss_all_ad_d
ess_all_ad_u<-ess_all_ad_l+ypy_i-mdss_all_ad_i
ess_all_ad_d<-(ess_all_ad_l+ess_all_ad_u)/2
ess_all_ad_i<-abs(ess_all_ad_u-ess_all_ad_l)/2
ess_all_un_l<-ypy_d-mdss_all_un_d
ess_all_un_u<-ess_all_un_l+ypy_i-mdss_all_un_i
ess_all_un_d<-(ess_all_un_l+ess_all_un_u)/2
ess_all_un_i<-abs(ess_all_un_u-ess_all_un_l)/2
esst_ad_l<-ypy_d-mdss_muebr_un_d
esst_ad_u<-esst_ad_l+ypy_i-mdss_muebr_un_i
esst_ad_d<-(esst_ad_u+esst_ad_l)/2
esst_ad_i<-abs(esst_ad_u-esst_ad_l)/2
esst_un_l<-ypy_d-mdss_muebr_ad_d
esst_un_u<-esst_un_l+ypy_i-mdss_muebr_ad_i
esst_un_d<-(esst_un_u+esst_un_l)/2
esst_un_i<-abs(esst_un_u-esst_un_l)/2
essb_ad_l<-ypy_d-mdss_tmuer_un_d
essb_ad_u<-essb_ad_l+ypy_i-mdss_tmuer_un_i
essb_ad_d<-(essb_ad_u+essb_ad_l)/2
essb_ad_i<-abs(essb_ad_u-essb_ad_l)/2
essb_un_l<-ypy_d-mdss_tmuer_ad_d
essb_un_u<-essb_un_l+ypy_i-mdss_tmuer_ad_i
essb_un_d<-(essb_un_u+essb_un_l)/2
essb_un_i<-abs(essb_un_u-essb_un_l)/2
essr_un_l<-ypy_d-mdss_tmueb_un_d
essr_un_u<-essr_un_l+ypy_i-mdss_tmueb_un_i
essr_un_d<-(essr_un_u+essr_un_l)/2
essr_un_i<-abs(essr_un_u-essr_un_l)/2
essr_ad_l<-ypy_d-mdss_tmueb_ad_d
essr_ad_u<-essr_ad_l+ypy_i-mdss_tmueb_ad_i
essr_ad_d<-(essr_ad_u+essr_ad_l)/2
essr_ad_i<-abs(essr_ad_u-essr_ad_l)/2
#####covariate adj trt ss adj to block
tss_ad_l<-esst_ad_d-ess_all_ad_d
tss_ad_u<-tss_ad_l+esst_ad_i-ess_all_ad_i
tss_ad_d<-(tss_ad_u+tss_ad_l)/2
tss_ad_i<-abs(tss_ad_u-tss_ad_l)/2
tss_un_l<-esst_un_d-ess_all_un_d
tss_un_u<-tss_un_l+esst_un_i-ess_all_un_i
tss_un_d<-(tss_un_u+tss_un_l)/2
tss_un_i<-abs(tss_un_u-tss_un_l)/2
bss_ad_l<-essb_ad_d-ess_all_un_d
bss_ad_u<-bss_ad_l+essb_ad_i-ess_all_un_i
bss_ad_d<-(bss_ad_u+bss_ad_l)/2
bss_ad_i<-abs(bss_ad_u-bss_ad_l)/2
bss_un_l<-essb_un_d-ess_all_ad_d
bss_un_u<-bss_un_l+essb_un_i-ess_all_ad_i
bss_un_d<-(bss_un_u+bss_un_l)/2
bss_un_i<-abs(bss_un_u-bss_un_l)/2
rss_un_l<-essr_un_d-ess_all_un_d
rss_un_u<-rss_un_l+essr_un_i-ess_all_un_i
rss_un_d<-(rss_un_u+rss_un_l)/2
rss_un_i<-abs(rss_un_u-rss_un_l)/2
rss_ad_l<-essr_ad_d-ess_all_ad_d
rss_ad_u<-rss_ad_l+essr_ad_i-ess_all_ad_i
rss_ad_d<-(rss_ad_u+rss_ad_l)/2
rss_ad_i<-abs(rss_ad_u-rss_ad_l)/2
###cov ss
gss_l<-gamma_1_ad_d*EXy_all_ad_d
gss_u<-gss_l+gamma_1_ad_d*EXy_all_ad_i+gamma_1_ad_i*EXy_all_ad_d+
  gamma_1_ad_i*EXy_all_ad_i
gss_d<-(gss_l+gss_u)/2
gss_i<-abs(gss_u-gss_l)/2
##totalss
toss_l<-t(dy)%*%(dy)
toss_u<-toss_l+2*t(dy)%*%iy+t(iy)%*%iy
toss_d<-(toss_l+toss_u)/2
toss_i<-abs(toss_u-toss_l)/2
# Example values (replace with your computed ones)
SV <- c("Covariate","Row","Treatment unadj", "Column adj", "Error",
         "Column unadj", "Treatment adj", "Error", "Total")

DF <- c(1,b-1, v-1, k-1, v*r-b-v-k+1, k-1,v-1, v*r-b-v-k+1, v*r-1)

SSN <- list(c(gss_l,gss_u),
            c(rss_ad_l,rss_ad_u),
            c(tss_un_l,tss_un_u),
            c(bss_ad_l,bss_ad_u),
            c(ess_all_un_l,ess_all_un_u),
            c(bss_un_l,bss_un_u),
            c(tss_ad_l,tss_ad_u),
            c(ess_all_ad_l,ess_all_ad_u),
            c(toss_l,toss_u))


MSSN <- list(
  c(gss_l/1,gss_u/1),
  c(rss_un_l/(b-1),rss_un_u/(b-1)),
  c(tss_un_l/(v-1),tss_un_u/(v-1)),
  c(bss_ad_l/(k-1),bss_ad_u/(k-1)),
  c(ess_all_un_l/((v*r-b-v-k+1)),ess_all_un_u/((v*r-b-v-k+1))),
  c(bss_un_l/(k-1),bss_un_u/(k-1)),
  c(tss_ad_l/(v-1),tss_ad_u/(v-1)),
  c(ess_all_ad_l/((v*r-b-v-k+1)),ess_all_ad_u/((v*r-b-v-k+1))),
  c(NA, NA)
)



F_C_L<-(gss_d/1)/(ess_all_ad_d/(v*r-b-v-k+1))
F_C_U<-F_C_L+((gss_d/1)*(ess_all_ad_i/(v*r-b-v-k+1))+
                (gss_i/1)*(ess_all_ad_d/(v*r-b-v-k+1))+
                (gss_i/1)*(ess_all_ad_i/(v*r-b-v-k+1)))/((ess_all_ad_d/(v*r-b-v-k+1))*(ess_all_ad_d/(v*r-b-v-k+1)))
F_r_un_L<-(rss_un_d/(b-1))/(ess_all_un_d/(v*r-b-v-k+1))
F_r_un_U<-F_r_un_L+(((rss_un_d/(b-1))*(ess_all_un_i/(v*r-b-v-k+1)))+
                      ((rss_un_i/(b-1))*(ess_all_un_d/(v*r-b-v-k+1)))+
                      ((rss_un_i/(b-1))*(ess_all_un_i/(v*r-b-v-k+1))))/((ess_all_un_d/(v*r-b-v-k+1))*(ess_all_un_d/(v*r-b-v-k+1)))
F_t_un_L<-(tss_un_d/(v-1))/(ess_all_un_d/(v*r-b-v-k+1))
F_t_un_U<-F_t_un_L+(((tss_un_d/(v-1))*(ess_all_un_i/(v*r-b-v-k+1)))+
                      ((tss_un_i/(v-1))*(ess_all_un_d/(v*r-b-v-k+1)))+
                      ((tss_un_i/(v-1))*(ess_all_un_i/(v*r-b-v-k+1))))/((ess_all_un_d/(v*r-b-v-k+1))*(ess_all_un_d/(v*r-b-v-k+1)))
F_b_ad_L<-(bss_ad_d/(k-1))/(ess_all_un_d/(v*r-b-v-k+1))
F_b_ad_U<-F_b_ad_L+(((bss_ad_d/(k-1))*(ess_all_un_i/(v*r-b-v-k+1)))+
                      ((bss_ad_i/(k-1))*(ess_all_un_d/(v*r-b-v-k+1)))+
                      ((bss_ad_i/(k-1))*(ess_all_un_i/(v*r-b-v-k+1))))/((ess_all_un_d/(v*r-b-v-k+1))*(ess_all_un_d/(v*r-b-v-k+1)))

F_t_ad_L<-(tss_ad_d/(v-1))/(ess_all_ad_d/(v*r-b-v-k+1))
F_t_ad_U<-F_t_ad_L+(((tss_ad_d/(v-1))*(ess_all_ad_i/(v*r-b-v-k+1)))+
                      ((tss_ad_i/(b-1))*(ess_all_ad_d/(v*r-b-v-k+1)))+
                      ((tss_ad_i/(b-1))*(ess_all_ad_i/(v*r-b-v-k+1))))/((ess_all_ad_d/(v*r-b-v-k+1))*(ess_all_ad_d/(v*r-b-v-k+1)))
F_b_un_L<-(bss_un_d/(k-1))/(ess_all_ad_d/(v*r-b-v-k+1))
F_b_un_U<-F_b_un_L+(((bss_un_d/(k-1))*(ess_all_ad_i/(v*r-b-v-k+1)))+
                      ((bss_un_i/(k-1))*(ess_all_ad_d/(v*r-b-v-k+1)))+
                      ((bss_un_i/(k-1))*(ess_all_ad_i/(v*r-b-v-k+1))))/((ess_all_ad_d/(v*r-b-v-k+1))*(ess_all_ad_d/(v*r-b-v-k+1)))

FN <- list(
  c(F_C_L,F_C_U),        # Covariate
  c(F_r_un_L,F_r_un_U),  #replication
  c(F_t_un_L,F_t_un_U),  # Treatment unadj
  c(F_b_ad_L,F_b_ad_U),  # block adj
  c(NA, NA),             # Error
  c(F_b_un_L,F_b_un_U),  # Block unadj
  c(F_t_ad_L,F_t_ad_U),  # Treatment adj
  c(NA, NA),             # Error adj
  c(NA, NA)              # Total
)
length(FN)
length(MSSN)


format_interval <- function(x) {
  if (any(is.na(x))) return("")
  paste0("[", round(x[1],2), ", ", round(x[2],2), "]")
}
SSN_fmt  <- sapply(SSN, format_interval)
MSSN_fmt <- sapply(MSSN, format_interval)
FN_fmt   <- sapply(FN, format_interval)
nancova_table <- data.frame(
  SV = SV,
  DF = DF,
  SSN = SSN_fmt,
  MSSN = MSSN_fmt,
  FN = FN_fmt,
  stringsAsFactors = FALSE
)
F_critical_5 <- qf(1-alpha,DF[3],DF[7])
F_critical_1 <- qf(1-alpha/5,DF[3],DF[7])
F_critical_001 <- qf(1-alpha/50,DF[3],DF[7])

significance <- c(

  ifelse(F_C_L > F_critical_001 & F_C_U > F_critical_001,"***",
         ifelse(F_C_L > F_critical_1 & F_C_U > F_critical_1,"**",
                ifelse(F_C_L > F_critical_5 & F_C_U > F_critical_5,"*",
                       ifelse(F_C_L < F_critical_5 & F_C_U < F_critical_5,"NS","ID")))),

  ifelse(F_r_un_L > F_critical_001 & F_r_un_U > F_critical_001,"***",
         ifelse(F_r_un_L > F_critical_1 & F_r_un_U > F_critical_1,"**",
                ifelse(F_r_un_L > F_critical_5 & F_r_un_U > F_critical_5,"*",
                       ifelse(F_r_un_L < F_critical_5 & F_r_un_U < F_critical_5,"NS","ID")))),

  ifelse(F_t_un_L > F_critical_001 & F_t_un_U > F_critical_001,"***",
         ifelse(F_t_un_L > F_critical_1 & F_t_un_U > F_critical_1,"**",
                ifelse(F_t_un_L > F_critical_5 & F_t_un_U > F_critical_5,"*",
                       ifelse(F_t_un_L < F_critical_5 & F_t_un_U < F_critical_5,"NS","ID")))),

  ifelse(F_b_ad_L > F_critical_001 & F_b_ad_U > F_critical_001,"***",
         ifelse(F_b_ad_L > F_critical_1 & F_b_ad_U > F_critical_1,"**",
                ifelse(F_b_ad_L > F_critical_5 & F_b_ad_U > F_critical_5,"*",
                       ifelse(F_b_ad_L < F_critical_5 & F_b_ad_U < F_critical_5,"NS","ID")))),

  "",

  ifelse(F_b_un_L > F_critical_001 & F_b_un_U > F_critical_001,"***",
         ifelse(F_b_un_L > F_critical_1 & F_b_un_U > F_critical_1,"**",
                ifelse(F_b_un_L > F_critical_5 & F_b_un_U > F_critical_5,"*",
                       ifelse(F_b_un_L < F_critical_5 & F_b_un_U < F_critical_5,"NS","ID")))),

  ifelse(F_t_ad_L > F_critical_001 & F_t_ad_U > F_critical_001,"***",
         ifelse(F_t_ad_L > F_critical_1 & F_t_ad_U > F_critical_1,"**",
                ifelse(F_t_ad_L > F_critical_5 & F_t_ad_U > F_critical_5,"*",
                       ifelse(F_t_ad_L < F_critical_5 & F_t_ad_U < F_critical_5,"NS","ID")))),

  "",
  ""
)

nancova_table <- data.frame(
  SV = SV,
  DF = DF,
  SSN = SSN_fmt,
  MSSN = MSSN_fmt,
  FN = FN_fmt,
  Significance = significance,
  stringsAsFactors = FALSE
)

if(verbose){
message("\nNeutrosophic Analysis of Covariance Table\n\n")
print(nancova_table,row.names=FALSE)

####LSD Test

mse_l <- ess_all_ad_l/(v*r-b-v-k+1)
mse_u <- ess_all_ad_u/(v*r-b-v-k+1)

if(significance[6] != "NS"){

  LSD_l <- qt(1-alpha/2,(v*r-b-v-k+1))*sqrt((2*mse_l)/r)

  LSD_u <- LSD_l+
    qt(1-alpha/2,(v*r-b-v-k+1))*sqrt((2*mse_u)/r)-
    qt(1-alpha/2,(v*r-b-v-k+1))*sqrt((2*mse_l)/r)

  message("\nTreatment effect is significant. Hence multiple comparison using LSD is performed.\n")

  message("\nLSD Interval : [",
      round(LSD_l,4),", ",
      round(LSD_u,4),"]\n")

  trt_means_l <- tapply(as.vector(Lower_y),as.vector(design),mean)

  trt_means_u <- tapply(as.vector(Upper_y),as.vector(design),mean)

  trt_means_d <- (trt_means_l+trt_means_u)/2

  trt_means_i <- abs(trt_means_u-trt_means_l)/2

  comparison <- data.frame(
    Treatment1 = character(),
    Treatment2 = character(),
    Difference = character(),
    LSD = character(),
    Result = character(),
    stringsAsFactors = FALSE
  )

  for(i in 1:(v-1)){

    for(j in (i+1):v){

      diff_l <- abs(trt_means_d[i]-trt_means_d[j])

      diff_u <- diff_l+trt_means_i[i]+trt_means_i[j]

      decision <- ifelse(diff_l > LSD_u & diff_u > LSD_u,"S",
                         ifelse(diff_l < LSD_l & diff_u < LSD_l,"NS","ID"))

      comparison <- rbind(
        comparison,
        data.frame(
          Treatment1 = i,
          Treatment2 = j,
          Difference = paste0("[",round(diff_l,2),", ",round(diff_u,2),"]"),
          LSD = paste0("[",round(LSD_l,2),", ",round(LSD_u,2),"]"),
          Result = decision,
          stringsAsFactors = FALSE
        )
      )

    }
  }

  message("\nTreatment Comparisons Using LSD\n")

  print(comparison,row.names=FALSE)

}else{

  message("\nTreatment effect is non significant. Hence multiple comparison is not performed.\n")

}

message("\nSignificance Codes:")
message("*** : Significant at p < 0.001")
message("**  : Significant at p < 0.01")
message("*   : Significant at p < 0.05")
message("NS  : Non Significant")
message("ID  : Indeterminate")
}
result <- list(
  nancova_table = nancova_table,
  comparison = if(exists("comparison")) comparison else NULL,
  LSD = if(exists("LSD_l"))
    c(Lower = LSD_l, Upper = LSD_u) else NULL
)

invisible(result)

}


