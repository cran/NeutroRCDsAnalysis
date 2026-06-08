#' Neutrosophic Analysis of Covariance for Latin Square Design
#'
#' Performs Neutrosophic Analysis of Covariance (NANCOVA) for Latin
#' square designs using interval-valued response and covariate
#' observations. The function computes neutrosophic sums of squares,
#' mean squares, interval-valued F-statistics, significance tests,
#' and Least Significant Difference (LSD)-based treatment
#' comparisons. For crisp data, enter identical lower and upper
#' values to obtain the corresponding classical ANCOVA results.
#'
#' @usage
#' LSDnsANCOVA(Lower_y, Upper_y, Lower_z, Upper_z, design, alpha = 0.05, verbose = FALSE)
#'
#' @param Lower_y Matrix containing lower bounds of response observations.
#' @param Upper_y Matrix containing upper bounds of response observations.
#' @param Lower_z Matrix containing lower bounds of covariate observations.
#' @param Upper_z Matrix containing upper bounds of covariate observations.
#' @param design Matrix representing Latin square treatment allocation.
#' @param alpha Significance level for the F-test and LSD test.Default is 0.05.
#' @param verbose Logical. If TRUE, displays the ANCOVA table,
#' LSD interval, treatment comparisons and significance codes.
#' Default is FALSE.
#'
#' @details
#' Input matrix structure:
#' \itemize{
#'   \item Rows represent blocks (or rows of the design).
#'   \item Columns represent treatment positions within each block.
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
#' \item \code{comparison}: LSD treatment comparisons, if applicable.
#' \item \code{LSD}: Lower and upper limits of the LSD interval.
#' }
#'
#' @examples
#' Lower_y <- matrix(c(
#' 18.86,15.49,18.62,16.37,12.04,
#' 18.23,26.80,26.70,15.89,16.96,
#' 26.26,12.16,27.15,18.30,22.93,
#' 25.31,26.71,18.36,9.98,28.34,
#' 29.78,21.20,21.99,18.97,18.71
#' ), nrow = 5, byrow = TRUE)
#'
#' Upper_y <- matrix(c(
#' 21.14,22.51,25.38,21.63,19.96,
#' 27.77,29.20,33.30,18.11,25.04,
#' 29.74,19.84,36.85,25.70,31.07,
#' 34.69,29.29,21.64,18.02,33.66,
#' 32.22,28.80,24.01,25.03,21.29
#' ), nrow = 5, byrow = TRUE)
#'
#' Lower_z <- matrix(c(
#' 9.11,7.39,7.09,4.10,5.92,
#' 6.91,4.26,13.10,5.86,12.54,
#' 8.89,4.98,4.32,8.80,9.77,
#' 6.40,7.28,3.29,3.81,4.45,
#' 6.42,2.51,6.06,3.25,5.73
#' ), nrow = 5, byrow = TRUE)
#'
#' Upper_z <- matrix(c(
#' 16.89,10.61,14.91,9.90,8.08,
#' 11.09,9.74,16.90,12.14,17.46,
#' 13.11,11.02,11.68,11.20,18.23,
#' 9.60,12.72,10.71,8.19,7.55,
#' 15.58,7.49,9.94,12.75,14.27
#' ), nrow = 5, byrow = TRUE)
#'
#' design <- matrix(c(
#' 2,1,5,4,3,
#' 1,4,3,5,2,
#' 3,2,4,1,5,
#' 5,3,1,2,4,
#' 4,5,2,3,1
#' ), nrow = 5, byrow = TRUE)
#'
#' LSDnsANCOVA(Lower_y, Upper_y, Lower_z, Upper_z, design, alpha = 0.05, verbose = TRUE)
#'
#' @importFrom MASS ginv
#' @importFrom stats qf qt
#' @export


LSDnsANCOVA <- function(Lower_y, Upper_y, Lower_z, Upper_z, design, alpha = 0.05,  verbose = FALSE){

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

  v <- max(design)
  b <- nrow(design)
  k <- ncol(design)

  rep_counts <- table(design)

  r <- as.numeric(unique(rep_counts))

  n=b*k

  ####X matrix

  trt_vec <- as.vector(t(design))

  Tmat<- matrix(0, n, v)

  for(i in 1:n){
    Tmat[i, trt_vec[i]] <- 1
  }

  mu <- matrix(1, n, 1)

  Bmat <- diag(v)[rep(1:v, times = v), ]

  rep_vec <- rep(1:r, each = b/r * k)

  Rmat<- matrix(0, n, r)

  for(i in 1:n){
    Rmat[i, rep_vec[i]] <- 1
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

  zpz_l<-t(dz)%*%dz

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

  ##Treatment unadjusted

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

  ay_t_ad_l <- ginv(Ct_adj) %*% Qy_t_ad_d

  ay_t_ad_u <- ay_t_ad_l+Qy_t_ad_i/v

  ay_t_ad_d<-(ay_t_ad_l +ay_t_ad_u)/2
  ay_t_ad_i<-abs(ay_t_ad_l -ay_t_ad_u)/2

  ##Column unadjusted

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

  ay_b_ad_l <- ginv(Cb_adj) %*% Qy_b_ad_d

  ay_b_ad_u <- ay_b_ad_l+Qy_b_ad_i/v

  ay_b_ad_d<-(ay_b_ad_l +ay_b_ad_u)/2
  ay_b_ad_i<-abs(ay_b_ad_l-ay_b_ad_u)/2

  ##Row unadjusted

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

  ay_r_u <- ay_r_l+Qy_r_un_i/v

  ay_r_d<-(ay_r_l +ay_r_u)/2
  ay_r_i<-abs(ay_r_l-ay_r_u)/2

  ##Mean

  Qy_mu_l<-sum(round(dy,6))

  Qy_mu_u<-Qy_mu_l+sum(round(iy,7))

  Qy_mu_d<-(Qy_mu_l+Qy_mu_u)/2
  Qy_mu_i<-abs(Qy_mu_l-Qy_mu_u)/2

  ay_mu_l<-Qy_mu_d*1/(b*k)

  ay_mu_u<-ay_mu_l+Qy_mu_d/(b*k)

  ay_mu_d<-(ay_mu_l+ay_mu_u)/2
  ay_mu_i<-abs(ay_mu_l-ay_mu_u)/2

  Qy_all_ad_d<-c(Qy_t_ad_d,Qy_mu_d,Qy_b_ad_d,Qy_r_un_d)
  Qy_all_ad_i<-c(Qy_t_ad_i,Qy_mu_i,Qy_b_ad_i,Qy_r_un_i)

  Qy_tmueb_ad_d<-c(Qy_t_ad_d,Qy_mu_d,Qy_b_ad_d)
  Qy_tmueb_ad_i<-c(Qy_t_ad_i,Qy_mu_i,Qy_b_ad_i)

  Qy_tmuer_ad_d<-c(Qy_t_ad_d,Qy_mu_d,Qy_r_un_d)
  Qy_tmuer_ad_i<-c(Qy_t_ad_i,Qy_mu_i,Qy_r_un_i)

  Qy_muebr_ad_d<-c(Qy_mu_d,Qy_b_ad_d,Qy_r_un_d)
  Qy_muebr_ad_i<-c(Qy_mu_i,Qy_b_ad_i,Qy_r_un_i)

  ay_all_ad_d<-c(ay_t_ad_d,ay_mu_d,ay_b_ad_d,ay_r_d)
  ay_all_ad_i<-c(ay_t_ad_i,ay_mu_i,ay_b_ad_i,ay_r_i)

  ay_tmueb_ad_d<-c(ay_t_ad_d,ay_mu_d,ay_b_ad_d)
  ay_tmueb_ad_i<-c(ay_t_ad_i,ay_mu_i,ay_b_ad_i)

  ay_tmuer_ad_d<-c(ay_t_ad_d,ay_mu_d,ay_r_d)
  ay_tmuer_ad_i<-c(ay_t_ad_i,ay_mu_i,ay_r_i)

  ay_muebr_ad_d<-c(ay_mu_d,ay_b_ad_d,ay_r_d)
  ay_muebr_ad_i<-c(ay_mu_i,ay_b_ad_i,ay_r_i)

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

  az_t_ad_l <- ginv(Ct_adj) %*% Qz_t_ad_d

  az_t_ad_u <- az_t_ad_l+Qz_t_ad_i/v

  az_t_ad_d<-(az_t_ad_l +az_t_ad_u)/2
  az_t_ad_i<-abs(az_t_ad_l -az_t_ad_u)/2

  ##col

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

  az_b_ad_u <- az_b_ad_l+Qz_b_ad_i/v

  az_b_ad_d<-(az_b_ad_l +az_b_ad_u)/2
  az_b_ad_i<-abs(az_b_ad_l-az_b_ad_u)/2

  ###row

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

  az_r_u <- az_r_l+Qz_r_un_i/v

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

  Qz_all_ad_d<-c(Qz_t_ad_d,Qz_mu_d,Qz_b_ad_d,Qz_r_un_d)
  Qz_all_ad_i<-c(Qz_t_ad_i,Qz_mu_i,Qz_b_ad_i,Qz_r_un_i)

  Qz_tmueb_ad_d<-c(Qz_t_ad_d,Qz_mu_d,Qz_b_ad_d)
  Qz_tmueb_ad_i<-c(Qz_t_ad_i,Qz_mu_i,Qz_b_ad_i)

  Qz_tmuer_ad_d<-c(Qz_t_ad_d,Qz_mu_d,Qz_r_un_d)
  Qz_tmuer_ad_i<-c(Qz_t_ad_i,Qz_mu_i,Qz_r_un_i)

  Qz_muebr_ad_d<-c(Qz_mu_d,Qz_b_ad_d,Qz_r_un_d)
  Qz_muebr_ad_i<-c(Qz_mu_i,Qz_b_ad_i,Qz_r_un_i)

  az_all_ad_d<-c(az_t_ad_d,az_mu_d,az_b_ad_d,az_r_d)
  az_all_ad_i<-c(az_t_ad_i,az_mu_i,az_b_ad_i,az_r_i)

  az_tmueb_ad_d<-c(az_t_ad_d,az_mu_d,az_b_ad_d)
  az_tmueb_ad_i<-c(az_t_ad_i,az_mu_i,az_b_ad_i)

  az_tmuer_ad_d<-c(az_t_ad_d,az_mu_d,az_r_d)
  az_tmuer_ad_i<-c(az_t_ad_i,az_mu_i,az_r_i)

  az_muebr_ad_d<-c(az_mu_d,az_b_ad_d,az_r_d)
  az_muebr_ad_i<-c(az_mu_i,az_b_ad_i,az_r_i)

  ##qz*az

  QZAZ_all_ad_l<-c(Qz_t_ad_d*az_t_ad_d,Qz_mu_d*az_mu_d,Qz_b_ad_d*az_b_ad_d,Qz_r_un_d*az_r_d)

  QZAZ_all_ad_u<-QZAZ_all_ad_l+
    c(Qz_t_ad_d*az_t_ad_i,Qz_mu_d*az_mu_i,Qz_b_ad_d*az_b_ad_i,Qz_r_un_d*az_r_i)+
    c(Qz_t_ad_i*az_t_ad_d,Qz_mu_i*az_mu_d,Qz_b_ad_i*az_b_ad_d,Qz_r_un_i*az_r_d)+
    c(Qz_t_ad_i*az_t_ad_i,Qz_mu_i*az_mu_i,Qz_b_ad_i*az_b_ad_i,Qz_r_un_i*az_r_i)

  QZAZ_all_ad_d<-sum((QZAZ_all_ad_l+QZAZ_all_ad_u)/2)
  QZAZ_all_ad_i<-sum(abs(QZAZ_all_ad_u-QZAZ_all_ad_l)/2)

  QZAZ_tmueb_ad_l<-c(Qz_t_ad_d*az_t_ad_d,Qz_mu_d*az_mu_d,Qz_b_ad_d*az_b_ad_d)

  QZAZ_tmueb_ad_u<-QZAZ_tmueb_ad_l+
    c(Qz_t_ad_d*az_t_ad_i,Qz_mu_d*az_mu_i,Qz_b_ad_d*az_b_ad_i)+
    c(Qz_t_ad_i*az_t_ad_d,Qz_mu_i*az_mu_d,Qz_b_ad_i*az_b_ad_d)+
    c(Qz_t_ad_i*az_t_ad_i,Qz_mu_i*az_mu_i,Qz_b_ad_i*az_b_ad_i)

  QZAZ_tmueb_ad_d<-sum((QZAZ_tmueb_ad_l+QZAZ_tmueb_ad_u)/2)
  QZAZ_tmueb_ad_i<-sum(abs((QZAZ_tmueb_ad_u-QZAZ_tmueb_ad_l)/2))

  QZAZ_tmuer_ad_l<-c(Qz_t_ad_d*az_t_ad_d,Qz_mu_d*az_mu_d,Qz_r_un_d*az_r_d)

  QZAZ_tmuer_ad_u<-QZAZ_tmuer_ad_l+
    c(Qz_t_ad_d*az_t_ad_i,Qz_mu_d*az_mu_i,Qz_r_un_d*az_r_i)+
    c(Qz_t_ad_i*az_t_ad_d,Qz_mu_i*az_mu_d,Qz_r_un_i*az_r_d)+
    c(Qz_t_ad_i*az_t_ad_i,Qz_mu_i*az_mu_i,Qz_r_un_i*az_r_i)

  QZAZ_tmuer_ad_d<-sum((QZAZ_tmuer_ad_l+QZAZ_tmuer_ad_u)/2)
  QZAZ_tmuer_ad_i<-sum(abs((QZAZ_tmuer_ad_u-QZAZ_tmuer_ad_l)/2))

  QZAZ_muebr_ad_l<-c(Qz_mu_d*az_mu_d,Qz_b_ad_d*az_b_ad_d,Qz_r_un_d*az_r_d)

  QZAZ_muebr_ad_u<-QZAZ_muebr_ad_l+
    c(Qz_mu_d*az_mu_i,Qz_b_ad_d*az_b_ad_i,Qz_r_un_d*az_r_i)+
    c(Qz_mu_i*az_mu_d,Qz_b_ad_i*az_b_ad_d,Qz_r_un_i*az_r_d)+
    c(Qz_mu_i*az_mu_i,Qz_b_ad_i*az_b_ad_i,Qz_r_un_i*az_r_i)

  QZAZ_muebr_ad_d<-sum((QZAZ_muebr_ad_l+QZAZ_muebr_ad_u)/2)
  QZAZ_muebr_ad_i<-sum(abs(QZAZ_muebr_ad_u-QZAZ_muebr_ad_l)/2)

  ##qz*ay

  QZAy_all_ad_l<-c(Qz_t_ad_d*ay_t_ad_d,Qz_mu_d*ay_mu_d,Qz_b_ad_d*ay_b_ad_d,Qz_r_un_d*ay_r_d)

  QZAy_all_ad_u<-QZAy_all_ad_l+
    c(Qz_t_ad_d*ay_t_ad_i,Qz_mu_d*ay_mu_i,Qz_b_ad_d*ay_b_ad_i,Qz_r_un_d*ay_r_i)+
    c(Qz_t_ad_i*ay_t_ad_d,Qz_mu_i*ay_mu_d,Qz_b_ad_i*ay_b_ad_d,Qz_r_un_i*ay_r_d)+
    c(Qz_t_ad_i*ay_t_ad_i,Qz_mu_i*ay_mu_i,Qz_b_ad_i*ay_b_ad_i,Qz_r_un_i*ay_r_i)

  QZAy_all_ad_d<-sum((QZAy_all_ad_l+QZAy_all_ad_u)/2)
  QZAy_all_ad_i<-sum(abs(QZAy_all_ad_u-QZAy_all_ad_l)/2)

  QZAy_tmueb_ad_l<-c(Qz_t_ad_d*ay_t_ad_d,Qz_mu_d*ay_mu_d,Qz_b_ad_d*ay_b_ad_d)

  QZAy_tmueb_ad_u<-QZAy_tmueb_ad_l+
    c(Qz_t_ad_d*ay_t_ad_i,Qz_mu_d*ay_mu_i,Qz_b_ad_d*ay_b_ad_i)+
    c(Qz_t_ad_i*ay_t_ad_d,Qz_mu_i*ay_mu_d,Qz_b_ad_i*ay_b_ad_d)+
    c(Qz_t_ad_i*ay_t_ad_i,Qz_mu_i*ay_mu_i,Qz_b_ad_i*ay_b_ad_i)

  QZAy_tmueb_ad_d<-sum((QZAy_tmueb_ad_l+QZAy_tmueb_ad_u)/2)
  QZAy_tmueb_ad_i<-sum(abs((QZAy_tmueb_ad_u-QZAy_tmueb_ad_l)/2))

  QZAy_tmuer_ad_l<-c(Qz_t_ad_d*ay_t_ad_d,Qz_mu_d*ay_mu_d,Qz_r_un_d*ay_r_d)

  QZAy_tmuer_ad_u<-QZAy_tmuer_ad_l+
    c(Qz_t_ad_d*ay_t_ad_i,Qz_mu_d*ay_mu_i,Qz_r_un_d*ay_r_i)+
    c(Qz_t_ad_i*ay_t_ad_d,Qz_mu_i*ay_mu_d,Qz_r_un_i*ay_r_d)+
    c(Qz_t_ad_i*ay_t_ad_i,Qz_mu_i*ay_mu_i,Qz_r_un_i*ay_r_i)

  QZAy_tmuer_ad_d<-sum((QZAy_tmuer_ad_l+QZAy_tmuer_ad_u)/2)
  QZAy_tmuer_ad_i<-sum(abs((QZAy_tmuer_ad_u-QZAy_tmuer_ad_l)/2))

  QZAy_muebr_ad_l<-c(Qz_mu_d*ay_mu_d,Qz_b_ad_d*ay_b_ad_d,Qz_r_un_d*ay_r_d)

  QZAy_muebr_ad_u<-QZAy_muebr_ad_l+
    c(Qz_mu_d*ay_mu_i,Qz_b_ad_d*ay_b_ad_i,Qz_r_un_d*ay_r_i)+
    c(Qz_mu_i*ay_mu_d,Qz_b_ad_i*ay_b_ad_d,Qz_r_un_i*ay_r_d)+
    c(Qz_mu_i*ay_mu_i,Qz_b_ad_i*ay_b_ad_i,Qz_r_un_i*ay_r_i)

  QZAy_muebr_ad_d<-sum((QZAy_muebr_ad_l+QZAy_muebr_ad_u)/2)
  QZAy_muebr_ad_i<-sum(abs(QZAy_muebr_ad_u-QZAy_muebr_ad_l)/2)

  ##GAMMA

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

  EXX_muebr_ad_l<-zpz_d-QZAZ_muebr_ad_d
  EXX_muebr_ad_u<-EXX_muebr_ad_l+zpz_i-QZAZ_muebr_ad_i

  EXX_muebr_ad_d<-(EXX_muebr_ad_l+EXX_muebr_ad_u)/2
  EXX_muebr_ad_i<-abs(EXX_muebr_ad_u-EXX_muebr_ad_l)/2



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

  EXy_muebr_ad_l<-zpy_d-QZAy_muebr_ad_d
  EXy_muebr_ad_u<-EXy_muebr_ad_l+zpy_i-QZAy_muebr_ad_i
  EXy_muebr_ad_d<-(EXy_muebr_ad_l+EXy_muebr_ad_u)/2
  EXy_muebr_ad_i<-abs(EXy_muebr_ad_u-EXy_muebr_ad_l)/2

  gamma_1_ad_l<-EXy_all_ad_d/EXX_all_ad_d
  gamma_1_ad_u<-gamma_1_ad_l+(EXy_all_ad_d*EXX_all_ad_i+
                                EXy_all_ad_i*EXX_all_ad_d+
                                EXy_all_ad_i*EXX_all_ad_i)/
    (EXX_all_ad_d*EXX_all_ad_d)
  gamma_1_ad_d<-(gamma_1_ad_l+gamma_1_ad_u)/2
  gamma_1_ad_i<-abs(gamma_1_ad_u-gamma_1_ad_l)/2

  gamma_2_ad_l<-EXy_tmueb_ad_d/EXX_tmueb_ad_d
  gamma_2_ad_u<-gamma_2_ad_l+(EXy_tmueb_ad_d*EXX_tmueb_ad_i+
                                EXy_tmueb_ad_i*EXX_tmueb_ad_d+
                                EXy_tmueb_ad_i*EXX_tmueb_ad_i)/
    (EXX_tmueb_ad_d*EXX_tmueb_ad_d)
  gamma_2_ad_d<-(gamma_2_ad_l+gamma_2_ad_u)/2
  gamma_2_ad_i<-abs(gamma_2_ad_u-gamma_2_ad_l)/2

  gamma_3_ad_l<-EXy_tmuer_ad_d/EXX_tmuer_ad_d
  gamma_3_ad_u<-gamma_3_ad_l+(EXy_tmuer_ad_d*EXX_tmuer_ad_i+
                                EXy_tmuer_ad_i*EXX_tmuer_ad_d+
                                EXy_tmuer_ad_i*EXX_tmuer_ad_i)/
    (EXX_tmuer_ad_d*EXX_tmuer_ad_d)
  gamma_3_ad_d<-(gamma_3_ad_l+gamma_3_ad_u)/2
  gamma_3_ad_i<-abs(gamma_3_ad_u-gamma_3_ad_l)/2

  gamma_4_ad_l<-EXy_muebr_ad_d/EXX_muebr_ad_d
  gamma_4_ad_u<-gamma_4_ad_l+(EXy_muebr_ad_d*EXX_muebr_ad_i+
                                EXy_muebr_ad_i*EXX_muebr_ad_d+
                                EXy_muebr_ad_i*EXX_muebr_ad_i)/
    (EXX_muebr_ad_d*EXX_muebr_ad_d)
  gamma_4_ad_d<-(gamma_4_ad_l+gamma_4_ad_u)/2
  gamma_4_ad_i<-abs(gamma_4_ad_u-gamma_4_ad_l)/2
  ############gamma az

  gaz_all_ad_l<-az_all_ad_d%*%gamma_1_ad_d

  gaz_all_ad_u<-gaz_all_ad_l+az_all_ad_d%*%gamma_1_ad_i+
    az_all_ad_i%*%gamma_1_ad_d+az_all_ad_i%*%gamma_1_ad_i

  gaz_all_ad_d<-(gaz_all_ad_l+gaz_all_ad_u)/2
  gaz_all_ad_i<-abs(gaz_all_ad_u-gaz_all_ad_l)/2

  acap_all_ad_l<-ay_all_ad_d-gaz_all_ad_d

  acap_all_ad_u<-acap_all_ad_l+ay_all_ad_i-gaz_all_ad_i

  acap_all_ad_d<-(acap_all_ad_l+acap_all_ad_u)/2
  acap_all_ad_i<-abs(acap_all_ad_u-acap_all_ad_l)/2

  gaz_tmueb_ad_l<-az_tmueb_ad_d%*%gamma_2_ad_d

  gaz_tmueb_ad_u<-gaz_tmueb_ad_l+az_tmueb_ad_d%*%gamma_2_ad_i+
    az_tmueb_ad_i%*%gamma_2_ad_d+az_tmueb_ad_i%*%gamma_2_ad_i

  gaz_tmueb_ad_d<-(gaz_tmueb_ad_l+gaz_tmueb_ad_u)/2
  gaz_tmueb_ad_i<-abs(gaz_tmueb_ad_u-gaz_tmueb_ad_l)/2

  acap_tmueb_ad_l<-ay_tmueb_ad_d-gaz_tmueb_ad_d

  acap_tmueb_ad_u<-acap_tmueb_ad_l+ay_tmueb_ad_i-gaz_tmueb_ad_i

  acap_tmueb_ad_d<-(acap_tmueb_ad_l+acap_tmueb_ad_u)/2
  acap_tmueb_ad_i<-abs(acap_tmueb_ad_u-acap_tmueb_ad_l)/2

  gaz_tmuer_ad_l<-az_tmuer_ad_d%*%gamma_3_ad_d

  gaz_tmuer_ad_u<-gaz_tmuer_ad_l+az_tmuer_ad_d%*%gamma_3_ad_i+
    az_tmuer_ad_i%*%gamma_3_ad_d+az_tmuer_ad_i%*%gamma_3_ad_i

  gaz_tmuer_ad_d<-(gaz_tmuer_ad_l+gaz_tmuer_ad_u)/2
  gaz_tmuer_ad_i<-abs(gaz_tmuer_ad_u-gaz_tmuer_ad_l)/2

  acap_tmuer_ad_l<-ay_tmuer_ad_d-gaz_tmuer_ad_d

  acap_tmuer_ad_u<-acap_tmuer_ad_l+ay_tmuer_ad_i-gaz_tmuer_ad_i

  acap_tmuer_ad_d<-(acap_tmuer_ad_l+acap_tmuer_ad_u)/2
  acap_tmuer_ad_i<-abs(acap_tmuer_ad_u-acap_tmuer_ad_l)/2

  gaz_muebr_ad_l<-az_muebr_ad_d%*%gamma_4_ad_d

  gaz_muebr_ad_u<-gaz_muebr_ad_l+az_muebr_ad_d%*%gamma_4_ad_i+
    az_muebr_ad_i%*%gamma_4_ad_d+az_muebr_ad_i%*%gamma_4_ad_i

  gaz_muebr_ad_d<-(gaz_muebr_ad_l+gaz_muebr_ad_u)/2
  gaz_muebr_ad_i<-abs(gaz_muebr_ad_u-gaz_muebr_ad_l)/2

  acap_muebr_ad_l<-ay_muebr_ad_d-gaz_muebr_ad_d

  acap_muebr_ad_u<-acap_muebr_ad_l+ay_muebr_ad_i-gaz_muebr_ad_i

  acap_muebr_ad_d<-(acap_muebr_ad_l+acap_muebr_ad_u)/2
  acap_muebr_ad_i<-abs(acap_muebr_ad_u-acap_muebr_ad_l)/2

  acapss_all_ad_l<-Qy_all_ad_d*acap_all_ad_d

  acapss_all_ad_u<-acapss_all_ad_l+Qy_all_ad_d*acap_all_ad_i+
    Qy_all_ad_i*acap_all_ad_d+Qy_all_ad_i*acap_all_ad_i

  acapss_all_ad_d<-(acapss_all_ad_u+acapss_all_ad_l)/2
  acapss_all_ad_i<-abs(acapss_all_ad_u-acapss_all_ad_l)/2

  acapss_tmueb_ad_l<-Qy_tmueb_ad_d*acap_tmueb_ad_d

  acapss_tmueb_ad_u<-acapss_tmueb_ad_l+Qy_tmueb_ad_d*acap_tmueb_ad_i+
    Qy_tmueb_ad_i*acap_tmueb_ad_d+Qy_tmueb_ad_i*acap_tmueb_ad_i

  acapss_tmueb_ad_d<-(acapss_tmueb_ad_u+acapss_tmueb_ad_l)/2
  acapss_tmueb_ad_i<-abs(acapss_tmueb_ad_u-acapss_tmueb_ad_l)/2

  acapss_tmuer_ad_l<-Qy_tmuer_ad_d*acap_tmuer_ad_d

  acapss_tmuer_ad_u<-acapss_tmuer_ad_l+Qy_tmuer_ad_d*acap_tmuer_ad_i+
    Qy_tmuer_ad_i*acap_tmuer_ad_d+Qy_tmuer_ad_i*acap_tmuer_ad_i

  acapss_tmuer_ad_d<-(acapss_tmuer_ad_u+acapss_tmuer_ad_l)/2
  acapss_tmuer_ad_i<-abs(acapss_tmuer_ad_u-acapss_tmuer_ad_l)/2

  acapss_muebr_ad_l<-Qy_muebr_ad_d*acap_muebr_ad_d

  acapss_muebr_ad_u<-acapss_muebr_ad_l+Qy_muebr_ad_d*acap_muebr_ad_i+
    Qy_muebr_ad_i*acap_muebr_ad_d+Qy_muebr_ad_i*acap_muebr_ad_i

  acapss_muebr_ad_d<-(acapss_muebr_ad_u+acapss_muebr_ad_l)/2
  acapss_muebr_ad_i<-abs(acapss_muebr_ad_u-acapss_muebr_ad_l)/2

  ss_all_ad_l<-sum(acapss_all_ad_l)
  ss_all_ad_u<-sum(acapss_all_ad_u)

  ss_all_ad_d<-(ss_all_ad_l+ss_all_ad_u)/2
  ss_all_ad_i<-sum(acapss_all_ad_i)

  ss_tmueb_ad_l<-sum(acapss_tmueb_ad_l)
  ss_tmueb_ad_u<-sum(acapss_tmueb_ad_u)

  ss_tmueb_ad_d<-(ss_tmueb_ad_l+ss_tmueb_ad_u)/2
  ss_tmueb_ad_i<-sum(acapss_tmueb_ad_i)

  ss_tmuer_ad_l<-sum(acapss_tmuer_ad_l)
  ss_tmuer_ad_u<-sum(acapss_tmuer_ad_u)

  ss_tmuer_ad_d<-(ss_tmuer_ad_l+ss_tmuer_ad_u)/2
  ss_tmuer_ad_i<-sum(acapss_tmuer_ad_i)

  ss_muebr_ad_l<-sum(acapss_muebr_ad_l)
  ss_muebr_ad_u<-sum(acapss_muebr_ad_u)

  ss_muebr_ad_d<-(ss_muebr_ad_l+ss_muebr_ad_u)/2
  ss_muebr_ad_i<-sum(acapss_muebr_ad_i)

  gzy_1_ad_l<-gamma_1_ad_d*zpy_d

  gzy_1_ad_u<-gzy_1_ad_l+gamma_1_ad_d*zpy_i+
    gamma_1_ad_i*zpy_d+gamma_1_ad_i*zpy_i

  gzy_1_ad_d<-(gzy_1_ad_l+gzy_1_ad_u)/2
  gzy_1_ad_i<-abs(gzy_1_ad_u-gzy_1_ad_l)/2

  gzy_2_ad_l<-gamma_2_ad_d*zpy_d

  gzy_2_ad_u<-gzy_2_ad_l+gamma_2_ad_d*zpy_i+
    gamma_2_ad_i*zpy_d+gamma_2_ad_i*zpy_i

  gzy_2_ad_d<-(gzy_2_ad_l+gzy_2_ad_u)/2
  gzy_2_ad_i<-abs(gzy_2_ad_u-gzy_2_ad_l)/2

  gzy_3_ad_l<-gamma_3_ad_d*zpy_d

  gzy_3_ad_u<-gzy_3_ad_l+gamma_3_ad_d*zpy_i+
    gamma_3_ad_i*zpy_d+gamma_3_ad_i*zpy_i

  gzy_3_ad_d<-(gzy_3_ad_l+gzy_3_ad_u)/2
  gzy_3_ad_i<-abs(gzy_3_ad_u-gzy_3_ad_l)/2

  gzy_4_ad_l<-gamma_4_ad_d*zpy_d

  gzy_4_ad_u<-gzy_4_ad_l+gamma_4_ad_d*zpy_i+
    gamma_4_ad_i*zpy_d+gamma_4_ad_i*zpy_i

  gzy_4_ad_d<-(gzy_4_ad_l+gzy_4_ad_u)/2
  gzy_4_ad_i<-abs(gzy_4_ad_u-gzy_4_ad_l)/2

  mdss_all_ad_l<-ss_all_ad_d+gzy_1_ad_d

  mdss_all_ad_u<-mdss_all_ad_l+ss_all_ad_i+gzy_1_ad_i

  mdss_all_ad_d<-(mdss_all_ad_u+mdss_all_ad_l)/2
  mdss_all_ad_i<-abs(mdss_all_ad_u-mdss_all_ad_l)/2

  mdss_tmueb_ad_l<-ss_tmueb_ad_d+gzy_2_ad_d

  mdss_tmueb_ad_u<-mdss_tmueb_ad_l+ss_tmueb_ad_i+gzy_2_ad_i

  mdss_tmueb_ad_d<-(mdss_tmueb_ad_u+mdss_tmueb_ad_l)/2
  mdss_tmueb_ad_i<-abs(mdss_tmueb_ad_u-mdss_tmueb_ad_l)/2

  mdss_tmuer_ad_l<-ss_tmuer_ad_d+gzy_3_ad_d

  mdss_tmuer_ad_u<-mdss_tmuer_ad_l+ss_tmuer_ad_i+gzy_3_ad_i

  mdss_tmuer_ad_d<-(mdss_tmuer_ad_u+mdss_tmuer_ad_l)/2
  mdss_tmuer_ad_i<-abs(mdss_tmuer_ad_u-mdss_tmuer_ad_l)/2

  mdss_muebr_ad_l<-ss_muebr_ad_d+gzy_4_ad_d

  mdss_muebr_ad_u<-mdss_muebr_ad_l+ss_muebr_ad_i+gzy_4_ad_i

  mdss_muebr_ad_d<-(mdss_muebr_ad_u+mdss_muebr_ad_l)/2
  mdss_muebr_ad_i<-abs(mdss_muebr_ad_u-mdss_muebr_ad_l)/2

  ####error sum of sq

  ess_all_ad_l<-ypy_d-mdss_all_ad_d

  ess_all_ad_u<-ess_all_ad_l+ypy_i-mdss_all_ad_i

  ess_all_ad_d<-(ess_all_ad_l+ess_all_ad_u)/2
  ess_all_ad_i<-abs(ess_all_ad_u-ess_all_ad_l)/2

  esst_ad_l<-ypy_d-mdss_muebr_ad_d

  esst_ad_u<-esst_ad_l+ypy_i-mdss_muebr_ad_i

  esst_ad_d<-(esst_ad_u+esst_ad_l)/2
  esst_ad_i<-abs(esst_ad_u-esst_ad_l)/2

  essb_ad_l<-ypy_d-mdss_tmuer_ad_d

  essb_ad_u<-essb_ad_l+ypy_i-mdss_tmuer_ad_i

  essb_ad_d<-(essb_ad_u+essb_ad_l)/2
  essb_ad_i<-abs(essb_ad_u-essb_ad_l)/2

  essr_ad_l<-ypy_d-mdss_tmueb_ad_d

  essr_ad_u<-essr_ad_l+ypy_i-mdss_tmueb_ad_i

  essr_ad_d<-(essr_ad_u+essr_ad_l)/2
  essr_ad_i<-abs(essr_ad_u-essr_ad_l)/2


  #####covariate adj trt ss adj to block

  tss_ad_l<-esst_ad_d-ess_all_ad_d
  tss_ad_u<-tss_ad_l+esst_ad_i-ess_all_ad_i
  tss_ad_d<-(tss_ad_u+tss_ad_l)/2
  tss_ad_i<-abs(tss_ad_u-tss_ad_l)/2

  bss_ad_l<-essb_ad_d-ess_all_ad_d
  bss_ad_u<-bss_ad_l+essb_ad_i-ess_all_ad_i
  bss_ad_d<-(bss_ad_u+bss_ad_l)/2
  bss_ad_i<-abs(bss_ad_u-bss_ad_l)/2

  rss_ad_l<-essr_ad_d-ess_all_ad_d
  rss_ad_u<-rss_ad_l+essr_ad_i-ess_all_ad_i
  rss_ad_d<-(rss_ad_u+rss_ad_l)/2
  rss_ad_i<-abs(rss_ad_u-rss_ad_l)/2

  ###cov ss

  gss_l<-gamma_1_ad_d*EXy_all_ad_d

  gss_u<-gss_l+gamma_1_ad_d*EXy_all_ad_i+
    gamma_1_ad_i*EXy_all_ad_d+
    gamma_1_ad_i*EXy_all_ad_i

  gss_d<-(gss_l+gss_u)/2
  gss_i<-abs(gss_u-gss_l)/2

  ##totalss

  toss_l<-t(dy)%*%(dy)

  toss_u<-toss_l+2*t(dy)%*%iy+t(iy)%*%iy

  toss_d<-(toss_l+toss_u)/2
  toss_i<-abs(toss_u-toss_l)/2

  ####Mean Squares

  msg_l<-gss_l/1
  msg_u<-gss_u/1

  msr_l<-rss_ad_l/(v-1)
  msr_u<-rss_ad_u/(v-1)

  mst_l<-tss_ad_l/(v-1)
  mst_u<-tss_ad_u/(v-1)

  msb_l<-bss_ad_l/(v-1)
  msb_u<-bss_ad_u/(v-1)

  mse_l<-ess_all_ad_l/((v-1)*(v-2)-1)
  mse_u<-ess_all_ad_u/((v-1)*(v-2)-1)

  ####F Statistics

  F_g_L<-(gss_d/1)/(ess_all_ad_d/((v-1)*(v-2)-1))

  F_g_U<-F_g_L+((gss_d/1)*(ess_all_ad_i/((v-1)*(v-2)-1))+
                  (gss_i/1)*(ess_all_ad_d/((v-1)*(v-2)-1))+
                  (gss_i/1)*(ess_all_ad_i/((v-1)*(v-2)-1)))/
    ((ess_all_ad_d/((v-1)*(v-2)-1))*(ess_all_ad_d/((v-1)*(v-2)-1)))

  F_b_ad_L<-((bss_ad_d/(v-1))/(ess_all_ad_d/((v-1)*(v-2)-1)))

  F_b_ad_U<-F_b_ad_L+((((bss_ad_d/(v-1))*(ess_all_ad_i/((v-1)*(v-2)-1)))+
                         ((bss_ad_i/(v-1))*(ess_all_ad_d/((v-1)*(v-2)-1)))+
                         ((bss_ad_i/(v-1))*(ess_all_ad_i/((v-1)*(v-2)-1))))/
                        ((ess_all_ad_d/((v-1)*(v-2)-1))*(ess_all_ad_d/((v-1)*(v-2)-1))))

  F_r_ad_L<-((rss_ad_d/(v-1))/(ess_all_ad_d/((v-1)*(v-2)-1)))

  F_r_ad_U<-F_r_ad_L+((((rss_ad_d/(v-1))*(ess_all_ad_i/((v-1)*(v-2)-1)))+
                         ((rss_ad_i/(v-1))*(ess_all_ad_d/((v-1)*(v-2)-1)))+
                         (rss_ad_i/(v-1))*(ess_all_ad_i/((v-1)*(v-2)-1)))/
                        ((ess_all_ad_d/((v-1)*(v-2)-1))*(ess_all_ad_d/((v-1)*(v-2)-1))))

  F_t_ad_L<-((tss_ad_d/(v-1))/(ess_all_ad_d/((v-1)*(v-2)-1)))

  F_t_ad_U<-F_t_ad_L+((((tss_ad_d/(v-1))*(ess_all_ad_i/((v-1)*(v-2)-1)))+
                         ((tss_ad_i/(v-1))*(ess_all_ad_d/((v-1)*(v-2)-1)))+
                         ((tss_ad_i/(v-1))*(ess_all_ad_i/((v-1)*(v-2)-1))))/
                        ((ess_all_ad_d/((v-1)*(v-2)-1))*(ess_all_ad_d/((v-1)*(v-2)-1))))

  ####Significance

  F_critical_5 <- qf(1-alpha,v-1,((v-1)*(v-2)-1))
  F_critical_1 <- qf(1-alpha/5,v-1,((v-1)*(v-2)-1))
  F_critical_001 <- qf(1-alpha/50,v-1,((v-1)*(v-2)-1))

  significance <- c(

    ifelse(F_g_L > F_critical_001 & F_g_U > F_critical_001,"***",
           ifelse(F_g_L > F_critical_1 & F_g_U > F_critical_1,"**",
                  ifelse(F_g_L > F_critical_5 & F_g_U > F_critical_5,"*",
                         ifelse(F_g_L < F_critical_5 & F_g_U < F_critical_5,"NS","ID")))),

    ifelse(F_r_ad_L > F_critical_001 & F_r_ad_U > F_critical_001,"***",
           ifelse(F_r_ad_L > F_critical_1 & F_r_ad_U > F_critical_1,"**",
                  ifelse(F_r_ad_L > F_critical_5 & F_r_ad_U > F_critical_5,"*",
                         ifelse(F_r_ad_L < F_critical_5 & F_r_ad_U < F_critical_5,"NS","ID")))),

    ifelse(F_t_ad_L > F_critical_001 & F_t_ad_U > F_critical_001,"***",
           ifelse(F_t_ad_L > F_critical_1 & F_t_ad_U > F_critical_1,"**",
                  ifelse(F_t_ad_L > F_critical_5 & F_t_ad_U > F_critical_5,"*",
                         ifelse(F_t_ad_L < F_critical_5 & F_t_ad_U < F_critical_5,"NS","ID")))),

    ifelse(F_b_ad_L > F_critical_001 & F_b_ad_U > F_critical_001,"***",
           ifelse(F_b_ad_L > F_critical_1 & F_b_ad_U > F_critical_1,"**",
                  ifelse(F_b_ad_L > F_critical_5 & F_b_ad_U > F_critical_5,"*",
                         ifelse(F_b_ad_L < F_critical_5 & F_b_ad_U < F_critical_5,"NS","ID")))),

    "",
    ""
  )

  ####Formatting

  format_interval <- function(x){
    if(any(is.na(x))) return("")
    paste0("[", round(x[1],2), ", ", round(x[2],2), "]")
  }

  ####ANCOVA Table

  SV <- c("Covariate","Row","Treatment","Column","Error","Total")

  DF <- c(1,v-1,v-1,v-1,(v-1)*(v-2)-1,v*v-1)

  SSN <- c(
    format_interval(c(gss_l,gss_u)),
    format_interval(c(rss_ad_l,rss_ad_u)),
    format_interval(c(tss_ad_l,tss_ad_u)),
    format_interval(c(bss_ad_l,bss_ad_u)),
    format_interval(c(ess_all_ad_l,ess_all_ad_u)),
    format_interval(c(toss_l,toss_u))
  )

  MSSN <- c(
    format_interval(c(msg_l,msg_u)),
    format_interval(c(msr_l,msr_u)),
    format_interval(c(mst_l,mst_u)),
    format_interval(c(msb_l,msb_u)),
    format_interval(c(mse_l,mse_u)),
    ""
  )

  FN <- c(
    format_interval(c(F_g_L,F_g_U)),
    format_interval(c(F_r_ad_L,F_r_ad_U)),
    format_interval(c(F_t_ad_L,F_t_ad_U)),
    format_interval(c(F_b_ad_L,F_b_ad_U)),
    "",
    ""
  )

  nancova_table <- data.frame(
    SV = SV,
    DF = DF,
    SSN = SSN,
    MSSN = MSSN,
    FN = FN,
    Significance = significance,
    stringsAsFactors = FALSE
  )


  ####LSD Test

  if(significance[3] != "NS"){

    LSD_l <- qt(1-alpha/2,((v-1)*(v-2)-1))*sqrt((2*mse_l)/r)

    LSD_u <- LSD_l+
      qt(1-alpha/2,((v-1)*(v-2)-1))*sqrt((2*mse_u)/r)-
      qt(1-alpha/2,((v-1)*(v-2)-1))*sqrt((2*mse_l)/r)

    trt_means_l <- apply(Lower_y,2,mean)
    trt_means_u <- apply(Upper_y,2,mean)

    trt_means_d <- (trt_means_l+trt_means_u)/2
    trt_means_i <- abs(trt_means_u-trt_means_l)/2

    comparison <- data.frame(
      Treatment1 = character(),
      Treatment2 = character(),
      Mean_Difference = character(),
      Decision = character(),
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
            Treatment1 = paste0("T",i),
            Treatment2 = paste0("T",j),
            Mean_Difference = paste0("[",round(diff_l,4),", ",round(diff_u,4),"]"),
            Decision = decision,
            stringsAsFactors = FALSE
          )
        )

      }
    }

    if(verbose){

      message("\nNeutrosophic Analysis of Covariance Table")
      print(nancova_table, row.names = FALSE)

      if(significance[3] != "NS"){

        message("\nTreatment effect is significant. Hence multiple comparison using LSD is performed.")

        message(
          "\nLSD Interval : [",
          round(LSD_l, 4), ", ",
          round(LSD_u, 4), "]"
        )

        message("\nTreatment Comparisons Using LSD")
        print(comparison, row.names = FALSE)

      } else {

        message("\nTreatment effect is non significant. Hence multiple comparison is not performed.")

      }

      message("\nSignificance Codes:")
      message("*** : Significant at p < 0.001")
      message("**  : Significant at p < 0.01")
      message("*   : Significant at p < 0.05")
      message("S  : Significant")
      message("NS : Non Significant")
      message("ID : Indeterminate")
    }
    invisible(list(
      nancova_table = nancova_table,
      comparison = if(significance[3] != "NS") comparison else NULL,
      LSD = if(significance[3] != "NS")
        c(Lower = LSD_l, Upper = LSD_u) else NULL
    ))
  }
}
