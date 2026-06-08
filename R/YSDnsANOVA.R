#' Neutrosophic Analysis of Variance for Youden Square Design
#'
#' Performs Neutrosophic Analysis of Variance (NANOVA) for Youden
#' square designs using interval-valued observations. The function
#' computes neutrosophic sums of squares, mean squares,
#' interval-valued F-statistics, significance tests, and Least
#' Significant Difference (LSD)-based treatment comparisons.
#' For crisp data, enter identical lower and upper values to obtain
#' the corresponding classical ANOVA results.
#'
#' @usage YSDnsANOVA(Lower_y, Upper_y, design, alpha = 0.05,  verbose = FALSE)
#'
#' @param Lower_y Matrix containing lower bounds of observations.
#' @param Upper_y Matrix containing upper bounds of observations.
#' @param design Matrix representing Youden square treatment allocation.
#' @param alpha Significance level for the F-test and LSD test. Default is 0.05.
#' @param verbose Logical. If TRUE, displays the NANOVA table,
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
#' }
#'
#' @return
#' A list containing:
#' \itemize{
#' \item \code{nanova_table}: Neutrosophic ANOVA table.
#' \item \code{comparison}: LSD treatment comparisons, if applicable.
#' \item \code{LSD}: Lower and upper limits of the LSD interval.
#' }
#'
#' @examples
#' Lower_y <- matrix(c(
#' 4.51,7.77,3.53,7.92,5.97,7.42,6.37,
#' 6.37,4.36,4.68,8.19,6.08,6.80,4.86,
#' 4.88,4.64,8.69,7.81,9.50,9.60,9.10,
#' 8.05,7.92,5.61,8.41,6.53,11.36,5.94
#' ), nrow = 4, byrow = TRUE)
#'
#' Upper_y <- matrix(c(
#' 8.12,11.06,6.44,11.16,8.30,10.17,9.81,
#' 10.31,7.50,8.01,10.82,8.49,10.69,7.77,
#' 7.76,7.24,10.82,11.04,13.20,11.63,11.59,
#' 11.66,10.45,8.92,12.03,8.95,14.64,9.14
#' ), nrow = 4, byrow = TRUE)
#'
#' design <- matrix(c(
#' 2,3,4,5,6,7,1,
#' 7,1,2,3,4,5,6,
#' 6,7,1,2,3,4,5,
#' 5,6,7,1,2,3,4
#' ), nrow = 4, byrow = TRUE)
#'
#' YSDnsANOVA(Lower_y, Upper_y, design, alpha = 0.05,  verbose = TRUE)
#'
#' @importFrom MASS ginv
#' @importFrom stats qf qt
#' @export

YSDnsANOVA <- function(Lower_y, Upper_y, design, alpha = 0.05,  verbose = FALSE){
  if(!is.matrix(Lower_y) || !is.matrix(Upper_y) || !is.matrix(design)){stop("Inputs must be matrices.")}
  if(any(dim(Lower_y) != dim(Upper_y))){stop("Lower_y and Upper_y must have same dimensions.")}
  if(any(dim(Lower_y) != dim(design))){stop("Design and response matrices must have same dimensions.")}
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
  Bmat <- diag(v)[rep(1:k,times=b),]
  rep_vec <- rep(1:r,each=b/r*k)
  Rmat <- matrix(0,n,r)
  for(i in 1:n){
    Rmat[i,rep_vec[i]] <- 1
  }
  X <- cbind(Tmat,mu,Bmat,Rmat)

  #### Response decomposition
  d_y <- (Lower_y+Upper_y)/2
  i_y <- abs(Upper_y-Lower_y)/2
  dy <- as.vector(t(d_y))
  iy <- as.vector(t(i_y))
  ybar_l <- mean(Lower_y)
  ybar_u <- mean(Upper_y)
  ybar_d <- mean(d_y)
  ybar_i <- mean(i_y)
  yminusbar_l <- d_y-ybar_d
  yminusbar_u <- yminusbar_l+i_y-ybar_i
  yminusbar_d <- (yminusbar_l+yminusbar_u)/2
  yminusbar_i <- abs(yminusbar_u-yminusbar_l)/2
  yminusbar_sq_l <- yminusbar_d*yminusbar_d
  yminusbar_sq_u <- yminusbar_sq_l+2*yminusbar_d*yminusbar_i+yminusbar_i*yminusbar_i
  yminusbar_sq_d <- (yminusbar_sq_l+yminusbar_sq_u)/2
  yminusbar_sq_i <- abs(yminusbar_sq_u-yminusbar_sq_l)/2

  ##totalss
  toss_l <- sum(yminusbar_sq_l)
  toss_u <- sum(yminusbar_sq_u)
  toss_d <- sum(yminusbar_sq_d)
  toss_i <- abs(toss_l-toss_u)/2

  ##Trt unadj
  Ct_unadj <- t(Tmat) %*% Tmat
  bcdy_t_un_d<-t(cbind(Tmat,Bmat,Rmat))%*%mu%*%ginv(t(mu)%*%mu)%*%t(mu)%*%dy
  bcdy_t_un_i<-t(cbind(Tmat,Bmat,Rmat))%*%mu%*%ginv(t(mu)%*%mu)%*%t(mu)%*%iy
  Qy_un_l<- (t(cbind(Tmat,Bmat,Rmat))%*%dy-bcdy_t_un_d)
  Qy_un_u<- Qy_un_l+t(cbind(Tmat,Bmat,Rmat))%*%iy-bcdy_t_un_i
  Qy_t_un_d<-(Qy_un_l[1:v,]+Qy_un_u[1:v,])/2
  Qy_t_un_i<-abs(Qy_un_l[1:v,]-Qy_un_u[1:v,])/2

  ##Trt adj
  Z <- cbind(mu,Bmat,Rmat)
  M <- diag(nrow(X)) - Z %*% ginv(t(Z) %*% Z) %*% t(Z)
  Ct_adj <- t(Tmat) %*% M %*% Tmat
  bcdy_t_ad_d<-t(Tmat)%*%cbind(mu,Bmat,Rmat)%*%ginv(t(cbind(mu,Bmat,Rmat))%*%cbind(mu,Bmat,Rmat))%*%t(cbind(mu,Bmat,Rmat))%*%dy
  bcdy_t_ad_i<-t(Tmat)%*%cbind(mu,Bmat,Rmat)%*%ginv(t(cbind(mu,Bmat,Rmat))%*%cbind(mu,Bmat,Rmat))%*%t(cbind(mu,Bmat,Rmat))%*%iy
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
  bcdy_b_un_d<-t(cbind(Bmat,Tmat,Rmat))%*%mu%*%ginv(t(mu)%*%mu)%*%t(mu)%*%dy
  bcdy_b_un_i<-t(cbind(Bmat,Tmat,Rmat))%*%mu%*%ginv(t(mu)%*%mu)%*%t(mu)%*%iy
  Qy_un_l<- (t(cbind(Bmat,Tmat,Rmat))%*%dy-bcdy_b_un_d)
  Qy_un_u<- Qy_un_l+t(cbind(Bmat,Tmat,Rmat))%*%iy-bcdy_b_un_i
  Qy_b_un_d<-(Qy_un_l[1:k,]+Qy_un_u[1:k,])/2
  Qy_b_un_i<-abs(Qy_un_l[1:k,]-Qy_un_u[1:k,])/2

  ##B adj y
  Z <- cbind(mu,Tmat,Rmat)
  M <- diag(nrow(X)) - Z %*% ginv(t(Z) %*% Z) %*% t(Z)
  Cb_adj <- t(Bmat) %*% M %*% Bmat
  bcdy_b_ad_d<-t(Bmat)%*%cbind(mu,Tmat,Rmat)%*%ginv(t(cbind(mu,Tmat,Rmat))%*%cbind(mu,Tmat,Rmat))%*%t(cbind(mu,Tmat,Rmat))%*%dy
  bcdy_b_ad_i<-t(Bmat)%*%cbind(mu,Tmat,Rmat)%*%ginv(t(cbind(mu,Tmat,Rmat))%*%cbind(mu,Tmat,Rmat))%*%t(cbind(mu,Tmat,Rmat))%*%iy
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
  bcdy_r_un_d<-t(cbind(Rmat,Bmat,Tmat))%*%mu%*%ginv(t(mu)%*%mu)%*%t(mu)%*%dy
  bcdy_r_un_i<-t(cbind(Rmat,Bmat,Tmat))%*%mu%*%ginv(t(mu)%*%mu)%*%t(mu)%*%iy
  Qy_un_l<- (t(cbind(Rmat,Bmat,Tmat))%*%dy-bcdy_r_un_d)
  Qy_un_u<- Qy_un_l+t(cbind(Rmat,Bmat,Tmat))%*%iy-bcdy_r_un_i
  Qy_r_un_d<-(Qy_un_l[1:r,]+Qy_un_u[1:r,])/2
  Qy_r_un_i<-abs(Qy_un_l[1:r,]-Qy_un_u[1:r,])/2
  ay_r_l <- ginv(Cr_unadj) %*% Qy_r_un_d
  ay_r_u <- ay_r_l+ginv(Cr_unadj)%*%Qy_r_un_i
  ay_r_d<-(ay_r_l +ay_r_u)/2
  ay_r_i<-abs(ay_r_l-ay_r_u)/2

  ####sum of square

  ss_t_ad_l<-sum(Qy_t_ad_d*ay_t_ad_d)
  ss_t_ad_u<-ss_t_ad_l+sum(Qy_t_ad_d*ay_t_ad_i+Qy_t_ad_i*ay_t_ad_d+Qy_t_ad_i*ay_t_ad_i)
  ss_t_ad_d<-(ss_t_ad_l+ss_t_ad_u)/2
  ss_t_ad_i<-abs(ss_t_ad_u-ss_t_ad_l)/2

  ss_t_un_l<-sum(Qy_t_un_d*ay_t_un_d)
  ss_t_un_u<-ss_t_un_l+sum(Qy_t_un_d*ay_t_un_i+Qy_t_un_i*ay_t_un_d+Qy_t_un_i*ay_t_un_i)
  ss_t_un_d<-(ss_t_un_l+ss_t_un_u)/2
  ss_t_un_i<-abs(ss_t_un_u-ss_t_un_l)/2

  ss_b_ad_l<-sum(Qy_b_ad_d*ay_b_ad_d)
  ss_b_ad_u<-ss_b_ad_l+sum(Qy_b_ad_d*ay_b_ad_i+Qy_b_ad_i*ay_b_ad_d+Qy_b_ad_i*ay_b_ad_i)
  ss_b_ad_d<-(ss_b_ad_l+ss_b_ad_u)/2
  ss_b_ad_i<-abs(ss_b_ad_u-ss_b_ad_l)/2

  ss_b_un_l<-sum(Qy_b_un_d*ay_b_un_d)
  ss_b_un_u<-ss_b_un_l+sum(Qy_b_un_d*ay_b_un_i+Qy_b_un_i*ay_b_un_d+Qy_b_un_i*ay_b_un_i)
  ss_b_un_d<-(ss_b_un_l+ss_b_un_u)/2
  ss_b_un_i<-abs(ss_b_un_u-ss_b_un_l)/2

  ss_r_un_l<-sum(Qy_r_un_d*ay_r_d)
  ss_r_un_u<-ss_r_un_l+sum(Qy_r_un_d*ay_r_i+Qy_r_un_i*ay_r_d+Qy_r_un_i*ay_r_i)
  ss_r_un_d<-(ss_r_un_l+ss_r_un_u)/2
  ss_r_un_i<-abs(ss_r_un_u-ss_r_un_l)/2

  mdss_all_ad_l<-ss_t_ad_d+ss_b_un_d+ss_r_un_d
  mdss_all_ad_u<-mdss_all_ad_l+ss_t_ad_i+ss_b_un_i+ss_r_un_i
  mdss_all_ad_d<-(mdss_all_ad_u+mdss_all_ad_l)/2
  mdss_all_ad_i<-abs(mdss_all_ad_u-mdss_all_ad_l)/2

  mdss_all_un_l<-ss_t_un_d+ss_b_ad_d+ss_r_un_d
  mdss_all_un_u<-mdss_all_un_l+ss_t_un_i+ss_b_ad_i+ss_r_un_i
  mdss_all_un_d<-(mdss_all_un_u+mdss_all_un_l)/2
  mdss_all_un_i<-abs(mdss_all_un_u-mdss_all_un_l)/2

  ess_all_ad_l<-toss_d-mdss_all_ad_d
  ess_all_ad_u<-ess_all_ad_l+toss_i-mdss_all_ad_i
  ess_all_ad_d<-(ess_all_ad_l+ess_all_ad_u)/2
  ess_all_ad_i<-abs(ess_all_ad_u-ess_all_ad_l)/2

  ess_all_un_l<-toss_d-mdss_all_un_d
  ess_all_un_u<-ess_all_un_l+toss_i-mdss_all_un_i
  ess_all_un_d<-(ess_all_un_l+ess_all_un_u)/2
  ess_all_un_i<-abs(ess_all_un_u-ess_all_un_l)/2

  SV <- c("Row","Treatment unadj","Column adj","Error",
          "Column unadj","Treatment adj","Error","Total")

  DF <- c(b-1,v-1,k-1,(v*r-b-v-k+2),
          k-1,v-1,(v*r-b-v-k+2),v*r-1)

  SSN <- list(
    c(ss_r_un_l,ss_r_un_u),
    c(ss_t_un_l,ss_t_un_u),
    c(ss_b_ad_l,ss_b_ad_u),
    c(ess_all_un_l,ess_all_un_u),
    c(ss_b_un_l,ss_b_un_u),
    c(ss_t_ad_l,ss_t_ad_u),
    c(ess_all_ad_l,ess_all_ad_u),
    c(toss_l,toss_u)
  )

  MSSN <- list(
    c(ss_r_un_l/(b-1),ss_r_un_u/(b-1)),
    c(ss_t_un_l/(v-1),ss_t_un_u/(v-1)),
    c(ss_b_ad_l/(k-1),ss_b_ad_u/(k-1)),
    c(ess_all_un_l/((v*r-b-v-k+2)),ess_all_un_u/((v*r-b-v-k+2))),
    c(ss_b_un_l/(k-1),ss_b_un_u/(k-1)),
    c(ss_t_ad_l/(v-1),ss_t_ad_u/(v-1)),
    c(ess_all_ad_l/((v*r-b-v-k+2)),ess_all_ad_u/((v*r-b-v-k+2))),
    c(NA, NA)
  )

  F_r_un_L<-(ss_r_un_d/(b-1))/(ess_all_un_d/(v*r-b-v-k+2))
  F_r_un_U<-F_r_un_L+(((ss_r_un_d/(b-1))*(ess_all_un_i/(v*r-b-v-k+2)))+
                        ((ss_r_un_i/(b-1))*(ess_all_un_d/(v*r-b-v-k+2)))+
                        ((ss_r_un_i/(b-1))*(ess_all_un_i/(v*r-b-v-k+2))))/
    ((ess_all_un_d/(v*r-b-v-k+2))*(ess_all_un_d/(v*r-b-v-k+2)))

  F_t_un_L<-(ss_t_un_d/(v-1))/(ess_all_un_d/(v*r-b-v-k+2))
  F_t_un_U<-F_t_un_L+(((ss_t_un_d/(v-1))*(ess_all_un_i/(v*r-b-v-k+2)))+
                        ((ss_t_un_i/(v-1))*(ess_all_un_d/(v*r-b-v-k+2)))+
                        ((ss_t_un_i/(v-1))*(ess_all_un_i/(v*r-b-v-k+2))))/
    ((ess_all_un_d/(v*r-b-v-k+2))*(ess_all_un_d/(v*r-b-v-k+2)))

  F_b_ad_L<-(ss_b_ad_d/(k-1))/(ess_all_un_d/(v*r-b-v-k+2))
  F_b_ad_U<-F_b_ad_L+(((ss_b_ad_d/(k-1))*(ess_all_un_i/(v*r-b-v-k+2)))+
                        ((ss_b_ad_i/(k-1))*(ess_all_un_d/(v*r-b-v-k+2)))+
                        ((ss_b_ad_i/(k-1))*(ess_all_un_i/(v*r-b-v-k+2))))/
    ((ess_all_un_d/(v*r-b-v-k+2))*(ess_all_un_d/(v*r-b-v-k+2)))

  F_t_ad_L<-(ss_t_ad_d/(v-1))/(ess_all_ad_d/(v*r-b-v-k+2))
  F_t_ad_U<-F_t_ad_L+(((ss_t_ad_d/(v-1))*(ess_all_ad_i/(v*r-b-v-k+2)))+
                        ((ss_t_ad_i/(v-1))*(ess_all_ad_d/(v*r-b-v-k+2)))+
                        ((ss_t_ad_i/(v-1))*(ess_all_ad_i/(v*r-b-v-k+2))))/
    ((ess_all_ad_d/(v*r-b-v-k+2))*(ess_all_ad_d/(v*r-b-v-k+2)))

  F_b_un_L<-(ss_b_un_d/(k-1))/(ess_all_ad_d/(v*r-b-v-k+2))
  F_b_un_U<-F_b_un_L+(((ss_b_un_d/(k-1))*(ess_all_ad_i/(v*r-b-v-k+2)))+
                        ((ss_b_un_i/(k-1))*(ess_all_ad_d/(v*r-b-v-k+2)))+
                        ((ss_b_un_i/(k-1))*(ess_all_ad_i/(v*r-b-v-k+2))))/
    ((ess_all_ad_d/(v*r-b-v-k+2))*(ess_all_ad_d/(v*r-b-v-k+2)))

  FN <- list(
    c(F_r_un_L,F_r_un_U),
    c(F_t_un_L,F_t_un_U),
    c(F_b_ad_L,F_b_ad_U),
    c(NA, NA),
    c(F_b_un_L,F_b_un_U),
    c(F_t_ad_L,F_t_ad_U),
    c(NA, NA),
    c(NA, NA)
  )

  format_interval <- function(x){
    if(any(is.na(x))) return("")
    paste0("[",round(x[1],2),", ",round(x[2],2),"]")
  }

  SSN_fmt <- sapply(SSN,format_interval)
  MSSN_fmt <- sapply(MSSN,format_interval)
  FN_fmt <- sapply(FN,format_interval)

  F_critical_5 <- qf(1-alpha,DF[2],DF[4])
  F_critical_1 <- qf(1-alpha/5,DF[2],DF[4])
  F_critical_001 <- qf(1-alpha/50,DF[2],DF[4])

  significance <- c(

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

  nanova_table <- data.frame(
    SV = SV,
    DF = DF,
    SSN = SSN_fmt,
    MSSN = MSSN_fmt,
    FN = FN_fmt,
    Significance = significance,
    stringsAsFactors = FALSE
  )

  if(verbose){
  message("\nNeutrosophic Analysis of Variance Table\n")
  print(nanova_table,row.names=FALSE)

  ####LSD Test
  mse_l <- ess_all_ad_l/(v*r-b-v-k+2)
  mse_u <- ess_all_ad_u/(v*r-b-v-k+2)

  if(F_t_ad_L < F_critical_5 & F_t_ad_U < F_critical_5){

    message("\nTreatment effect is non significant. Hence multiple comparison is not performed.\n")

  }else{

    LSD_l <- qt(1-alpha/2,(v*r-b-v-k+2))*sqrt((2*mse_l)/r)

    LSD_u <- LSD_l+
      qt(1-alpha/2,(v*r-b-v-k+2))*sqrt((2*mse_u)/r)-
      qt(1-alpha/2,(v*r-b-v-k+2))*sqrt((2*mse_l)/r)

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

    message("\nTreatment effect is significant. Hence multiple comparison using LSD is performed.\n")

    message("\nLSD Interval : [",
        round(LSD_l,4),", ",
        round(LSD_u,4),"]\n")

    message("\nTreatment Comparisons Using LSD\n")

    print(comparison,row.names=FALSE)

  }
  message("\nSignificance Codes:")
  message("*** : Significant at p < 0.001")
  message("**  : Significant at p < 0.01")
  message("*   : Significant at p < 0.05")
  message("NS  : Non Significant")
  message("ID  : Indeterminate")
  }
  result <- list(
    nanova_table = nanova_table,
    comparison = if(exists("comparison")) comparison else NULL,
    LSD = if(exists("LSD_l"))
      c(Lower = LSD_l, Upper = LSD_u) else NULL
  )

  invisible(result)
}






