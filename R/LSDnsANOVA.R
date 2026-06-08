#' Neutrosophic Analysis of Variance for Latin Square Design
#'
#' Performs Neutrosophic Analysis of Variance (NANOVA) for Latin square
#' designs using interval-valued observations. The function computes
#' neutrosophic sums of squares, mean squares, interval-valued F-statistics,
#' significance tests, and Least Significant Difference (LSD)-based
#' treatment comparisons. When the data are crisp, identical lower and
#' upper values may be entered to obtain the corresponding classical
#' ANOVA results.
#'
#' @usage LSDnsANOVA(Lower_y, Upper_y, design, alpha = 0.05, verbose = FALSE)
#'
#' @param Lower_y Matrix containing lower bounds of observations.
#' @param Upper_y Matrix containing upper bounds of observations.
#' @param design Matrix representing Latin square treatment allocation.
#' @param alpha Significance level (default is 0.05.) for the F-test and LSD test.
#' @param verbose Logical. If TRUE, displays the NANOVA table,
#' LSD interval, treatment comparisons, and significance codes.
#' Default is FALSE.
#'
#' @details
#' Input matrix structure:
#' \itemize{
#'   \item Rows represent blocks (or rows of the design).
#'   \item Columns represent treatment positions within each block.
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
#' 243.28,233.45,232.41,239.27,
#' 237.38,223.30,240.86,240.33,
#' 232.12,234.67,230.45,220.81,
#' 219.26,231.43,248.44,227.88
#' ), nrow = 4, byrow = TRUE)
#'
#' Upper_y <- matrix(c(
#' 246.72,238.55,237.59,242.73,
#' 242.62,226.70,245.14,245.67,
#' 235.88,239.33,233.55,223.19,
#' 224.74,236.57,253.57,230.12
#' ), nrow = 4, byrow = TRUE)
#'
#' design <- matrix(c(
#' 3,4,2,1,
#' 1,2,4,3,
#' 4,3,1,2,
#' 2,1,3,4
#' ), nrow = 4, byrow = TRUE)
#'
#' LSDnsANOVA(Lower_y, Upper_y, design, alpha = 0.05, verbose = TRUE)
#'
#' @importFrom MASS ginv
#' @importFrom stats qf qt
#' @export


LSDnsANOVA <- function(Lower_y, Upper_y, design, alpha = 0.05, verbose = FALSE){

  if(!is.matrix(Lower_y) || !is.matrix(Upper_y) || !is.matrix(design)){
    stop("Inputs must be matrices.")
  }

  if(any(dim(Lower_y) != dim(Upper_y))){
    stop("Lower_y and Upper_y must have same dimensions.")
  }

  if(any(dim(Lower_y) != dim(design))){
    stop("Design and response matrices must have same dimensions.")
  }

  v <- max(design)
  b <- nrow(design)
  k <- ncol(design)

  rep_counts <- table(design)
  r <- unique(rep_counts)

  interval <- function(L, U){
    list(
      d = (L+U)/2,
      i = abs(U-L)/2
    )
  }

  n=b*k

  trt_vec <- as.vector(t(design))
  Tmat<- matrix(0, n, v)

  for(i in 1:n){
    Tmat[i, trt_vec[i]] <- 1
  }

  mu <- matrix(1, n, 1)

  blk <- rep(1:v, each = v)
  Bmat <- diag(v)[rep(1:v, times = v), ]

  rep_vec <- rep(1:r, each = b/r * k)
  Rmat<- matrix(0, n, r)

  for(i in 1:n){
    Rmat[i, rep_vec[i]] <- 1
  }

  X <- cbind(Tmat,mu,Bmat,Rmat)

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

  toss_l <- sum(yminusbar_sq_l)
  toss_u <- sum(yminusbar_sq_u)
  toss_d <- sum(yminusbar_sq_d)
  toss_i <- abs(toss_l-toss_u)/2

  Z <- cbind(mu, Bmat,Rmat)

  M <- diag(nrow(X)) - Z %*% ginv(t(Z) %*% Z) %*% t(Z)

  Ct_adj <- t(Tmat) %*% M %*% Tmat

  bcdy_t_ad_d <- t(Tmat)%*%cbind(mu,Bmat,Rmat)%*%ginv(t(cbind(mu,Bmat,Rmat))%*%cbind(mu,Bmat,Rmat))%*%t(cbind(mu,Bmat,Rmat))%*%dy

  bcdy_t_ad_i <- t(Tmat)%*%cbind(mu,Bmat,Rmat)%*%ginv(t(cbind(mu,Bmat,Rmat))%*%cbind(mu,Bmat,Rmat))%*%t(cbind(mu,Bmat,Rmat))%*%iy

  Qy_ad_l <- t(Tmat)%*%dy-bcdy_t_ad_d
  Qy_ad_u <- Qy_ad_l[(1:v),]+t(Tmat)%*%iy-bcdy_t_ad_i

  Qy_t_ad_d <- (Qy_ad_l+Qy_ad_u)/2
  Qy_t_ad_i <- abs(Qy_ad_l-Qy_ad_u)/2

  ay_t_ad_l <- ginv(Ct_adj) %*% Qy_t_ad_d
  ay_t_ad_u <- ay_t_ad_l+Qy_t_ad_i/v

  ay_t_ad_d <- (ay_t_ad_l+ay_t_ad_u)/2
  ay_t_ad_i <- abs(ay_t_ad_l-ay_t_ad_u)/2

  Z <- cbind(mu, Tmat,Rmat)

  M <- diag(nrow(X)) - Z %*% ginv(t(Z) %*% Z) %*% t(Z)

  Cb_adj <- t(Bmat) %*% M %*% Bmat

  bcdy_b_ad_d <- t(Bmat)%*%cbind(mu,Tmat,Rmat)%*%ginv(t(cbind(mu,Tmat,Rmat))%*%cbind(mu,Tmat,Rmat))%*%t(cbind(mu,Tmat,Rmat))%*%dy

  bcdy_b_ad_i <- t(Bmat)%*%cbind(mu,Tmat,Rmat)%*%ginv(t(cbind(mu,Tmat,Rmat))%*%cbind(mu,Tmat,Rmat))%*%t(cbind(mu,Tmat,Rmat))%*%iy

  Qy_ad_l <- t(Bmat)%*%dy-bcdy_b_ad_d
  Qy_ad_u <- Qy_ad_l[(1:v),]+t(Bmat)%*%iy-bcdy_b_ad_i

  Qy_b_ad_d <- (Qy_ad_l+Qy_ad_u)/2
  Qy_b_ad_i <- abs(Qy_ad_l-Qy_ad_u)/2

  ay_b_ad_l <- ginv(Cb_adj) %*% Qy_b_ad_d
  ay_b_ad_u <- ay_b_ad_l+Qy_b_ad_i/v

  ay_b_ad_d <- (ay_b_ad_l+ay_b_ad_u)/2
  ay_b_ad_i <- abs(ay_b_ad_l-ay_b_ad_u)/2

  Cr_unadj <- t(Rmat) %*% Rmat

  bcdy_r_ad_d <- t(Rmat)%*%cbind(mu,Tmat,Bmat)%*%ginv(t(cbind(mu,Tmat,Bmat))%*%cbind(mu,Tmat,Bmat))%*%t(cbind(mu,Tmat,Bmat))%*%dy

  bcdy_r_ad_i <- t(Rmat)%*%cbind(mu,Tmat,Bmat)%*%ginv(t(cbind(mu,Tmat,Bmat))%*%cbind(mu,Tmat,Bmat))%*%t(cbind(mu,Tmat,Bmat))%*%iy

  Qy_ad_l <- t(Rmat)%*%dy-bcdy_r_ad_d
  Qy_ad_u <- Qy_ad_l[(1:v),]+t(Rmat)%*%iy-bcdy_r_ad_i

  Qy_r_ad_d <- (Qy_ad_l[1:r,]+Qy_ad_u[1:r,])/2
  Qy_r_ad_i <- abs(Qy_ad_l[1:r,]-Qy_ad_u[1:r,])/2

  ay_r_ad_l <- ginv(Cr_unadj) %*% Qy_r_ad_d
  ay_r_ad_u <- ay_r_ad_l+Qy_r_ad_i/v

  ay_r_ad_d <- (ay_r_ad_l+ay_r_ad_u)/2
  ay_r_ad_i <- abs(ay_r_ad_l-ay_r_ad_u)/2

  ss_t_ad_l <- sum(Qy_t_ad_d*ay_t_ad_d)
  ss_t_ad_u <- ss_t_ad_l+sum(Qy_t_ad_d*ay_t_ad_i+Qy_t_ad_i*ay_t_ad_d+Qy_t_ad_i*ay_t_ad_i)

  ss_t_ad_d <- (ss_t_ad_l+ss_t_ad_u)/2
  ss_t_ad_i <- abs(ss_t_ad_u-ss_t_ad_l)/2

  ss_b_ad_l <- sum(Qy_b_ad_d*ay_b_ad_d)
  ss_b_ad_u <- ss_b_ad_l+sum(Qy_b_ad_d*ay_b_ad_i+Qy_b_ad_i*ay_b_ad_d+Qy_b_ad_i*ay_b_ad_i)

  ss_b_ad_d <- (ss_b_ad_l+ss_b_ad_u)/2
  ss_b_ad_i <- abs(ss_b_ad_u-ss_b_ad_l)/2

  ss_r_ad_l <- sum(Qy_r_ad_d*ay_r_ad_d)
  ss_r_ad_u <- ss_r_ad_l+sum(Qy_r_ad_d*ay_r_ad_i+Qy_r_ad_i*ay_r_ad_d+Qy_r_ad_i*ay_r_ad_i)

  ss_r_ad_d <- (ss_r_ad_l+ss_r_ad_u)/2
  ss_r_ad_i <- abs(ss_r_ad_u-ss_r_ad_l)/2

  mdss_all_ad_l <- ss_t_ad_d+ss_b_ad_d+ss_r_ad_d
  mdss_all_ad_u <- mdss_all_ad_l+ss_t_ad_i+ss_b_ad_i+ss_r_ad_i

  mdss_all_ad_d <- (mdss_all_ad_u+mdss_all_ad_l)/2
  mdss_all_ad_i <- abs(mdss_all_ad_u-mdss_all_ad_l)/2

  ess_all_ad_l <- toss_d-mdss_all_ad_d
  ess_all_ad_u <- ess_all_ad_l+toss_i-mdss_all_ad_i

  ess_all_ad_d <- (ess_all_ad_l+ess_all_ad_u)/2
  ess_all_ad_i <- abs(ess_all_ad_u-ess_all_ad_l)/2

  SV <- c("Treatment","Row","column","Error","Total")

  DF <- c(v-1,v-1,v-1,(v-1)*(v-2),v*v-1)

  SSN <- list(
    c(ss_t_ad_l,ss_t_ad_u),
    c(ss_r_ad_l,ss_r_ad_u),
    c(ss_b_ad_l,ss_b_ad_u),
    c(ess_all_ad_l,ess_all_ad_u),
    c(toss_l,toss_u)
  )

  MSSN <- list(
    c(ss_t_ad_l/(v-1),ss_t_ad_u/(v-1)),
    c(ss_r_ad_l/(v-1),ss_r_ad_u/(v-1)),
    c(ss_b_ad_l/(v-1),ss_b_ad_u/(v-1)),
    c(ess_all_ad_l/((v-1)*(v-2)),ess_all_ad_u/((v-1)*(v-2))),
    c(NA,NA)
  )

  F_t_ad_L <- (ss_t_ad_d/(v-1))/(ess_all_ad_d/((v-1)*(v-2)))
  F_t_ad_U <- F_t_ad_L+(((ss_t_ad_d/(v-1))*(ess_all_ad_i/(((v-1)*(v-2)))))+((ss_t_ad_i/(v-1))*(ess_all_ad_d/(((v-1)*(v-2)))))+((ss_t_ad_i/(v-1))*(ess_all_ad_i/(((v-1)*(v-2))))))/((ess_all_ad_d/(((v-1)*(v-2))))*(ess_all_ad_d/((v-1)*(v-2))))
  F_b_ad_L <- (ss_b_ad_d/(b-1))/(ess_all_ad_d/((v-1)*(v-2)))
  F_b_ad_U <- F_b_ad_L+((((ss_b_ad_d/(b-1))*(ess_all_ad_i/((v-1)*(v-2))))+((ss_b_ad_i/(b-1))*(ess_all_ad_d/((v-1)*(v-2))))+((ss_b_ad_i/(b-1))*(ess_all_ad_i/((v-1)*(v-2)))))/((ess_all_ad_d/((v-1)*(v-2)))*(ess_all_ad_d/((v-1)*(v-2)))))
  F_r_ad_L <- (ss_r_ad_d/(v-1))/(ess_all_ad_d/((v-1)*(v-2)))
  F_r_ad_U <- F_r_ad_L+((((ss_r_ad_d/(v-1))*(ess_all_ad_i/((v-1)*(v-2))))+((ss_r_ad_i/(v-1))*(ess_all_ad_d/((v-1)*(v-2))))+((ss_r_ad_i/(v-1))*(ess_all_ad_i/((v-1)*(v-2)))))/((ess_all_ad_d/((v-1)*(v-2)))*(ess_all_ad_d/((v-1)*(v-2)))))

  F_critical_5 <- qf(1-alpha,v-1,(v-1)*(v-2))
  F_critical_1 <- qf(1-alpha/5,v-1,(v-1)*(v-2))
  F_critical_001 <- qf(1-alpha/50,v-1,(v-1)*(v-2))

  significance <- c(

    ifelse(F_t_ad_L > F_critical_001 & F_t_ad_U > F_critical_001,"***",
           ifelse(F_t_ad_L > F_critical_1 & F_t_ad_U > F_critical_1,"**",
                  ifelse(F_t_ad_L > F_critical_5 & F_t_ad_U > F_critical_5,"*",
                         ifelse(F_t_ad_L < F_critical_5 & F_t_ad_U < F_critical_5,"NS","ID")))),

    ifelse(F_r_ad_L > F_critical_001 & F_r_ad_U > F_critical_001,"***",
           ifelse(F_r_ad_L > F_critical_1 & F_r_ad_U > F_critical_1,"**",
                  ifelse(F_r_ad_L > F_critical_5 & F_r_ad_U > F_critical_5,"*",
                         ifelse(F_r_ad_L < F_critical_5 & F_r_ad_U < F_critical_5,"NS","ID")))),

    ifelse(F_b_ad_L > F_critical_001 & F_b_ad_U > F_critical_001,"***",
           ifelse(F_b_ad_L > F_critical_1 & F_b_ad_U > F_critical_1,"**",
                  ifelse(F_b_ad_L > F_critical_5 & F_b_ad_U > F_critical_5,"*",
                         ifelse(F_b_ad_L < F_critical_5 & F_b_ad_U < F_critical_5,"NS","ID")))),

    "",
    ""
  )

  FN <- list(
    c(F_t_ad_L,F_t_ad_U),
    c(F_r_ad_L,F_r_ad_U),
    c(F_b_ad_L,F_b_ad_U),
    c(NA,NA),
    c(NA,NA)
  )

  format_interval <- function(x){
    if(any(is.na(x))) return("")
    paste0("[", round(x[1],2), ", ", round(x[2],2), "]")
  }

  SSN_fmt <- sapply(SSN, format_interval)
  MSSN_fmt <- sapply(MSSN, format_interval)
  FN_fmt <- sapply(FN, format_interval)

  nanova_table <- data.frame(
    SV = SV,
    DF = DF,
    SSN = SSN_fmt,
    MSSN = MSSN_fmt,
    FN = FN_fmt,
    Significance = significance,
    stringsAsFactors = FALSE
  )

  trt_means_l <- tapply(as.vector(Lower_y),as.vector(design),mean)
  trt_means_u <- tapply(as.vector(Upper_y),as.vector(design),mean)

  MSE_l <- ess_all_ad_l/((v-1)*(v-2))
  MSE_u <- ess_all_ad_u/((v-1)*(v-2))

  LSD_l <- qt(1-alpha/2,((v-1)*(v-2)))*sqrt((2*MSE_l)/r)
  LSD_u <- qt(1-alpha/2,((v-1)*(v-2)))*sqrt((2*MSE_u)/r)

  comparison <- data.frame()

  for(i in 1:(v-1)){
    for(j in (i+1):v){

      diff_l <- abs(trt_means_l[i]-trt_means_l[j])
      diff_u <- abs(trt_means_u[i]-trt_means_u[j])

      result <- ifelse(diff_l > LSD_u,"S",
                       ifelse(diff_u > LSD_l,"S","NS"))

      comparison <- rbind(comparison,
                          data.frame(
                            Treatment1=i,
                            Treatment2=j,
                            Difference=paste0("[",round(diff_l,2),", ",round(diff_u,2),"]"),
                            LSD=paste0("[",round(LSD_l,2),", ",round(LSD_u,2),"]"),
                            Result=result
                          ))
    }
  }

  if(verbose){

    message("\nNeutrosophic Analysis of Variance Table")
    print(nanova_table, row.names = FALSE)

    if(significance[1] != "NS"){

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
    message("NS  : Non Significant")
    message("ID  : Indeterminate")
  }
    invisible(list(
      nanova_table = nanova_table,
      comparison = if(significance[1] != "NS") comparison else NULL,
      LSD = c(Lower = LSD_l, Upper = LSD_u)
    ))
  }

