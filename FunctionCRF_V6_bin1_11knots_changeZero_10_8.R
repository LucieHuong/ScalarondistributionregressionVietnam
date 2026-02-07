#------------------------Today 22 SEPT 2022
#====18 basis function and order = 4

#===============================================================================
#===================0.1============================================================
#-----------------Trapezoidal formula for numerical integration

trapezoidF <- function(grid_value, f)
{
  #---code of a function to approximate integral, modify from trapzc()
  # The function will calculate with give a same interval, i.e length(grid_value) == 1
  #or different intervals, i.e length(grid_value) > 1, for the approximate function.
  # grid_value: The length of evaluated interval or the list of points.
  # f: grid evaluation of densityf
  # In case length(grid_value) > 1, f if give at each grid_value point
  
  if (length(grid_value) == 1) {
    #h: The same interval length for f.
    h = grid_value
    n <- length(f)
    fxdx <- rep(0, n-1)
    for (i in 1 : (n-1))
    {
      fxdx[i] <- (h / 2) * (f[i] + f[i+1])
    }
    fxdxSUM  <- sum(fxdx)
    return(fxdxSUM )
  } else {
    n <- length(f)
    fxdx <- rep(0, n-1)
    for (i in 1 : (n-1))
    {
      h <- grid_value[i+1] - grid_value[i]
      fxdx[i] <- (h / 2) * (f[i] + f[i+1])
    }
    
    fxdxSUM  <- sum(fxdx)
    return(fxdxSUM )
  }
  
}


#===============================================================================
#===============================================================================
#============ Function mapping from B^2 into L^2

fcenLRF <- function (z, density) 
{
  
  #---check if Inf value of log(density)
  
  InfN <- which(is.infinite(log(density)))
  if (length(InfN) > 0) 
    {
    density <- density[-(InfN)] # exclude the 0 density
    z <- z[-(InfN)]
  }
  out = log(density) - trapezoidF (z, log(density))/(max(z) -   min(z))
  outALL <- data.frame(gridpoints= z, estimateclrf = out)
  #gridpoints  is grid values to evaluate the function
  # estimateclrf: Estimated clr(f), here in L2
  return( outALL)
}

#===============================================================================
#===============================================================================
#============ Function mapping from #L^2 to B^2


fcenLRinvF <- function (z_stepvalue, fcenLRvalue, k = 1) 
{
  #----Modify function fcenLRin in robCompositions package
  # compared to   fcenLRinv(), no need the 1st argument
  
  out = (exp(fcenLRvalue)/trapezoidF(z_stepvalue, exp(fcenLRvalue))) * k
  
  outALL <- data.frame(knots = z_stepvalue, estimate = out)
  return( outALL)
  
}


#=================Function to create a ZB basis based on the given knots and order

ZBsplineDesign <- function(knots_beta, order = 4)
{
  k = order # Order of B-splines
  r_beta = length(knots_beta)
  Total_knots_beta = 2 * (k - 1) + r_beta 
  # Count total of knots_beta, i.e lambda_beta: (r_beta-2) inside knots_beta, 2 border knots_beta and (2*(k - 1)) extra knots_beta.
  
  #-------------- The following codes create knots_beta lambda_beta, including inside, border and extra knots_beta
  lambda_beta = c()
  for (i in 1:(Total_knots_beta)) {
    if (i <= k - 1) {
      lambda_beta[i] = knots_beta[1] # border left knots_beta and extra left knots_beta
    }
    if ((i > k - 1) && (i <= r_beta + k - 1)) {
      lambda_beta[i] = knots_beta[i - (k - 1)] # inside knots_beta
    }
    if (i > r_beta + k - 1) {
      lambda_beta[i] = knots_beta[r_beta] # border right knots_beta and extra right knots_beta
    }
  }
  
  
  Partition = seq(min(lambda_beta), max(lambda_beta), length = 1401) 
  lambda_beta_index = c(0:(r_beta - 1))
  g_beta = lambda_beta_index[length(lambda_beta_index) - 1] # g: the position of the last inside knot lambda_beta.
  N_Beta = g_beta + (k - 1) + 1 # k is order of the spline 
  
  #-------------- The following codes create B-spline basic BB, which are used later to draw ZB-splines based on lambda_beta knots_beta.
  
  BB_beta  = array(0, c(length(Partition), N_Beta)) 
  
  l = c() # Store cumulative knots_beta lambda_beta
  
  for (i in (1:N_Beta)) {
    for (j in 1:(k + 1)) {
      l[j] = lambda_beta[i + j - 1]
    }
    BB_beta [, i] = splineDesign(l, Partition, k, outer.ok = TRUE)
    # B-spline basics, based on knots_beta l, evaluation points Partition and k order
  }
  
  
  #=====Construct matrix D and K
  difference_beta = lambda_beta[(1 + k):(r_beta + 2 * (k - 1))] - lambda_beta[(1:(r_beta + k - 2))]
  D_beta = (k) * diag(1/difference_beta) 
  K_beta = array(0, c(N_Beta, N_Beta - 1))
  K_beta[1, 1] = 1
  K_beta[N_Beta, N_Beta - 1] = -1
  for (j in (2:(N_Beta - 1))) {
    K_beta[j, j - 1] = (-1)
    K_beta[j, j] = 1
  }
  
  Zbasis_finner_Beta = BB_beta %*% D_beta %*% K_beta # Dimension 1401x(N_Beta-1), ZB-Splines
  
  return(list(Partition = Partition, #values of the Objective function at grid points
              Zbasis_finner= Zbasis_finner_Beta, # Coefficients of ZB-basis
              knots = knots_beta#, D = D_beta, K = K_beta
              ))
}

#===============================================================================
#===============================================================================
#============================compositional Spline function


compositionalSplineF <- function (t, clrf, knots, 
                                  knots_beta0,
                                  w, order, der,
                                  alpha, spline.plot = FALSE, 
                                  basis.plot = FALSE,
                                  namemain = NULL, tw = NULL) 
{
  # t: class midpoints, i.e 
  # clrf  clr transformed values at class midpoints, i.e., fcenLR(f(t))
  # knots	 sequence of knots for histogram of X
  # knots_beta0 of the beta curve solutions
  # w	 weights
  # order	 order of the spline (i.e., degree + 1)
  # der	  lth derivation, for penalized term in the minimize functions of smoothing splines, usually der = 2.
  # alpha	  smoothing parameter
  # spline.plot	 if TRUE, the resulting spline is plotted
  # basis.plot	  if TRUE, the ZB-spline basis system is plotted
  
  if (alpha <= 0 | alpha > 1) 
    stop("parameter alpha must be from interval (0,1]")
  k = order # Order of B-splines
  # k: (k-1) is the degree of the B- spline function, a cubic spline has k = 4.
  # Note: in "Compositional splines for representation of density functions, 
  #a cubic spline has k = 3.
  r = length(knots)
  Total_knots = 2 * (k - 1) + r 
  # Count total of knots, i.e lambda: (r-2) inside knots, 2 border knots and (2*(k - 1)) extra knots.
  
  #-------------- The following codes create knots lambda, including inside, border and extra knots
  lambda = c()
  for (i in 1:(Total_knots)) {
    if (i <= k - 1) {
      lambda[i] = knots[1] # border left knots and extra left knots
    }
    if ((i > k - 1) && (i <= r + k - 1)) {
      lambda[i] = knots[i - (k - 1)] # inside knots
    }
    if (i > r + k - 1) {
      lambda[i] = knots[r] # border right knots and extra right knots
    }
  }
  
  #-------------- The following codes create collocation matrix B
  B = splineDesign(lambda, t, ord = k, outer.ok = TRUE) 
  # lambda:a numeric vector of knot positions;
  # t: a numeric vector of values at which to evaluate the B-spline functions
  # k: (k-1) is the order of the spline function, a cubic spline has k = 4.
  # B is the Collocation matrix
  # B has dimension length(t) x(length(lambda)-ord), or length(t) x( 2 * (k - 1) + r - k)
  # or length(t) x( k  + r - 2), or  length(t) x( k  + g), g = r-2
  
  W = diag(w) # weights of data points
  
  Partition = seq(min(lambda), max(lambda), length = 1401) 
  # partite interval [a,b]  = [min(lambda), max(lambda] to prepare for the matrix BB
  lambda_index = c(0:(r - 1))
  g = lambda_index[length(lambda_index) - 1] # g: the position of the last inside knot lambda.
  N_B = g + (k - 1) + 1 # k is order of the spline 
  # N_B is also the number of B-splines basic, by splineDesign() function, given (g = r-2) insides knots, 2 border knots and order (k-1)
  
  #-------------- The following codes create B-spline basic BB, which are used later to draw ZB-splines based on lambda knots.

  BB = array(0, c(length(Partition), N_B)) 
  # BB matrix has dimensition length(Partition) row X N_B columns, or 1401xN
  # length(Partition) = 1401 as above
  
  l = c() # Store cumulative knots lambda
  
  for (i in (1:N_B)) {
    for (j in 1:(k + 1)) {
      l[j] = lambda[i + j - 1]
    }
    BB[, i] = splineDesign(l, Partition, k, outer.ok = TRUE)
    # B-spline basics, based on knots l, evaluation points Partition and k order
  }
  
  if (length(t) <= N_B) 
    stop("length(t) must be higher then Dimension(space of splines)")
  if (qr(B)$rank != N_B) 
    stop("Collocation matrix does not have full column rank.")
  
  #-------------- The following codes create matrix S, one of main components in the optimal function.
  
  S = array(0)
  if (der == 0) {
    S = diag(1, c(N_B, N_B))
  }
  
  if (der >= (k - 1)) {
    stop("Error. The lth derivation is not from the set {1,2,...,order-2}.")
  } else {
    S_pom = diag(1, N_B, N_B)
    # S is upper triangular matrix,  N = g + (k - 1) + 1 = g+k
    for (j in 1:der) {
      #  der	  lth derivation, for penalized term in the minimize functions of smoothing splines, usually der = 2.
      D_mat = array(0)
      Diff_Lambda = lambda[(1 + k):(N_B + k - j)] - lambda[(1 +  j):(N_B)] # denominator for D_mat, different of two lambda knots.
      D_mat = (k - j) * diag(1/Diff_Lambda) # D_mat is a diagonal matrix
      # given j, D_mat has dimension (g+k-j)x(g+k-j) = (N-j)x(N-j)
      L_mat = array(0, c(N_B - j, N_B - j + 1)) # L_mat has dimension (N_B - j)x(N_B - j + 1)
      # Matrix L, dimension 
      for (J in (1:(N_B - j))) {
        L_mat[J, J] = (-1)
        L_mat[J, J + 1] = 1
      }
      S_pom = D_mat %*% L_mat %*% S_pom 
    }
    S = S_pom  # Dimension (g+k-der)x(g+k) or (N-der)xN
  }
  
  #-------------- The following codes create matrix M (based on BBB, a B-splines to make scalar product in L^2[min(lambda), max(lambda]))
  #---------------It also creates matrix D, K and then G: some components in the optimal function.
  
  kk = k - der # order of B-splines in the scalar product in L^2, i.e matrix M.
  Total_knots = 2 * (kk - 1) + r 
  # Re-Count total of knots, i.e Lambda: (r-2) inside knots, 2 border knots and (2*(kk - 1)) extra knots.
  # Here, knots Lambda with L capital, above knots, lambda!!
  Lambda = c()
  for (i in 1:Total_knots) {
    if (i <= (kk - 1)) {
      Lambda[i] = knots[1] # border left knots and extra left knots
    }
    if ((i > kk - 1) && (i <= r + kk - 1)) {
      Lambda[i] = knots[i - (kk - 1)]  # Inside knots
    }
    if (i > (r + (kk - 1))) {
      Lambda[i] = knots[r] # border right knots and extra right knots
    }
  }
   
  Partition = seq(min(Lambda), max(Lambda), length = 10000)  # Partition interval [a,b]  = [min(Lambda), max(Lambda] to prepare for the matrix BBB
  Lambda_index = c(0:(r - 1))
  G = Lambda_index[length(Lambda_index) - 1] # the position of the last inside knot Lambda, G = r-2
  NN = G + (kk - 1) + 1 # kk is order of the spline 
  # NN is also the number of B-splines basic, by splineDesign() function, given (G = r-2) inside knots, 2 border knots and order (kk-1)
 
   BBB = splineDesign(Lambda, Partition, kk, outer.ok = TRUE)
  
   step = diff(Partition[1:2]) # Partition is an equidistant list. Step for trapezoidF() for matrix M
  M = array(0, c(NN, NN))
  for (i in 1:NN) {
    for (j in 1:NN) {
      nenulove = c()
      soucin = BBB[, i] * BBB[, j]  # scalar product of 2 B-splines in L^2
      for (m in 1:length(Partition)) {
        if (soucin[m] != 0) {
          nenulove[m] = soucin[m] # scalar product of 2 B-splines for m_ij
        }
      }
      M[i, j] = trapezoidF(step, soucin)# [, 2] # trapzc in
      # M has dimension NNxNN, NN = G + (kk - 1) + 1 # kk is order of the spline, 
      #G inside knots. 
      # In addition, G = r-2 and kk = k - der
      # NN = r-2 + k- der -1+1 = g+k-der, g = r-2,
      #  N = g + (k - 1) + 1 then M has dimension (N-der)x(N-der)
    }
  }
  
  difference = lambda[(1 + k):(r + 2 * (k - 1))] - lambda[(1:(r + k - 2))]
  D = (k) * diag(1/difference) # Matrix D, dimension N_BxN_B, note that g = r-2
  K = array(0, c(N_B, N_B - 1)) # Dimension N_Bx(N_B-1), i.e N_B = g + degree + 1
  K[1, 1] = 1
  K[N_B, N_B - 1] = -1
  for (j in (2:(N_B - 1))) {
    K[j, j - 1] = (-1)
    K[j, j] = 1
  }
  U = D %*% K # Dimension N_Bx(N_B-1), i.e N_B = g + (k - 1) + 1 = g+k, g = r - 2
  GG = t(U) %*% ((1 - alpha) * t(S) %*% M %*% S + alpha * t(B) %*%   W %*% B) %*% U 
  # GG has dimension (N_B-1)xN_B and inside the (...) has dimension N_BxN_B # check by Huong on 8Feb 2022!!
  # In detail:
  #B has dimension: length(t) x( k  + g), or length(t)xN_B
  # U has dimension N_Bx(N_B-1), S has (N_B-der)xN_B size; M has dimension (N_B-der)x(N_B-der)
  # W has dimension length(t)xlength(t)
  gg = alpha * t(K) %*% t(D) %*% t(B) %*% W %*% clrf
  # D has dimension NxN, K has dimension Nx(N-1), clrf has dimension length(t)x1
  # gg has dimension (N-1)x1
  
  #-------------The following codes show the optimal z
  
  z = solve(GG) %*% gg # length = (k-1) + g = N_B-1  # Optimal z, which is the solution of minimizing problem
  #  eq (24) in Machalova et al. 2021
  
  Zbasis = B %*% D %*% K # Dimension length(t)x(N_B-1)
  # ZB spline, evaluated at the midpoints for optimization (4NOV) 
  
  Zbasis_finner = BB %*% D %*% K # Dimension 1401x(N_B-1), ZB-Splines
  # BB has dimension: 1401xN_B, length(Partition) = 1401
  #Zbasis_finner is a ZB-spline basis function, evaluated at the Partition for plotting (4NOV) 
  #ZB spline and ZB spline_finner are the SAME, just different in term of evaluation points
 
   #------Calculate matrix M_finner
  
  M_finner = array(0, c(N_B, N_B))
  Partition = seq(min(lambda), max(lambda), length = 1401) 
  
  for (i in 1:N_B) {
    for (j in 1:N_B) {
      nenulove = c()
      soucin = BB[, i] * BB[, j]  # scalar product of 2 B-splines in L^2
      for (m in 1:length(Partition)) {
        if (soucin[m] != 0) {
          nenulove[m] = soucin[m] # scalar product of 2 B-splines for m_ij
        }
      }
      M_finner[i, j] = trapezoidF(step, soucin)# [, 2] # trapzc in
      # dimension N_BxN_B, N_B = g + (k - 1) + 1 = g+k
    }
  }
  
  
  #-------------The following codes show some figures, draw ZB-splines 
  if (basis.plot == TRUE) {
    matplot(Partition, Zbasis_finner, type = "l", lty = 1, 
            las = 1, xlab = "Temperature", ylab = "Clr(density)", 
            col = rainbow(dim(Zbasis_finner)[2]), 
            main = "ZB-spline basis")
    abline(v = knots, col = "gray", lty = 2)
  }
  
  
  spline0 = Zbasis_finner %*% z # Smoothing splines
  if (spline.plot == TRUE) {
    matplot(Partition, spline0, type = "l", las = 1, 
            xlab = "Temperature", ylab = "Clr(density)", 
            col = "darkblue", lwd = 2, cex.lab = 2, 
            cex.axis = 2, 
            ylim = c(min(c(min(clrf), min(spline0))),
                     max(c(max(clrf), 
                           max(spline0)))))
    #, 
     #       main = paste("Compositional spline of ", tw, " at province ",   namemain))
    matpoints(t, clrf, pch = 8, col = "darkblue", cex = 1.3)
    abline(h = 0, col = "red", lty = 2, lwd = 1)
    abline(v = knots, col = "gray", lty = 2, lwd = 1)
  }
  
  clrf = as.matrix(clrf)
  Hmat = (B %*% D %*% K) %*% solve(GG) %*% (alpha * t(K) %*% 
                                              t(D) %*% t(B) %*% W)
  # hat matrix H, in CV formula 
  clrfhat = (B %*% D %*% K) %*% z # Zbasis%*% z: smoothing spline based on the optimal z
  
  
  #----------The following codes shows the CV and GCV to find optimal smoothing splines, i.e alpha
  reziduals = (clrf - clrfhat) # part of the numerator in GCV
  Hmat_diag = c()
  for (i in 1:length(clrf)) Hmat_diag[i] = Hmat[i, i]
  Hmat_diag_mean = (1/length(clrf)) * sum(Hmat_diag) #Trace of hat matrix H
  CV = (1/length(clrf)) * sum((reziduals/(rep(1, length(Hmat_diag)) - 
                                            Hmat_diag))^2)
  GCV = (1/length(clrf)) * (sum((reziduals)^2))/((1 - Hmat_diag_mean)^2)
  
  #---------The following code show the evaluated value of Objective Function J
  
  J = (1 - alpha) * t(z) %*% t(U) %*% t(S) %*% M %*% S %*% U %*% z + alpha * t(clrf - B %*% D %*% K %*% z) %*% W %*% 
    (clrf - B %*% D %*% K %*% z)
  
  #==========================================================================
  #-------------Add Sigma matrix: the inner products of all pair ZB spline basis
  #M_finner has dimension NxN 
  # D has dimension NxN, K has dimension N_Bx(N_B-1)
  
  knots_beta <- knots_beta0 #c(12.00, 21.64, 24.99, 27.03, 28.49, 
                 # 29.57, 30.40, 31.08, 31.71, 32.35, 33.12, 
                #  34.23, 40.00 )
  Zbasis_beta <- ZBsplineDesign(knots_beta, order = 4)
  
  
  #========================Get matrix M_Z (below equation (40), Talska et al. 2021 - Regression paper)
  #Dimension of ZB basis for histogram
  N_ZB_h <- dim(Zbasis_finner)[2] 
  #Dimension of ZB basis for beta curve
  N_ZB_beta <- dim(Zbasis_beta$Zbasis_finner)[2]
 # Zbasis_finner is a ZB basis for histogram
  ZB_beta <- Zbasis_beta$Zbasis_finner # ZB- spline basis for beta curve evaluated at Partition
  
  MZ = array(0, c(N_ZB_h, N_ZB_beta))
  
  for (i in 1:N_ZB_h) {
    for (j in 1:N_ZB_beta) {
      Temp_Z = c()
      Inner_Z = Zbasis_finner[, i] * ZB_beta [, j]  # scalar product of 2 ZB-splines in L^2
      for (m in 1:length(Partition)) {
        if (Inner_Z[m] != 0) {
          Temp_Z[m] = Inner_Z[m] # scalar product of 2 ZB-splines for m_ij
        }
      }
      MZ[i, j] = trapezoidF(step, Inner_Z)
    }
  } #
  
  Sigma =  MZ  # today 8NOV 2022: It is matrix M_Z in Taslka et al 2021, regression paper
  
  #----- C*Sigma
  C_Sigma = t(z)%*%Sigma # 1x(N_B-1)
  
  
  #-----------Output of the function
  return(list(J = J, #values of the Objective function at grid points
              ZB_coef = z, # Coefficients of ZB-basis
              CV = CV, # alpha by CV
              GCV = GCV, # alpha by GCV
              Partition = Partition, # Evaluation grid, 1401 points
              spline0 = spline0, # clr(f) curve
              Zbasis_finner = Zbasis_finner, # ZB-Basis
              Zbasis_finner_beta = Zbasis_beta$Zbasis_finner,
              Sigma = Sigma,  #the inner products of all pair ZB spline basis==R, equation below (40) in Talska (2021)
              C_Sigma =  C_Sigma, # clr(f) spline at grid points
              clrf = clrf, # clr(f) 
              t =t, #  Grid values, equal distance = 0.1
              knots=  knots, knots_beta = knots_beta,
              clrfhat = clrfhat))
}


#================================================================================
#=================================Applying compositional splines for climate data
#-------------In this V2 version, there is only one function for tmax and tmin-------------------------
#-------------Keep the same B-splines for tmax, similar for tmin.

#------------- Function for density: density by frequency, then input zero
#--------------and then add more step to smooth histogram

histF <- function(data.t, bin_t)
{
  #data.t: weather data, whether tmin or tmax
  # bin_t: ALL the left point of bin and the last point
  # Output of this function is input of fcenLRF()
  #--------histogram by frequency
  
  
  #--------------First, hisogram by frenquency
  n0 <- length(bin_t) - 1 # Take only the number of leftpoints
  histempt <- data.frame(density = 0, leftpoints = bin_t[1:n0])
  
  for (b0 in 2: length(histempt$density) )
  {
    i0 <- bin_t[b0-1]
    i1 <- bin_t[b0]
    
    b1 <- data.t[i0 <= data.t & data.t < i1]
    histempt$density[b0-1] <- length(b1)/length(data.t) #correct 27OCT 2022
  }
  
  #=====The last left point
  
  i0 <- bin_t[n0]
  
  b1 <- data.t[i0 <= data.t ]
  histempt$density[n0] <- length(b1)/length(data.t) #correct 27OCT 2022
  
  #-------Second, Zero replacement
  
  n_zero <- length( histempt $density[histempt $density==0 ] )
  histempt $density[histempt $density==0 ]  <- 10^(-8)
  
  histempt $density <- histempt $density*(1/(1+ n_zero*10^(-8)))
  
  #====third, adjustment for bin width
  
  histempt$h <- NA
  for (i in 1:(n0-1))
  {
    histempt$h[i] <- histempt$leftpoint[i+1] - histempt$leftpoint[i]
  }
  
  histempt$h[n0] <- max(bin_t) - histempt$leftpoint[n0]
  
  for (i in 1:n0)
  {
    histempt$H[i] <- histempt$density[i]/histempt$h[i]
  }
  
  
  return(list(leftpoints = histempt$leftpoints,  
              density = histempt$H, # density with bandwidth adjustment
              densityRaw = histempt$density,
              maxpoint = max(bin_t))# The maximum value of the right knots
  )
}


#=============Today 27 OCT 2022
refinehistF <- function(density0)
{
  #histF0 is the output from histF() function
  
  gridpoints <- seq(min(density0$leftpoints),
                    density0$maxpoint, length.out = 1401)
  
  point.bin <- c(density0$leftpoints, density0$maxpoint) # all the points: all leftpoints and the maximum points
  
  histempt2 <- data.frame(density = NA, 
                          gridpoints = gridpoints)# 400 gridpoints
  for (i in 1: length(histempt2$density) )
  {
    
    j0 <- findInterval(histempt2$gridpoints[i], point.bin) #the last gridpoints do not enter in any Interval
    histempt2$density[i] <- density0$density[j0] 
  }
  
  # for the last gridpoints
  histempt2$density[length(histempt2$density)] <- histempt2$density[length(histempt2$density)-1] 
  return(list(gridpoints = histempt2$gridpoints,  
              density = histempt2$density,
              densityRaw = density0$density))
  
}



funCS <- function(data.tweather, type.weather)
{
  # type.weather has two option: tmin or tmax
  if (type.weather == "tmax")
  {
    bins <- c(12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 
              24, 25, 26, 27, 28, 29, 30, 31, 
              32, 33, 34, 35, 36, 37, 38, 39, 40)
    knots0 <- c(12.00, 22.51, 25.90, 27.96, 29.38, 30.40, 31.21, 31.96, 32.79, 33.96, 40.00 )
  }
  if (type.weather == "tmin")
  {
    bins <- c( 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
               18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29) # g = 13 inside bins, 
    knots0 <- c(7.00, 14.84, 17.95, 19.79, 21.16, 22.41, 23.52, 24.43, 25.27, 26.25, 29.00 )
  }
  
  #--------histogram by frequency in percentage
  histF0 <-  histF(data.tweather,   bins)
  

  #========Transform from B2 to L2
  midpoint <- bins[1:(length(bins)-1)]+1
  clrhisttempt <- fcenLRF(midpoint, 
                          histF0$density ) 
  
  #---define t, clrf and w for compositionalSplineF() later
  
  clrf <-  clrhisttempt$estimateclrf
  t  <- clrhisttempt$gridpoints
  
  #bins defined by quantile of all tmax, 14JAN2022
  w <- rep( 1, length(t))# rep( 1/length(t), length(t))
  
  return(list(clrf = clrf, # Estimated clr value of f in L2
              midpoint = midpoint, # gridpoints, at which we evaluate clr(f)
              bins = bins, # the leftpoints to calculate histogram frequency
              w = w,
              knots = knots0,
               # ,
             # densityRaw = histF0$densityRaw, # estimated histogram frequency
              density = histF0$density
  ))
}


#======================================================================================
#--------------------The following code to find optimal smoothing parameter alpha, by GCV, given the weather data

Listalphaln0 <- data.frame(alphal_m =  seq(-10, 10, by = 0.2)) # M = 5 and mu = alphaln in this step.

Listalphaln0$alpha_ld <- exp(Listalphaln0$alphal_m) # here: lamnda=exp(mu)

Listalphaln0 <- Listalphaln0 %>% mutate (alpha = 1/(alpha_ld+1) )


Listalphaln0 <- Listalphaln0[Listalphaln0$alpha <= 1,  ]
Listalphaln0$CV <- NA
Listalphaln0$GCV <- NA

Funlambda <- function(data.tweather, type.weather)
{
  
  alphaALL <- Listalphaln0
  
  # type.weather has two option: tmin or tmax
  if (type.weather == "tmax")
  {
    don <- funCS(data.tweather, "tmax" )
  }
  
  if (type.weather == "tmin")
  {
    don <- funCS(data.tweather, "tmin" )
  }

  for(i in 1: dim(Listalphaln0)[1] )
  {
    test1 <- NULL
    test1 <- compositionalSplineF(t = don$midpoint, clrf = don$clrf, 
                                  knots = don$knots, 
                                  knots_beta0 = don$knots, 
                                  w = don$w, order, der, 
                                  alpha = alphaALL$alpha[i], 
                                  spline.plot = FALSE, 
                                  basis.plot = FALSE)
    alphaALL[i, "CV"] <- test1$CV
    alphaALL[i, "GCV"] <- test1$GCV
  }
  
  alphaopti <-  alphaALL[which.min(alphaALL[, "GCV"]), "alpha"]
  return(alphaopti)
}


#==========today 18 May 2022

#===========================function: func_clrinv_all
# The aim of this function is to merge several steps in estimation: 
# By given a province and typeweathername, we apply funCS() function and 
# compositionalSplineF() and then back to B^2
func_clrinv_all <- function(knots_beta0, yrname, provincename, 
                            typeweathername)
{
  # type.weather has two option: tmin or tmax
  don_weather0 <- NULL
  yrname <- as.character(yrname)
  if (typeweathername == "tmax")
  {
    weather0 <- LtmaxALL[[yrname]][[provincename]]$tmaxALL 
    don_weather0 <- funCS(weather0, "tmax")
  }
  
  if (typeweathername == "tmin")
  {
    weather0 <- LtminALL[[yrname]][[provincename]]$tminALL
    don_weather0 <- funCS(weather0, "tmin")
  }
  
  
  sfun_weather0 <- stepfun(x = don_weather0$bins[-c(1, length(don_weather0$bins))],
                           y = don_weather0$density)
  
  order = 4
  der = 2
  
  if (typeweathername == "tmin")
  {
    alpha0 <- AlphaALL_tmin[AlphaALL_tmin$PRO == provincename & AlphaALL_tmin$year == yrname, "alpha"]
  }
  
  if (typeweathername == "tmax")
  {
    alpha0 <- AlphaALL_tmax[AlphaALL_tmax$PRO == provincename & AlphaALL_tmax$year == yrname, "alpha"]
  }
  
  
  CFR <-  compositionalSplineF(t = don_weather0$midpoint,
                               clrf = don_weather0$clrf, 
                               knots = don_weather0$knots,
                               knots_beta0 = knots_beta0,
                               w = don_weather0$w, order,
                               der, alpha = alpha0, 
                               spline.plot = FALSE, 
                               basis.plot = FALSE, tw = "tmax",
                               namemain = paste (yrname, provincename))
  f.fcenLRinv <- fcenLRinvF(CFR$t, CFR$clrfhat) # Back B^2
  
  #-------------------Today 4 NOV 2022
  #====Refine the results of CFR for plotting
  
  CFR.hist.spline <- CFR$Zbasis_finner %*% CFR$ZB_coef
  f.fcenLRinv_Refine <- fcenLRinvF(CFR$Partition, CFR.hist.spline) 
  
  
  return(list(sfun_weather0 = sfun_weather0, 
              knots = CFR$knots,
              ZB_coef = CFR$ ZB_coef,
              density.estimate = f.fcenLRinv_Refine$estimate,
              clrfhat = CFR.hist.spline, t = CFR$t, C_Sigma = CFR$C_Sigma,
              Partition = CFR$Partition))
}



#===============Today 3 NOV 2022
#=================Calculate histogram for the classical histogram 

#------------- Function for density: density by frequency, then input zero

histF.ilr <- function(data.t, type.weather)
{
  
  #data.t: weather data, whether tmin or tmax
  # interval.t: the max and minimun of weather data
  # type.weather has two option: tmin or tmax
  # Output of this function is to apply ILR to regression
  #--------histogram by frequency
  
  # type.weather has two option: tmin or tmax
  # bin_t is the left point of bin
  
  if (type.weather == "tmax")
  {
    bin_t <- c(12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 
               24, 25, 26, 27, 28, 29, 30, 31, 
               32, 33, 34, 35, 36, 37, 38, 39, 40)
    
  }
  
  if (type.weather == "tmin")
  {
    bin_t <- c( 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
               18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29)
  }
  
  #--------------First, hisogram by frenquency
  nb0 <- length(bin_t)
  histempt <- data.frame(density = 0, leftpoints = bin_t[-nb0])
  
  for (b0 in 2: length(histempt$density) )
  {
    i0 <- bin_t[b0-1]
    i1 <- bin_t[b0]
    
    b1 <- data.t[i0 <= data.t & data.t < i1]
    histempt$density[b0-1] <- length(b1)/length(data.t) #correct 22 August 2022
  }
  
  #=====The last left point
  
  i0 <- bin_t[(nb0-1)]
  
  b1 <- data.t[i0 <= data.t ]
  histempt$density[(nb0-1)] <- length(b1)/length(data.t) #correct 27 SEPT 2022
  
  
  #-------Second, Modify histog2 DEC 2021: replace zero for a small more
  
  n0 <- length( histempt $density[histempt $density==0 ] )
  histempt $density[histempt $density==0 ]  <- 10^(-8)
  
  histempt $density <- histempt $density*(1/(1+n0*10^(-8)))
  
  
  return(histempt)
}


#=============Today 17 NOV 2022
refinebetacurve <- function(hist0, bins, k) 
{
  #hist0 is one discrete of the smooth beta curve
  # bins is the list of bins
  # k is the number of sub-bin from each interval in bins
  # Today, 17NOV2022, we then discretize using the sub-bin then we do not need this funtions
  
  k0 <- diff(c(bins[1], bins[2]))/k # the bin width of sub-bin
  bins_rf <- seq(from = bins[1], to = bins[length(bins)], by = k0)  # List bins after refine
  
  histemptrf <- data.frame(density = NA, 
                           gridpoints = bins_rf ) 
  for (i in 1: length( histemptrf$density) )
  {
    
    j0 <- findInterval( histemptrf$gridpoints[i], bins) #the last gridpoints do not enter in any Interval
    histemptrf$density[i] <- hist0[j0] 
  }
  
  # for the last gridpoints
  
  histemptrf$density[length(histemptrf$density)] <-  histemptrf$density[length( histemptrf$density)-1] 
  return(list(gridpoints =  histemptrf$gridpoints,  
              density =  histemptrf$density))
  
}


#=============Today 17 NOV 2022
refinecurve <- function(curve0, bins, k)
{
  #curve0 is one curve which will be refined
  # bins is the list of bins
  # k is the number of sub-bin from each interval in bins
  
  k0 <- diff(c(bins[1], bins[2]))/k # the bin width of sub-bin
  bins_rf <- seq(from = bins[1], to = bins[length(bins)], by = k0)  # List bins after refine
  
  temptrf <- data.frame(valuerf= NA, 
                        gridpoints = bins_rf ) 
  for (i in 1: length(temptrf$valuerf) )
  {
    
    j0 <- findInterval(temptrf$gridpoints[i], bins) #the last gridpoints do not enter in any Interval
    temptrf$valuerf[i] <- curve0[j0] 
  }
  
  # for the last gridpoints
  temptrf$valuerf[length(temptrf$valuerf)] <- temptrf$valuerf[length(temptrf$valuerf)-1] 
  return(list(gridpoints = temptrf$gridpoints,  
              valuerf= temptrf$valuerf))
  
}

