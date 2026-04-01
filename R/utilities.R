################Data with NA#########################
completed_data <- function(x, y, defined_domain=c(0,100)) {
  
  if (length(defined_domain) != 2) {
    stop("Domain must have a form like c(starting point,ending point)")
  }
  start_domain <- defined_domain[1]
  end_domain <- defined_domain[2]
  
  # Perform linear interpolation using approx() function
  xnew <- seq(start_domain, end_domain, by = 1)  # New x values with steps of 1
  interpolated <- approx(na.omit(x), na.omit(y), xout = xnew)  # Interpolate y values
  
  # Return new y values
  return(interpolated$y)
}


###################Utilities Functions###########################

sigma_to_fwhm <- function(sigma){
  return(sigma*(2 * sqrt(2 * log(2))))
}

fwhm_to_sigma <- function(fwhm){
  return(fwhm / (2 * sqrt(2 * log(2))))
}



#unsmoothed guassian noise curve
noise_guassian_curve <- function(number_of_curves, continuum_size){ 
  data <- matrix(rnorm(number_of_curves*continuum_size),
                 nrow = number_of_curves, ncol = continuum_size)
  return(data)
  
}


#apply this scaling factor to any smoothed data to ensure it has unit variance after smoothing
# from _set_scale function in https://github.com/0todd0000/power1d/blob/master/src/power1d/random.py#L35

set_scale <- function(nNodes, SD) {
  # Define a small epsilon to prevent division by zero
  eps <- .Machine$double.eps
  
  # Step 1: Define a Gaussian kernel
  t <- seq(-0.5 * (nNodes - 1), 0.5 * (nNodes - 1), length.out = nNodes)
  gf <- exp(-(t^2) / (2 * SD^2 + eps))
  
  # Step 2: Normalize the Gaussian kernel
  gf <- gf / sum(gf)
  
  # Step 3: Calculate the expected variance of the smoothed data
  # Perform FFT and compute power spectrum
  AG <- fft(gf)
  Pag <- Mod(AG)^2  # Equivalent to AG * Conj(AG)
  
  # Calculate the autocovariance by inverse FFT
  COV <- Re(fft(Pag, inverse = TRUE)) / length(Pag)
  svar <- COV[1]  # Variance of the smoothed field
  
  # Step 4: Calculate the scaling factor
  SCALE <- sqrt(1.0 / svar)
  
  return(SCALE)
}


# Computes a 1-D Gaussian (or its derivative) kernel.
# sigma: standard deviation
# order: 0 for the standard Gaussian, positive integers for derivatives
# radius: half-width of the kernel; if not given, use truncate * sigma.
# The same function as scipy.ndimage$gaussian_filter1d(x, SD, mode='wrap')#$tolist()
gaussian_kernel1d <- function(sigma, order = 0, radius) {
  if (order < 0)
    stop("order must be non-negative")
  
  # Create grid: from -radius to radius.
  x <- seq(-radius, radius)
  
  # Compute the basic Gaussian (unnormalized)
  phi_x <- exp(-0.5 * (x / sigma)^2)
  phi_x <- phi_x / sum(phi_x)  # normalize so that sum(phi_x) == 1
  
  if (order == 0) {
    return(phi_x)
  } else {
    # For derivatives we need to multiply the Gaussian by a polynomial.
    # The SciPy code computes the polynomial coefficients via a matrix method.
    # Here we mimic that process.
    expo <- 0:order  # exponents 0,1,...,order
    q <- numeric(order + 1)
    q[1] <- 1  # q[0] = 1  (R indexing: first element corresponds to exponent 0)
    
    # Build the “differentiation” matrices D (superdiagonal) and P (subdiagonal)
    D <- matrix(0, nrow = order + 1, ncol = order + 1)
    for (i in 1:order) {
      D[i, i + 1] <- expo[i + 1]
    }
    P <- matrix(0, nrow = order + 1, ncol = order + 1)
    for (i in 2:(order + 1)) {
      P[i, i - 1] <- 1 / (-sigma^2)
    }
    Q_deriv <- D + P
    # Apply the operator repeatedly (order times)
    for (i in 1:order) {
      q <- Q_deriv %*% q
    }
    # For each x value, evaluate the polynomial:
    # Compute sum_{j=0}^{order} q[j] * x^j.
    # (outer(x, expo, `^`) builds a matrix whose (i,j) entry is x[i]^(expo[j]).)
    poly_val <- as.vector(outer(x, expo, `^`) %*% q)
    return(poly_val * phi_x)
  }
}


# Performs a 1-D correlation using periodic ("wrap") boundary conditions.
# input: a numeric vector.
# weights: the 1-D kernel (assumed to have odd length).
correlate1d_wrap <- function(input, weights) {
  n <- length(input)
  k <- length(weights)
  # Assume kernel size is odd so that it is symmetric around its center:
  half <- (k - 1) / 2
  result <- numeric(n)
  
  # For each element of the output, sum over the kernel with periodic indexing.
  for (i in seq_len(n)) {
    # Offsets: from -half to +half.
    offsets <- (-half):half
    # Compute wrapped indices: in R indices run 1..n.
    idx <- ((i - 1 + offsets) %% n) + 1
    result[i] <- sum(input[idx] * weights)
  }
  return(result)
}

# The main function: a 1-D Gaussian filter with periodic (wrap) boundary conditions.
#
# Arguments:
#  - input: numeric vector to filter.
#  - sigma: standard deviation of the Gaussian.
#  - order: the order of the derivative (default 0 means no derivative).
#  - truncate: how many sigmas to include in the kernel (ignored if radius is provided).
#  - radius: if provided, the half-length of the kernel; otherwise computed as round(truncate * sigma).
#
# Note: In the SciPy version the kernel is reversed before calling the correlate1d.
# We do the same here.
gaussian_filter1d <- function(input, fwhm, order = 0, truncate = 4.0, radius = NULL) {
  sigma <- fwhm_to_sigma(fwhm)
  if (!is.numeric(input)){
    stop("input must be numeric")
  }
  
  # Determine kernel half-width (radius)
  if (is.null(radius)) {
    radius <- round(truncate * sigma)
  } else if (radius < 0 || radius != as.integer(radius)) {
    stop("radius must be a nonnegative integer")
  }
  
  # Create the kernel.
  weights <- gaussian_kernel1d(sigma, order, radius)
  # In SciPy, the kernel is reversed (because of the use of correlation rather than convolution).
  weights <- rev(weights)
  
  # Apply the filter using periodic (wrap) boundary conditions.
  output <- t(apply(input, 1, function(x) correlate1d_wrap(x, weights)))
  return(output)
}


smoothed_gussian_curves <- function(data, mu, sig, fwhm) {
  
  # Step 1: Smooth each curve in the data
  smoothed_data <- gaussian_filter1d(data, fwhm)
  
  
  # Step 2: Normalize the smoothed data to have unit variance
  nNodes <- ncol(smoothed_data)
  SD <- fwhm_to_sigma(fwhm)
  scale_factor <- set_scale(nNodes, SD)
  
  # Step 3: Scale the smoothed data
  smoothed_data_scaled <- smoothed_data * scale_factor
  
  # Step 4: Transform to have mean = mu and standard deviation = sig
  # if mu and sig was a single number
  smoothed_data_final <- smoothed_data_scaled *
    matrix(sig, nrow = dim(smoothed_data_scaled)[1],
           ncol = dim(smoothed_data_scaled)[2], byrow = TRUE) +
    matrix(mu,  nrow = dim(smoothed_data_scaled)[1],
           ncol = dim(smoothed_data_scaled)[2], byrow = TRUE)
  
  return(t(smoothed_data_final))
}





Noise_generator <- function(Sample_size, Continuum_size, Noise_mu, Noise_sig, Noise_fwhm){
  
  noise1 <- noise_guassian_curve(number_of_curves = Sample_size,
                                 continuum_size = Continuum_size)
  noise1 <- smoothed_gussian_curves(noise1, Noise_mu, Noise_sig, Noise_fwhm)
  
  sd_noise <- apply(noise1, 1, sd)
  mean_noise <- apply(noise1, 1, mean)
  
  return(list(noise = noise1, SD = sd_noise, Mean = mean_noise))
}

# Generating data (Baseline+noise+signal) or (Baseline+noise) or (Baseline+signal)
data_generator <- function(data,signal,noise) {
  sample_size <- dim(noise)[2]
  
  if(missing(data) & missing(signal)){
    data_out <- noise
  } else if (missing(data)){
    signal_baseline <- matrix(rep(signal,time=sample_size),
                              ncol = sample_size)
    data_out <- signal_baseline + noise
  }
  else if (missing(signal)) {
    data_baseline <- matrix(rep(data,time=sample_size),
                            ncol = sample_size)
    data_out <- data_baseline + noise
  } else if (missing(noise)) {
    data_baseline <- matrix(rep(data,time=sample_size),
                            ncol = sample_size)
    signal_baseline <- matrix(rep(signal,time=sample_size),
                              ncol = sample_size)
    data_out <- data_baseline + signal_baseline
  } else {
    data_baseline <- matrix(rep(data,time=sample_size),
                            ncol = sample_size)
    signal_baseline <- matrix(rep(signal,time=sample_size),
                              ncol = sample_size)
    data_out <- data_baseline + signal_baseline + noise 
  }
  return(data_out)
}



gaussian_pulse <- function(center, fwhm, continuum_size) {
  
  sigma = fwhm_to_sigma((fwhm/100)*continuum_size)
  x_values = seq(0, continuum_size-1, by = 1)
  dens <- dnorm(x_values, mean = center, sd = sigma)
  return(list(density_val=dens, x_values=x_values))
}

amplitude_pulse <- function(data, amp){
  scaling_factor = amp / max(data)
  y_values = scaling_factor * data
  return(y_values)
}


square_pulse <- function(start_end_pulse, start_height,
                         pulse_height, continuum_size=101) {
  
  if (length(start_end_pulse) != 2 ||
      start_end_pulse[2]>(continuum_size-1) ||
      start_end_pulse[1] > start_end_pulse[2] ) {
    stop("start_end_pulse must have a form like c(start, end) and be within the continuum_size")
  }
  start = start_end_pulse[1]+1
  end = start_end_pulse[2]+1
  pulse = rep(start_height, continuum_size)
  pulse[start:end] = pulse_height
  return(pulse)
}

dense_data <- function(data1, data2, eval_numbers){
  fdObject_data1 = fda::Data2fd(argvals = seq(0,dim(data1)[1]-1,by=1), y = data1) #data must be a matrix with each column as a curve
  fdObject_data2 = fda::Data2fd(argvals = seq(0,dim(data2)[1]-1,by=1), y = data2) #data must be a matrix with each column as a curve
  rangeval1 <- fdObject_data1$basis$rangeval
  rangeval2 <- fdObject_data2$basis$rangeval
  if (sum(rangeval1 == rangeval2) != 2) {
    stop(
      "The range of values of {.arg data1} must be the same as the range of
        values of {.arg data2}."
    )
  }
  abscissa <- seq(rangeval1[1], rangeval1[2], length.out = eval_numbers)
  coeff1 <- t(fda::eval.fd(fdobj = fdObject_data1, evalarg = abscissa))
  coeff2 <- t(fda::eval.fd(fdobj = fdObject_data2, evalarg = abscissa))
  
  return(list(dense_data1 = coeff1, dense_data2 = coeff2))
}


col_diff <- function(data){
  return(data[,2]-data[,1])
}

##### Power Estimation Functions #####

check_true_positives_rate <- function(sig_matrix, D1) {
  apply(sig_matrix, 2, function(col) {
    sum(col[D1]) / length(D1)
  })
  
}

omnibus_power <- function(test_outputs, alpha = 0.05){
  # test_outputs: a matrix contains p-values/reject-notreject for each iteration (columns) and each point (rows).
  # alpha: significance level for determining if a p-value indicates a significant result.
  
  if(is.logical(test_outputs)){
    # If test_outputs is already a logical matrix (TRUE for significant, FALSE for not), we can directly calculate power.
    significant_iterations <- apply(test_outputs, 2, any)
  } else if(is.numeric(test_outputs)){
    # If test_outputs is a numeric matrix of p-values, we need to determine significance based on alpha.
    # Determine if any point is significant in each iteration
    significant_iterations <- apply(test_outputs < alpha, 2, any)
  } else {
    stop("test_outputs must be either a logical matrix or a numeric matrix of p-values.")
  }
  # Calculate power as the proportion of iterations with at least one significant node
  power_results <- mean(significant_iterations)
  return(power_results)
}

sensitivity <- function(test_outputs, D1, alpha = 0.05){
  # test_outputs: a matrix contains p-values/reject-notreject for each iteration (columns) and each point (rows).
  # D1: a numeric vector indicating which points are truly different.
  # alpha: significance level for determining if a p-value indicates a significant result.
  
  if(is.logical(test_outputs)){
    # If test_outputs is already a logical matrix (TRUE for significant, FALSE for not), we can directly calculate sensitivity.
    true_positives_ratio <- check_true_positives_rate(test_outputs, D1)
  } else if(is.numeric(test_outputs)){
    # If test_outputs is a numeric matrix of p-values, we need to determine significance based on alpha.
    significant_iterations <- test_outputs < alpha
    true_positives_ratio <- check_true_positives_rate(significant_iterations, D1)
  } else {
    stop("test_outputs must be either a logical matrix or a numeric matrix of p-values.")
  }
  sensitivity_results <- mean(true_positives_ratio)
  return(sensitivity_results)
  
}



ROI_power <- function(test_outputs, R, alpha = 0.05){
  # test_outputs: a matrix contains p-values/reject-notreject for each iteration (columns) and each point (rows).
  # R: a numeric vector indicating which region is of interest.
  # alpha: significance level for determining if a p-value indicates a significant result.
  
  if(is.logical(test_outputs)){
    # If test_outputs is already a logical matrix (TRUE for significant, FALSE for not), we can directly calculate sensitivity.
    ROI <- apply(test_outputs[R,], 2, any)
  } else if(is.numeric(test_outputs)){
    # If test_outputs is a numeric matrix of p-values, we need to determine significance based on alpha.
    significant_iterations <- test_outputs < alpha
    ROI <- apply(significant_iterations[R,], 2, any)
  } else {
    stop("test_outputs must be either a logical matrix or a numeric matrix of p-values.")
  }
  ROI_power <- mean(ROI)
  return(ROI_power)
}



SAE_function <- function(m1, m2, sd1, sd2, n1, n2){
  pooled_sd <- sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
  effect <- abs(m1 - m2) / pooled_sd
  return(effect)
}

compute_D1p <- function(m1, m2, sd1, sd2, n1, n2, taus) {
  
  if (any(taus <= 0)) {
    stop("taus must be a positive number(s)")
  }
  
  effect <- SAE_function(m1, m2, sd1, sd2, n1, n2)
  
  results <- list()
  
  for (i in seq_along(taus)) {
    tau <- taus[i]
    idx <- which(effect > tau)
    
    if (length(idx) > 0) {
      runs <- split(idx, cumsum(c(1, diff(idx) != 1))) #if the effect exceeds the 
      #threshold at indices 5, 6, 7, and 10, 11, it splits them into two distinct groups (runs).
      
      for (r in runs) {
        #append
        results[[length(results)+1]] <-
          data.frame(
            xmin = min(r),
            xmax = max(r),
            tau_value = tau
          )
      }
    }
  }
  
  return(bind_rows(results))
}


make_indices <- function(df) {
  split(df, df$tau_value) |>
    lapply(function(subdf) {
      unlist(mapply(seq, subdf$xmin, subdf$xmax, SIMPLIFY = FALSE))
    })
}

##### Inferential Methods #####

TWT <- function(data1,data2){
  TWT2=TWT2_new(data1,data2,mu=0,B=1000,paired=FALSE,dx=NULL)
  pvalue_adj_TWT2=TWT2$adjusted_pval
  return(pvalue_adj_TWT2)
}


IATSE <- function(data1,data2) {
  
  # t() because in Pvalue_calculator we t(data1) and t(data2)
  #since the other methods works with t() of data
  data1= t(data1)
  data2= t(data2)
  
  d <- dim(data1)[1]
  n_data <-dim(data1)[2]
  c_s <- create_curve_set(list(r=0:(d-1),obs=cbind(data1,data2)))
  
  group <- factor(rep(c(1,2),each=n_data))
  
  test <- frank.fanova(
    
    nsim = 1000,
    curve_set = c_s,
    groups = group,
    variances = "equal", #unequal variances can be used, but here we assume equal variances due to the fact that method is based on permutation
    test.equality = "mean",
    typeone = "fdr",
    algorithm = "IATSE",
    contrasts =TRUE,
    alpha = 0.05
  )
  test_result <- (test$obs>test$hi) | (test$obs<test$lo) 
  return(test_result)
  
}




centered_ranges <- function(percentages, domain = c(0, 100)) {
  midpoint <- mean(domain)
  range_widths <- (percentages / 100) * diff(domain)
  
  starts <- midpoint - range_widths / 2
  ends <- midpoint + range_widths / 2
  
  result <- data.frame(
    Percentage = percentages,
    Start = starts,
    End = ends
  )
  
  return(result)
}




##########Estimate parameters of the noise and data##########

estimate_fwhm <- function(R) {
  
  # Sum of squares across rows for each column
  ssq <- colSums(R^2)
  
  # Machine epsilon for numerical stability
  eps <- .Machine$double.eps
  
  n_rows <- nrow(R)
  n_cols <- ncol(R)
  
  # Initialize matrix to store approximated derivatives (gradient) along columns
  dx <- matrix(NA, nrow = n_rows, ncol = n_cols)
  
  # Compute gradient along columns for each row:
  for (i in 1:n_rows) {
    # Forward difference for the first column
    dx[i, 1] <- R[i, 2] - R[i, 1]
    # Backward difference for the last column
    dx[i, n_cols] <- R[i, n_cols] - R[i, n_cols - 1]
    # Central differences for interior columns (if available)
    if (n_cols > 2) {
      dx[i, 2:(n_cols - 1)] <- (R[i, 3:n_cols] - R[i, 1:(n_cols - 2)]) / 2
    }
  }
  
  # Sum the squared gradients across rows for each column
  v <- colSums(dx^2)
  
  # Normalize by the column-wise sum of squares (plus eps for stability)
  v <- v / (ssq + eps)
  
  # Remove any NaN values (in case some columns had zero variance)
  v <- v[!is.na(v)]
  
  # Compute resels per node (using the relation with FWHM)
  reselsPerNode <- sqrt(v / (4 * log(2)))
  
  # The global FWHM estimate is the reciprocal of the average resels per node
  FWHM <- 1 / mean(reselsPerNode)
  
  return(FWHM)
}


residuals_data <- function(data1, data2) {
  # Subtract the column-wise means of each group
  r1 <- data1 - matrix(colMeans(data1), nrow = nrow(data1), ncol = ncol(data1), byrow = TRUE)
  r2 <- data2 - matrix(colMeans(data2), nrow = nrow(data2), ncol = ncol(data2), byrow = TRUE)
  return(rbind(r1, r2))
}

pooled_sd <- function(data1, data2) {
  n1 <- NCOL(data1)
  n2 <- ncol(data2)
  var1 <- apply(data1, 1, var)
  var2 <- apply(data2, 1, var)
  pooled_var <- ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2)
  return(sqrt(pooled_var))
}


