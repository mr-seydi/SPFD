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


###################Functions###########################

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



##############Methods functions###########

Initialize_method_list <- function(Methods, Conti_size=101, Iter_number=100){
  method_list <- list()
  for (M in Methods) {
    method_list[[M]] <- matrix(,nrow = Conti_size, ncol = 0)
  }
  return(method_list)
}

Power_data_generator <- function(Sample_size, Data,
                                 Signal, Conti_size = 101,
                                 Noise_mu, Noise_sig, Noise_fwhm,
                                 n_evaluation_points = NULL){
  
  noise1 <- noise_guassian_curve(number_of_curves = Sample_size,
                                 continuum_size = Conti_size)
  noise2 <- noise_guassian_curve(number_of_curves = Sample_size,
                                 continuum_size = Conti_size)
  noise1 <- smoothed_gussian_curves(noise1, Noise_mu, Noise_sig, Noise_fwhm)
  noise2 <- smoothed_gussian_curves(noise2, Noise_mu, Noise_sig, Noise_fwhm)
  
  if (is.null(Data)){
    data1 <- data_generator(signal = Signal, noise = noise1)
    data2 <- data_generator(noise = noise2)    
  } else if (is.null(Signal)) {
    data1 <- data_generator(data = Data[,1], noise = noise1)
    data2 <- data_generator(data = Data[,2], noise = noise2)
  } else {
    data1 <- data_generator(data = Data, signal = Signal, noise = noise1)
    data2 <- data_generator(data = Data, noise = noise2)      
  }
  
  if (!is.null(n_evaluation_points)) {
    dense_data <- dense_data(data1=data1, data2=data2, eval_numbers = n_evaluation_points)
    data1 <- t(dense_data$dense_data1)
    data2 <- t(dense_data$dense_data2)
  }
  
  return(list(data1 = data1, data2 = data2))
}

Pvalue_calculator <- function(method_list, data1, data2){
  
  Methods <- names(method_list)
  
  for (M in Methods) {
    Pvalues <- Pval_method(sampel1 = t(data1), sample2 = t(data2),
                           method = M)
    #p_values dimension is continuum_size*Iter_number
    method_list[[M]] <- cbind(method_list[[M]], Pvalues)
  }
  return(method_list) #Filled in method list
}


Power_calculator <- function(Pvalues, Iter_number=NULL, Alpha){
  
  Methods <- names(Pvalues)
  
  for (M in Methods) {
    
    if (is.null(Iter_number)) {
      Iter_number <- ncol(Pvalues[[M]])
    }
    
    # Check if the method is either "ERL" or "IATSE"
    #Because those methods do not return p-values
    if (M %in% c("ERL", "IATSE")) {
      pvalue_less_alpha <- Pvalues[[M]]
      # For "ERL" and "IATSE", only run this code
      power <- sum(colSums(pvalue_less_alpha) > 0) / Iter_number
    } else {
      pvalue_less_alpha <- Pvalues[[M]] < Alpha
      power <- sum(colSums(pvalue_less_alpha) > 0) / Iter_number
    }
    
    Pvalues[[M]] <- power
  }
  
  return(Pvalues)
  
}

check_false_positives <- function(sig_matrix, D1) {
  apply(sig_matrix, 2, function(col) {
    sum(col[D1]) / length(D1)
  })
  
}

Sensetivity_calculator <- function(Pvalues, D1, Alpha){
  
  Methods <- names(Pvalues)
  
  for (M in Methods) {
    
    
    # Check if the method is either "ERL" or "IATSE"
    #Because those methods do not return p-values
    if (M %in% c("ERL", "IATSE")) {
      pvalue_less_alpha <- Pvalues[[M]]
      
      #check false positives
      sensetivity <- mean(check_false_positives(pvalue_less_alpha, D1))
      
    } else {
      pvalue_less_alpha <- Pvalues[[M]] < Alpha
      sensetivity <-  mean(check_false_positives(pvalue_less_alpha, D1))
    }
    
    Pvalues[[M]] <- sensetivity
  }
  
  return(Pvalues)
  
}


Pval_method <- function(sampel1,sample2,method) {
  if (method=="IWT") {
    pval <- IWT(sampel1,sample2)
  } else if (method=="TWT"){
    pval <- TWT(sampel1,sample2)
  } else if (method=="Parametric_SPM"){
    pval <- p_spm(sampel1,sample2)
  } else if (method=="Nonparametric_SPM"){
    pval <- p_snpm(sampel1,sample2)
  } else if (method=="ERL"){
    pval <- ERL(sampel1,sample2)
  } else if (method=="IATSE"){
    pval <- IATSE(sampel1,sample2)
  } else {
    stop("Choose a method between options")
  }
  return(pval)
}


p_spm <- function(data1, data2){
  # spm  <- spm1d$stats$ttest2(data1, data2, equal_var=FALSE)
  # p_val <- spm1d$rft1d$f$sf((spm$z)^2, spm$df, spm$Q, spm$fwhm, withBonf=TRUE)
  p_val <- SPM(data1, data2, variance.equal = FALSE, Clevel = 0.95)
  return(p_val)
}

p_snpm <- function(data1, data2, B = 1000){
  
  n1 = dim(data1)[1]
  n2 = dim(data2)[1]
  
  
  group12 = factor(c(rep(1,n1),rep(2,n2)))
  data_group12 <- rbind(data1,data2)
  
  # Create a data frame that includes both data_group12 and group12
  combined_data <- data.frame(data_group12, group12)
  
  # Pass the formula with the combined data to Fmax
  Fmax_pval <- Fmax(data_group12 ~ group12, DATA = combined_data)
  return(Fmax_pval$adjusted_pval_F)
  
}







############## Parallel processing functions ####################
# Custom combine function to accumulate matrices across iterations
# parallel_combine_function <- function(x, y) {
#   for (M in names(x)) {
#     # Combine matrices from the same method in both x and y
#     x[[M]] <- cbind(x[[M]], y[[M]])
#   }
#   return(x)
# }

# 2) A generic combine that dispatches on matrix vs numeric
# parallel_combine_function <- function(x, y) {
#   for (nm in names(x)) {
#     xi <- x[[nm]]
#     yi <- y[[nm]]
#     
#     #–– both sub‐lists?  Recurse:
#     if (is.list(xi) && is.list(yi)) {
#       x[[nm]] <- parallel_combine_function(xi, yi)
#       
#       #–– either matrices or data.frames?  cbind them:
#     } else if ((is.matrix(xi)    || is.data.frame(xi)) &&
#                (is.matrix(yi)    || is.data.frame(yi))) {
#       # normalize to data.frame so columns line up
#       x[[nm]] <- cbind(as.data.frame(xi), as.data.frame(yi))
#       
#       #–– both numeric vectors?  concatenate:
#     } else if (is.numeric(xi) && is.numeric(yi)) {
#       x[[nm]] <- c(xi, yi)
#       
#     } else {
#       stop(sprintf(
#         "Cannot combine element '%s' classes %s vs %s",
#         nm, class(xi)[1], class(yi)[1]
#       ))
#     }
#   }
#   x
# }

parallel_combine_function <- function(x, y) {
  # if one side is NULL (e.g. first pass), just take the other
  if (is.null(x)) return(y)
  if (is.null(y)) return(x)
  
  # 1) stitch together the two method_list’s
  #    each is itself a list of two matrices of the same shape
  ml_combined <- Map(function(mx, my) {
    # both mx & my are matrices, so cbind them
    cbind(mx, my)
  }, x$method_list, y$method_list)
  
  # 2) concatenate the FWHM estimates
  noise_fwhm_combined <- c(x$est_noise_fwhm, y$est_noise_fwhm)
  data_fwhm_combined  <- c(x$est_data_fwhm,  y$est_data_fwhm)
  
  list(
    method_list    = ml_combined,
    est_noise_fwhm = noise_fwhm_combined,
    est_data_fwhm  = data_fwhm_combined
  )
}





## Replace Excel output with fast binary .fst files
# Install fst if needed:
# install.packages("fst")

# Handle the loop result list and save each element as a separate .fst file
write_results_to_fst <- function(loop, base_path = "results") {
  
  # Iterate over each element in the list
  for (method_name in names(loop)) {
    # Extract the data frame/matrix
    df <- as.data.frame(loop[[method_name]])
    
    out_file <- file.path(paste0(base_path,"_",method_name, ".fst"))
    
    # Write using high-speed fst format
    write_fst(df, out_file, compress = 50)  # compress level 0-100; 50 is a good balance
  }
  
}

write_estimated_noisefwhm_to_fst <- function(fwhm_est, base_path = "noisefwhm") {
  # Convert to data frame
  fwhm_df <- as.data.frame(fwhm_est)
  
  # Write using high-speed fst format
  write_fst(fwhm_df, file.path(paste0(base_path,"_noisefwhm_est.fst")), compress = 50)
}

write_estimated_datafwhm_to_fst <- function(fwhm_est, base_path = "datafwhm") {
  # Convert to data frame
  fwhm_df <- as.data.frame(fwhm_est)
  
  # Write using high-speed fst format
  write_fst(fwhm_df, file.path(paste0(base_path,"_datafwhm_est.fst")), compress = 50)
}



############simulation function####################
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

######dense data#####

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


