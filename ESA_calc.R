# Katharine Sink
# Development of Evaporation State Angle

library(tidyverse)
library(sf)
library(data.table)
library(lme4)

###############################################
### TIME SERIES DATA ###
###############################################
# use csv files for each individual basin 
directory = "F:/ESA/MACH_data"
all_files = list.files(directory, pattern = "\\.csv$", full.names = TRUE)

# get character vector of site numbers needed
gauges = site_info %>% filter(!SITENO %in% gauges)

# find files according to filename and basin id number
file_site_numbers = sub(".*basin_(\\d{8})_MACH\\.csv", "\\1", basename(all_files))
matching_files <- all_files[file_site_numbers %in% gauges]

# read in as one data frame 
combined_data <- map_dfr(matching_files, function(file) {
  # Read each CSV, adding the site number as a column
  df <- fread(file)
  df[, SITENO := sub(".*basin_(\\d{8})_MACH\\.csv", "\\1", basename(file))]
  return(df)
})

# select columns for ESA calculation
combined_data = combined_data %>% select(SITENO, DATE, PET, PRCP, AET)

# folder to save updated csv files
dest_dir = "F:/ESA/ESA_2026"

# create updated csv files for selected basins and columns 
for (file in matching_files) {
  # read the CSV
  df <- fread(file)
  
  # subset to specific columns (assuming they exist)
  df_subset <- df[, .(SITENO = sub(".*basin_(\\d{8})_MACH\\.csv", "\\1", basename(file)),
                      DATE, PRCP, PET, AET)]
  
  # define new filename (same as original)
  new_filename <- basename(file)
  new_path <- file.path(dest_dir, new_filename)
  
  # write the subsetted data to a new CSV
  fwrite(df_subset, new_path, row.names = FALSE)
}

################################################
### ESA ###
################################################
# adapted MATLAB code (Ferguson) to R
# input file PRCP, ET, PET columns 

# daily data files 
dir_path = "F:/ESA/ESA_2026"

# read MACH files 
  rdMACH = function(filename) {
    data = read_csv(filename, col_types = cols(
      SITENO = col_character(),
      DATE = col_date(),
      PRCP = col_double(),
      PET = col_double(),
      AET = col_double()))
    
    # read in columns as matrix with days as rows and 3 columns
    Y = as.matrix(data[, c("PRCP", "PET", "AET")])  
    DT = as.Date(data$DATE) # date vector
    SN = unique(data$SITENO) # site number vector
    
    if (length(SN) > 1) stop("Multiple site numbers detected")
    if (length(SN) == 0) stop("No valid site number found")
    
    return(list(Y = Y, DT = DT, SN = SN[1]))
  }

#########################################
# QEP function to calculate ESA statistic
  
QEP = function(A) {
  n = nrow(A)  # matrix A with rows as time steps
  m = ncol(A) # PRCP, PET, AET

# Budyko analysis, define indices 
  EP = A[,3] / A[,1]  # AET/PRCP
  EpP = A[,2] / A[,1]  # PET/PRCP

# evaporation state angle    
  delta = sin(135/2 * pi/180) / sqrt(2)  # constant based on 67.5 degree angle for scaling
  esa = numeric(n)
  
 for (i in 1:n) {
   # create unit vector  
   nhat = c(EpP[i] - 1, EP[i] - 1) / sqrt(sum((c(EpP[i] - 1, EP[i] - 1))^2))  
   # calculate ESA converting radians to degrees, shift by 67.5 and scale by 67.5 
    esa[i] = (atan2(nhat[2], nhat[1]) * 180/pi + 67.5) / 67.5  
  }
  return(list(esa = esa, delta = delta, EP = EP, EpP = EpP, Ynew = A))  
}

#########################################
# firsmth function
# lowpass FIR filter to smooth time series data
# removes high frequency noise (daily/seasonal fluctuations) and focuses on long-term trends

# y = matrix of data to filter 
# dt = sampling interval (1 day)
# fc = lowpass filter cutoff frequency (1/a), i.e. 0.5/365.242 is 1 cycle per 2 years, removes fluctuations shorter than 2 years
# dftb = transition bandwidth in cycles per day
  
firsmth = function(y, dt, fc, dftb) {
  y = as.matrix(y)  # PRCP, PET, AET data
  nsam = nrow(y) # number of samples/rows (16,071 days)
  nchn = ncol(y) # number of channels (columns)
  
  # determine filter length
  nfir = 4 / (dt * dftb) # based on transition bandwidth
  nfir = min(nfir, nsam)  # caps nfir at number of samples 
  nfir = floor(nfir) # rounds to an integer (16,000)
  if (nfir %% 2 == 0) nfir <- nfir - 1  # ensure filter length is odd (16,000 becomes 15,999) for symmetry, centered at zero
  
  n = 0:(nfir - 1) # sequence from 0 to 15,998 (length 15,999)
  nfiro2 = (nfir - 1) / 2  # half the filter length (15,999 - 1)/2 = 7,999 used for centering
  n_centered = n - nfiro2  # shifts indices to be symmetric around zero from -7999 to 7999, ensures filter is zero-phase (no time shift after convolution)
  
  # define Hamming window
  # reduces side lobes in frequency response, minimizing ripple effects in the filtered signal
  # normalizes indices from 0 to 1 over the filter length
  # values range from 0.07 (edges) to 1 (center) tapering the filter coefficients
  hamming_window = 0.53 - 0.46 * cos(2 * pi * n / (nfir - 1))
  
  # define sinc filter (not defined in R)
  # ideal low pass filter in the time domain 
  sinc_func = function(x) {
    y = rep(0, length(x))
    y[x == 0] = 1  # at x = 0, sinc(0) = 1
    y[x != 0] = sin(pi * x[x != 0]) / (pi * x[x != 0])  # for x not equal 0, oscillates and decays
    return(y)
  }
  
  # establish filter coefficients (15,999 values)
  # sinc function applies to 2*fc*dt*n_centered, creates core lowpass filter shape, passing frequencies below fc
  # 2*fc*dt normalizes cutoff frequency for sampling rate
  # fc = 0.5/365.242 ~= 0.00137 cycles/day, dt = 1, so 2*fc*dt ~= 0.00274
  # flp is vector of 15,999 filter coefficients, defining lowpass filter that attenuates frequencies above fc
  flp = 2 * fc * dt * sinc_func(2 * fc * dt * n_centered) * hamming_window
  
  # initialize output, smoothed matrix with same dimensions as y
  # contains smoothed time series for each channel
  s = matrix(0, nrow = nsam, ncol = nchn)
  
  # convolve each channel and apply zero-phase shift
  for (i in 1:nchn) {
    stmp = convolve(y[, i], flp, type = "open")  # full convolution, applies filter to each column of y
    start_idx = nfiro2 + 1 # trims the start of the convolution output (7999 + 1 = 8000)
    end_idx = length(stmp) - nfiro2 # trims the end (32,069 - 7,999 = 24,070)
    s[, i] = stmp[start_idx:end_idx] # ensures output s has same length as y (16,071) and is centered (no phase shift)
  }
  
  return(list(s = s, flp = flp))
}

#########################################    
# QEPanal matlab function 
# process single site's data through filtering, decimation, and ESA calculation

# stacod = site code
# Y = data matrix and DT dates from rdMACH
# head = column labels
# fc = cutoff frequency (0.5 cycles/year)
  
QEPanal = function(stacod, Y, DT, head, fc) {
  tyr = 365.242   # tropical year (days)
  dt = 1   # sample interval (days)
  fNfac = 2 # Nyquist frequency factor for decimation, 2 ensures the new sampling rate avoids aliasing
  
  nd = nrow(Y) # number of days
  T = 0:(nd-1) # time index (0 to 16,070) for detrending
  
# set headers, rename P = PRCP, Ep = PET, E = AET
  head = c("P", "Ep", "E")
  
# transition bandwidth calculation
  df = 1/(nd*dt) # frequency resolution = 1/(16,071 * 1) = 0.0000622 cycles/day
  dftb = 4*df # transition bandwidth in cycles per day, 4/16071 ~=0.00025 cycles/days, defines sharpness of cutoff

# detrending each column to make zero mean, focus on variability rather than long term trends
# fits linear model for each column against time T, storing coefficients in plin (intercept, slope)  
  Ytmp = matrix(0, nd, 3) # create temporary matrix to store data
  plin = matrix(0, 3, 2) # detrend the data by removing linear trends
  for (j in 1:3) {
    plin[j,] = lm(Y[,j] ~ T)$coefficients
    Ytmp[,j] = Y[,j] - plin[j,1] - plin[j,2]*T  # subtract linear trend, making data zero mean and trend free
  }

# low pass filtering    
# apply lowpass filter using firsmth function 
  filt = firsmth(Ytmp, dt, fc/tyr, dftb)  # with fc/tyr = 0.5/365.242 = 0.00137 cycles/day, retaining fluctuations longer than 2 years  
  Ytmp = filt$s # updated with smoothed data (16,071 rows and 3 columns)
  flp = filt$flp # filter coefficients (length 15,999)
  
# equivalent width of filter
# determines EQ of filter, used later to scale decimated data
# ratio of the sum of filter coefficients to maximum coefficient, representing effective smoothing window size  
  EW = sum(flp)/max(flp)

# decimation    
  fNnew = fNfac*fc/tyr   # Set new Nyquist frequency (2 * (0.5/365.242) = 0.00274 cycles/day)
  dtnew = 1/(2*fNnew) # new sampling interval (1/2 * fNnew = 182.5 days) meaning data will be sampled every 6 months
  
# decimation factor
ndec = max(1, ceiling(dtnew/dt)) # every 183rd sample retained
ndd = max(1, floor((nd - 1)/ndec) + 1) # number of decimated samples (16071-1)/183)+1 = 88, 88 biannual samples

# decimation
# subsamples smoothed data (Ytmp) every 183 days, reducing from 16071 daily samples to 88 biannual samples
Ynew = matrix(0, ndd, 3)
Tnew = numeric(ndd) # updates corresponding time indices
DTnew = numeric(ndd) # updates corresponding dates

j = 1
for (i in seq(1, nd, by=ndec)) {
  if (i > nd || j > ndd) break
  Tnew[j] = T[i]
  DTnew[j] = DT[i] + 719529 # Convert Unix epoch (seconds) to MATLAB datenum (days)
  Ynew[j,] = Ytmp[i,]
  j = j + 1
}

# if no samples are selected, uses first sample to ensure Ynew has valid data (no empty output)
if (j == 1) {
  Ynew[1,] = Ytmp[1,]
  Tnew[1] = T[1]
  DTnew[1] = DT[1]  + 719529
}

# restore trends and mean to ensure long term trend is preserved
  for (j in 1:3) {
    Ynew[,j] = Ynew[,j] + plin[j,1] + plin[j,2]*Tnew
  }

# scale by equivalent width  
  Ynew <- EW * Ynew # EW = 15999, converts filtered means to window sums, adjusting for filter's smoothing effect

# prepare for ESA calculation    
  A = Ynew  # final decimated and scaled data (88 rows, 3 columns)
  DTa = DTnew  # corresponding biannual dates
  
# Budyko analysis
  result = QEP(A)

# esaa = ESA values
# EPa = evaporative index values, EpPa = aridity index values 
# Ynew = decimated PRCP, PET, AET values   
  return(list(esaa=result$esa, DTa=DTa, delta=result$delta, 
              EPa=result$EP, EpPa=result$EpP, Ynew = A))
}

#########################################
# runQEPanal function

# create output file     
runQEPanal = function(fc, dir_path, output_file = "qep_output_march.csv") {
  library(readr)
  
  # validate directory
  path = dir_path
  if (!dir.exists(path)) stop("Directory does not exist")
  
  # compile list of *.csv files
  filelst = list.files(path, pattern = "\\.csv$", full.names = TRUE)
  nfile = length(filelst)
  
  if (nfile == 0) stop("No *.csv files were found")
  
  # initialize results data frame
  results = data.frame()
  
  for (i in 1:nfile) {
    filename = basename(filelst[i])
    stacod = sub("\\.csv$", "", filename)
    
    # read data
    dat = rdMACH(filelst[i])
    X = dat$Y
    DT = dat$DT + 719529 # Convert R days to MATLAB datenum
    head = c("P", "Ep", "E")
    
    # run QEP analysis
    result = QEPanal(stacod, X, DT, head, fc)
    
    # combine results
    file_results = data.frame(
      SITENO = rep(stacod, length(result$esaa)),
      ESA = result$esaa,  # filtered/decimated ESA
      P = result$Ynew[,1],  # filtered/decimated PRCP
      Ep = result$Ynew[,2],  # filtered/decimated PET
      E = result$Ynew[,3],  # filtered/decimated AET
      EP = result$EPa,  # filtered/decimated AET/PRCP
      EpP = result$EpPa  # filtered/decimated PET/PRCP
    )
    results = rbind(results, file_results)
  }
  
  # write to CSV
  write_csv(results, output_file)
  
  return(results)
}

###################################
# run all code for directory of files with established fc value
# will save in directory path listed above 
    
runQEPanal(fc = 0.5, dir_path)