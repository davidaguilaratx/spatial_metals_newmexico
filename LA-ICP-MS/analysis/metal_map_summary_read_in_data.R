library(openxlsx)
library(tidyverse)
U = read.xlsx('c:/Users/david/Documents/Research/PhD/Spatial_Omics/LA-ICP-MS/data/bc_samples/121124_sample15_curve5__20um_PbUVWK/sample15 U238_ppm matrix.xlsx',
              colNames = FALSE,skipEmptyCols = F, skipEmptyRows = F,
              na.strings=c("","NA"))

V = read.xlsx('c:/Users/david/Documents/Research/PhD/Spatial_Omics/LA-ICP-MS/data/bc_samples/121124_sample15_curve5__20um_PbUVWK/sample15 V51_ppm matrix.xlsx',
              colNames = FALSE,skipEmptyCols = F, skipEmptyRows = F,
              na.strings=c("","NA"))

Pb = read.xlsx('c:/Users/david/Documents/Research/PhD/Spatial_Omics/LA-ICP-MS/data/bc_samples/121124_sample15_curve5__20um_PbUVWK/sample15 Pb208_ppm matrix.xlsx',
               colNames = FALSE,skipEmptyCols = F, skipEmptyRows = F,
               na.strings=c("","NA"))

W = read.xlsx('c:/Users/david/Documents/Research/PhD/Spatial_Omics/LA-ICP-MS/data/bc_samples/121124_sample15_curve5__20um_PbUVWK/sample15 W182_ppm matrix.xlsx',
              colNames = FALSE,skipEmptyCols = F, skipEmptyRows = F,
              na.strings=c("","NA"))


sub_values = function(df, low_cal=0.01) {
  neg_vals = (df < 0) & !is.na(df)
  df[neg_vals] <- 0
  
  below_vals = (df > 0) & (df < low_cal) & !is.na(df)
  df[below_vals] <- low_cal / 2
  
  return(df)
}

U = sub_values(U)
V = sub_values(V)
Pb = sub_values(Pb)
W = sub_values(W)


test = unlist(U)
summary(test)
quantile(test, probs = seq(0, 1, 1/10), na.rm=T)
quantile(test, probs = c(0.9,0.95,0.99, 0.999, 0.9995, 0.9999), na.rm=T)


library(terra)
# --- Step 1: Convert your dataframes to SpatRaster objects ---
# We assume 'U', 'V', 'Pb', and 'W' are m x n dataframes or matrices.
# The rast() function can take a matrix directly.
U_r <- rast(as.matrix(U))
V_r <- rast(as.matrix(V))
Pb_r <- rast(as.matrix(Pb))
W_r <- rast(as.matrix(W))


# --- Step 2: Combine all layers into a single SpatRaster object ---
# This is your "tensor".
metal_stack <- c(U_r, V_r, Pb_r, W_r)
names(metal_stack) <- c("U", "V", "Pb", "W")

# Assign the real-world dimensions based on your 20um spot size.
# We'll set the origin (bottom-left) to (0, 0).
m <- nrow(U)
n <- ncol(U)

# Extent is defined as: (xmin, xmax, ymin, ymax)
# We multiply the number of columns/rows by the spot size (20 um).
xmin <- 0
xmax <- n * 20
ymin <- 0
ymax <- m * 20

ext(metal_stack) <- c(xmin, xmax, ymin, ymax)

# Optional: You can also set a dummy Coordinate Reference System (CRS)
# This just tells R to treat the units as "micrometers" (or whatever)
# and not latitude/longitude.
crs(metal_stack) <- "local"

# Now, explore your object:
print(metal_stack)
# class       : SpatRaster
# dimensions  : m, n, 4  (rows, cols, nlyr)
# resolution  : 20, 20  (x, y)
# extent      : 0, [xmax], 0, [ymax]  (xmin, xmax, ymin, ymax)
# coord. ref. : local
# names       : U, V, Pb, W

# You can also plot it to check
plot(metal_stack)


# Calculate the pixel-wise correlation matrix
# This returns a 4x4 matrix
cor_matrix <- layerCor(metal_stack, fun = "pearson")



print(cor_matrix)
