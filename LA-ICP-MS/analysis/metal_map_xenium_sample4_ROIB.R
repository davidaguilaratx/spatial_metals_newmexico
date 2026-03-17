library(tidyverse)
U = read.csv('c:/Users/david/Documents/Research/PhD/Spatial_Omics/LA-ICP-MS/data/bc_samples/with_xenium/250610_Aguilar_005-22-09_ROI_B_PbUVWAsCdFeNiZnK/sample9b U238_ppm matrix.csv',
              blank.lines.skip = FALSE, header=F,
              na.strings=c("","NA"))

V = read.csv('c:/Users/david/Documents/Research/PhD/Spatial_Omics/LA-ICP-MS/data/bc_samples/with_xenium/250610_Aguilar_005-22-09_ROI_B_PbUVWAsCdFeNiZnK/sample9b V51_ppm matrix.csv',
             blank.lines.skip = FALSE, header=F,
             na.strings=c("","NA"))

Pb = read.csv('c:/Users/david/Documents/Research/PhD/Spatial_Omics/LA-ICP-MS/data/bc_samples/with_xenium/250610_Aguilar_005-22-09_ROI_B_PbUVWAsCdFeNiZnK/sample9b Pb208_ppm matrix.csv',
              blank.lines.skip = FALSE, header=F,
              na.strings=c("","NA"))

W = read.csv('c:/Users/david/Documents/Research/PhD/Spatial_Omics/LA-ICP-MS/data/bc_samples/with_xenium/250610_Aguilar_005-22-09_ROI_B_PbUVWAsCdFeNiZnK/sample9b W182_ppm matrix.csv',
             blank.lines.skip = FALSE, header=F,
             na.strings=c("","NA"))

Ni = read.csv('c:/Users/david/Documents/Research/PhD/Spatial_Omics/LA-ICP-MS/data/bc_samples/with_xenium/250610_Aguilar_005-22-09_ROI_B_PbUVWAsCdFeNiZnK/sample9b Ni60_ppm matrix.csv',
              blank.lines.skip = FALSE, header=F,
              na.strings=c("","NA"))

Cd = read.csv('c:/Users/david/Documents/Research/PhD/Spatial_Omics/LA-ICP-MS/data/bc_samples/with_xenium/250610_Aguilar_005-22-09_ROI_B_PbUVWAsCdFeNiZnK/sample9b Cd111_ppm matrix.csv',
              blank.lines.skip = FALSE, header=F,
              na.strings=c("","NA"))

As = read.csv('c:/Users/david/Documents/Research/PhD/Spatial_Omics/LA-ICP-MS/data/bc_samples/with_xenium/250610_Aguilar_005-22-09_ROI_B_PbUVWAsCdFeNiZnK/sample9b As75_ppm matrix.csv',
              blank.lines.skip = FALSE, header=F,
              na.strings=c("","NA"))

P = read.csv('c:/Users/david/Documents/Research/PhD/Spatial_Omics/LA-ICP-MS/data/bc_samples/with_xenium/250610_Aguilar_005-22-09_ROI_B_PbUVWAsCdFeNiZnK/sample9b P31_CPS matrix.csv',
             blank.lines.skip = FALSE, header=F,
             na.strings=c("","NA"))
#asd
test = unlist(P)
summary(test)
quantile(test, probs = seq(0, 1, 1/10), na.rm=T)
quantile(test, probs = c(0.9,0.95,0.99, 0.999, 0.9995, 0.9999), na.rm=T)

sub_values = function(df, low_cal=0.01) {
  neg_vals = ((df < 0) & !is.na(df))
  df[neg_vals] <- 0
  
  below_vals = ((df > 0) & (df < low_cal) & !is.na(df))
  df[below_vals] <- low_cal / 2
  
  return(df)
}

U = sub_values(U)
V = sub_values(V)
Pb = sub_values(Pb)
W = sub_values(W)
Ni = sub_values(Ni,low_cal = 0.2)
Cd = sub_values(Cd)
As = sub_values(As, low_cal = 0.02)
P = sub_values(P, low_cal = 8000)

test = unlist(Cd)
summary(test)
quantile(test, probs = seq(0, 1, 1/10), na.rm=T)
quantile(test, probs = c(0.9,0.95,0.99, 0.999, 0.9995, 0.9999, 0.9999), na.rm=T)






library(terra)
library(spatialEco)
# --- Step 1: Convert your dataframes to SpatRaster objects ---
# We assume 'U', 'V', 'Pb', and 'W' are m x n dataframes or matrices.
# The rast() function can take a matrix directly.
U_r <- rast(as.matrix(U))
V_r <- rast(as.matrix(V))
Pb_r <- rast(as.matrix(Pb))
W_r <- rast(as.matrix(W))
Ni_r <- rast(as.matrix(Ni))
Cd_r <- rast(as.matrix(Cd))
As_r <- rast(as.matrix(As))
P_r <- rast(as.matrix(P))

# --- Step 2: Combine all layers into a single SpatRaster object ---
# This is your "tensor".
metal_stack <- c(U_r, V_r, Pb_r, W_r, Ni_r, Cd_r, As_r, P_r)
names(metal_stack) <- c("U", "V", "Pb", "W", "Ni", "Cd", "As", "P")

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


metal_names = names(metal_stack)


# Calculate the pixel-wise correlation matrix
# This returns a 4x4 matrix
cor_matrix <- layerCor(metal_stack, fun = "cor")
print(cor_matrix)




vals <- values(metal_stack, mat=TRUE)

# Find pixels where ANY metal has a value (i.e., tissue was measured)
tissue_mask <- rowSums(!is.na(vals[, metal_names])) > 0
tissue_idx <- which(tissue_mask)




cat("Total pixels:", nrow(vals), "\n")
cat("Tissue pixels:", sum(tissue_mask), "\n")
cat("NA pixels (no measurement):", sum(!tissue_mask), "\n")
# 
# Step 2: Keep only tissue pixels
vals_tissue <- vals[tissue_mask, ]

# 
# # Step 3: Subsample from tissue pixels only
set.seed(281330800)
sample_size <- 50000
if(nrow(vals_tissue) > sample_size) {
  sample_idx <- sample(nrow(vals_tissue), sample_size)
  vals_sample <- vals_tissue[sample_idx, ]
} else {
  vals_sample <- vals_tissue
}




cor_matrix <- cor(vals_sample, use="pairwise.complete.obs")


# Calculate p-values for each pair
n_layers <- nlyr(metal_stack)
p_values <- matrix(NA, n_layers, n_layers)
corr_est = matrix(NA, n_layers, n_layers)

list_of_corr.tests = replicate(64, NA, simplify = F)
dim(list_of_corr.tests) = c(8,8)
rownames(list_of_corr.tests) = names(metal_stack)
colnames(list_of_corr.tests) = names(metal_stack)

for(i in 1:n_layers) {
  for(j in 1:n_layers) {
    if(i != j) {
      test <- cor.test(vals_sample[,i], vals_sample[,j], method = 'spearman')
      p_values[i,j] <- test$p.value
      list_of_corr.tests[[i,j]] = test
    }
  }
}

list_of_corr.tests


extract_corr_value <- function(list_of_tests, value = "p.value") {
  p_values <- matrix(sapply(list_of_tests, function(x) if(is.list(x)) x[[value]] else NA), nrow=8, ncol=8, byrow=TRUE)
  dimnames(p_values) <- list(metal_names, metal_names)
  return(p_values)
}


p_values <- extract_corr_value(list_of_corr.tests, "p.value")
estimates <- extract_corr_value(list_of_corr.tests, "estimate")
p_values
#               U             V            Pb             W            Ni            Cd           As             P
# U             NA  3.705571e-76  0.000000e+00  5.659391e-16 1.274978e-112 1.584197e-299 2.694600e-14  0.000000e+00
# V   3.705571e-76            NA 1.012824e-159  5.103529e-96  0.000000e+00  6.900768e-02 5.222137e-12  3.115725e-92
# Pb  0.000000e+00 1.012824e-159            NA 2.772725e-137 1.298073e-173  5.096056e-23 5.684031e-03 3.856238e-258
# W   5.659391e-16  5.103529e-96 2.772725e-137            NA  5.287500e-40  7.372916e-08 4.917597e-01  2.402168e-11
# Ni 1.274978e-112  0.000000e+00 1.298073e-173  5.287500e-40            NA  3.081173e-06 2.888285e-13 1.190257e-232
# Cd 1.584197e-299  6.900768e-02  5.096056e-23  7.372916e-08  3.081173e-06            NA 2.419371e-04  0.000000e+00
# As  2.694600e-14  5.222137e-12  5.684031e-03  4.917597e-01  2.888285e-13  2.419371e-04           NA  1.278359e-25
# P   0.000000e+00  3.115725e-92 3.856238e-258  2.402168e-11 1.190257e-232  0.000000e+00 1.278359e-25            NA

estimates
#             U           V         Pb           W         Ni           Cd         As           P
# U          NA 0.082454827 0.21421193  0.03619739 0.10060349  0.164308470 0.03403380  0.44190918
# V  0.08245483          NA 0.11999653  0.09278675 0.23255206  0.008132064 0.03084850  0.09090897
# Pb 0.21421193 0.119996534         NA  0.11118984 0.12513139  0.044163127 0.01236741  0.15259101
# W  0.03619739 0.092786747 0.11118984          NA 0.05915219 -0.024065427 0.00307472 -0.02986447
# Ni 0.10060349 0.232552058 0.12513139  0.05915219         NA  0.020862256 0.03263650  0.14488014
# Cd 0.16430847 0.008132064 0.04416313 -0.02406543 0.02086226           NA 0.01641494  0.26382430
# As 0.03403380 0.030848499 0.01236741  0.00307472 0.03263650  0.016414938         NA  0.04676723
# P  0.44190918 0.090908969 0.15259101 -0.02986447 0.14488014  0.263824300 0.04676723          NA





# Method 2: Local Spatial Correlation (Moving Window)
# This is a true spatial method. It answers: 
# "How does the correlation between U and V change across the tissue?"
w_size <- 3

# Get the names of your layers
metal_names <- names(metal_stack)

# Create an empty list to store the 21 output rasters
local_cor_results <- list()

# --- Run the Loop ---

# Loop through all unique pairs
# nlyr(metal_stack) is now 7
for (i in 1:(nlyr(metal_stack) - 1)) {
  for (j in (i + 1):nlyr(metal_stack)) {
    
    # Get the layer names
    name_i <- metal_names[i]
    name_j <- metal_names[j]
    
    # Create a clear name for the output, e.g., "cor_U_V"
    list_name <- paste0("cor_", name_i, "_", name_j)
    
    # Let us know what's happening
    cat("Calculating local correlation for:", list_name, "\n")
    
    # Run the terra function on the two layers
    cor_map <- rasterCorrelation(metal_stack[[i]], 
                                 metal_stack[[j]], 
                                 s = w_size)
    
    # Store the resulting map in the list
    local_cor_results[[list_name]] <- cor_map
  }
}

cat("--- Calculation Complete: 21 pairs processed --- \n")



print(local_cor_results$cor_U_Pb)
# $cor_U_V# $cor_U_V# $cor_U_V
# class       : SpatRaster
# ...
# $cor_U_Pb
# class       : SpatRaster
# ...
# [and 4 more]

# To access and plot a specific map:
# For example, the correlation between Lead (Pb) and Tungsten (W)
plot(local_cor_results$cor_U_Pb, 
     main = "Local Correlation (r) between Pb and U (3x3 window)")






autocor(metal_stack$P, w=matrix(c(1,1,1,1,0,1,1,1,1),nrow=3), method="moran", global=T)



# bivaraite morans I

# library(spdep)
# library(terra)
# 
# Step 1: Identify tissue pixels (where at least one metal is measured)
metal_df <- as.data.frame(metal_stack, xy = TRUE, na.rm = FALSE)
coords <- as.matrix(metal_df[, c("x", "y")])
# 
# Find pixels where ANY metal has a value (i.e., tissue was measured)
tissue_mask <- rowSums(!is.na(metal_df[, metal_names])) > 0
tissue_idx <- which(tissue_mask)
# # 
cat("Total pixels:", nrow(metal_df), "\n")
cat("Tissue pixels:", sum(tissue_mask), "\n")
cat("NA pixels (no measurement):", sum(!tissue_mask), "\n")
# 
# Step 2: Keep only tissue pixels
metal_df_tissue <- metal_df[tissue_mask, ]
coords_tissue <- coords[tissue_mask, ]
# 
# # Step 3: Subsample from tissue pixels only
set.seed(281330800)
sample_size <- 50000
if(nrow(metal_df_tissue) > sample_size) {
  sample_idx <- sample(nrow(metal_df_tissue), sample_size)
  metal_df_sample <- metal_df_tissue[sample_idx, ]
  coords_sample <- coords_tissue[sample_idx, ]
} else {
  metal_df_sample <- metal_df_tissue
  coords_sample <- coords_tissue
}

cat("Sampled pixels:", nrow(metal_df_sample), "\n")
# 
# # Step 4: Create spatial weights ONLY for tissue pixels
# # This ensures neighbors are only within actual tissue
# knn <- knearneigh(coords_sample, k=8)
# nb <- knn2nb(knn)
# lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
# 
# # Step 5: Bivariate Moran's I function (handles pairwise NAs)
# bivariate_morans_I <- function(var1, var2, coords_all, listw_all) {
#   # Find pixels where BOTH variables have values
#   complete <- complete.cases(var1, var2)
# 
#   if(sum(complete) < 100) {
#     return(list(statistic = NA, p.value = NA, n = sum(complete)))
#   }
# 
#   # Standardize
#   var1_std <- scale(var1_clean)[,1]
#   var2_std <- scale(var2_clean)[,1]
# 
#   # Calculate Bivariate Moran's I
#   n <- length(var1_std)
#   W_mat <- listw2mat(lw)
#   I <- (t(var1_std) %*% (W_mat %*% var2_std)) / n
# 
#   # Permutation test
#   n_perm <- 499
#   I_perm <- replicate(n_perm, {
#     var2_perm <- sample(var2_std)
#     (t(var1_std) %*% (W_mat %*% var2_perm)) / n
#   })
# 
#   p_val <- (sum(abs(I_perm) >= abs(I)) + 1) / (n_perm + 1)
# 
#   list(statistic = as.numeric(I), p.value = p_val, n = sum(complete))
# }
# 
# # Step 6: Calculate for all pairs
# morans_I_matrix <- matrix(NA, 8, 8)
# morans_p_matrix <- matrix(NA, 8, 8)
# morans_n_matrix <- matrix(NA, 8, 8)
# dimnames(morans_I_matrix) <- dimnames(morans_p_matrix) <- list(metal_names, metal_names)
# 
# 
# var1 <- metal_df_sample[[metal_names[1]]]
# var2 <- metal_df_sample[[metal_names[8]]]
# 
# 
# result <- moran_bv(var1, var2, lw, nsim = 499, scale = TRUE)
# 
# test = xapp(x=metal_stack$U, y=metal_stack$Pb, fun=cor)
# 
# cat("\nCalculating Bivariate Moran's I for all pairs...\n")
# for(i in 1:8) {
#   for(j in 1:8) {
#     if(i != j) {
#       result <- bivariate_morans_I(
#         metal_df_sample[[metal_names[i]]],
#         metal_df_sample[[metal_names[j]]],
#         coords_sample,
#         lw
#       )
#       morans_I_matrix[i,j] <- result$statistic
#       morans_p_matrix[i,j] <- result$p.value
#       morans_n_matrix[i,j] <- result$n
#     }
#   }
#   cat("Completed row", i, "of 8\n")
# }
# 
# # Step 7: View results
# print("\nBivariate Moran's I values:")
# print(round(morans_I_matrix, 3))
# 
# print("\nP-values:")
# print(round(morans_p_matrix, 3))
# 
# print("\nSample sizes (pixels with both measurements):")
# print(morans_n_matrix)
# 
# # Step 8: Visualize
# library(corrplot)
# corrplot(morans_I_matrix, method = "color", type = "upper",
#          p.mat = morans_p_matrix, sig.level = 0.05, insig = "pch",
#          tl.col = "black", tl.srt = 45,
#          title = "Bivariate Moran's I\n(* = p < 0.05)",
#          mar = c(0,0,2,0),
#          is.corr = FALSE,
#          col = colorRampPalette(c("blue", "white", "red"))(200),
#          na.label = " ")
# 
# 
# saveRDS(morans_I_matrix, 'c:/users/david/documents/research/PhD/spatial_omics/LA-ICP-MS/Data/morans_i_matrix.rds')
# saveRDS(morans_pvalue_matrix, 'c:/users/david/documents/research/PhD/spatial_omics/LA-ICP-MS/Data/morans_pvalue_matrix.rds')
# saveRDS(morans_n_matrix, 'c:/users/david/documents/research/PhD/spatial_omics/LA-ICP-MS/Data/morans_n_matrix.rds')




library(spatialEco)
# 
# spatial_cor <- raster.modified.ttest(
#   x = metal_stack$U,
#   y = metal_stack[[3]],
#   d = "auto",              # Auto-determine neighborhood distance
#   sample = "hexagonal",    # Recommended sampling strategy
#   p = 0.06              # Sample 10% of pixels
# )


spatial_cor <- raster.modified.ttest(
  x = rast(as.matrix(metal_df_sample$U)),
  y = rast(as.matrix(metal_df_sample$Pb)),
  d = "auto",              # Auto-determine neighborhood distance
  sample = "hexagonal",    # Recommended sampling strategy
  p = 0.1              # Sample 10% of pixels
)




# library(SpatialPack)

# result <- modified.ttest(
#   x = metal_df_sample$U,
#   y = metal_df_sample$Pb,
#   coords = coords_sample,
#   nclass = 13 # default
# )

