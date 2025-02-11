## Analysis of fluorescence line profiles of each cell 

# Load necessary libraries
library(readxl)      # For reading Excel files
library(ggplot2)     # For plotting
library(dplyr)       # For data manipulation
library(pracma)      # For interpolation
library(tidyr)       # For data reshaping
library(writexl)     # For writing to multiple sheets in an Excel file
library(viridis)     # Colorblind-friendly palette

# ----------------------------
# SETUP: File Paths (Modify as needed)
# ----------------------------
file_path <- "raw_data_fluorescence_m_smeg.xlsx"  # Input file
output_excel_file <- "interpolated_fluorescence_profiles.xlsx"  # Processed output file
output_csv_file <- "Interpolated_Fluorescence_Profiles.csv"  # Processed output CSV
output_plot_all <- "fluorescence_profiles_all_cells.png"  # Plot of all cell profiles
output_plot_mean <- "fluorescence_mean.png"  # Plot of mean fluorescence

# ----------------------------
# STEP 1: Load and Preprocess Data
# ----------------------------
# Load first sheet
raw_data <- read_excel(file_path, sheet = "Sheet1", col_names = FALSE)

# Extract cell numbers from the first row
cell_numbers <- as.character(raw_data[1, ])
cell_numbers <- cell_numbers[seq(1, length(cell_numbers), by = 2)]

# Generate new column names
new_colnames <- c()
for (i in 1:length(cell_numbers)) {
  new_colnames <- c(new_colnames, paste0(cell_numbers[i], "_x"), paste0(cell_numbers[i], "_y"))
}
colnames(raw_data) <- new_colnames  # Apply new column names

# Remove first two rows (headers) and convert to numeric
raw_data <- raw_data[-c(1, 2), ] %>% mutate_all(as.numeric)

# ----------------------------
# STEP 2: Normalize Data
# ----------------------------
# Normalize x values to fractional cell length (0-1)
for (i in seq(1, length(colnames(raw_data)), by = 2)) {
  x_col <- colnames(raw_data)[i]  
  raw_data[[x_col]] <- (raw_data[[x_col]] - min(raw_data[[x_col]], na.rm = TRUE)) /
    (max(raw_data[[x_col]], na.rm = TRUE) - min(raw_data[[x_col]], na.rm = TRUE))
}

# Normalize fluorescence values (y) using min-max scaling
for (i in seq(1, length(colnames(raw_data)), by = 2)) {
  y_col <- colnames(raw_data)[i+1]  
  raw_data[[y_col]] <- (raw_data[[y_col]] - min(raw_data[[y_col]], na.rm = TRUE)) /  
    (max(raw_data[[y_col]], na.rm = TRUE) - min(raw_data[[y_col]], na.rm = TRUE))
}

# ----------------------------
# STEP 3: Interpolate Data
# ----------------------------
num_points <- 100   # Define interpolation points
interpolated_profiles <- list()

for (i in seq(1, length(colnames(raw_data)), by = 2)) {  
  x_col <- colnames(raw_data)[i]    
  y_col <- colnames(raw_data)[i+1]  
  
  x <- as.numeric(na.omit(raw_data[[x_col]]))
  y <- as.numeric(na.omit(raw_data[[y_col]]))
  
  x_norm <- (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  new_x <- seq(0, 1, length.out = num_points)   
  y_interp <- approx(x_norm, y, xout = new_x, method = "linear", rule = 2)$y  
  
  interpolated_profiles[[x_col]] <- y_interp
}

interpolated_matrix <- do.call(rbind, interpolated_profiles)
interpolated_df <- as.data.frame(t(interpolated_matrix))
interpolated_df$Fractional_Length <- seq(0, 1, length.out = num_points)

# Convert to long format for plotting
interpolated_long <- pivot_longer(interpolated_df, cols = -Fractional_Length, 
                                  names_to = "Cell_ID", values_to = "Fluorescence")
interpolated_long$Cell_ID <- gsub("_x", "", interpolated_long$Cell_ID)

# ----------------------------
# STEP 4: Save Processed Data
# ----------------------------
write.csv(interpolated_long, output_csv_file, row.names = FALSE)

cell_data_list <- lapply(split(interpolated_long, interpolated_long$Cell_ID), function(df) {
  df <- df[, c("Fractional_Length", "Fluorescence")]
  return(df)
})
write_xlsx(cell_data_list, output_excel_file)

# ----------------------------
# STEP 5: Visualization
# ----------------------------
p1 <- ggplot(interpolated_long, aes(x = Fractional_Length, y = Fluorescence, group = Cell_ID, color = Cell_ID)) +
  geom_line(size = 1, alpha = 0.8) +  
  scale_color_viridis_d(option = "plasma") +  
  labs(x = "Normalized Cell Length", y = "Normalized Fluorescence", title = "Fluorescence Profiles for All Cells") +  
  theme_classic(base_size = 14) +  
  theme(legend.position = "top", axis.text.x = element_blank(), axis.ticks.x = element_blank())  

# Compute mean & standard deviation profile
mean_profile <- interpolated_df %>%
  select(-Fractional_Length) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(everything(), names_to = "Fractional_Length", values_to = "Mean_Fluorescence")
mean_profile$Fractional_Length <- seq(0, 1, length.out = num_points)

std_profile <- interpolated_df %>%
  select(-Fractional_Length) %>%
  summarise(across(everything(), sd, na.rm = TRUE)) %>%
  pivot_longer(everything(), names_to = "Fractional_Length", values_to = "SD_Fluorescence")
std_profile$Fractional_Length <- seq(0, 1, length.out = num_points)

mean_profile <- left_join(mean_profile, std_profile, by = "Fractional_Length")

p2 <- p1 +
  geom_line(data = mean_profile, aes(x = Fractional_Length, y = Mean_Fluorescence), color = "black", size = 1.2) +
  geom_ribbon(data = mean_profile, aes(x = Fractional_Length, ymin = Mean_Fluorescence - SD_Fluorescence, 
                                       ymax = Mean_Fluorescence + SD_Fluorescence), fill = "gray", alpha = 0.3)

# Save plots
ggsave(output_plot_all, p1, width = 8, height = 6, dpi = 300)
ggsave(output_plot_mean, p2, width = 8, height = 6, dpi = 300)

print(p1)
print(p2)
