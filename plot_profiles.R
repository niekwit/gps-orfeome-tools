### plot_profiles.R
# This script will plot the barcode profiles of a given set of genes in a grid.

#Input: positional arguments
#1. Path to _barcode.summary.csv file (GPSW output)
#2. Path to text file containing ORF names (one per line)
#3. String indicating dimensions of the plot (e.g., "2x3")
#4. String indicating whether to include barcodes with
#   twin peaks (twinpeaks=TRUE/FALSE)

# Redirect R output to log
log <- file(
  paste0(
    "plot_profiles_",
    Sys.Date(),
    "_",
    format(Sys.time(), "%H-%M-%S"),
    ".log"
  ),
  open = "wt"
)

sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(cowplot)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 4) {
  stop(paste0(
    "Error: Four arguments are required: ",
    "<barcode_summary.csv> <orf_names.txt> <dimensions> <twinpeaks>"
  ))
}

# Check if valid GPSW output file is provided
csv_file <- args[1]
if (!file.exists(csv_file)) {
  stop(paste("Error: The file", csv_file, "does not exist."))
} else if (!grepl("_barcode\\.summary\\.csv$", csv_file)) {
  # Check if it ends with _barcode.summary.csv
  stop(paste(
    "Error: The file",
    csv_file,
    "is not a valid GPSW output file (should end with '_barcode.summary.csv')."
  ))
} else {
  csv <- read_csv(csv_file, show_col_types = FALSE)
}

# Check if ORF names file is provided
orf_file <- args[2]
if (!file.exists(orf_file)) {
  stop(paste("Error: The file", orf_file, "does not exist."))
} else {
  orfs_df <- read_csv(orf_file, col_names = FALSE, show_col_types = FALSE) %>%
    # rename columns to orf_id, category
    rename(orf_id = X1, category = X2)
  orfs <- orfs_df[["orf_id"]]
}

# Check if dimensions are provided
dimensions <- args[3]
if (is.na(dimensions) || !grepl("^[0-9]+x[0-9]+$", dimensions)) {
  stop(
    paste0(
      "Error: Invalid dimensions format. ",
      "Please provide dimensions in the format 'NxM' ",
      "(e.g., '2x3')."
    )
  )
} else {
  dimensions <- str_split(dimensions, "x")[[1]]
  nrow <- as.numeric(dimensions[1])
  ncol <- as.numeric(dimensions[2])
}

# Check if twin peaks option is provided
twinpeaks <- args[4]
if (
  is.na(twinpeaks) ||
    !twinpeaks %in% c("keeptwinpeaks=TRUE", "keeptwinpeaks=FALSE")
) {
  stop(paste0(
    "Error: Invalid twin peaks option. ",
    "Please provide 'keeptwinpeaks=TRUE' or 'keeptwinpeaks=FALSE'."
  ))
}

if (twinpeaks == "keeptwinpeaks=TRUE") {
  twinpeaks <- TRUE
} else {
  twinpeaks <- FALSE
}

# Extract reference and test sample names from the CSV file name
file_name <- basename(csv_file)
file_name <- sub("_barcode\\.summary\\.csv$", "", file_name)
ref_sample <- str_split(file_name, "_vs_")[[1]][1]
test_sample <- str_split(file_name, "_vs_")[[1]][2]

# Filter out twin peaks if specified
if (twinpeaks) {
  # Only keep barcodes without twin peaks
  df <- csv %>%
    filter(twin_peaks == "FALSE")
} else {
  # Keep all barcodes
  df <- csv
}

# Subset csv for the specified genes
# Also create columns for colour gradients
df <- df %>%
  mutate(gene.id = paste0(gene, "_", orf_id)) %>%
  filter(orf_id %in% orfs) %>%
  left_join(orfs_df, by = "orf_id", multiple = "all") %>%
  dplyr::select(
    gene.id,
    category,
    starts_with(ref_sample),
    starts_with(test_sample),
    twin_peaks,
    delta_PSI_mean,
    delta_PSI_SD,
  ) %>%
  group_by(gene.id) %>%
  mutate(barcode_number = row_number(), barcode_sum = n()) %>%
  ungroup() %>%
  pivot_longer(
    cols = c(starts_with(ref_sample), starts_with(test_sample)),
    names_to = "sample",
    values_to = "proportion"
  ) %>%
  separate(sample, into = c("condition", "bin"), sep = "\\_") %>%
  group_by(gene.id, condition, barcode_number) %>%
  mutate(
    colour = ifelse(
      condition == ref_sample,
      scales::seq_gradient_pal("black", "grey", "Lab")(seq(
        0,
        1,
        length.out = barcode_sum[1]
      ))[barcode_number],
      scales::seq_gradient_pal("red4", "red", "Lab")(seq(
        0,
        1,
        length.out = barcode_sum[1]
      ))[barcode_number]
    )
  )

# Create a small dummy data frame to define legend entries
legend_df <- tibble::tibble(
  bin = 1,
  proportion = 0,
  colour = c("#8B0000", "#000000"),
  gene.id = df$gene.id[1] # or any of the valid gene IDs
)

# Create a label df
label_df <- df %>%
  group_by(gene.id) %>%
  summarise(
    delta_PSI_mean = first(delta_PSI_mean),
    delta_PSI_SD = first(delta_PSI_SD)
  ) %>%
  ungroup() %>%
  mutate(
    proportion_max = max(df$proportion, na.rm = TRUE) * 0.95
  )

# Create the plot
p <- ggplot(
  df,
  aes(
    x = bin,
    y = proportion,
    group = interaction(condition, barcode_number),
    colour = colour
  )
) +
  facet_wrap(~gene.id, nrow = nrow, ncol = ncol) +
  theme_cowplot(8) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  geom_point(
    data = legend_df,
    aes(x = bin, y = proportion, colour = colour),
    size = 0, # invisible, just to trigger the legend
    show.legend = TRUE,
    inherit.aes = FALSE
  ) +
  scale_colour_identity(
    guide = "legend",
    labels = c(ref_sample, test_sample),
    breaks = c("#000000", "#8B0000"),
    name = "Condition"
  ) +
  labs(
    x = "Bin",
    y = "Proportion of reads",
  ) +
  geom_text(
    data = label_df,
    aes(
      x = 2,
      y = proportion_max,
      label = paste0(
        "dPSI: ",
        round(delta_PSI_mean, 2),
        " ± ",
        round(delta_PSI_SD, 2)
      )
    ),
    colour = "black",
    size = 2,
    hjust = 0.5,
    vjust = 0.5,
    inherit.aes = FALSE
  )

# Save plot to file
output_file <- paste0(
  "barcode_profiles_",
  Sys.Date(),
  "_",
  format(Sys.time(), "%H-%M-%S"),
  ".pdf"
)
ggsave(
  filename = output_file,
  plot = p,
  width = 2.2 * ncol,
  height = 2 * nrow,
)
