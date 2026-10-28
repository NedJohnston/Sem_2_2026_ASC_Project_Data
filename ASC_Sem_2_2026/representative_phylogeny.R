
###### Set Up #######

getwd()
setwd("C:/Users/u6867511/Documents/Phylogenetics")


install.packages("pacman")

# Use pacman to install (if required) and load libraries

pacman::p_load(
  rio,             # import/export
  here,            # relative file paths
  tidyverse,       # general data management and visualization
  ape,             # to import and export phylogenetic files
  ggtree,          # to visualize phylogenetic files
  treeio,          # to visualize phylogenetic files
  ggnewscale,      # to add additional layers of color schemes
  msa,
  readxl,
  Biostrings,
  DECIPHER,
  dplyr,
  ggplot2
)

BiocManager::install("phangorn")
BiocManager::install("msa")

library(phangorn)
library(msa)

####### Begin analysis #######

# Import metadata
data <- read_csv("representative_phylogeny.csv")
View(data)
nrow(data)

# Take IDs of QC-filtered sequences
data_ids <- data$`Specimen ID`

# Import sequence data (FASTA format)
seq <- read.dna("Data/NED_Plate1to6.fasta", format = "fasta")
seq

# Keep only QC-filtered sequences
seq_filt <- seq[names(seq) %in% data_ids]


# Save filtered file
write.dna(seq_filt,
          file = "NED_Plate1to6_QC.fasta",
          format = "fasta",
          colsep = "")

# Orient sequences the same way using DECIPHER package
# This step probably not necessary
seq_ss <- readDNAStringSet("NED_Plate1to6_QC.fasta")
seq_or <- OrientNucleotides(seq_ss)

# Convert to dataframe
seq_or_df <- data.frame(
  ID = names(seq_or),
  Sequence = as.character(seq_or)
)

View(seq_or_df)

# Save oriented sequences
write_csv(seq_or_df, "NED_oriented.csv")

# Align sequences
aln <- msa(seq_or, method = "Muscle")
aligned <- as(aln, "DNAStringSet")

# Convert to dataframe
aln_df <- data.frame(
  ID = names(aligned),
  Sequence = as.character(aligned)
)

View(aln_df)

# Save aligned sequences
write_csv(aln_df, "NED_aligned.csv")


# Conclusion -- I think they were already oriented the same way
# Aligning only fixes one sample (NED0000198) which is missing a base pair at the beginning
# But probably good to do it anyway, to double check


###### Build distance matrix #######

# Choose a substitution model:
# Raw p-distance (raw)
# Kimura 2-parameter (K2P, K81)
# Jukes-Cantor (JC)
# Felsenstein's 1981 model (F81)

# Test different models using modelTest

seq_pd <- as.phyDat(aln)
seq_pd

mt <- modelTest(seq_pd, model = c("JC", "F81", "K81", "raw"))

# Better model fit = Lower AIC and BIC values, logLik closer to 0

# Sort results by AIC values
mt[order(mt$AIC), ][1:12, ]

# Best model = K81+G(4)+I
# Kimura 3-parameter model
# +G allows gamma-distributed rate variation
# +I allows invariant sites

# Convert aligned sequences back to DNA bin object
seq_db <- as.DNAbin(aligned)

# View model options in dist.dna function
?dist.dna

# Make distance matrix using K81 model
dist <- dist.dna(seq_db, model = "K81")


###### Build phylogenetic tree from the distance matrix ########

# Read more about different tree types, pros and cons of each, etc.

# Options:

# Neighbour Joining tree (NJ)

tree_nj <- nj(dist)

# Save as nwk file
# Can be later imported into iTOL
write.tree(tree_nj, "subset_NED_K81_NJ_tree.nwk")


# Maximum Likelihood tree (ML)

# Using phagorn package

# Convert to phyDat object (phagorn requires this)
phy <- phyDat(as.matrix(seq_db), type = "DNA")

# Make ML distance matrix
dm <- dist.ml(phy)

# Make NJ tree to start with
tree_nj_ml <- NJ(dm)

# Fit ML model
fit <- pml(tree_nj_ml, data = phy)

# Optimise model
# Chose Kimura 2-parameter tree, allowing gamma-dist and invariant sites
# Based on previous modelTest results

fit_K81 <- optim.pml(
  fit,
  model = "K81",
  optInv = TRUE,
  optGamma = TRUE,
  rearrangement = "stochastic"
)
fit_K81

# Save tree
write.tree(
  fit_K81$tree,
  file = "NED_ML_tree.nwk"
)

####### iTOL annotation: labels with MOTU + colourstrip from "colour" column ########


# Map colour names to hex codes
colour_map <- c(
  black = "#000000",
  brown = "#8B4513",   # saddle brown
  pale  = "#F0DA97"    # pale brown (tan)
)

# Check that every colour in the csv has a hex code
stopifnot(all(data$colour %in% names(colour_map)))

# Data frame with Specimen ID, label (e.g. NED0000001 (MOTU 1)) and hex colour
itol_labels <- data.frame(
  specimen_id = data$`Specimen ID`,
  label = paste0(data$`Specimen ID`, " (MOTU ", data$motu, ")"),
  colour_name = data$colour,
  hex = unname(colour_map[data$colour])
)
head(itol_labels)

# 1 ADDS COLOUR/SYMBOL FOR POLLINATOR
key_motus <- c(1, 2, 4, 5, 7, 8, 10, 17, 51, 52)

# Check all key MOTUs have a representative in the csv
setdiff(key_motus, data$motu)   # should return numeric(0)

itol_labels$key <- itol_labels$specimen_id %in%
  data$`Specimen ID`[data$motu %in% key_motus]
sum(itol_labels$key)   # should equal length(key_motus)

key_rows <- itol_labels[itol_labels$key, ]

# Colour for pollinators
key_hex <- "#D62728"   # red; change if you like

writeLines(
  c(
    "TREE_COLORS",
    "SEPARATOR COMMA",
    "DATA",
    paste(key_rows$specimen_id, "label", key_hex, "bold", "1.2", sep = ",")
  ),
  "subset_NED_key_MOTU_text.txt"
)

# Star for pollinators

writeLines(
  c(
    "DATASET_BINARY",
    "SEPARATOR COMMA",
    "DATASET_LABEL,Key MOTU",
    paste0("COLOR,", key_hex),
    "FIELD_SHAPES,3",            # 3 = star
    "FIELD_LABELS,Key MOTU",
    paste0("FIELD_COLORS,", key_hex),
    "DATA",
    paste(key_rows$specimen_id, "1", sep = ",")
  ),
  "subset_NED_key_MOTU_star.txt"
)

# 1. LABELS file: changes the text shown at each tip
writeLines(
  c(
    "LABELS",
    "SEPARATOR COMMA",
    "DATA",
    paste(itol_labels$specimen_id, itol_labels$label, sep = ",")
  ),
  "subset_NED_MOTU_labels.txt"
)

# 2. COLORSTRIP file: coloured strip beside each tip
# Format: ID,colour,label (the label shows when hovering over the strip)
filename <- "subset_NED_colourstrip.txt"

writeLines(
  c(
    "DATASET_COLORSTRIP",
    "SEPARATOR COMMA",
    "DATASET_LABEL,Colour",
    "COLOR,#000000",
    "LEGEND_TITLE,Colour",
    "LEGEND_SHAPES,1,1,1",
    paste0("LEGEND_COLORS,", paste(colour_map, collapse = ",")),
    paste0("LEGEND_LABELS,", paste(names(colour_map), collapse = ",")),
    "DATA"
  ),
  filename
)

write.table(
  itol_labels[, c("specimen_id", "hex", "colour_name")],
  file = filename,
  sep = ",",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  append = TRUE
)

# Upload both files to iTOL (drag and drop onto the tree)