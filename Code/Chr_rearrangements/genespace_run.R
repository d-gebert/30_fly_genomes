# In Terminal: open -na rstudio

# Load genespace
library(GENESPACE)

# Set working directory and MCScanX path
wd <- "/Users/dgebert/genespace/all_species"
path2mcscanx <- "/usr/local/MCScanX"

# Organise species subgroups
subgroup1 <- c("Dmau", "Dsim", "Dsec", "Dmel", "Dtei", "Dyak", "Dere", "Dtak")
subgroup2 <- c("Dboc", "Djam", "Dkik", "Druf", "Dtri")
subgroup3 <- c("Dmal_pal", "Dmal_mal", "Dbip", "Dpar", "Dpse_pse", "Dana")
subgroup4 <- c("Dper", "Dpse", "Dsub")
subgroup5 <- c("Dequ", "Dpau", "Dwil", "Dtro", "Dins")
subgroup6 <- c("Dvir", "Dame", "Dlit")
all_specs <- c(subgroup1, subgroup2, subgroup3, subgroup4, subgroup5, subgroup6)

# Declare which species to run (all)
genomes2run <- all_specs

# Checking that everything is in order to run genespace proper
gpar <- init_genespace(wd = wd, path2mcscanx = path2mcscanx)

# Run genespace
out <- run_genespace(gpar, overwrite = TRUE)

# General appearance: white background
ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

# General appearance: chromosome colors
customPal <- colorRampPalette(
  c("#C00000", "#07A84E", "#003FA2", "#F6B200", "#7030A0", "#EE7500"))

# Riparian plots: all species
# Invert chromosomes with wrong orientation
invchr <- data.frame(
  genome = c("Dpse_pse", "Dana", "Dper", "Dpse"),
  chr = c("F", "F", "F", "F"))

# Plot riparian
ripDat <- plot_riparian(
  gsParam = out, 
  refGenome = "Dmel",
  genomeIDs = rev(all_specs),
  forceRecalcBlocks = FALSE,
  useOrder = FALSE,
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  minChrLen2plot = 1000000,
  invertTheseChrs = invchr,
  gapProp = 0.01, 
  scaleBraidGap = 1,
  scaleGapSize = .75,
  howSquare = 10,
  chrExpand = 0.5,
  customRefChrOrder = c("A","B","C","D","E","F"))

# Riparian plots: subgroup1
ripDat <- plot_riparian(
  gsParam = out, 
  refGenome = "Dmel",
  genomeIDs = rev(subgroup1),
  forceRecalcBlocks = FALSE,
  useOrder = FALSE,
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  minChrLen2plot = 1000000,
  gapProp = 0.01, 
  scaleBraidGap = 1,
  scaleGapSize = .75,
  howSquare = 10,
  chrExpand = 0.5,
  customRefChrOrder = c("A","B","C","D","E","F"))

# Riparian plots: subgroup2
ripDat <- plot_riparian(
  gsParam = out, 
  refGenome = "Djam",
  genomeIDs = rev(subgroup2),
  forceRecalcBlocks = FALSE,
  useOrder = FALSE,
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  minChrLen2plot = 1000000,
  gapProp = 0.01, 
  scaleBraidGap = 1,
  scaleGapSize = .75,
  howSquare = 10,
  chrExpand = 0.5,
  customRefChrOrder = c("A","B","C","D","E","F"))

# Riparian plots: subgroup3
invchr <- data.frame(
  genome = c("Dpse_pse", "Dana"),
  chr = c("F", "F"))
# Plot
ripDat <- plot_riparian(
  gsParam = out, 
  refGenome = "Dbip",
  genomeIDs = rev(subgroup3),
  forceRecalcBlocks = FALSE,
  useOrder = FALSE,
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  minChrLen2plot = 1000000,
  invertTheseChrs = invchr,
  gapProp = 0.01,
  scaleBraidGap = 1,
  scaleGapSize = .75,
  howSquare = 10,
  chrExpand = 0.5,
  customRefChrOrder = c("A","B","C","D","E","F"))

# Riparian plots: subgroup4
invchr <- data.frame(
  genome = c("Dper", "Dpse"), 
  chr = c("F", "F"))
# Plot
ripDat <- plot_riparian(
  gsParam = out, 
  refGenome = "Dsub",
  genomeIDs = rev(subgroup4),
  forceRecalcBlocks = FALSE,
  useOrder = FALSE,
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  minChrLen2plot = 1000000,
  invertTheseChrs = invchr,
  gapProp = 0.01, 
  scaleBraidGap = 1,
  scaleGapSize = .75,
  howSquare = 10,
  chrExpand = 0.5,
  customRefChrOrder = c("A","B","C","D","E","F"))

# General appearance: chromosome colors
customPal <- colorRampPalette(
  c("#C00000", "#F6B200", "#07A84E", "#003FA2", "#7030A0", "#EE7500")
)
# Riparian plots: subgroup4 A-D
invchr <- data.frame(
  genome = c("Dere", "Dper", "Dpse", "Dsub"),
  chr = c("D", "D", "D", "D")
)
# Plot
pdf("/Users/dgebert/genespace/FigS.pdf", width = 6, height = 6)
ripDat <- plot_riparian(
  gsParam = out,
  refGenome = "Dsub",
  genomeIDs = rev(c("Dere", "Dper", "Dpse", "Dsub")),
  forceRecalcBlocks = FALSE,
  useOrder = FALSE,
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  minChrLen2plot = 1000000,
  invertTheseChrs = invchr,
  gapProp = 0.01,
  scaleBraidGap = 1,
  scaleGapSize = .75,
  howSquare = 10,
  chrExpand = 0.5,
  customRefChrOrder = c("A", "D", "B", "C", "E", "F")
)
dev.off()

# General appearance: chromosome colors
customPal <- colorRampPalette(
  c("#C00000", "#07A84E", "#003FA2", "#F6B200", "#7030A0", "#EE7500")
)

# Riparian plots: subgroup5
ripDat <- plot_riparian(
  gsParam = out, 
  refGenome = "Dpau",
  genomeIDs = rev(subgroup5),
  forceRecalcBlocks = FALSE,
  useOrder = FALSE,
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  minChrLen2plot = 1000000,
  gapProp = 0.01, 
  scaleBraidGap = 1,
  scaleGapSize = .75,
  howSquare = 10,
  chrExpand = 0.5,
  customRefChrOrder = c("A","B","C","D","E","F"))

# Riparian plots: subgroup6
ripDat <- plot_riparian(
  gsParam = out, 
  refGenome = "Dvir",
  genomeIDs = rev(subgroup6),
  forceRecalcBlocks = FALSE,
  useOrder = FALSE,
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  minChrLen2plot = 1000000,
  gapProp = 0.01, 
  scaleBraidGap = 1,
  scaleGapSize = .75,
  howSquare = 10,
  chrExpand = 0.5,
  customRefChrOrder = c("A","B","C","D","E","F"))

# Region of interest: Muller A of Dmel, Dpib, Dvir
roi <- data.frame(
  genome = c("Dmel", "Dbip", "Dvir"),
  chr = c("A", "A", "A"), 
  color = c("#C00000", "#C00000", "#C00000"))

ripDat <- plot_riparian(
  gsParam = out, 
  highlightBed = roi,
  backgroundColor = NULL,
  refGenome = "Dmel",
  genomeIDs = rev(c("Dmel", "Dbip", "Dvir")),
  forceRecalcBlocks = FALSE,
  useOrder = FALSE,
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  minChrLen2plot = 1000000,
  gapProp = 0.01,
  scaleBraidGap = 1,
  scaleGapSize = .75,
  howSquare = 10,
  chrExpand = 0.5)

# Region of interest: Muller E of Dmel, Dpib, Dvir
roi <- data.frame(
  genome = c("Dsub", "Dlit", "Dvir"),
  chr = c("E", "E", "E"),
  color = c("#7030A0", "#7030A0", "#7030A0")
)

ripDat <- plot_riparian(
  gsParam = out,
  highlightBed = roi,
  backgroundColor = NULL,
  refGenome = "Dsub",
  genomeIDs = rev(c("Dsub", "Dlit", "Dvir")),
  forceRecalcBlocks = FALSE,
  useOrder = FALSE,
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  minChrLen2plot = 1000000,
  gapProp = 0.01,
  scaleBraidGap = 1,
  scaleGapSize = .75,
  howSquare = 10,
  chrExpand = 0.5
)


roi <- data.frame(
  genome = c("Dmel", "Dmel"),
  chr = c("E", "E"),
  start = c(7322683, 17301055),
  end = c(7570501, 17708850)
)
roibed <- roi[, c("genome", "chr", "start", "end")]
roibed$color <- c("#7030A0", "#7030C0")

ripd <- plot_riparian(
  refGenome = "Dmel",
  genomeIDs = rev(all_specs),
  gsParam = out,
  useOrder = FALSE,
  useRegions = TRUE,
  forceRecalcBlocks = FALSE,
  minChrLen2plot = 10000000,
  highlightBed = roibed,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  gapProp = 0.01,
  scaleBraidGap = 1,
  scaleGapSize = .75,
  howSquare = 10,
  chrExpand = 0.5,
  backgroundColor = "#F2F2F2",
  customRefChrOrder = rev(c("A", "B", "C", "D", "E"))
)