#####
library(ggplot2)
library(Seurat)
library(dplyr)
library(magrittr)
library(BiocParallel)
library(progressr)
library(spatstat)
library(sf)
library(ggsci)
library(pheatmap)
mol.type <- "microns"
coord.space <- "micron"
z.stack <- 3L

fig_dir = "/Users/ruoqing/Projects/aged_PLX3/figure/"
#####
allsection = c("/Volumes/T7/202308281442_Janos-aged-PLX3_VMSC03901/region_0/",
               "/Volumes/T7/202308281442_Janos-aged-PLX3_VMSC03901/region_1/",
               "/Volumes/T7/202308281442_Janos-aged-PLX3_VMSC03901/region_2/",
               "/Volumes/T7/202308281442_Janos-aged-PLX3_VMSC03901/region_3/",
               "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_0/",
               "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_1/",
               "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_2/",
               "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_3/",
               "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_4/",
               "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_5/",
               "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_6/",
               "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_7/",
               "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_8/",
               "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_9/",
               "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_0/",
               "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_1/",
               "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_2/",
               "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_3/",
               "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_4/",
               "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_5/",
               "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_6/",
               "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_7/",
               "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_8/",
               "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_9/")
sectionames = c("B3S0","B3S1","B3S2","B3S3",
                "B4S0","B4S1","B4S2","B4S3","B4S4","B4S5","B4S6","B4S7","B4S8","B4S9","B4S10",
                "B5S0","B5S1","B5S2","B5S3","B5S4","B5S5","B5S6","B5S7","B5S8","B5S9")
#####
B3S0 = LoadVizgen(data.dir = "/Volumes/T7/202308281442_Janos-aged-PLX3_VMSC03901/region_0/",
                  fov = "B3S0",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B3S0@meta.data$section = "B3S0"
B3S0@meta.data$cell = rownames(B3S0@meta.data)
B3S0@meta.data$keep = "no"
p1 = ImageDimPlot(B3S0, fov = "B3S0", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B3S0@meta.data$keep[B3S0@meta.data$cell %in% B3S0[["B3S0"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B3S0 = subset(B3S0, keep == "yes")

B3S1 = LoadVizgen(data.dir = "/Volumes/T7/202308281442_Janos-aged-PLX3_VMSC03901/region_1/",
                  fov = "B3S1",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B3S1@meta.data$section = "B3S1"
B3S1@meta.data$cell = rownames(B3S1@meta.data)
B3S1@meta.data$keep = "no"
p1 = ImageDimPlot(B3S1, fov = "B3S1", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B3S1@meta.data$keep[B3S1@meta.data$cell %in% B3S1[["B3S1"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B3S1 = subset(B3S1, keep == "yes")

B3S2 = LoadVizgen(data.dir = "/Volumes/T7/202308281442_Janos-aged-PLX3_VMSC03901/region_2/",
                  fov = "B3S2",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B3S2@meta.data$section = "B3S2"
B3S2@meta.data$cell = rownames(B3S2@meta.data)
B3S2@meta.data$keep = "no"
p1 = ImageDimPlot(B3S2, fov = "B3S2", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B3S2@meta.data$keep[B3S2@meta.data$cell %in% B3S2[["B3S2"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B3S2 = subset(B3S2, keep == "yes")

B3S3 = LoadVizgen(data.dir = "/Volumes/T7/202308281442_Janos-aged-PLX3_VMSC03901/region_3/",
                  fov = "B3S3",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B3S3@meta.data$section = "B3S3"
B3S3@meta.data$cell = rownames(B3S3@meta.data)
B3S3@meta.data$keep = "no"
p1 = ImageDimPlot(B3S3, fov = "B3S3", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B3S3@meta.data$keep[B3S3@meta.data$cell %in% B3S3[["B3S3"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B3S3 = subset(B3S3, keep == "yes")

B4S0 = LoadVizgen(data.dir = "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_0/",
                  fov = "B4S0",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B4S0@meta.data$section = "B4S0"
B4S0@meta.data$cell = rownames(B4S0@meta.data)
B4S0@meta.data$keep = "no"
p1 = ImageDimPlot(B4S0, fov = "B4S0", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B4S0@meta.data$keep[B4S0@meta.data$cell %in% B4S0[["B4S0"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B4S0 = subset(B4S0, keep == "yes")

B4S1 = LoadVizgen(data.dir = "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_1/",
                  fov = "B4S1",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B4S1@meta.data$section = "B4S1"
B4S1@meta.data$cell = rownames(B4S1@meta.data)
B4S1@meta.data$keep = "no"
p1 = ImageDimPlot(B4S1, fov = "B4S1", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B4S1@meta.data$keep[B4S1@meta.data$cell %in% B4S1[["B4S1"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B4S1 = subset(B4S1, keep == "yes")

B4S2 = LoadVizgen(data.dir = "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_2/",
                  fov = "B4S2",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B4S2@meta.data$section = "B4S2"
B4S2@meta.data$cell = rownames(B4S2@meta.data)
B4S2@meta.data$keep = "no"
p1 = ImageDimPlot(B4S2, fov = "B4S2", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B4S2@meta.data$keep[B4S2@meta.data$cell %in% B4S2[["B4S2"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B4S2 = subset(B4S2, keep == "yes")

B4S3 = LoadVizgen(data.dir = "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_3/",
                  fov = "B4S3",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B4S3@meta.data$section = "B4S3"
B4S3@meta.data$cell = rownames(B4S3@meta.data)
B4S3@meta.data$keep = "no"
p1 = ImageDimPlot(B4S3, fov = "B4S3", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B4S3@meta.data$keep[B4S3@meta.data$cell %in% B4S3[["B4S3"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B4S3 = subset(B4S3, keep == "yes")

B4S4 = LoadVizgen(data.dir = "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_4/",
                  fov = "B4S4",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B4S4@meta.data$section = "B4S4"
B4S4@meta.data$cell = rownames(B4S4@meta.data)
B4S4@meta.data$keep = "no"
p1 = ImageDimPlot(B4S4, fov = "B4S4", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B4S4@meta.data$keep[B4S4@meta.data$cell %in% B4S4[["B4S4"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B4S4 = subset(B4S4, keep == "yes")

B4S5 = LoadVizgen(data.dir = "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_5/",
                  fov = "B4S5",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B4S5@meta.data$section = "B4S5"
B4S5@meta.data$cell = rownames(B4S5@meta.data)
B4S5@meta.data$keep = "no"
p1 = ImageDimPlot(B4S5, fov = "B4S5", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B4S5@meta.data$keep[B4S5@meta.data$cell %in% B4S5[["B4S5"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B4S5 = subset(B4S5, keep == "yes")

B4S5 = LoadVizgen(data.dir = "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_5/",
                  fov = "B4S5",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B4S5@meta.data$section = "B4S5"
B4S5@meta.data$cell = rownames(B4S5@meta.data)
B4S5@meta.data$keep = "no"
p1 = ImageDimPlot(B4S5, fov = "B4S5", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B4S5@meta.data$keep[B4S5@meta.data$cell %in% B4S5[["B4S5"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B4S5 = subset(B4S5, keep == "yes")

B4S6 = LoadVizgen(data.dir = "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_6/",
                  fov = "B4S6",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B4S6@meta.data$section = "B4S6"
B4S6@meta.data$cell = rownames(B4S6@meta.data)
B4S6@meta.data$keep = "no"
p1 = ImageDimPlot(B4S6, fov = "B4S6", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B4S6@meta.data$keep[B4S6@meta.data$cell %in% B4S6[["B4S6"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B4S6 = subset(B4S6, keep == "yes")

B4S7 = LoadVizgen(data.dir = "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_7/",
                  fov = "B4S7",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B4S7@meta.data$section = "B4S7"
B4S7@meta.data$cell = rownames(B4S7@meta.data)
B4S7@meta.data$keep = "no"
p1 = ImageDimPlot(B4S7, fov = "B4S7", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B4S7@meta.data$keep[B4S7@meta.data$cell %in% B4S7[["B4S7"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B4S7 = subset(B4S7, keep == "yes")

B4S8 = LoadVizgen(data.dir = "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_7/",
                  fov = "B4S8",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B4S8@meta.data$section = "B4S8"
B4S8@meta.data$cell = rownames(B4S8@meta.data)
B4S8@meta.data$keep = "no"
p1 = ImageDimPlot(B4S8, fov = "B4S8", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B4S8@meta.data$keep[B4S8@meta.data$cell %in% B4S8[["B4S8"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B4S8 = subset(B4S8, keep == "yes")

B4S9 = LoadVizgen(data.dir = "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_8/",
                  fov = "B4S9",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B4S9@meta.data$section = "B4S9"
B4S9@meta.data$cell = rownames(B4S9@meta.data)
B4S9@meta.data$keep = "no"
p1 = ImageDimPlot(B4S9, fov = "B4S9", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B4S9@meta.data$keep[B4S9@meta.data$cell %in% B4S9[["B4S9"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B4S9 = subset(B4S9, keep == "yes")

B4S9 = LoadVizgen(data.dir = "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_8/",
                  fov = "B4S9",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B4S9@meta.data$section = "B4S9"
B4S9@meta.data$cell = rownames(B4S9@meta.data)
B4S9@meta.data$keep = "no"
p1 = ImageDimPlot(B4S9, fov = "B4S9", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B4S9@meta.data$keep[B4S9@meta.data$cell %in% B4S9[["B4S9"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B4S9 = subset(B4S9, keep == "yes")

B4S10 = LoadVizgen(data.dir = "/Volumes/T7/202309221159_Janos-agedPLX-4CP_VMSC03901/region_9/",
                  fov = "B4S10",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B4S10@meta.data$section = "B4S10"
B4S10@meta.data$cell = rownames(B4S10@meta.data)
B4S10@meta.data$keep = "no"
p1 = ImageDimPlot(B4S10, fov = "B4S10", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B4S10@meta.data$keep[B4S10@meta.data$cell %in% B4S10[["B4S10"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B4S10 = subset(B4S10, keep == "yes")

B5S0 = LoadVizgen(data.dir = "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_0/",
                   fov = "B5S0",
                   assay = "Vizgen",
                   metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                   mol.type = mol.type, 
                   type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                   z = z.stack,
                   add.zIndex = TRUE, # add z slice section to a cell
                   update.object = TRUE,
                   use.BiocParallel = TRUE,
                   workers.MulticoreParam = 14, # for `BiocParallel` processing
                   verbose = T)
B5S0@meta.data$section = "B5S0"
B5S0@meta.data$cell = rownames(B5S0@meta.data)
B5S0@meta.data$keep = "no"
p1 = ImageDimPlot(B5S0, fov = "B5S0", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B5S0@meta.data$keep[B5S0@meta.data$cell %in% B5S0[["B5S0"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B5S0 = subset(B5S0, keep == "yes")

B5S1 = LoadVizgen(data.dir = "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_1/",
                  fov = "B5S1",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B5S1@meta.data$section = "B5S1"
B5S1@meta.data$cell = rownames(B5S1@meta.data)
B5S1@meta.data$keep = "no"
p1 = ImageDimPlot(B5S1, fov = "B5S1", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B5S1@meta.data$keep[B5S1@meta.data$cell %in% B5S1[["B5S1"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B5S1 = subset(B5S1, keep == "yes")

B5S2 = LoadVizgen(data.dir = "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_2/",
                  fov = "B5S2",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B5S2@meta.data$section = "B5S2"
B5S2@meta.data$cell = rownames(B5S2@meta.data)
B5S2@meta.data$keep = "no"
p1 = ImageDimPlot(B5S2, fov = "B5S2", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B5S2@meta.data$keep[B5S2@meta.data$cell %in% B5S2[["B5S2"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B5S2 = subset(B5S2, keep == "yes")

B5S3 = LoadVizgen(data.dir = "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_3/",
                  fov = "B5S3",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B5S3@meta.data$section = "B5S3"
B5S3@meta.data$cell = rownames(B5S3@meta.data)
B5S3@meta.data$keep = "no"
p1 = ImageDimPlot(B5S3, fov = "B5S3", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B5S3@meta.data$keep[B5S3@meta.data$cell %in% B5S3[["B5S3"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B5S3 = subset(B5S3, keep == "yes")

B5S3 = LoadVizgen(data.dir = "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_3/",
                  fov = "B5S3",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B5S3@meta.data$section = "B5S3"
B5S3@meta.data$cell = rownames(B5S3@meta.data)
B5S3@meta.data$keep = "no"
p1 = ImageDimPlot(B5S3, fov = "B5S3", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B5S3@meta.data$keep[B5S3@meta.data$cell %in% B5S3[["B5S3"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B5S3 = subset(B5S3, keep == "yes")

B5S4 = LoadVizgen(data.dir = "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_4/",
                  fov = "B5S4",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B5S4@meta.data$section = "B5S4"
B5S4@meta.data$cell = rownames(B5S4@meta.data)
B5S4@meta.data$keep = "no"
p1 = ImageDimPlot(B5S4, fov = "B5S4", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B5S4@meta.data$keep[B5S4@meta.data$cell %in% B5S4[["B5S4"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B5S4 = subset(B5S4, keep == "yes")

B5S5 = LoadVizgen(data.dir = "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_5/",
                  fov = "B5S5",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B5S5@meta.data$section = "B5S5"
B5S5@meta.data$cell = rownames(B5S5@meta.data)
B5S5@meta.data$keep = "no"
p1 = ImageDimPlot(B5S5, fov = "B5S5", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B5S5@meta.data$keep[B5S5@meta.data$cell %in% B5S5[["B5S5"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B5S5 = subset(B5S5, keep == "yes")

B5S6 = LoadVizgen(data.dir = "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_6/",
                  fov = "B5S6",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B5S6@meta.data$section = "B5S6"
B5S6@meta.data$cell = rownames(B5S6@meta.data)
B5S6@meta.data$keep = "no"
p1 = ImageDimPlot(B5S6, fov = "B5S6", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B5S6@meta.data$keep[B5S6@meta.data$cell %in% B5S6[["B5S6"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B5S6 = subset(B5S6, keep == "yes")

B5S7 = LoadVizgen(data.dir = "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_7/",
                  fov = "B5S7",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B5S7@meta.data$section = "B5S7"
B5S7@meta.data$cell = rownames(B5S7@meta.data)
B5S7@meta.data$keep = "no"
p1 = ImageDimPlot(B5S7, fov = "B5S7", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B5S7@meta.data$keep[B5S7@meta.data$cell %in% B5S7[["B5S7"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B5S7 = subset(B5S7, keep == "yes")

B5S8 = LoadVizgen(data.dir = "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_8/",
                  fov = "B5S8",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B5S8@meta.data$section = "B5S8"
B5S8@meta.data$cell = rownames(B5S8@meta.data)
B5S8@meta.data$keep = "no"
p1 = ImageDimPlot(B5S8, fov = "B5S8", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B5S8@meta.data$keep[B5S8@meta.data$cell %in% B5S8[["B5S8"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B5S8 = subset(B5S8, keep == "yes")

B5S9 = LoadVizgen(data.dir = "/Volumes/T7/202309291111_Janos-aged-PLX5CP_VMSC03901/region_9/",
                  fov = "B5S9",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B5S9@meta.data$section = "B5S9"
B5S9@meta.data$cell = rownames(B5S9@meta.data)
B5S9@meta.data$keep = "no"
p1 = ImageDimPlot(B5S9, fov = "B5S9", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B5S9@meta.data$keep[B5S9@meta.data$cell %in% B5S9[["B5S9"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B5S9 = subset(B5S9, keep == "yes")

B6S0 = LoadVizgen(data.dir = "/Volumes/T7/202311151132_Janos-aged-PLX-longi-1cp_VMSC03901/region_0/",
                  fov = "B6S0",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B6S0@meta.data$section = "B6S0"
B6S0@meta.data$cell = rownames(B6S0@meta.data)
p1 = ImageDimPlot(B6S0, fov = "B6S0", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B6S0@meta.data$keep = "no"
B6S0@meta.data$keep[B6S0@meta.data$cell %in% B6S0[["B6S0"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B6S0@meta.data$keep[B6S0@meta.data$cell %in% B6S0[["B6S0"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "no"
B6S0 = subset(B6S0,keep == "yes")


B6S1 = LoadVizgen(data.dir = "/Volumes/T7/202311151132_Janos-aged-PLX-longi-1cp_VMSC03901/region_1/",
                  fov = "B6S1",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B6S1@meta.data$section = "B6S1"
B6S1@meta.data$cell = rownames(B6S1@meta.data)
B6S1 = subset(B6S1,nFeature_Vizgen > 1 & nCount_Vizgen > 1)
p1 = ImageDimPlot(B6S1, fov = "B6S1", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B6S1@meta.data$keep = "no"
B6S1@meta.data$keep[B6S1@meta.data$cell %in% B6S1[["B6S1"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B6S1@meta.data$keep[B6S1@meta.data$cell %in% B6S1[["B6S1"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "no"
B6S1 = subset(B6S1,keep == "yes")

B6S2 = LoadVizgen(data.dir = "/Volumes/T7/202311151132_Janos-aged-PLX-longi-1cp_VMSC03901/region_2/",
                  fov = "B6S2",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B6S2@meta.data$section = "B6S2"
B6S2@meta.data$cell = rownames(B6S2@meta.data)
B6S2 = subset(B6S2,nFeature_Vizgen > 1 & nCount_Vizgen > 1)
p1 = ImageDimPlot(B6S2, fov = "B6S2", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B6S2@meta.data$keep = "no"
B6S2@meta.data$keep[B6S2@meta.data$cell %in% B6S2[["B6S2"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B6S2@meta.data$keep[B6S2@meta.data$cell %in% B6S2[["B6S2"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "no"
B6S2 = subset(B6S2,keep == "yes")

B6S3 = LoadVizgen(data.dir = "/Volumes/T7/202311151132_Janos-aged-PLX-longi-1cp_VMSC03901/region_3/",
                  fov = "B6S3",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B6S3@meta.data$section = "B6S3"
B6S3@meta.data$cell = rownames(B6S3@meta.data)
B6S3 = subset(B6S3,nFeature_Vizgen > 1 & nCount_Vizgen > 1)
p1 = ImageDimPlot(B6S3, fov = "B6S3", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B6S3@meta.data$keep = "no"
B6S3@meta.data$keep[B6S3@meta.data$cell %in% B6S3[["B6S3"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B6S3 = subset(B6S3,keep == "yes")

B6S4 = LoadVizgen(data.dir = "/Volumes/T7/202311151132_Janos-aged-PLX-longi-1cp_VMSC03901/region_4/",
                  fov = "B6S4",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B6S4@meta.data$section = "B6S4"
B6S4@meta.data$cell = rownames(B6S4@meta.data)
B6S4 = subset(B6S4,nFeature_Vizgen > 1 & nCount_Vizgen > 1)
p1 = ImageDimPlot(B6S4, fov = "B6S4", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B6S4@meta.data$keep = "no"
B6S4@meta.data$keep[B6S4@meta.data$cell %in% B6S4[["B6S4"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B6S4@meta.data$keep[B6S4@meta.data$cell %in% B6S4[["B6S4"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "no"
B6S4 = subset(B6S4,keep == "yes")

B6S5 = LoadVizgen(data.dir = "/Volumes/T7/202311151132_Janos-aged-PLX-longi-1cp_VMSC03901/region_5/",
                  fov = "B6S5",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B6S5@meta.data$section = "B6S5"
B6S5@meta.data$cell = rownames(B6S5@meta.data)
B6S5 = subset(B6S5,nFeature_Vizgen > 1 & nCount_Vizgen > 1)
p1 = ImageDimPlot(B6S5, fov = "B6S5", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B6S5@meta.data$keep = "no"
B6S5@meta.data$keep[B6S5@meta.data$cell %in% B6S5[["B6S5"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B6S5@meta.data$keep[B6S5@meta.data$cell %in% B6S5[["B6S5"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "no"
B6S5 = subset(B6S5,keep == "yes")

B6S6 = LoadVizgen(data.dir = "/Volumes/T7/202311151132_Janos-aged-PLX-longi-1cp_VMSC03901/region_6/",
                  fov = "B6S6",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B6S6@meta.data$section = "B6S6"
B6S6@meta.data$cell = rownames(B6S6@meta.data)
B6S6 = subset(B6S6,nFeature_Vizgen > 1 & nCount_Vizgen > 1)
p1 = ImageDimPlot(B6S6, fov = "B6S6", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B6S6@meta.data$keep = "no"
B6S6@meta.data$keep[B6S6@meta.data$cell %in% B6S6[["B6S6"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B6S6@meta.data$keep[B6S6@meta.data$cell %in% B6S6[["B6S6"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "no"
B6S6 = subset(B6S6,keep == "yes")


B6S7 = LoadVizgen(data.dir = "/Volumes/T7/202311151132_Janos-aged-PLX-longi-1cp_VMSC03901/region_7/",
                  fov = "B6S7",
                  assay = "Vizgen",
                  metadata = c("volume", "fov","solidity","permimeter_area_ratio","DAPI_high_pass"," PolyT_high_pass"), 
                  mol.type = mol.type, 
                  type = c("segmentations", "centroids"), # type of cell spatial coord matrices
                  z = z.stack,
                  add.zIndex = TRUE, # add z slice section to a cell
                  update.object = TRUE,
                  use.BiocParallel = TRUE,
                  workers.MulticoreParam = 14, # for `BiocParallel` processing
                  verbose = T)
B6S7@meta.data$section = "B6S7"
B6S7@meta.data$cell = rownames(B6S7@meta.data)
B6S7 = subset(B6S7,nFeature_Vizgen > 1 & nCount_Vizgen > 1)
p1 = ImageDimPlot(B6S7, fov = "B6S7", size = 2, cols = "black",coord.fixed = TRUE,dark.background = F)
B6S7@meta.data$keep = "no"
B6S7@meta.data$keep[B6S7@meta.data$cell %in% B6S7[["B6S7"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "yes"
B6S7@meta.data$keep[B6S7@meta.data$cell %in% B6S7[["B6S7"]]$centroids@cells[as.numeric(CellSelector(p1))]] = "no"
B6S7 = subset(B6S7,keep == "yes")
#####
sections.seurat = merge(B3S1,
                        y = c(B3S2,B3S3,B4S2,B4S3,B4S5,B4S6,B4S9,B5S2,B5S7,B5S9,B6S0,B6S1,B6S2,B6S3,B6S4,B6S5,B6S6,B6S7),
                        add.cell.ids = c("B3S1","B3S2","B3S3","B4S2","B4S3","B4S5","B4S6","B4S9",
                                         "B5S2","B5S7","B5S9","B6S0","B6S1","B6S2","B6S3","B6S4","B6S5","B6S6","B6S7"),project = "opticnerve")
Idents(sections.seurat) = "section"
sectionames = c("B3S1","B3S2","B3S3","B4S2","B4S3","B4S5","B4S6","B4S9",
                "B5S2","B5S7","B5S9","B6S0","B6S1","B6S2","B6S3","B6S4","B6S5","B6S6","B6S7")
for (i in sectionames){
  DefaultBoundary(sections.seurat[[i]]) <- "centroids"
}
p1 = ImageDimPlot(sections.seurat, fov = sectionames, size = 2, cols = "black",border.size = 0,border.color = "black",
                  coord.fixed = TRUE,dark.background = F,crop = F)
ggsave("/Users/ruoqing/Projects/aged_PLX3/analysis345/alltogether2.pdf",p1,width = 15,height = 15,dpi = 300)

Idents(sections.seurat) = "orig.ident"
p1 = VlnPlot(sections.seurat,features = c("nCount_Vizgen", "nFeature_Vizgen","volume"))
ggsave("/Users/ruoqing/Projects/aged_PLX3/analysis345/VlnQC.pdf",p1,width = 6,height = 6,dpi = 300)

p1 = ggplot(data = sections.seurat@meta.data,aes(x = volume)) + 
  geom_histogram(color="black",fill= "white")+
  geom_vline(aes(xintercept=100),color="red", linetype="dashed", size=1)+
  geom_vline(aes(xintercept=2200),color="red", linetype="dashed", size=1)+
  theme_classic()
ggsave("/Users/ruoqing/Projects/aged_PLX3/analysis345/volume_cutoff.pdf",p1,width = 4,height = 3,dpi = 300)

sections.seurat.clean10 = subset(sections.seurat,nFeature_Vizgen > 10  & volume >= 100 & nCount_Vizgen <= 750)
p1 = VlnPlot(sections.seurat.clean10,features = c("nCount_Vizgen", "nFeature_Vizgen","volume"))
ggsave("/Users/ruoqing/Projects/aged_PLX3/analysis345/VlnQC10.pdf",p1,width = 6,height = 6,dpi = 300)

p1 = ImageDimPlot(sections.seurat.clean10, fov = sectionames, size = 2, cols = "black",border.size = 0,border.color = "black",
             coord.fixed = TRUE,dark.background = F,crop = F)
ggsave("/Users/ruoqing/Projects/aged_PLX3/analysis345/alltogethercutoff10.pdf",p1,width = 15,height = 15,dpi = 300)
for (i in sectionames){
  DefaultBoundary(sections.seurat.clean10[[i]]) <- "centroids"
}
p1 = ImageDimPlot(sections.seurat.clean10, fov = sectionames, size = 2, cols = "black",border.size = 0,border.color = "black",
                  coord.fixed = TRUE,dark.background = F,crop = F)
ggsave("/Users/ruoqing/Projects/aged_PLX3/analysis345/alltogethercutoff102.pdf",p1,width = 15,height = 15,dpi = 300)

p1 <- ImageFeaturePlot(sections.seurat, fov = sectionames,features = "nCount_Vizgen",size = 2,
                       border.color = "black",coord.fixed = TRUE,dark.background = F,crop = F)
ggsave("/Users/ruoqing/Projects/aged_PLX3/analysis345/alltogetherfeatnum.pdf",p1,width = 15,height = 15,dpi = 300)


sections.seurat.clean5 = subset(sections.seurat,nFeature_Vizgen >= 5 & volume >= 100 & nCount_Vizgen <= 750)
p1 = VlnPlot(sections.seurat.clean5,features = c("nCount_Vizgen", "nFeature_Vizgen","volume"))
ggsave("/Users/ruoqing/Projects/aged_PLX3/analysis345/VlnQC5.pdf",p1,width = 6,height = 6,dpi = 300)
Idents(sections.seurat.clean5) = "section"
p1 = ImageDimPlot(sections.seurat.clean5, fov = sectionames, size = 2, cols = "black",border.size = 0,border.color = "black",
                  coord.fixed = TRUE,dark.background = F,crop = F)
ggsave("/Users/ruoqing/Projects/aged_PLX3/analysis345/alltogethercutoff5.pdf",p1,width = 15,height = 15,dpi = 300)
for (i in c("B3S1","B3S2","B3S3",
            "B4S2","B4S3","B4S5","B4S6","B4S9",
            "B5S2","B5S7","B5S9")){
  DefaultBoundary(sections.seurat.clean5[[i]]) <- "segmentation"
}
p1 = ImageDimPlot(sections.seurat.clean5, fov = sectionames, size = 2, cols = "black",border.size = 0,border.color = "black",
                  coord.fixed = TRUE,dark.background = F,crop = F)
ggsave("/Users/ruoqing/Projects/aged_PLX3/analysis345/alltogethercutoff52.pdf",p1,width = 15,height = 15,dpi = 300)

#####

sections.seurat.clean5 <- SCTransform(sections.seurat.clean5, assay = "Vizgen", clip.range = c(-10, 10),
                                      variable.features.n = 350)
VariableFeaturePlot(sections.seurat.clean5)
sections.seurat.clean5 <- RunPCA(sections.seurat.clean5, npcs = 50, features = rownames(sections.seurat.clean5))
ElbowPlot(sections.seurat.clean5,ndims = 50)
sections.seurat.clean5 <- RunUMAP(sections.seurat.clean5, dims = 1:30)
sections.seurat.clean5 <- FindNeighbors(sections.seurat.clean5, reduction = "pca", dims = 1:30)
sections.seurat.clean5 <- FindClusters(sections.seurat.clean5, resolution = 8)

sections.seurat.clean5@meta.data$tmp = as.numeric(sections.seurat.clean5@meta.data$SCT_snn_res.0.8) - 1
sections.seurat.clean5@meta.data$cluster = "Oligo"
sections.seurat.clean5@meta.data$cluster[sections.seurat.clean5@meta.data$tmp %in% c(4,6,7,10,12,16)] = "Fibro/Mural"
sections.seurat.clean5@meta.data$cluster[sections.seurat.clean5@meta.data$tmp %in% c(0,5,15)] = "Astro"
sections.seurat.clean5@meta.data$cluster[sections.seurat.clean5@meta.data$tmp %in% c(8,14)] = "Micro"
sections.seurat.clean5@meta.data$cluster[sections.seurat.clean5@meta.data$tmp == 11] = "BAM"
sections.seurat.clean5@meta.data$cluster[sections.seurat.clean5@meta.data$tmp == 9] = "Endo"
sections.seurat.clean5@meta.data$cluster[sections.seurat.clean5@meta.data$tmp == 13] = "OPC"

sections.seurat.clean5@meta.data$group[sections.seurat.clean5@meta.data$section == "B4S2"] = "adult"
sections.seurat.clean5@meta.data$group[sections.seurat.clean5@meta.data$section == "B4S3"] = "aged"



sections.seurat.clean5@meta.data$tmp = as.numeric(sections.seurat.clean5@meta.data$SCT_snn_res.8) - 1
sections.seurat.clean5@meta.data$cluster[sections.seurat.clean5@meta.data$tmp == 60] = "T"
Idents(sections.seurat.clean5) = "cluster"

sections.seurat.clean5@meta.data$cd8 = as.numeric(sections.seurat.clean5@assays$SCT@data["Cd8a",])
sections.seurat.clean5@meta.data$cd3 = as.numeric(sections.seurat.clean5@assays$SCT@data["Cd3e",])

sections.seurat.clean5@meta.data$cluster[sections.seurat.clean5@meta.data$cd8 > 0 | sections.seurat.clean5@meta.data$cd3 > 0] = "T"
Idents(sections.seurat.clean5) = "cluster"

p1 = DimPlot(sections.seurat.clean5,label = T,cols = color)
ggsave(paste0(fig_dir,"dimplot.pdf"),p1,width = 6,height = 5,dpi = 300)


markers = FindAllMarkers(sections.seurat.clean5,logfc.threshold = 0)

top10 <- markers %>% group_by(cluster) %>% top_n(n = 6, wt = avg_log2FC)
use.gene = top10$gene

tmp = AverageExpression(sections.seurat.clean5,assays = "SCT",features = use.gene,return.seurat = T)
DoHeatmap(tmp,features = use.gene)

tmp = tmp@assays$SCT@scale.data

bk = seq(-2,2,0.1)
p1 = pheatmap(tmp,breaks = bk,
              color = colorRampPalette(c("navy","white", "firebrick3"))(length(bk)),
              cluster_rows = F,cluster_cols = F,border_color = NA)
ggsave(paste0(fig_dir,"marker_heatmap.pdf"),p1,width = 3,height = 8,dpi = 300)


sections.seurat.clean5@meta.data$group = "aged-PLX"
sections.seurat.clean5@meta.data$group[sections.seurat.clean5@meta.data$section %in% c("B4S3","B4S6","B5S2","B6S1","B6S2")] = "aged"
sections.seurat.clean5@meta.data$group[sections.seurat.clean5@meta.data$section %in% c("B4S2","B4S5","B5S7")] = "adult"

sections.seurat.clean5@meta.data$method = "cross"
sections.seurat.clean5@meta.data$method[sections.seurat.clean5@meta.data$section %in% c("B6S0","B6S1","B6S2","B6S3","B6S4","B6S5","B6S6","B6S7")] = "long"

color = ggsci::pal_npg()(8)
names(color) = c("Endo","OPC","BAM","T","Astro","Fibro/Mural","Micro","Oligo")


group.color = randomcoloR::distinctColorPalette(3)
names(group.color) = c("aged-PLX","aged","adult")

p1 = ggplot(sections.seurat.clean5@meta.data,aes(x = group, fill = cluster))+
  geom_bar(stat = "count",position = "fill")+
  theme_classic()+
  scale_y_continuous(expand=c(0,0))+
  xlab("")+ylab("")+
  scale_fill_manual(values = color, breaks = names(color))
ggsave(paste0(fig_dir,"barplot_group.pdf"),p1,width = 4,height = 4,dpi = 300)

p1 = ggplot(sections.seurat.clean5@meta.data[sections.seurat.clean5@meta.data$method == "cross",],aes(x = group, fill = cluster))+
  geom_bar(stat = "count",position = "fill")+
  theme_classic()+
  scale_y_continuous(expand=c(0,0))+
  xlab("")+ylab("")+
  scale_fill_manual(values = color, breaks = names(color))
ggsave(paste0(fig_dir,"barplot_group_cross_section.pdf"),p1,width = 4,height = 4,dpi = 300)

score = read.table("/Users/ruoqing/Projects/aged_PLX3/score.txt",header = F)
score$V2 = paste0(score$V2,"_",score$V3)
score = score[,c(1,2)]
colnames(score) = c("gene","module")
exclude_gene = c("Lpl","Fabp5","Timp1","Psmb8","Cd63","Cxcl11","Ccl6","Cxcl5","Il1a","Il1b","Tlr2","Serping1")
score = score[score$gene %in% rownames(sections.seurat.clean5@assays$SCT@counts),]
score = score[score$gene %in% setdiff(score$gene,exclude_gene),]
score = score[score$module %in% setdiff(unique(score$module),"Cytokine_response"),]
score$tmp = paste0(score$gene,"_",score$module)
score = score[score$tmp %in% setdiff(unique(score$tmp),"H2-D1_Final_Astrocytes"),]
score = score[,c(1,2)]
score$module[score$gene == "Ccl3"] = "Final_microglia"
score$module[score$gene %in% c("Cxcl10","Cxcl14")] = "Final_Astrocytes"

score$gene[score$gene == "H2-D1" & score$module == "Final_microglia"] = "H2-K1"
score$gene[score$gene == "Cd74" & score$module == "Final_microglia"] = "Cdkn1a"
score$gene[score$gene == "Lyz2" & score$module == "Final_microglia"] = "Cd68"
score$module[score$gene == "Clec7a" & score$module == "Final_microglia"] = "remove"

score = score[score$module != "remove",]

score$gene[score$gene == "Hspb1" & score$module == "Final_Astrocytes"] = "Hif1a"
tmp = data.frame(gene = c("Stat3","Isg15","H2-K1"),module = rep("Final_Astrocytes",3))
score = rbind(score,tmp)

score$gene[score$gene == "Lyz2" & score$module == "Final_Oligo"] = "H2-K1"

geneSets <- lapply(unique(score$module), function(x){print(x);score$gene[score$module == x]})
names(geneSets) <- unique(score$module)

tmp = AddModuleScore(sections.seurat.clean5,features = geneSets,ctrl = 10,assay = "SCT",name = names(geneSets))

tmp = tmp@meta.data[,c("cluster","group","Final_microglia1","Final_Astrocytes2","Final_Oligo3")]
colnames(tmp) = c("cluster","group","Micro Activation","Astro Activation","Oligo Activation")
tmp = tmp[tmp$cluster %in% c("Astro","Oligo","Micro"),]
tmp = reshape2::melt(tmp)
tmp$new = paste0(tmp$cluster,"_",tmp$variable)
tmp = tmp[tmp$new %in% c("Micro_Micro Activation","Astro_Astro Activation","Oligo_Oligo Activation"),]
tmp$variable = factor(tmp$variable,levels = c("Oligo Activation","Micro Activation","Astro Activation"),ordered = T)
p1 = ggplot(tmp,aes(x = group, y = value,fill = cluster))+
  geom_violin(linetype="blank",scale = "width")+
  facet_grid(variable~. ,scales = 'free_y')+
  scale_fill_manual(values = as.character(color),breaks = names(color))+
  theme_classic()+xlab("")+ylab("")+
  theme(panel.grid = element_blank(),
        panel.background = element_blank())
ggsave(paste0(fig_dir,"activescore_update.pdf"),p1,width = 5,height = 4,dpi = 300)

### if sig
tmp = AddModuleScore(sections.seurat.clean5,features = geneSets,ctrl = 10,assay = "SCT",name = names(geneSets))

tmp = tmp@meta.data[,c("cluster","group","Final_microglia1","Final_Astrocytes2","Final_Oligo3")]
colnames(tmp) = c("cluster","group","Micro Activation","Astro Activation","Oligo Activation")
tmp = tmp[tmp$cluster %in% c("Astro","Oligo","Micro"),]
tmp = reshape2::melt(tmp)
tmp$new = paste0(tmp$cluster,"_",tmp$variable)
tmp = tmp[tmp$new %in% c("Micro_Micro Activation","Astro_Astro Activation","Oligo_Oligo Activation"),]

library(ggpubr)
p1 = ggplot(tmp[tmp$cluster == "Oligo",],aes(x = group, y = value, fill = group))+
  geom_boxplot()+
  stat_compare_means(comparisons = list(c("adult","aged"),c("adult","aged-PLX"),c("aged-PLX","aged")),method = "wilcox.test",label="p.signif")+
  theme_classic()+theme(legend.position = "none")+ggtitle("Oligo")
p2 = ggplot(tmp[tmp$cluster == "Oligo",],aes(x = group, y = value, fill = group))+
  geom_boxplot()+
  stat_compare_means(comparisons = list(c("adult","aged"),c("adult","aged-PLX"),c("aged-PLX","aged")),method = "wilcox.test")+
  theme_classic()+theme(legend.position = "none")+ggtitle("Oligo")
p1+p2

tmp = sections.seurat.clean5@meta.data[,c("cluster","group","Cytokine_response")]
tmp = tmp[tmp$cluster %in% c("Astro","Oligo","Micro","OPC","T" ),]
tmp = reshape2::melt(tmp)
p1 = ggplot(tmp[tmp$value > 0,],aes(x = group, y = value,fill = cluster))+
  geom_violin(linetype="blank",scale = "width")+
  facet_grid(cluster~. ,scales = 'free_y')+
  scale_fill_manual(values = as.character(color),breaks = names(color))+
  theme_classic()+xlab("")+ylab("Cytokine_response")+
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = 'none')
ggsave(paste0(fig_dir,"Cytokine_response_score.pdf"),p1,width = 4,height = 4,dpi = 300)

for (i in sectionames){
  DefaultBoundary(sections.seurat.clean5[[i]]) <- "segmentation"
}

Idents(sections.seurat.clean5) = "cluster"
for (i in unique(sections.seurat.clean5@meta.data$section[sections.seurat.clean5@meta.data$method == "cross"])){
  p1 = ImageDimPlot(sections.seurat.clean5, fov = i, size = 2, border.size = 0,border.color = "grey66",
                    coord.fixed = TRUE,dark.background = F,crop = F,cols = color)
  ggsave(paste0(fig_dir,"sections/",i,"_",unique(sections.seurat.clean5@meta.data$group[sections.seurat.clean5@meta.data$section == i]),".pdf"),p1,width = 8,height = 7,dpi = 300)
}

for (i in unique(sections.seurat.clean5@meta.data$section[sections.seurat.clean5@meta.data$method == "long"])){
  p1 = ImageDimPlot(sections.seurat.clean5, fov = i, size = 2, border.size = 0,border.color = "grey66",
                    coord.fixed = TRUE,dark.background = F,crop = F,cols = color)
  ggsave(paste0(fig_dir,"sections/",i,"_",unique(sections.seurat.clean5@meta.data$group[sections.seurat.clean5@meta.data$section == i]),".pdf"),p1,width = 10,height = 6,dpi = 300)
}

for (i in c("B6S1","B6S2")){
  p1 = ImageDimPlot(sections.seurat.clean5, fov = i, size = 2, border.size = 0,border.color = "grey66",
                    coord.fixed = TRUE,dark.background = F,crop = F,cols = color)
  ggsave(paste0(fig_dir,"sections/",i,"_",unique(sections.seurat.clean5@meta.data$group[sections.seurat.clean5@meta.data$section == i]),".pdf"),p1,width = 6,height = 10,dpi = 300)
}

color2 = c(c("grey30","grey60","grey90"))
names(color2) = c("aged-PLX","aged","adult")

tmp = matrix(NA,ncol = length(unique(sections.seurat.clean5@meta.data$cluster)),nrow = length(unique(sections.seurat.clean5@meta.data$section[sections.seurat.clean5@meta.data$method == "cross"])))
colnames(tmp) = unique(sections.seurat.clean5@meta.data$cluster)
rownames(tmp) = unique(unique(sections.seurat.clean5@meta.data$section[sections.seurat.clean5@meta.data$method == "cross"]))

for (i in colnames(tmp)){
  for (j in rownames(tmp)){
    tmp[j,i] = nrow(sections.seurat.clean5@meta.data[sections.seurat.clean5@meta.data$cluster == i & sections.seurat.clean5@meta.data$section == j,])/nrow(sections.seurat.clean5@meta.data[sections.seurat.clean5@meta.data$section == j,])
  }
}
tmp = as.data.frame(tmp)
tmp$section = rownames(tmp)
write.table(tmp,paste0(fig_dir,"clu_composition.txt"),quote = F,row.names = F)

tmp$group = "aged-PLX"
tmp$group[tmp$section %in% c("B4S3","B4S6","B5S2","B6S1","B6S2")] = "aged"
tmp$group[tmp$section %in% c("B4S2","B4S5","B5S7")] = "adult"
tmp = reshape2::melt(tmp)
p1 = ggplot(tmp)+
  geom_boxplot(aes(x = variable,y = value, fill = group), position = position_dodge(1),outlier.shape = NA)+
  scale_fill_manual(values = color2,breaks = names(color2))+
  theme_classic()+
  xlab("")+ylab("")
ggsave(paste0(fig_dir,"boxplot_cross_section.pdf"),p1,width = 8,height = 5,dpi = 300)


p1 = ImageFeaturePlot(sections.seurat.clean5,fov = "B3S3",features = c("Cxcl10"),dark.background = F)
ggsave(paste0(fig_dir,"b3s3_gene_expcxcl10.pdf"),p1,width = 6,height = 5,dpi = 300)
p1 = ImageFeaturePlot(sections.seurat.clean5,fov = "B3S3",features = c("Lgals3"),dark.background = F)
ggsave(paste0(fig_dir,"b3s3_gene_expLgals3.pdf"),p1,width = 6,height = 5,dpi = 300)

p1 = ImageFeaturePlot(sections.seurat.clean5,fov = "B3S1",features = c("Cxcl10"),dark.background = F)
ggsave(paste0(fig_dir,"b3s1_gene_expcxcl10.pdf"),p1,width = 6,height = 5,dpi = 300)
p1 = ImageFeaturePlot(sections.seurat.clean5,fov = "B3S1",features = c("Lgals3"),dark.background = F)
ggsave(paste0(fig_dir,"b3s1_gene_expLgals3.pdf"),p1,width = 6,height = 5,dpi = 300)

#####
##neighbor 
t_section = unique(sections.seurat.clean5@meta.data$section[sections.seurat.clean5@meta.data$cluster == "T"])
t_neighbors = c()

for (i in t_section) {
  tmp = sections.seurat.clean5@images[[i]]$centroids@coords
  tmp = as.data.frame(tmp)
  rownames(tmp) = sections.seurat.clean5@images[[i]]$centroids@cells
  
  tmp1 = BiocNeighbors::findNeighbors(tmp,threshold = 50)
  tmp$cell = rownames(tmp)
  tmp$cluster = "s"
  
  for (j in tmp$cell) {
    tmp$cluster[tmp$cell == j] = sections.seurat.clean5@meta.data$cluster[sections.seurat.clean5@meta.data$cell == j]
  }
  
  rownames(tmp) = NULL
  
  tmp$neighbor = "non-neighbor"
  tmp$neighbor[unlist(tmp1$index[as.numeric(rownames(tmp[tmp$cluster == "T",]))])] = "neighbor"
  tmp$neighbor[tmp$cluster == "T"] = "non-neighbor"
  
  t_neighbors = c(t_neighbors,tmp$cell[tmp$neighbor == "neighbor"])
}

sections.seurat.clean5@meta.data$t_neighbor = "non-neighbor"
sections.seurat.clean5@meta.data$t_neighbor[sections.seurat.clean5@meta.data$cell %in% t_neighbors] = "neighbor"
sections.seurat.clean5@meta.data$t_neighbor[sections.seurat.clean5@meta.data$cluster == "T"] = "T"

Idents(sections.seurat.clean5) = "t_neighbor"

tmp = CreateSeuratObject(sections.seurat.clean5@assays$SCT@counts,meta.data = sections.seurat.clean5@meta.data)
tmp <- NormalizeData(tmp)
Idents(tmp) = "t_neighbor"



s = FindMarkers(tmp,ident.1 = "neighbor",ident.2 = "non-neighbor",min.cells.group = 1, 
                min.cells.feature = 1,
                min.pct = 0,
                logfc.threshold = 0,
                only.pos = FALSE)
s$sig = "non"
s$sig[s$avg_log2FC >= 0.3 & s$p_val_adj <= 0.1] = "neighbor"
s$sig[s$avg_log2FC <= -0.3 & s$p_val_adj <= 0.1] = "non_neighbor"
s$gene = rownames(s)
ggplot(s,aes(x = avg_log2FC,y = -log10(p_val_adj),color = sig,label = gene))+
  geom_point()+
  geom_text_repel(data = s[s$sig %in% c("neighbor","non_neighbor"),])+
  theme_bw()+
  theme(legend.position = "none")



microglia_section = unique(sections.seurat.clean5@meta.data$section[sections.seurat.clean5@meta.data$cluster == "Micro"])
micro_neighbors = c()

for (i in microglia_section) {
  tmp = sections.seurat.clean5@images[[i]]$centroids@coords
  tmp = as.data.frame(tmp)
  rownames(tmp) = sections.seurat.clean5@images[[i]]$centroids@cells
  
  tmp1 = BiocNeighbors::findNeighbors(tmp,threshold = 50)
  tmp$cell = rownames(tmp)
  tmp$cluster = "s"
  
  for (j in tmp$cell) {
    tmp$cluster[tmp$cell == j] = sections.seurat.clean5@meta.data$cluster[sections.seurat.clean5@meta.data$cell == j]
  }
  
  rownames(tmp) = NULL
  
  tmp$neighbor = "non-neighbor"
  tmp$neighbor[unlist(tmp1$index[as.numeric(rownames(tmp[tmp$cluster == "T",]))])] = "neighbor"
  tmp$neighbor[tmp$cluster == "T"] = "non-neighbor"
  
  micro_neighbors = c(micro_neighbors,tmp$cell[tmp$neighbor == "neighbor"])
}

sections.seurat.clean5@meta.data$micro_neighbors = "non-neighbor"
sections.seurat.clean5@meta.data$micro_neighbors[sections.seurat.clean5@meta.data$cell %in% micro_neighbors] = "neighbor"
sections.seurat.clean5@meta.data$micro_neighbors[sections.seurat.clean5@meta.data$cluster == "Micro"] = "Micro"

Idents(sections.seurat.clean5) = "micro_neighbors"

tmp = CreateSeuratObject(sections.seurat.clean5@assays$SCT@counts,meta.data = sections.seurat.clean5@meta.data)
tmp <- NormalizeData(tmp)
Idents(tmp) = "micro_neighbors"

micro.deg = FindMarkers(tmp,ident.1 = "neighbor",ident.2 = "non-neighbor",min.cells.group = 1, 
                        min.cells.feature = 1,
                        min.pct = 0,
                        logfc.threshold = 0,
                        only.pos = FALSE)

micro.deg$sig = "non"
micro.deg$sig[micro.deg$avg_log2FC >= 0.3 & micro.deg$p_val_adj <= 0.1] = "neighbor"
micro.deg$sig[micro.deg$avg_log2FC <= -0.3 & micro.deg$p_val_adj <= 0.1] = "non_neighbor"
micro.deg$gene = rownames(micro.deg)
pmac = ggplot(micro.deg,aes(x = avg_log2FC,y = -log10(p_val_adj),color = sig,label = gene))+
  geom_point()+
  geom_text_repel(data = micro.deg[micro.deg$sig %in% c("neighbor","non_neighbor"),])+
  scale_color_manual(breaks = c("non","neighbor","non_neighbor"),values = c("grey66","#C53300","navy"))+
  geom_vline(xintercept = 0.3,linetype = 2)+
  geom_vline(xintercept = -0.3,linetype = 2)+
  geom_hline(yintercept = -log10(0.1),linetype = 2)+
  annotate(geom="text", x= 1, y=11, label="Micro neighbor",  fontface="bold",colour='#C53300', size=5)+
  annotate(geom="text", x= -1, y=11, label="non-neighbor",  fontface="bold",colour='navy', size=5)+
  theme_bw()+theme_classic()+
  theme(legend.position = "none")
ggsave(paste0(fig_dir,"micro_neighbor.pdf"),pmac,width = 6,height = 6,dpi = 300)

s = FindMarkers(tmp,ident.1 = "neighbor",ident.2 = "non-neighbor",min.cells.group = 1, 
                min.cells.feature = 1,
                min.pct = 0,
                logfc.threshold = 0,
                only.pos = FALSE)
s$sig = "non"
s$sig[s$avg_log2FC >= 0.3 & s$p_val_adj <= 0.1] = "neighbor"
s$sig[s$avg_log2FC <= -0.3 & s$p_val_adj <= 0.1] = "non_neighbor"
s$gene = rownames(s)
pt = ggplot(s,aes(x = avg_log2FC,y = -log10(p_val_adj),color = sig,label = gene))+
  geom_point()+
  geom_text_repel(data = s[s$sig %in% c("neighbor","non_neighbor"),])+
  scale_color_manual(breaks = c("non","neighbor","non_neighbor"),values = c("grey66","#C53300","navy"))+
  geom_vline(xintercept = 0.3,linetype = 2)+
  geom_vline(xintercept = -0.3,linetype = 2)+
  geom_hline(yintercept = -log10(0.1),linetype = 2)+
  annotate(geom="text", x= 1, y=11, label="T neighbor",  fontface="bold",colour='#C53300', size=5)+
  annotate(geom="text", x= -1, y=11, label="non-neighbor",  fontface="bold",colour='navy', size=5)+
  theme_bw()+theme_classic()+
  theme(legend.position = "none")
ggsave(paste0(fig_dir,"t_neighbor.pdf"),pt,width = 6,height = 6,dpi = 300)


####sccoda

tmp = matrix(NA,ncol = length(unique(sections.seurat.clean5@meta.data$cluster)),nrow = length(unique(sections.seurat.clean5@meta.data$section[sections.seurat.clean5@meta.data$method == "cross"])))
colnames(tmp) = unique(sections.seurat.clean5@meta.data$cluster)
rownames(tmp) = unique(unique(sections.seurat.clean5@meta.data$section[sections.seurat.clean5@meta.data$method == "cross"]))

for (i in colnames(tmp)){
  for (j in rownames(tmp)){
    tmp[j,i] = nrow(sections.seurat.clean5@meta.data[sections.seurat.clean5@meta.data$cluster == i & sections.seurat.clean5@meta.data$section == j,])
  }
}

write.table(tmp,"/Users/ruoqing/Projects/aged_PLX3/analysis345/sccoda.txt",quote = F)
