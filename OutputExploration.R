library(tidync)
library(sp)
library(rgdal)
library(rgeos)
library(rbgm)

# Install latest atlantisom
devtools::install_github("r4atlantis/atlantisom")
library(atlantisom)

# Try loading with atlantisom functions/functionality

# code from load_fgs()
result <- read.table(file = "//storage1.smast.umassd.edu/lab_fay/rwildermuth/ATLANTIS/scenarios/Oct23_08d/functionalGroups.prm", 
                     sep = " ",
                     header = TRUE, stringsAsFactors = FALSE)

load_bps(dir = "//storage1.smast.umassd.edu/lab_fay/rwildermuth/ATLANTIS/scenarios/Oct23_08d/", 
         fgs = "functionalGroups.prm", 
         file_init = "inneus_2007.nc")

boxes <- get_boundary(load_box(dir = "//storage1.smast.umassd.edu/lab_fay/rwildermuth/ATLANTIS/scenarios/Oct23_08d/", 
                               file_bgm = "neus30_2006.bgm"))

load_nc_physics(dir = "//storage1.smast.umassd.edu/lab_fay/rwildermuth/ATLANTIS/scenarios/Oct23_08d/", 
                file_nc = "Salt_neus.nc",
                physic_variables = "salt",
                aggregate_layers = TRUE,
                bboxes = c(5, 7, 8, 9, 10, 11, 19))
                
# Try from Sarah's notes
# taken from existing load_nc

file.nc <- "//storage1.smast.umassd.edu/lab_fay/rwildermuth/ATLANTIS/scenarios/Oct23_08d/neusDynEffort_Oct23_08d_.nc"

# Load ATLANTIS output!
at_out <- RNetCDF::open.nc(con = file.nc)

# Get info from netcdf file! (Filestructure and all variable names)
var_names_ncdf <- sapply(seq_len(RNetCDF::file.inq.nc(at_out)$nvars - 1),
                         function(x) RNetCDF::var.inq.nc(at_out, x)$name)
n_timesteps <- RNetCDF::dim.inq.nc(at_out, 0)$length
n_boxes     <- RNetCDF::dim.inq.nc(at_out, 1)$length
n_layers    <- RNetCDF::dim.inq.nc(at_out, 2)$length

RNetCDF::close.nc(at_out)

# check that needed groups are vars in output
varNames <- c("Diatom",
              "DinoFlag",
              "PicoPhytopl",
              "Lab_Det",
              "Ref_Det",
              "Pisciv_S_Fish",
              "Demersal_DC_Fish",
              "Macrobenth_Shallow",
              "Benthic_grazer",
              "Benthic_Carniv",
              "Deposit_Feeder",
              "Zoo",
              "MicroZoo",
              "Gelat_Zoo",
              "Pinniped",
              "Whale_Baleen",
              "Pisciv_D_Fish",
              "Pisciv_B_Fish",
              "Demersal_D_Fish",
              "Demersal_S_Fish",
              "Demersal_B_Fish",
              "Demersal_DC_Fish",
              "Demersal_O_Fish",
              "Demersal_F_Fish",
              "Shark_B",
              "Shark_D",
              "SkateRay",
              "Filter_Shallow",
              "Filter_Other",
              "Megazoobenthos",
              "Prawn",
              "Planktiv_L_Fish",
              "Planktiv_S_Fish",
              "Demersal_E_Fish",
              "Cephalopod",
              "Benthopel_Fish")

var_names_ncdf[grepl("Lab_Det", var_names_ncdf, fixed = TRUE)]

sapply(varNames, 1, FUN = grep, fixed = TRUE, value = TRUE, x = var_names_ncdf)
#---------------------------------------------------------
prodConnect <- tidync("//storage1.smast.umassd.edu/lab_fay/rwildermuth/ATLANTIS/scenarios/Oct23_08d/neusDynEffort_Oct23_08d_PROD.nc")

prodConnect %>% hyper_filter()
prodConnect %>% activate(DepositFeedProdn) %>% hyper_tibble()
prodConnect %>% activate(nominal_dz)
prodConnect %>% activate(topk)

totConnect <- tidync("//storage1.smast.umassd.edu/lab_fay/rwildermuth/ATLANTIS/scenarios/Oct23_08d/neusDynEffort_Oct23_08d_TOT.nc")

totConnect %>% hyper_filter()
test2 <- totConnect %>% activate("D1,D0") %>% hyper_vars()
test2$name

test5 <- totConnect %>% activate(Gelat_Zoo_N)
test5 %>% hyper_vars()
test5 %>% hyper_transforms()
test5 %>% hyper_dims()
test6 <- test5 %>% hyper_array()
str(test6)
test6$Gelat_Zoo_N

catchConnect <- tidync("//storage1.smast.umassd.edu/lab_fay/rwildermuth/ATLANTIS/scenarios/Oct23_08d/neusDynEffort_Oct23_08d_CATCH.nc")

catchConnect %>% hyper_filter()
test3 <- catchConnect %>% activate("D1,D0") %>% hyper_vars()
test3$name

saltConnect <- tidync("//storage1.smast.umassd.edu/lab_fay/rwildermuth/ATLANTIS/scenarios/Oct23_08d/Salt_neus.nc")
saltConnect %>% hyper_filter()
saltConnect %>% hyper_vars()
# See if I can extract info for density calcs for Box 10 (center of GB)
saltConnect %>% hyper_filter(b = index == 10, z = index < 3) %>%
  hyper_tibble() # long-form table with salinities at 730 time points, not sure what time units are

# Find bottom salinity vals
test7 <- saltConnect %>% hyper_filter(b = dplyr::between(index, 5, 11), 
                             z = index == 5) %>% hyper_tibble()
any(is.na(test7$salinity))


initConnect <- tidync("//storage1.smast.umassd.edu/lab_fay/rwildermuth/ATLANTIS/scenarios/Oct23_08d/inneus_2007.nc")
initConnect %>% hyper_filter()
test8 <- initConnect %>% activate(numlayers) %>% hyper_tibble()
test8$numlayers
# look at initialization specifications for benthic/epifauna
# Macroalgae_N, Seagrass_N, Filter_Shallow_N, Filter_Deep_N, Filter_Other_N, Macrobenth_Shallow_N, 
# Macrobenth_Deep_N, Megazoobenthos_N, Benthic_grazer_N, Macroalgae_Cover, Seagrass_Cover, 
# MicroPB_Cover, Filter_Shallow_Cover, Filter_Deep_Cover, Filter_Other_Cover
initConnect %>% activate(Seagrass_Cover) %>% hyper_tibble()
initConnect %>% activate(Macroalgae_Cover) %>% hyper_tibble()
initConnect %>% activate(MicroPB_Cover) %>% hyper_tibble()
initConnect %>% activate(Filter_Shallow_Cover) %>% hyper_tibble() # has stuff!
initConnect %>% activate(Filter_Deep_Cover) %>% hyper_tibble()
initConnect %>% activate(Filter_Other_Cover) %>% hyper_tibble() # has stuff!

# look at other benthic groups
initConnect %>% activate(Benthic_grazer_N) %>% hyper_tibble()
initConnect %>% activate(Megazoobenthos_N) %>% hyper_tibble()
initConnect %>% activate(Macrobenth_Shallow_N) %>% hyper_tibble()
initConnect %>% activate(Deposit_Feeder_N) %>% hyper_tibble()
initConnect %>% activate(Benthic_Carniv_N) %>% hyper_tibble()

test9 <- totConnect %>% activate(Planktiv_L_Fish_N) %>% hyper_tibble()
test9 %>% filter(Planktiv_L_Fish_N > 0)      
                   
# get areas of each box in m^2
# First modified the functionalGroups.csv from new model base to work with old output files
file_fgs = "NeusGroups_modforusewNEUS1.0.csv"
# function args
scenario = "neusDynEffort_Oct23_08d_"
dir = "//storage1.smast.umassd.edu/lab_fay/rwildermuth/ATLANTIS/scenarios/Oct23_08d"
file_bgm = "neus30_2006.bgm" # geography bgm for box definitions
load_boxarea(dir = dir, file_bgm = file_bgm)

#   polygon        area
# 6        5 11085441644
# 8        7 10363076990
# 9        8  6455559428
# 10       9 17308896405
# 11      10 11219585917
# 12      11 15998216129
# 20      19 17736767097

# Check overlap of lease areas with Atlantis boxes and GB EPU
shpEPU <- readOGR("C:/Users/rwildermuth/Dropbox/PhD_UMass/Drafts and Outlines/ImagesShapefiles/Extended_EPU")
shpGB <- subset(shpEPU, shpEPU$EPU %in% "GB")
windLeases <- readOGR("C:/Users/rwildermuth/Dropbox/PhD_UMass/Drafts and Outlines/ImagesShapefiles/BOEM-Renewable-Energy-Shapefiles",
                      layer = "BOEM_Lease_Areas_04_10_2017")
proj4string(windLeases) # should be: proj4string = CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0")

plot(windLeases, col = 4)
plot(shpGB, add = TRUE)

neusBoxes <- bgmfile(x = paste(dir, file_bgm, sep = "/")) # bgm package doesn't work for old codebase?

# 8 current (10/21/2020 from northeastoceandataportal.org) active leases in Boxes 5 and 7:
# OCS-A 0486 - DWW Rev I, LLC
# OCS-A 0487 - Deepwater Wind New England LLC
# OCS-A 0500 - Bay State Wind LLC
# OCS-A 0501 - Vineyard Wind LLC
# OCS-A 0517 - Deepwater Wind South Fork, LLC
# OCS-A 0520 - Equinor Wind US LLC
# OCS-A 0521 - Mayflower Wind Energy LLC
# OCS-A 0522 - Vineyard Wind LLC
# (plus Block Island turbines)
n <- 8

# sq km of nearshore (Box 5) habitat affected by wind 
nearAffected <- (17.8 * 0.2 + 0.5 + 0.1 + 10.5)*n 

# Assume 21.9 sq km of deep seafloor habitat affected
seafloorAffected <- 21.9

# area of box 5
area5 <- 11085441644 * 1e-6
# Are of box 7
area7 <- 10363076990 * 1e-6

# proportion of Box 5/nearshore affected
nearAffected/area5
# proportion of Box 7 affected
seafloorAffected/area7
