# Testing modifications to {atlantisom} for reading in old NEUS codebase
# Created: 10/21/2020, Robert Wildermuth

# use NEFSC's version for now
library(atlantisom)

source("C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/load_neus_v1_runprm.R")

# First modified the functionalGroups.csv from new model base to work with old output files
file_fgs = "NeusGroups_modforusewNEUS1.0.csv"

#---------------------------
# Test run_truth() with age info commented out

# function args
scenario = "neusDynEffort_Oct23_08d_"
dir = "//storage1.smast.umassd.edu/lab_fay/rwildermuth/ATLANTIS/scenarios/Oct23_08d"
file_bgm = "neus30_2006.bgm" # geography bgm for box definitions
select_groups = c("Diatom", # the group 'Name's in the 'file_fgs' file
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
                  "Shark_D","SkateRay",
                  "Filter_Shallow",
                  "Filter_Other",
                  "Megazoobenthos",
                  "Prawn",
                  "Planktiv_L_Fish",
                  "Planktiv_S_Fish",
                  "Demersal_E_Fish",
                  "Cephalopod",
                  "Benthopel_Fish")
# try shorter group list
select_groups = c("Pinniped",
                  "Whale_Baleen",
                  "Pisciv_D_Fish",
                  "Pisciv_B_Fish",
                  "Demersal_D_Fish",
                  "Demersal_S_Fish")
file_init = "inneus_2007.nc" # initial condition .nc file
file_biolprm = "at_biol_neus_Oct10B_Jason_FE.prm" # -b flag in the batch file
file_runprm = "Run_settings.xml" # needed for catch correction calculation, how to replace? -r flag in batch file
#file_fish = "Fisheries.csv" # used in age-specific results, comment out
verbose = FALSE
save = FALSE
annage = FALSE

# Try using modified run_truth()
test1 <- run_truth(scenario = scenario,
                   dir = dir,
                   file_fgs = "NeusGroups_modforusewNEUS1.0.csv",
                            file_bgm = "neus30_2006.bgm",
                            select_groups = select_groups,
                            file_init = "inneus_2007.nc",
                            file_biolprm = "at_biol_neus_Oct10B_Jason_FE.prm",
                            file_runprm = "NEUSv1",
                            save = FALSE)
# get 'catch_all' in terms of group names
fxnlGroups <- load_fgs(dir, file_fgs = "NeusGroups_modforusewNEUS1.0.csv")

test2 <- merge(test1$catch_all, fxnlGroups[, c("Code", "Name")], by.x = "species", by.y = "Code")

test1$catch_all$species <- test2$Name

# subset to groups in 'select_groups'
test1$catch_all <- test1$catch_all %>% dplyr::filter(species %in% select_groups)

compare_catchB <-ggplot() +
  geom_line(data=test1$bio_catch, aes(x=time,y=atoutput, #color="catch atlantisom"
                                      ), 
            alpha = 10/10) +
  geom_point(data=test1$catch_all, aes(x=time,y=atoutput, #color="txt output catch bio"
                                       ), # has all spp in Atlantis codes, not Names
             alpha = 1/10) + 
  theme_minimal() +
  theme(legend.position = "top") 
  #labs(colour=scenario.name)

compare_catchB + 
  facet_wrap(~species, ncol=3, nrow = 2, scales="free")

# need to sum among 'bio_catch' values between time intervals in 'catch_all'
test1$bio_catch$timeAgg <- 0

for(i in 2:nrow(test1$catch_all)){
  for(j in 1:nrow(test1$bio_catch)){
    if(test1$bio_catch$time[j] > test1$catch_all$time[i-1] & 
       test1$bio_catch$time[j] <= test1$catch_all$time[i]){
      test1$bio_catch$timeAgg[j] <- test1$catch_all$time[i]
    }
  }
}

test3 <- test1$bio_catch %>% group_by(species, timeAgg) %>%
            dplyr::summarize(catchAgg = sum(atoutput),
                             n = n())

check <- merge(test1$catch_all, test3,
                 by.x = c("species", "time"), by.y = c("species", "timeAgg"))
check$check <- with(check, atoutput / catchAgg)

check %>% group_by(time) %>%
  dplyr::summarize(minScale = min(check),
                   maxScale = max(check))
check %>% pivot_wider(id_cols = time, names_from = species, values_from = check)

# Read in the biomass pools
bps <- load_bps(dir = dir, fgs = file_fgs, file_init = file_init)
#[1] "Filter_Shallow"     "Filter_Other"       "Filter_Deep"        "Benthic_grazer"     "Macrobenth_Deep"   
#[6] "Megazoobenthos"     "Macrobenth_Shallow" "Macroalgae"         "Seagrass"

# Figure out time calculation in load_nc()
file.nc <- file.path(dir, paste0(scenario, "CATCH.nc"))
at_out <- RNetCDF::open.nc(con = file.nc)
n_timesteps <- RNetCDF::dim.inq.nc(at_out, 0)$length # load_nc() constructs time series based on length of time dimension
# for NEUS v1 this length is 203, but catch.txt only has 52

# look at example files in package
d <- system.file("extdata", "SETAS_Example", package = "atlantisom")
file.nc <- file.path(dir = d, file_nc = "outputsCATCH.nc")
at_out <- RNetCDF::open.nc(con = file.nc)
n_timesteps <- RNetCDF::dim.inq.nc(at_out, 0)$length

# Run example and look at time lengths
d <- system.file("extdata", "SETAS_Example", package = "atlantisom")
groups <- load_fgs(dir = d, "Functional_groups.csv")
truth <- run_truth(scenario = "outputs",
                   dir = d,
                   file_fgs = "Functional_groups.csv",
                   file_bgm = "Geography.bgm",
                   select_groups = groups[groups$IsTurnedOn > 0, "Name"],
                   file_init = "Initial_condition.nc",
                   file_biolprm = "Biology.prm",
                   file_runprm = "Run_settings.xml",
                   save = FALSE)
str(truth)
unique(truth$catch$time)
unique(truth$catch_all$time)
