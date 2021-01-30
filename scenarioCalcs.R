# Figure out a daily mortality rate for simulating effect of wind farms

library(tidync)
library(tidyverse)


# annual mortality factor
annMort <- 0.05

# daily mortalty
annMort/365

# Atlantis fishing mortality probability based on Sec 15.3.1 in Atlantis manual (pg 14)
mFC <- 1-exp(-annMort/365)


# Check that it's working

# load biomasses for the base scenario
biomassConnect <- tidync("C:/Users/rwildermuth/Documents/AtlantisTesting/5yrTest/NeusScenario3aa000_10daybase/neusDynEffort_Oct23_08d_.nc")

biomassConnect %>% hyper_filter()
baseDepFeedN <- biomassConnect %>% 
                    activate(what = "D2,D1,D0") %>% 
                    hyper_tibble() %>%
                    dplyr::select(t, z, b, Deposit_Feeder_N, Benthic_Carniv_N)
unique(baseDepFeedN$t)

baseFiltersN <- biomassConnect %>% 
                  activate(what = "D1,D0") %>% 
                  hyper_tibble() %>%
                  dplyr::select(t, b, Filter_Shallow_N, Filter_Other_N, 
                         Macrobenth_Shallow_N, Megazoobenthos_N, Benthic_grazer_N)

# load biomasses for the mortality scenario
biomassConnect <- tidync("C:/Users/rwildermuth/Documents/AtlantisTesting/5yrTest/NeusScenario3aa000_10daymort/neusDynEffort_Oct23_08d_.nc")

biomassConnect %>% hyper_filter()
mortDepFeedN <- biomassConnect %>% 
  activate(what = "D2,D1,D0") %>% 
  hyper_tibble() %>%
  dplyr::select(t, z, b, Deposit_Feeder_N, Benthic_Carniv_N)
unique(mortDepFeedN$t)

mortFiltersN <- biomassConnect %>% 
  activate(what = "D1,D0") %>% 
  hyper_tibble() %>%
  dplyr::select(t, b, Filter_Shallow_N, Filter_Other_N, 
         Macrobenth_Shallow_N, Megazoobenthos_N, Benthic_grazer_N)

test1 <- merge(baseDepFeedN, mortDepFeedN, by = c("t", "z", "b"))
test1 %>% filter(b %in% c(5,7,8,9,10,11,19), z == 5) %>% arrange(t)
all(test1$Benthic_Carniv_N.x == test1$Benthic_Carniv_N.y)


test2 <- merge(baseFiltersN, mortFiltersN, by = c("t", "b"))
test3 <- test2 %>% filter(b == 5) %>% arrange(t) %>% select(t, b, Filter_Shallow_N.x, Filter_Shallow_N.y)

test3$Filter_Shallow_N.x == test3$Filter_Shallow_N.y

#--------------------------------------------------
# Check that 1994 scenarios initialized at same level of 1994 projected from historical time series
# load biomasses for the base scenario
initBiom1994Connect <- tidync("C:/Users/rwildermuth/Documents/AtlantisTesting/FishingScenarios/NeusScenario3aa000/new_inneus_1994start_fixed_effort.nc")
initBiom1994N <- initBiom1994Connect %>% 
  activate(what = "D2,D1,D0") %>% 
  hyper_tibble() #%>%
#  dplyr::select(t, z, b)

initTime <- unique(initBiom1994N$t)[1]

biomassConnect <- tidync("C:/Users/rwildermuth/Documents/AtlantisTesting/FishingScenarios/NeusScenario3aaHISTORICAL/neusDynEffort_Oct23_08d_.nc")
baseHistN <- biomassConnect %>% 
  activate(what = "D2,D1,D0") %>% 
  hyper_tibble() 

baseHistN %>% dplyr::select(t, z, b)

histTest <- baseHistN %>% filter(t == initTime)
initTest <- initBiom1994N %>% filter(t == initTime) %>% dplyr::select(names(histTest))
dim(merge(histTest, initTest, all.x = TRUE))
