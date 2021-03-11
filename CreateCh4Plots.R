# Plot output from Atlantis NEUS and QNM, BDN MSE
# Created 3/1/2021, Robert Wildermuth


library(tidyverse)


# Time series plots of NEUS output ----------------------------------------

# load historical and scenario outputs
load(file = "C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/omNEUS/neusHistoricData.RData")
load(file = "C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/omNEUS/neusWind.RData")
load(file = "C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/omNEUS/neusIncFish.RData")
load(file = "C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/omNEUS/neusDecFish.RData")
load(file = "C:/Users/rwildermuth/Dropbox/PhD_UMass/ATLANTIS/omNEUS/neusBaseline.RData")

neusHistoricData <- neusHistoricData %>% mutate(year = 1965:(1965+49),
                                                scenario = "Historical") 
neusBaseline <- neusBaseline %>% mutate(year = 1995:(1995+69),
                                        scenario = "Baseline") %>%
                  filter(year >= 2014)
neusDecFish <- neusDecFish %>% mutate(year = 1995:(1995+69),
                                      scenario = "DecFish") %>%
                  filter(year >= 2014)
neusIncFish <- neusIncFish %>% mutate(year = 1995:(1995+69),
                                      scenario = "IncFish") %>%
                  filter(year >= 2014)
neusWind <- neusWind %>% mutate(year = 1995:(1995+69),
                                scenario = "Energy") %>%
                  filter(year >= 2014)
neusScenarioData <- rbind(neusHistoricData, neusBaseline, neusDecFish, neusIncFish, neusWind)

# plot the time series for a subset of indices
neusScenarioData %>% select(-PS, -time) %>% #select(year, scenario, forage, inverts, ground, recParticip, habDem, habPel, strat, SST) %>%
  pivot_longer(cols = c(-year, -scenario), names_to = "node", values_to = "vals") %>%
  ggplot() + geom_line(aes(x=year,y=vals, color = scenario)) +
  facet_wrap(~node, nrow = 7, scales = "free") +
  theme_minimal()
