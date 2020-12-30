#gets run parameter list object for neus v1 run based on historical period
load_neus_v1_runprm <- function() {
  run <- list(
    toutstart = 1,
    toutinc = 365, #90, RW: changed to match NeusScenario3aa runs
    toutfinc = 365, #90,
    tstop = 54760, #18260,
    nyears = 50,
    timestep = 12,
    timestepunit = "Hours",
    outputstep = 365, #90,
    outputstepunit = "days",
    hemisphere = "northern",
    nspp = 67,
    nfleet = 33
  )
  return(run)
}

#run_prm <- load_neus_v1_runprm()