#'
#' branchChangeSim
#' Simulate the expected change along a branch of length branchtime.
#' @export
#' @name branchChangeSim

branchChangeSim <- function(branchtime = 100, direction = "non", startval = 0, var = 0.2, increment = 1) { 
  # Create the output array, starting with the starting value
  simout <- array(dim = c(1,2), data = c(0, startval))
  varsimout <- NULL

  # Identify the time intervals within which we wish to simulate
  timeinc <- seq(from = 0, to = branchtime, by = increment) # divide the time up into increments
  if(max(timeinc) < branchtime) {
    timeinc <- c(timeinc, branchtime) # add on however long is left at the end. 
  } 

  # Loop through each of these time increments and simulate
  for (inc in 2:length(timeinc)) {
    multiplier <- timeinc[inc] - timeinc[inc - 1] # how much time happens along this increment?
    varsim  <- rnorm(n = 1, mean = 0, sqrt(var)) * multiplier # draw a random value from the distribution defined by var 
    # multiply by the amount of time that has occurred
    # Get the simulated value...
    if(direction == "neg") {
      simval <- simout[inc-1,2] - abs(varsim) 
    } else if(direction == "pos") {
      simval <- simout[inc-1,2] + abs(varsim) 
    } else if(direction == "non") {
      simval <- simout[inc-1,2] + varsim
    }

    # Add it to the table
    simout <- rbind(simout,c(timeinc[inc], simval))
  }
  return(simout)
}
