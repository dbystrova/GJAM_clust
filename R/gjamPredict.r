

gjamPredict <- function(output, newdata = NULL, y2plot = NULL, ylim = NULL,
                        FULL = FALSE){
  
  # output   - gjamGibbs object
  # newdata  - list with xdata and y needed for prediction
  # y2plot   - names of y for plotting
  # if y is included has only columns to condition on 
  
  PLOT <- T
  if(is.null(y2plot))PLOT <- F
  
  invisible( .gjamPrediction( output, newdata, y2plot, PLOT, ylim, FULL ) )
  
}
