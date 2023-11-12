#' #' Title
#' #'
#' #' @param
#' #'
#  #' @importFrom tidyr gather
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' plot_parameter <- function(Parmat,Mean = NULL,titlename){
#'   if(is.na(Mean)) Mean = TRUE
#'   # Number of seizures
#'   nsei = dim(Parmat)[1]
#'   # Number of parameters
#'   npar = dim(Parmat)[2]
#'   ### construct data frame
#'   df = as.data.frame(Parmat)
#'   df$seizure <- 1:nsei
#'   DF <- gather(df, parameter, value, - seizure)
#'   DF$parameter = factor(DF$parameter,levels = colnames(df))
#'   # plots of estimated parameters
#'   col.rng =  rainbow_hcl(nsei,c = 80, l = 50)
#'   # x axis labels
#'   if(Mean) {
#'     lag = npar -2
#'     greek_labels <- c(sapply(1:lag, function(i) substitute(alpha[i], list(i = i))),
#'                       expression(beta), expression(sigma))
#'   } else {
#'     lag = npar - 1
#'     greek_labels <- c(sapply(1:lag, function(i) substitute(alpha[i], list(i = i))),
#'                       expression(sigma))
#'   }
#'   ### coefficient plot
#'   fig = ggplot(DF, aes(x = parameter, y = value, group = seizure, color = as.factor(seizure))) +
#'     geom_point(size = 3) +
#'     geom_line(data = subset(DF, as.integer(parameter) <= lag), size = 1) +
#'     labs(x = "Parameter", y = "Value", color = "Seizure") +
#'     scale_color_manual(values = col.rng) +
#'     scale_x_discrete(labels = greek_labels) +
#'     ggtitle(titlename)
#'   return(fig)
#' }
