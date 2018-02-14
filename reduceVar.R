reduceVar <- function (x, breaks) 
{
        if (is.numeric(x) | is.integer(x)) {
                if (breaks=="quantile") {
                        breaks <- quantile(x)
                }
        }
              
        if (length(breaks) > 1) 
                x <- cut(x, breaks = breaks, include.lowest = TRUE, 
                         labels = FALSE)
        return(list(x = x, breaks = breaks))
}