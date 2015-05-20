#' lmEq
#'
#' A function to generate the linear model equation and r-squared for labelling plots of 
#' two variables against each other - e.g. branch lengths recoverd from rj bayes traits.
#' @param df data frame with x and y as the x and y variables.
#' @return a string that can be used as a label for ggplot.
#' @export
#' @examples
#' lmEq(data.frame(dat, x = og, y = rj))

lmEq = function(df){
    m = lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(coef(m)[1], digits = 2), 
              b = format(coef(m)[2], digits = 2), 
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
}
