
#' @title Bangdiwala Chart
#'
#' @description A function for creating a Bangdiwala Chart.
#'    Constructed using "The agreement chart" by Shrikant I Bangdiwala.
#'    To learn more about Bangdiwala Chart, see
#'    \url{http://www.stat.ncsu.edu/information/library/mimeo.archive/ISMS_1988_1859.pdf}, and
#'    book: \url{http://www.utdallas.edu/~pkc022000/agreement_book/} chapter 12: Categorical Data.
#'
#' @param x data as a table object in R.
#'    note: columns of x will turn into x-axis,
#'    rows of x will turn into y-axis,
#'    to reverse, use \code{t()} function (i.e., transpose).
#' @param main title of plot.
#' @param xlab title of x-axis (default: extracted from x).
#' @param ylab title of y-axis (default: extracted from x).
#' @param cat.label category labels (default: extracted from x).
#' @param las y-axis category labels rotated.
#' @param add.legend `TRUE/FALSE`, add a legend to the plot.
#' @param location.legend either `bottomright` or `topleft` (only relevant if `add.legend=TRUE`).
#' @param add.B.hat `TRUE/FALSE`, add Bangdiwala agreement measure to plot,
#'     (otherwise it is only given as the function's output),
#'     (agreement measure without weights).
#' @param location.B.hat  either `bottomright` or `topleft` (only relevant if `add.B.hat=TRUE`).
#' @param add.x.margin `TRUE/FALSE`, add marginal values for columns (x) to plot.
#' @param add.y.margin `TRUE/FALSE`, add marginal values for columns (y) to plot.
#' @param ordinal.data `TRUE/FALSE`, if `TRUE`, partial agreement squares are added.
#' @param col.line color of line (default: `firebrick`).
#' @param lwd.line width of line (default: `2`).
#' @param lty.line type of line (`solid`, `dotted`, etc.; default: `solid`).
#'
#' @import grDevices
#' @import graphics
#' @return A chart, and Bangdiwala agreement measure
#' @export
#' @references Bangdiwala, S. I. (1988). The Agreement Chart. Department of Biostatistics, University of North Carolina at Chapel Hill, Institute of Statistics Mimeo Series No. 1859, http://www.stat.ncsu.edu/information/library/mimeo.archive/ISMS_1988_1859.pdf
#'
#'      David Meyer, Achim Zeileis, and Kurt Hornik (2017). vcd: Visualizing Categorical Data. R package version 1.4-4.
#'
#'      Choudhary, P. K. and Nagaraja, H. N. (2017). Measuring Agreement: Models, Methods, and Applications, Wiley.
#'
#'      Westlund, K. B. and Kurland, L. T. (1953). Studies on multiple sclerosis in Winnipeg, Manitoba, and New Orleans, Louisiana I. Prevalence; comparison between the patient groups in Winnipeg and New Orleans. American Journal of Hygiene 57, 380–396.
bangdiwala_chart <- function(x, main="",
                             xlab=names(dimnames(x))[2],
                             ylab= names(dimnames(x))[1],
                             cat.label=colnames(x),
                             las=TRUE,
                             add.legend=TRUE,
                             location.legend="bottomright",
                             add.B.hat=TRUE	,
                             location.B.hat="topleft",
                             add.x.margin=TRUE,
                             add.y.margin=TRUE,
                             ordinal.data=TRUE,
                             col.line="grey"
                             ,
                             lwd.line=2,
                             lty.line=1){

  #####
  # 0. preliminary checks and calculations

  # ensure data is of the correct type:
  if(any(!is.numeric(x))) return("ERROR: some cells have non-numeric counts")
  if(any(is.na(x))) return("ERROR: some cells have missing counts")
  if(any(x < 0))  return("ERROR: some cells have negative counts")

  # ensure data is in the form of a table:
  if(!is.table(x)) x <- as.table(x)

  n <- sum(x)
  k <- length(cat.label)

  #####
  # 1. draw n x n square

  plot(1:n, 1:n, col="white", main=main, xlab=xlab, ylab=ylab,
       xlim=c(0,n), ylim=c(0,n), xaxs="i", yaxs="i", xaxt="n", yaxt="n")

  #####
  # 2. add axis labels
  x.margin <- c(0, cumsum(as.vector(apply(x, 2, sum))))
  y.margin <- c(0, cumsum(as.vector(apply(x, 1, sum))))
  # axis labels should be at the center of each rectangle in step 3
  axis(side=1, at=(x.margin[1:k] + x.margin[2:(k+1)])/2, labels=cat.label)
  axis(side=2, at=(y.margin[1:k] + y.margin[2:(k+1)])/2, labels=cat.label, las=las)

  #####
  # 3. k rectangles with width = column marginals, height = row marginals
  for(i in 1:k){
    rect(xleft=x.margin[i], ybottom=y.margin[i], xright=x.margin[i+1], ytop=y.margin[i+1])
  } # end for loop
  rm(i)

  #####
  # 4. partial agreement squares Using Bangdiwala partial agreement definition (1988, p 11)
  #    -- add only if your rating categories are ordinal

  if(ordinal.data){

    # darker gray --> closer to full agreement
    gray.vector <- gray.colors(n=k-2, start=0.3, end=0.6)

    for(j in rev(1:(k-2))){ # loop through colors

      # within +/- j of correct rating (make larger rectangles first)
      for(i in 1:k){ # loop through rating values

        # summation vector
        temp <- sort((i-j):(i+j), decreasing=FALSE)
        temp <- temp[ (temp > 0) & (temp <= k) ]

        partial.agreement.x <- sum(as.vector(x[temp,i]))
        partial.agreement.y <- sum(as.vector(x[i,temp]))

        if(temp[1]==1){
          add.x <- 0
          add.y <- 0
        }else{
          add.x <- sum(as.vector(x[1:(temp[1]-1),i]))
          add.y <- sum(as.vector(x[i,1:(temp[1]-1)]))
        } # end if/else

        rect(xleft=x.margin[i]+add.x, ybottom=y.margin[i]+add.y,
             xright=x.margin[i]+add.x+partial.agreement.x,
             ytop=y.margin[i]+add.y+partial.agreement.y,
             col=gray.vector[j])

        rm(temp, partial.agreement.x, partial.agreement.y, add.x, add.y)
      } # end inner for loop

    } # end outer for loop
  } # end if for ordinal.data

  #####
  # 5. full agreement squares (in black)
  full.agreement <- as.vector(diag(x))
  for(i in 1:k){
    if(i==1){
      add.x <- 0
      add.y <- 0
    }else{
      add.x <- sum(as.vector(x[1:(i-1),i]))
      add.y <- sum(as.vector(x[i,1:(i-1)]))
    } # end if/else

    rect(xleft=x.margin[i]+add.x, ybottom=y.margin[i]+add.y,
         xright=x.margin[i]+add.x+full.agreement[i], ytop=y.margin[i]+add.y+full.agreement[i],
         col="black")
    rm(add.x, add.y)
  } # end for loop
  rm(i)

  ####
  # 6. draw 45 degree diagonal line
  abline(a=0, b=1, col=col.line, lwd=lwd.line, lty=lty.line)

  ####
  # 7. add marginal counts to plot?
  if(add.x.margin){
    temp <- as.vector(apply(x, 2, sum))
    # axis labels should be at the center of each rectangle in step 3
    axis(side=3, at=(x.margin[1:k] + x.margin[2:(k+1)])/2, labels=temp,
         cex.axis=0.8, line=-0.8, tick=FALSE)
    rm(temp)
  } # end if

  if(add.y.margin){
    temp <- as.vector(apply(x, 1, sum))
    # axis labels should be at the center of each rectangle in step 3
    axis(side=4, at=(y.margin[1:k] + y.margin[2:(k+1)])/2, labels=temp, las=las,
         cex.axis=0.8, line=-0.8, tick=FALSE)
    rm(temp)
  } # end if

  ####
  # 8. add legend?

  if(add.legend){

    legend(location.legend,
           legend=c("full", paste("partial (+/- ",1:(k-1), ")", sep=""),
                    # PKC:                    expression(paste(45*degree, " line", sep=""))),
                    paste(45, " degree line", sep="")),
           col=c("black", gray.vector, "black", col.line),
           pch=c(rep(15, times=k-1),22, NA),
           pt.cex=2,
           lty=c(rep(NA, times=k), lty.line),
           lwd=c(rep(NA, times=k), lwd.line), bty="n")

  } # end if statement

  ####
  # 9. compute Bangdiwala agreement measure
  B.hat <- sum(full.agreement^2)/sum((as.vector(apply(x, 2, sum))*as.vector(apply(x, 1, sum))))
  # note: don't have to get rid of the first element of x.margin or y.margin because it is 0

  # add B.hat to plot?
  if(add.B.hat){
    legend(location.B.hat, legend=bquote(hat(B)==.(round(B.hat, digits=3))),
           fill=NA, bty="n", border="white")
  } # end if

  # output Bangdiwala agreement measure
  return(B.hat)

}

