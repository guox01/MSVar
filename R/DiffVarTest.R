#' Identify differentially variable proteins
#'
#' Conduct hypothesis test for differentially variable proteins. All
#' \code{\link{proObj}} objects should be associated with the same scaling of
#' estimated \emph{v}. Thus, they must have the same scaling results of
#' \emph{v}.
#'
#' Proteins which are located in the lower part of \emph{v} are obviously
#' non-DVP, their \emph{p}-values are set as \code{NA} directly.
#'
#' @param x,y \code{x} is any R object for which a \code{diffTest} method has
#'   been defined. For the method for class \code{\link{proObj}}, \code{x} and
#'   \code{y} are two \code{proObj} objects to be compared. They must be
#'   associated with estimated \emph{v}, Fisher information and scaled by the
#'   same baseline \eqn{v_0}.
#' @param scaleV.raw.FC A specified average basic fluctuation around mean
#'   expression of a protein. The default setting of \code{scaleV.raw.FC} is
#'   10\%, which amounts to an average fluctuation of 10\% around mean expression.
#'   Accordingly, the baseline fluctuation for each protein is defined as
#'   \eqn{log(log_2(1+10\%)^2)}. If \code{scaleV.raw.FC} is 0, the differential
#'   variability test would be conducted without scaling of estimated \emph{v}.
#' @param alternative The alternative hypothesis. Must be one of \code{"upper"}
#'   (default), \code{"lower"}, and \code{"two.sided"}. Can be abbreviated.
#' @param p.adjust.method The method for adjusting \emph{p}-values for multiple
#'   testing. See \code{\link[stats]{p.adjust}} for all options.
#'
#' @return An object of the \code{\link[base]{class}} \code{c("diffVar",
#'   "data.frame")}, which records the differential variability test results for
#'   each protein by each row.
#' @export
#'
DiffVar <- function(x, y, scaleV.raw.FC = 0.1,
                    alternative = c("upper", "lower", "two.sided"),
                    p.adjust.method = "BH") {
  if (!is(x, "proObj") | !is(y, "proObj")) {
    stop("x and y must be of the class \"proObj\".")
  }

  if (scaleV.raw.FC == 0) {
    if (is.null(x$resBio$est.bio.res) | is.null(y$resBio$est.bio.res)) {
      stop(paste("The biological variance have not been estimated yet.",
                 "Try applying first the biological-variance-fitting process.",
                 sep = "\n"))
    } else {
      est.bio.res1 <- x$resBio$est.bio.res
      est.bio.res2 <- y$resBio$est.bio.res
    }
  } else {
    if (is.null(x$resBio$est.bio.scale.res[[paste0("rawFCscale", scaleV.raw.FC)]]) |
        is.null(y$resBio$est.bio.scale.res[[paste0("rawFCscale", scaleV.raw.FC)]])) {
      stop(paste("The biological variance have not been scaled yet.",
                 "Try applying first the biological-variance-scaling process.",
                 sep = "\n"))
    }
    est.bio.res1 <- x$resBio$est.bio.scale.res[[paste0("rawFCscale", scaleV.raw.FC)]]
    est.bio.res2 <- y$resBio$est.bio.scale.res[[paste0("rawFCscale", scaleV.raw.FC)]]
  }

  raw.est.bio.res <- x$resBio$est.bio.res
  group1 <- attr(raw.est.bio.res, "group1")

  pro_name <- intersect(rownames(est.bio.res1), rownames(est.bio.res2))
  est.bio.res1 <- est.bio.res1[pro_name, ]
  est.bio.res2 <- est.bio.res2[pro_name, ]
  
  highv <- x$pro.name[group1]
  lowv <- setdiff(pro_name, highv)

  # Derive p-values.
  v1 <- est.bio.res1$v
  fisher.v1 <- est.bio.res1$fisher.v
  v2 <- est.bio.res2$v
  fisher.v2 <- est.bio.res2$fisher.v
  Sta <- (v1 - v2) / sqrt(1 / fisher.v1 + 1 / fisher.v2)
  names(Sta) <- pro_name

  alternative <- match.arg(alternative)
  if (alternative == "upper") {
    pval <- pnorm(Sta, 0, 1, lower.tail = F)
  } else if (alternative == "lower") {
    pval <- pnorm(Sta, 0, 1, lower.tail = T)
  } else if (alternative == "two.sided") {
    pval <- 2 * pnorm(-abs(Sta), 0, 1, lower.tail = T)
  }
  pval[lowv] <- NA
  padj <- p.adjust(pval, method = p.adjust.method)

  # Return.
  res <- data.frame(Sta = Sta,
                    v1 = v1, fisher.v1 = fisher.v1,
                    v2 = v2, fisher.v2 = fisher.v2,
                    pval = pval, padj = padj,
                    row.names = pro_name)

  attr(res, "scaleV.raw.FC") <- scaleV.raw.FC
  attr(res, "p.adjust.method") <- p.adjust.method
  attr(res, "alternative") <- alternative

  class(res) <- c("diffVar", "data.frame")
  res
}

#' Plot a \code{diffVar} Object
#'
#' Given a \code{diffVar} object, which records the results of calling
#' HVRs and/or LVRs from a \code{\link{proObj}}, this method creates a scatter
#' plot of \code{(ref.mean, v)} pairs from all observations, marking
#' specifically the ones that are associated with statistical significance.
#' Besides, the theoretical \eqn{v_0} associated with this round of test of
#' \eqn{v} is added to the plot, serving as a baseline to which each sample
#' variance can be compared.
#'
#' @param res A \code{diffVar} object.
#' @param ref.mean A numeric vector of mean among all reference samples.
#' @param P.cutoff A scalar between 0 and 1. Proteins with \emph{p}-values below
#'   this cutoff are considered statistically significant under the null
#'   hypothesis and are correspondingly highlighted in a distinct color in the
#'   plot.
#' @param p.type A character which indicates the type of \emph{p}-value to show
#'   in the plot. Must be one of \code{"pval"} (default), and \code{"padj"}. Can
#'   be abbreviated.
#' @param Nsample The sample size in this group.
#' @param sample.type1 The sample type1 of differential variability test.
#' @param sample.type2 The sample type2 of differential variability test.
#' @param outdir If \code{PLOT} is TRUE, the out put direction.
#' @param outadd If \code{PLOT} is TRUE, the additions to the out put documents.
#' @param args.points A named list of graphical parameters for plotting the
#'   point estimates (see \code{\link[graphics]{points}} for commonly used
#'   graphical parameters).
#' @param args.abline A named list of graphical parameters for plotting the null
#'   \emph{v_0} (see \code{\link[graphics]{abline}} for commonly used graphical
#'   parameters).
#' @param args.legend A named list of graphical parameters for legend (see
#'   \code{\link[graphics]{legend}} for commonly used graphical parameters).
#' @param xlab Label for the X axis.
#' @param ylab Label for the Y axis.
#' @param xlim Range of X coordinate.
#' @param ylim1 Range of Y coordinate in subplot1.
#' @param ylim2 Range of Y coordinate in subplot2.
#' @param ... Further arguments to be passed to \code{\link[graphics]{plot}}.
#'
#' @return This function returns \code{NULL} invisibly.
#' @export
#'
#' @seealso \code{\link{proObj}} for creating a \code{proObj} object;
#'   \code{\link{estTechVar}} for fitting technical variance given a
#'   \code{proObj} objects; \code{\link{estBioVar}} for estimating biological
#'   mean \emph{m} and \emph{log}-transformed biological variance \emph{v};
#'   \code{diffVar} for a test of differential-variability.
diffVarPlot <- function(res, ref.mean,
                        P.cutoff = 0.05,
                        p.type = c("pval", "padj"),
                        Nsample, sample.type1, sample.type2,
                        args.points = list(pch = 21, cex = 1.2),
                        args.abline = list(col = "darkorange", lwd = 2.5,
                                           lty = 5),
                        args.legend = list(x = "topright"),
                        xlab = "mean", ylab = "v = log(Var)",
                        xlim = NULL, ylim1 = NULL, ylim2 = NULL,
                        outdir = "./", outadd = "", ...) {
  if (!is(res, "diffVar")) {
    stop("res must be of the class \"diffVar\".")
  }

  alternative <- attr(res, "alternative")
  scaleV.raw.FC <- attr(res, "scaleV.raw.FC")

  v1 <- res$v1
  v2 <- res$v2
  ref.mean <- ref.mean[rownames(res)]
  scale.cut <- log(log2(1+scaleV.raw.FC)^2)

  if(p.type == "pval") {
    PVal <- res$pval
    PLeg <- "P val"
  } else if (p.type == "padj") {
    PVal <- res$padj
    PLeg <- "BH P"
  }

  main1 <- substitute(expression(atop("Mean-LogBioVar," ~ sample.type1 ~ alternative,
                                      Npro ~ "proteins," ~ Nsample ~ sample.type1 ~ "samples, scale v by" ~ log ~ (log[2]^2 ~ scaleV.raw.FC))),
                      list(alternative = alternative, scaleV.raw.FC = 1+scaleV.raw.FC,
                           Nsample = Nsample,
                           sample.type1 = sample.type1, Npro = length(ref.mean)))
  main2 <- substitute(expression(atop("Mean-LogBioVar," ~ sample.type1 ~ alternative,
                                      Npro ~ "proteins," ~ Nsample ~ sample.type2 ~ "samples, scale v by" ~ log ~ (log[2]^2 ~ scaleV.raw.FC))),
                      list(alternative = alternative, scaleV.raw.FC = 1+scaleV.raw.FC,
                           Nsample = Nsample,
                           sample.type2 = sample.type2, Npro = length(ref.mean)))

  if (is.null(ylim1)) {
    ylim1 <- c(max(min(v1), -22) - 0.2, max(v1) + 0.2)
  }
  if (is.null(ylim2)) {
    ylim2 <- c(max(min(v2), -22) - 0.2, max(v2) + 0.2)
  }
  if (is.null(xlim)) {
    xlim <- c(min(ref.mean) - 0.01, max(ref.mean) + 0.01)
  }

  non_sig <- which(PVal >= P.cutoff)
  sig <- which(PVal < P.cutoff)

  args.plot <- list(xlab = xlab, ylab = ylab,
                    xlim = xlim,
                    cex.lab = 1.4, cex.axis = 1.3, cex.main = 1.4, ...)

  temp <- list(legend = c("mean v scatter",
                          "baseline v",
                          paste0("DVP, ", PLeg, " < ", P.cutoff, ", ", length(sig))))
  args.legend <- c(args.legend, list(pch = c(21, NA, 21), lty = c(NA, 5, NA),
                                     col = c(alpha("blue", 0.2), "darkorange", alpha("red", 0.2)), pt.cex = 1.2))

  par(mfrow = c(1, 2))

  # subplot1
  # Plot the observed mean-variance pairs.
  do.call(plot, c(list(ref.mean[non_sig], v1[non_sig], type = "p"),
                  col = alpha("blue", 0.2),
                  main = main1, ylim = ylim1,
                  args.points, args.plot))
  do.call(abline, c(list(h = scale.cut), args.abline))
  do.call(points, c(list(ref.mean[sig], v1[sig]),
                    col = alpha("red", 0.2),
                    args.points))

  # Add the legend.
  do.call(legend, c(temp, args.legend))

  # subplot2
  # Plot the observed mean-variance pairs.
  do.call(plot, c(list(ref.mean[non_sig], v2[non_sig], type = "p"),
                  col = alpha("blue", 0.2),
                  main = main2, ylim = ylim2,
                  args.points, args.plot))

  do.call(abline, c(list(h = scale.cut), args.abline))
  do.call(points, c(list(ref.mean[sig], v2[sig]),
                    col = alpha("red", 0.2), args.points))

  # Add the legend
  do.call(legend, c(temp, args.legend))

  invisible()
}
