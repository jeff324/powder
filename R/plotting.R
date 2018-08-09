#' Trace plot of MCMC chains
#'
#' Displays a plot of iterations as a function of sampled values with a separate plot per variable.
#' @import ggplot2
#' @importFrom graphics plot
#' @param x a Powder.Individual or Powder.Hierarchical object
#' @param pars a character vector of parameters to plot
#' @param subject a numeric vector of subject numbers to plot. If not supplied, the group-level parameters will be plotted.
#' @param options a list that can supply any of the following components
#' * \code{with.burnin} a \code{\link{logical}} indicating whether to include burnin with plot. Defaults to false.
#' * \code{ggplot_theme} a ggplot theme
#' * \code{return.plot} a \code{\link{logical}} indicating whether to return the ggplot object. A plot will be displayed regardless. Defaults to false.
#' @md
#' @param ... additional arguments
#' @return A plot of the samples as a function of iteration, parameter, and if supplied, subject
#' @export
traceplot = function(x, pars=NULL, options=list(), ...) UseMethod('traceplot')

#' Parameter density plots
#'
#' \code{densplot}
#' @import ggplot2
#' @importFrom graphics plot
#' @param x a Powder.Individual or Powder.Hierarchical object
#' @param pars a character vector of parameters to plot
#' @param subject a numeric vector of subject numbers to plot. If not supplied, the group-level parameters will be plotted.
#' @param options a list that can supply any of the following components
#' @param ... additional arguments
#' \describe{
#' \item{\code{with.burnin}}{a \code{\link{logical}} indicating whether to include burnin with plot. Defaults to false.}
#' \item{\code{ggplot_theme}}{a ggplot theme}
#' \item{\code{return.plot}}{a \code{\link{logical}} indicating whether to return the ggplot object. A plot will be displayed regardless. Defaults to false.}
#' }
#' @return A density plot
#' @export
densplot = function(x,pars=NULL,options=list(),...) UseMethod('densplot')

#' @export
traceplot.Powder.Individual = function(x,pars=NULL,options=list(),...){

     pow.out = x

     if(is.null(options$with.burnin)) {
          options$with.burnin = FALSE
     }
     if (is.null(options$ggplot.theme)) {
          options$ggplot.theme = NULL
     }
     if(is.null(options$return.plot)) {
          options$return.plot = FALSE
     }
     if (is.null(pars)) {
          pars = colnames(pow.out$theta)
     }

     traceplot.individual(pow.out,pars,options)

}

#' @export
traceplot.Powder.Hierarchical = function(x, pars=NULL, options=list(), subject=NULL, ...){

     pow.out = x

     if(is.null(options$with.burnin)) {
          options$with.burnin = FALSE
     }
     if (is.null(options$ggplot.theme)) {
          options$ggplot.theme = NULL
     }
     if(is.null(options$return.plot)) {
          options$return.plot = FALSE
     }

     if (is.null(subject)) {
          traceplot.group(pow.out,pars,options)
     } else {
          traceplot.subject(pow.out,pars,subject,options)
     }

}

#' @export
densplot.Powder.Individual = function(x,pars=NULL,options=list(),...){

     pow.out = x

     if(is.null(options$with.burnin)) {
          options$with.burnin = FALSE
     }
     if (is.null(options$ggplot.theme)) {
          options$ggplot.theme = NULL
     }
     if(is.null(options$return.plot)) {
          options$return.plot = FALSE
     }
     if (is.null(pars)) {
          pars = colnames(pow.out$theta)
     }

     densplot.individual(pow.out,pars,options)

}

#' @export
densplot.Powder.Hierarchical = function(x, pars=NULL, options=list(), subject=NULL, ...){

     pow.out = x

     if(is.null(options$with.burnin)) {
          options$with.burnin = FALSE
     }
     if (is.null(options$ggplot.theme)) {
          options$ggplot.theme = NULL
     }
     if(is.null(options$return.plot)) {
          options$return.plot = FALSE
     }

     if (is.null(subject)) {
          densplot.group(pow.out,pars,options)
     } else {
          densplot.subject(pow.out,pars,subject,options)
     }

}

densplot.individual = function(pow.out,pars,options){

     theta = reshape.theta(pow.out,pars,options$with.burnin)

     theta = theta[,c('parameter','value')]
     p1 = ggplot2::ggplot(theta,aes(value)) +
          stat_density(geom='line') +
          facet_grid(.~parameter,scales='free_x') +
          theme(legend.position = 'none')
     plot(p1)
     if (options$return.plot) {
          return(p1)
     }

}

traceplot.individual = function(pow.out,pars,options){

     theta = reshape.theta(pow.out,pars,options$with.burnin)

     p1 = ggplot2::ggplot(theta,aes(x=iteration,y=value,group=chain,color=chain)) +
          geom_line() +
          facet_grid(parameter~.,scales='free_y') +
          theme(legend.position = 'none')
     plot(p1)


     if (options$return.plot) {
          return(p1)
     }
}


densplot.subject = function(pow.out,pars,subject,options){

     theta = reshape.theta(pow.out,pars,options$with.burnin,subject)
     theta = theta[,c('parameter','value','subject')]
     p1 = ggplot2::ggplot(theta,aes(value)) +
          stat_density(geom='line') +
          facet_grid(parameter~subject,scales='free_y') +
          theme(legend.position = 'none')
     plot(p1)
     if (options$return.plot) {
          return(p1)
     }
}

densplot.group = function(pow.out,pars,options){

     phi = reshape.phi(pow.out,pars,options$with.burnin)
     phi = phi[,c('parameter','value')]
     p1 = ggplot2::ggplot(phi,aes(value)) +
          stat_density(geom='line') +
          facet_grid(.~parameter,scales='free_x') +
          theme(legend.position = 'none')
     plot(p1)
     if (options$return.plot) {
          return(p1)
     }

}

traceplot.subject = function(pow.out,pars,subject,options){

     theta = reshape.theta(pow.out,pars,options$with.burnin,subject)

     p1 = ggplot2::ggplot(theta,aes(x=iteration,y=value,group=chain,color=chain)) +
          geom_line() +
          facet_grid(parameter~subject,scales='free_y') +
          theme(legend.position = 'none')
     plot(p1)


     if (options$return.plot) {
          return(p1)
     }
}

traceplot.group = function(pow.out,pars,options){

     phi = reshape.phi(pow.out,pars,options$with.burnin)

     p1 = ggplot2::ggplot(phi,aes(x=iteration,y=value,group=chain,color=chain)) +
          geom_line() +
          facet_grid(parameter~.,scales='free_y') +
          theme(legend.position = 'none')
     plot(p1)


     if (options$return.plot) {
          return(p1)
     }
}




