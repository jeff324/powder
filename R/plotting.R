#' Plots the samples by iteration of each chain
#'
#' \code{traceplot}
#' @import ggplot2
#' @importFrom graphics plot
#' @param x a Powder.Individual or Powder.Hierarchical object
#' @param pars a character vector of parameters to plot
#' @param subject a numeric vector of subject numbers to plot. If not supplied, the group-level parameters will be plotted.
#' @param options a list that can supply any of the following components
#' \describe{
#' \item{\code{with.burnin}}{a \code{\link{logical}} indicating whether to include burnin with plot. Defaults to false.}
#' \item{\code{ggplot_theme}}{a ggplot theme}
#' \item{\code{return.plot}}{a \code{\link{logical}} indicating whether to return the ggplot object. A plot will be displayed regardless. Defaults to false.}
#' }
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

     #use the group plotting function since theta has
     #same structure as phi
     pow.out$phi = pow.out$theta
     traceplot.group(pow.out,pars,options)

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

     #use the group plotting function since theta has
     #same structure as phi
     pow.out$phi = pow.out$theta
     densplot.group(pow.out,pars,options)

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


densplot.subject = function(pow.out,pars,subject,options){

     theta = reshape.theta(pow.out,pars,subject,options)
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

     phi = reshape.phi(pow.out,pars,options)
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

     theta = reshape.theta(pow.out,pars,subject,options)

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

     phi = reshape.phi(pow.out,pars,options)

     p1 = ggplot2::ggplot(phi,aes(x=iteration,y=value,group=chain,color=chain)) +
          geom_line() +
          facet_grid(parameter~.,scales='free_y') +
          theme(legend.position = 'none')
     plot(p1)


     if (options$return.plot) {
          return(p1)
     }
}


reshape.phi = function(pow.out,pars,options){
     if (is.null(pars)) {
          pars = colnames(pow.out$phi)
     }

     burnin = pow.out$options$burnin

     if (options$with.burnin) {
          phi = pow.out$phi
     } else {
          phi = pow.out$phi[,,burnin:length(pow.out$phi[1,1,])]
     }

     n.iter = length(phi[1,1,])
     n.chains = length(phi[,1,1])

     #reshape for ggplot
     phi.list = lapply(1:length(pars), function(x) cbind('iteration'=1:n.iter,'parameter'=pars[x],t(phi[,pars[x],])))
     phi = data.frame(do.call(rbind,phi.list),stringsAsFactors = FALSE)
     phi.chains = phi[,3:ncol(phi)]
     phi.chains[] = lapply(phi.chains,as.numeric)
     phi[,3:ncol(phi)] = phi.chains
     colnames(phi) = c('iteration','parameter',1:n.chains)
     phi$iteration = as.numeric(phi$iteration)
     phi = reshape2::melt(phi,id.vars = c('parameter','iteration'), value.name='value', variable.name = 'chain')
     return(phi)
}



reshape.theta = function(pow.out,pars,subject,options){
     if (is.null(pars)) {
          pars = colnames(pow.out$theta)
     }

     burnin = pow.out$options$burnin

     if (options$with.burnin) {
          theta = pow.out$theta
     } else {
          theta = pow.out$theta[,,,burnin:length(pow.out$theta[1,1,1,])]
     }

     n.iter = length(theta[1,1,1,])
     n.chains = length(theta[,1,1,1])

     theta.list = lapply(1:length(subject),
                         function(s) do.call(rbind,lapply(1:length(pars),
                                                          function(p) cbind('iteration'=1:n.iter,'parameter'=pars[p],t(theta[,pars[p],subject[s],])) )))
     theta.list = Map(cbind,subject=subject,theta.list)
     theta = data.frame(do.call(rbind,theta.list),stringsAsFactors = FALSE)
     theta.chains = theta[,4:ncol(theta)]
     theta.chains[] = lapply(theta.chains,as.numeric)
     theta[,4:ncol(theta)] = theta.chains
     colnames(theta) = c('subject','iteration','parameter',1:n.chains)
     theta$iteration = as.numeric(theta$iteration)
     theta = reshape2::melt(theta,id.vars = c('parameter','subject','iteration'),value.name='value',variable.name = 'chain')
     return(theta)
}


