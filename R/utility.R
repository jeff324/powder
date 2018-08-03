#' Load an RData file
#'
#' \code{load_rdata} loads an RData file into a specified variable
#'
#' @param  file_name a character string giving the name of the file to load
#' @return The objects in the RData file. If there is more than one object, a list will be returned.
#'
#' @examples
#' \dontrun{
#' x = pi # to ensure there is some data
#' x_list = list(pi = pi, pi_sq = pi^2)
#' save(x, file= "pi.RData")
#' save(x_list, file= "pi_list.RData")
#' rm(x)
#' rm(x_list)
#' x = load_rdata('pi.RData')
#' x_list = load_rdata('pi_list.RData')
#' }
#' @keywords internal
load_rdata = function(file_name){
     load(file_name)
     obj_names = ls()
     obj_names = obj_names[obj_names != 'file_name']
     if(length(obj_names) > 1){
          #return a list of objects
          dat = sapply(obj_names, function(x)get(x), simplify=FALSE, USE.NAMES=TRUE)
     }else{
          #return a single object
          dat = get(obj_names)
     }
     return(dat)
}

#' Compute Bayes factor from powder objects
#'
#' \code{bayes.factor} computes the Bayes factor from the powder objects
#'
#' @param m1 a Powder.Hierarchical or Powder.Individual object
#' @param m2 a Powder.Hierarchical or Powder.Individual object
#' @param log When TRUE, the natural log of the Bayes factor will be returned.
#' @return The natural log of the Bayes factor in terms of m1 / m2.
#' @export
bayes.factor = function(m1,m2,log=TRUE){
     bf = data.frame(Method = c('TI','Log Steppingstone','Harmonic Mean'),
                     Log.Bayes.factor= c(m1$Value[m1$Method == 'TI'] - m2$Value[m1$Method == 'TI'],
                                 m1$Value[m1$Method == 'Log Steppingstone'] - m2$Value[m1$Method == 'Log Steppingstone'],
                                 m1$Value[m1$Method == 'Harmonic Mean'] - m2$Value[m1$Method == 'Harmonic Mean']))
     if (log==FALSE) {
          bf$Log.Bayes.factor = exp(bf$Log.Bayes.factor)
          colnames(bf) = c('Method','Bayes factor')
     }
     return(bf)
}


#' Plots the samples by iteration of each chain
#'
#' \code{plot.Powder.Hierarchical}
#' @import ggplot2
#' @importFrom graphics plot
#' @param x a Powder.Hierarchical object
#' @param pars a character vector of parameters to plot
#' @param subject a numeric vector of subject numbers to plot. If not supplied, the group-level parameters will be plotted.
#' @param type a character vector indicating the type of plot. Two plots are supported: `density` and `iteration`.
#' @param options a list that can supply any of the following components
#' \describe{
#' \item{\code{with.burnin}}{a \code{\link{logical}} indicating whether to include burnin with plot. Defaults to false.}
#' \item{\code{ggplot_theme}}{a ggplot theme}
#' \item{\code{return.plot}}{a \code{\link{logical}} indicating whether to return the ggplot object. A plot will be displayed regardless. Defaults to false.}
#' }
#' @return A plot of the samples as a function of iteration, parameter, and if supplied, subject
#' @export
plot.Powder.Hierarchical = function(x, pars=NULL,subject=NULL,type='iteration',
                                    options=list(),...){

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
          plot.group(pow.out,pars,type,options)
     } else {
          plot.subject(pow.out,pars,subject,type,options)
     }

}

#' Plots the samples by iteration of each chain
#'
#' \code{plot.Powder.Individual}
#' @import ggplot2
#' @importFrom graphics plot
#' @param x a Powder.Individual object
#' @param pars a character vector of parameters to plot
#' @param type a character vector indicating the type of plot. Two plots are supported: `density` and `iteration`.
#' @param options a list that can supply any of the following components
#' \describe{
#' \item{\code{with.burnin}}{a \code{\link{logical}} indicating whether to include burnin with plot. Defaults to false.}
#' \item{\code{ggplot_theme}}{a ggplot theme}
#' \item{\code{return.plot}}{a \code{\link{logical}} indicating whether to return the ggplot object. A plot will be displayed regardless. Defaults to false.}
#' }
#' @return A plot of the samples as a function of iteration and parameter
#' @export
plot.Powder.Individual = function(x,pars=NULL,type='iteration',options=list(),...){

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
     plot.group(pow.out,pars,type,options)

}

plot.group = function(pow.out,pars,type,options){

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

     if (type=='iteration') {
          p1 = ggplot2::ggplot(phi,aes(x=iteration,y=value,group=chain,color=chain)) +
               geom_line() +
               facet_grid(parameter~.,scales='free_y') +
               theme(legend.position = 'none')
           plot(p1)
     }

     if (type=='density') {
          phi = phi[,c('parameter','value')]
          p1 = ggplot2::ggplot(phi,aes(value)) +
               stat_density(geom='line') +
               facet_grid(.~parameter,scales='free_x') +
               theme(legend.position = 'none')
          plot(p1)
     }

     if (options$return.plot) {
          return(p1)
     }
}

plot.subject = function(pow.out,pars,subject,type,options){
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

     if(type=='iteration'){
          p1 = ggplot2::ggplot(theta,aes(x=iteration,y=value,group=chain,color=chain)) +
               geom_line() +
               facet_grid(parameter~subject,scales='free_y') +
               theme(legend.position = 'none')
          plot(p1)
     }

     if (type=='density') {
          theta = theta[,c('parameter','value','subject')]
          p1 = ggplot2::ggplot(theta,aes(value)) +
               stat_density(geom='line') +
               facet_grid(parameter~subject,scales='free_y') +
               theme(legend.position = 'none')
          plot(p1)
     }

     if (options$return.plot) {
          return(p1)
     }
}
