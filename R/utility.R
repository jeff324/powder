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

#' Convert theta array to dataframe
#'
#' Reshapes the subject-level parameters in Powder.Individual or Powder.Hierarchical objects to a dataframe for easier analysis.
#'
#' @param pow.out Powder.Individual or Powder.Hierarchical object
#' @param pars subset of parameters to reshape. Defaults to all parameters.
#' @param subject subset of subjects to reshape. Defaults to all subjects.
#' @param with.burnin if false, burnin is discarded
#' @param ... additional arguments
#' @return A dataframe with columns corresponding to parameter, iteration, chain, value,
#' and subject if Powder.Hierarchical object.
#' @examples
#' \dontrun{
#' data('null')
#' lba = LBA$new()
#' #low number of samples for illustration purposes
#' out = powder(model=lba,data=null,burnin=10,n.samples=10,method='posterior')
#' phi_df = reshape.phi(out)
#' }
#' @export
reshape.theta = function(pow.out,pars=NULL,with.burnin=FALSE,...) UseMethod('reshape.theta')

#' @export
reshape.theta.Powder.Individual = function(pow.out,pars=NULL,with.burnin=FALSE,...){
     pow.out$phi = pow.out$theta
     theta = reshape.phi(pow.out,pars,with.burnin)
     return(theta)
}

#' @export
reshape.theta.Powder.Hierarchical = function(pow.out,pars=NULL,with.burnin=FALSE,subject=NULL, ...){

     if (is.null(pars)) {
          pars = colnames(pow.out$theta)
     }

     if (is.null(subject)) {
          subject = 1:length(pow.out$theta[1,1,,1])
     }

     burnin = pow.out$options$burnin
     n.iter = length(pow.out$theta[1,1,1,])
     n.chains = length(pow.out$theta[,1,1,1])

     if (with.burnin) {
          theta = pow.out$theta
          iter = 1:n.iter
     } else {
          theta = pow.out$theta[,,,burnin:n.iter]
          iter = burnin:n.iter
     }

     theta.list = lapply(1:length(subject),
                         function(s) do.call(rbind,lapply(1:length(pars),
                                                          function(p) cbind('iteration'=iter,'parameter'=pars[p],t(theta[,pars[p],subject[s],])) )))
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

#' Convert phi array to dataframe
#'
#'  Reshapes the group-level parameters in Powder.Individual or Powder.Hierarchical objects to a dataframe for easier analysis.
#'
#' @param pow.out Powder.Hierarchical object
#' @param pars subset of parameters to reshape. Defaults to all parameters.
#' @return A dataframe with columns corresponding to parameter, iteration, chain, and value
#' @param with.burnin if false, burnin is discarded
#' @examples
#' \dontrun{
#' data('null')
#' lba = LBA$new()
#' #low number of samples for illustration purposes
#' out = powder(model=lba,data=null,burnin=10,n.samples=10,method='posterior')
#' phi_df = reshape.phi(out)
#' }
#' @export
reshape.phi = function(pow.out,pars=NULL,with.burnin=FALSE){

     if (is.null(pars)) {
          pars = colnames(pow.out$phi)
     }

     burnin = pow.out$options$burnin
     n.iter = length(pow.out$phi[1,1,])
     n.chains = length(pow.out$phi[,1,1])

     if (with.burnin) {
          phi = pow.out$phi
          iter = 1:n.iter
     } else {
          phi = pow.out$phi[,,burnin:n.iter]
          iter = burnin:n.iter
     }

     phi.list = lapply(1:length(pars), function(x) cbind('iteration'=iter,'parameter'=pars[x],t(phi[,pars[x],])))
     phi = data.frame(do.call(rbind,phi.list),stringsAsFactors = FALSE)
     phi.chains = phi[,3:ncol(phi)]
     phi.chains[] = lapply(phi.chains,as.numeric)
     phi[,3:ncol(phi)] = phi.chains
     colnames(phi) = c('iteration','parameter',1:n.chains)
     phi$iteration = as.numeric(phi$iteration)
     phi = reshape2::melt(phi,id.vars = c('parameter','iteration'), value.name='value', variable.name = 'chain')
     return(phi)
}

