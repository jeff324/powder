#' Power posterior sampling using differential evolution
#'
#' \code{powder} Runs power posterior sampling using differential evolution markov chain monte carlo
#'
#' @param data A list of data where each element is the data for a specific subject
#' @param model See \code{\link{LBA}} for more information and an example.
#' @param num.temps Number of temperatures (i.e. number of power posteriors to sample from)
#' @param alpha controls the temperature schedule. 0.3 is recommended.
#' @param high.temps.first If true, then the posterior will be sampled from first, followed by lower temperature posteriors.
#' If false, then the prior will be sampled from first, followed by higher temperature posteriors.
#' @param n.sequences This is useful for running parallel power posteriors. If \code{n.sequences = 1}, then all power posteriors
#' will be run sequentially. If \code{n.sequences > 1} then \code{num.temps / n.sequences}
#' power posteriors will be run sequentially.
#' @param current.sequence If n.sequences > 1, then \code{num.temps / n.sequences} power posteriors will be run sequentially.
#' \code{current.sequence} is the index of the current sequence. For example, if \code{num.temps = 10}, \code{n.sequences = 2}, and
#' \code{current.sequence = 1} then temperatures 1 through 5 will be run.
#' If \code{current.sequence = 2} then temperatures 6 through 10 will be run.
#' @param n.samples The number of samples to draw from each power posterior
#' @param n.chains The number of chains to use. By default, \code{n.chains} is the number of subject-level parameters times 3.
#' @param burnin The number of samples to discard when computing the marginal likelihood.
#' These samples are included in the raw output.
#' @param meltin The power posteriors are sampled sequentially. When moving to the next power posterior, the sampling process
#' takes some time to adapt to the new power posterior. \code{meltin} is the number of samples to discard after moving to the next
#' power posterior. While \code{burnin} discards samples from the first power posterior, \code{meltin} discard samples for all
#' other power posteriors.
#' @param de_params A list containing the following options for DE-MCMC. See Turner et al. (2013) for details.
#' * \code{b} The parameter for the uniformly distributed noise term for the DE proposal. Default is .001.
#' * \code{migration} Circulates the states of the chains to remedy problem of outlier chains.
#' * \code{migration.freq} Number of iterations to wait between each migration step
#' * \code{migration.start} When to start migrating. This should be after chains are burned in.
#' * \code{migration.end} When to stop migrating. Migration should stop well before sampling is finished.
#' @md
#' @param sample.posterior If true, only the posterior will be sampled. This causes other arguments, like num.temps to be ignored.
#' @param return.samples If true, return subject and group-level samples.
#' Otherwise, return only the log likelihoods from each power posterior.
#' @param verbose Display progress
#' @param update Number of iterations between progress display updates
#' @return A list with the following elements.
#' @return \code{log.like.list} A list containing the log likelihoods for each temperature
#' @return \code{theta} An array containing the subject-level samples
#' @return \code{phi} An array containing the group-level samples
#' @return Note, if \code{return.samples = FALSE} then only \code{log.like.list} will be returned.
#' @examples
#' \dontrun{
#' model = LBA$new(b=T)
#' data('null',package='powder')
#' num.temps = 30
#' pow.out = powder(data=null, model=model, num.temps=num.temps)
#' est = marginal_likelihood(pow.out)
#' print(est)
#' }
#' @export
powder <- function(model,data,num.temps,alpha,high.temps.first,n.sequences,current.sequence,
                   n.samples,n.chains,burnin,meltin,de_params=list(),sample.posterior,
                   return.samples,verbose,update) UseMethod("powder")

#' @export
powder.Model.Individual = function(model,data,num.temps=NULL,alpha=.3,high.temps.first=FALSE,n.sequences=1,current.sequence=1,
                                     n.samples=1000,n.chains=NULL,burnin=500,meltin=250,de_params=list(),
                                   sample.posterior=FALSE,return.samples=TRUE,verbose=TRUE,update=100){

     theta.names = model$theta.names
     n.pars = length(theta.names)

     #check and set default parameters
     check_pars()

     #build the model structure
     prior = model$prior
     theta = array(NA,c(n.chains, n.pars, meltin*length(temperatures) + burnin + n.samples*length(temperatures)))
     weight = array(-Inf,c(meltin*length(temperatures) + burnin + n.samples*length(temperatures),n.chains))

     log.like.list = list()

     colnames(theta) = theta.names

     #options list for output
     opt = list(num.temps=num.temps,alpha=alpha,high.temps.first=high.temps.first,
                n.sequences=n.sequences,current.sequence=current.sequence,n.pars=n.pars,
                n.samples=n.samples,n.chains=n.chains,burnin=burnin,meltin=meltin,
                de_params=de_params,sample.posterior=sample.posterior,return.samples=return.samples,
                temperatures=temperatures)

     #run sampling
     idx = 2
     for (t in 1:length(temperatures)) {
          if(verbose){
               cat('\n')
               if(sample.posterior){
                    cat(paste(message))
               }else{
                    cat(paste(message, round(temperatures[t]),':',t,'/', length(temperatures)))
               }
          }
          if (t==1) {
               n.iter = burnin + n.samples
               for (i in 1:n.chains) {
                    while (weight[1, i]==-Inf) {
                         theta[i, ,1] = model$theta.init()
                         weight[1, i] = model$log.dens.like(theta[i, ,1],data=data,par.names=theta.names)
                    }
               }
          }else{
               n.iter = meltin + n.samples
               log.like[1, ] = sapply(1:n.chains, function(x) model$log.dens.like(theta[x, ,idx-1], data=data, par.names=theta.names))
               weight[idx-1, ] = log.like[1, ]
          }

          log.like = array(NA,c(n.iter, n.chains))
          for (i in 2:n.iter) {
               if(verbose){
                    if (i%%update==0)cat("\n ",i,'/',n.iter)
               }
               temp = t(sapply(1:n.chains, crossover_ind, n.chains=n.chains, b=b, pars=1:n.pars, use.theta=theta[ , ,idx-1], use.like=weight[idx-1, ],
                             data=data, hyper=prior, par.names=theta.names, temperature=temperatures[t], model=model))
               weight[idx, ] = temp[ ,1]
               theta[ , ,idx] = temp[ ,2:(n.pars+1)]
               log.like[i, ] = temp[ ,1]
               idx = idx + 1
          }
          log.like.list[[t]] = log.like
     }
     if(verbose){
          cat('\n')
          cat('Sampling Completed')
          cat('\n')
     }
     if (return.samples) {
          out = list(log.like.list=log.like.list, theta=theta, options=opt)
          class(out) = 'Powder.Individual'
          return(out)
     } else {
          out = list(log.like.list=log.like.list, options=opt)
          class(out) = 'Powder.Individual'
          return()
     }
}

#' @export
powder.Model.Hierarchical = function(model,data,num.temps=NULL,alpha=.3,high.temps.first=FALSE,n.sequences=1,current.sequence=1,
                  n.samples=1000,n.chains=NULL,burnin=500,meltin=250,
                  de_params=list(),sample.posterior=FALSE,return.samples=TRUE,verbose=TRUE,update=10){


     theta.names = model$theta.names
     phi.names = model$phi.names
     prior = model$prior

     n.subj = length(data)
     n.pars = length(theta.names)
     n.hpars = length(phi.names)

     check_pars()

     if (is.null(de_params$migration)) {
          de_params$migration = FALSE
          migration.freq = 0
          migration.start = -1
          migration.end = -1
     } else {
          if (is.null(de_params$migration.freq)) {
               warning('migration.freq not specified. Defaulting to 20',call. = FALSE,immediate. = TRUE)
               de_params$migration.freq = 20
               migration.freq = de_params$migration.freq
          }

          if (is.null(de_params$migration.start)) {
               warning('migration.start not specified. Defaulting to burnin',call. = FALSE,immediate. = TRUE)
               de_params$migration.start = burnin
               migration.start = burnin
          }

          if (is.null(de_params$migration.end)) {
               warning('migration.end not specified. Defaulting to migration.start + round(n.samples/5)',
                       call. = FALSE,immediate. = TRUE)
               de_params$migration.end = de_params$migration.start + round(n.samples/5)
               migration.end = de_params$migration.start + round(n.samples/5)
          }
     }


     theta=array(NA, c(n.chains, n.pars, n.subj, meltin*length(temperatures) + burnin + n.samples*length(temperatures)))
     phi=array(NA, c(n.chains, n.hpars, meltin*length(temperatures) + burnin + n.samples*length(temperatures)))
     weight=array(-Inf, c(meltin*length(temperatures) + burnin + n.samples*length(temperatures),n.chains,n.subj))
     log.like.list = list()

     colnames(theta) = theta.names
     colnames(phi) = phi.names

     for (i in 1:n.chains) {
          phi[i, ,1] = model$phi.init()
     }

     opt = list(num.temps=num.temps, alpha=alpha, high.temps.first=high.temps.first,
                n.sequences=n.sequences, current.sequence=current.sequence,
                n.subj=n.subj,n.pars=n.pars,n.hpars=n.hpars,
                n.samples=n.samples, n.chains=n.chains,burnin=burnin, meltin=meltin,
                de_params=de_params, sample.posterior=sample.posterior, return.samples=return.samples,
                temperatures=temperatures)

     idx = 2
     for (t in 1:length(temperatures)) {
          if(verbose){
               cat('\n')
               if(sample.posterior){
                    cat(paste(message))
               }else{
                    cat(paste(message, round(temperatures[t]),':',t,'/', length(temperatures)))
               }
          }
          if (t==1) {
               n.iter = burnin + n.samples

          } else {
               n.iter = meltin + n.samples
          }
          log_like = array(NA,dim=c(n.iter, n.chains, n.subj))
          if (t==1) {
               for (i in 1:n.chains) {
                    for (j in 1:n.subj) {
                         while (weight[1,i,j]==-Inf) {
                              theta[i, ,j,1] = model$theta.init()
                              weight[1,i,j] = model$log.dens.like(theta[i, ,j,1],data=data[[j]],par.names=theta.names)
                         }
                    }
               }
          } else {
               for (j in 1:n.subj) {
                    log_like[1, ,j] = sapply(1:n.chains,function(x)model$log.dens.like(theta[x, ,j,idx-1],data=data[[j]],par.names=theta.names))
                    weight[idx-1, ,j] = log_like[1, ,j]
               }
          }

          for (i in 2:n.iter) {
               if(verbose){
                    if (i%%update == 0)cat("\n ", i, '/', n.iter, ' ')
               }
               phi[,,idx] = phi[,,idx-1]
               rand.samp = sample(1:n.chains,n.chains)
               for (p in theta.names) {
                    which.theta = match(x=p, table=theta.names)
                    which.phi = grep(paste0('^',p),phi.names)
                    if (idx %% migration.freq == 0 & idx > migration.start & idx < migration.end) {
                         phi[ , ,idx] = migration.crossover_hyper(pars=which.phi, n.chains=n.chains, use.theta=theta[rand.samp,which.theta,,idx-1],
                                                              use.phi=phi[ , ,idx], prior=prior[[p]], model)
                    } else {
                         phi[ , ,idx] = t(sapply(1:n.chains, crossover_hyper, pars=which.phi, n.chains=n.chains, b=b,
                                             use.theta=theta[rand.samp, which.theta, ,idx-1], use.phi=phi[ , ,idx], prior=prior[[p]], model))
                    }
               }
               rand.samp = sample(1:n.chains, n.chains)
               hyper = phi[rand.samp, , idx]

               for (j in 1:n.subj) {
                    #cat(j)
                    if (idx %% migration.freq == 0 & idx > migration.start & idx < migration.end) {
                         temp = migration.crossover(pars=1:n.pars, n.chains=n.chains, use.theta=theta[,,j,idx-1], use.like=weight[idx-1,,j],
                                                  data=data[[j]],hyper=hyper, par.names=theta.names, model)
                    } else {
                         temp = t(sapply(1:n.chains, crossover, pars=1:n.pars, n.chains=n.chains, b=b, use.theta=theta[ , ,j,idx-1],
                                       use.like=weight[idx-1, ,j], data=data[[j]], hyper=hyper,par.names=theta.names,temperatures[t],model))
                    }
                    weight[idx, ,j] = temp[ ,1]
                    theta[ , ,j,idx] = temp[ ,2:(n.pars+1)]
                    log_like[i, ,j] = temp[ ,1]
               }
               idx = idx + 1
          }
          log.like.list[[t]] = log_like

     }
     if(verbose){
         cat('\n')
          cat('Sampling Completed')
          cat('\n')
     }
     if (return.samples) {
          out = list(log.like.list=log.like.list, theta=theta, phi=phi, options=opt)
          class(out) = 'Powder.Hierarchical'
          return(out)
     } else {
          out = list(log.like.list=log.like.list, options=opt)
          class(out) = 'Powder.Hierarchical'
          return(out)
     }
}


get_temperatures = function(num.temps,alpha,high.temps.first,n.sequences,current.sequence){

     temperatures = (0:(num.temps-1)/(num.temps-1))^(1/alpha)

     if(high.temps.first){
          temperatures = rev(temperatures)
     }

     temperatures = split(temperatures, ceiling(seq_along(temperatures)/(num.temps/n.sequences)))[[current.sequence]]

     return(temperatures)
}

check_pars = function(){

     de_params = get('de_params',parent.frame())

     if (is.null(de_params$b)) {
          de_params$b = .001
          b = de_params$b
          assign('de_params', value = de_params, envir = parent.frame())
          assign('b', value = b, envir = parent.frame())
     }

     if(is.null(get('num.temps',parent.frame()))){
          assign('num.temps', 30, envir = parent.frame())
     }

     if(!get('sample.posterior',parent.frame())){
          temperatures = get_temperatures(get('num.temps',parent.frame()),
                                          get('alpha',parent.frame()),
                                          get('high.temps.first',parent.frame()),
                                          get('n.sequences',parent.frame()),
                                          get('current.sequence',parent.frame()))

          assign('temperatures',temperatures,envir = parent.frame())
          assign(message, 'Sampling power posterior @ temperature', envir = parent.frame())
     }else{
          assign('temperatures', 1, envir = parent.frame())
          assign('meltin', 0, envir = parent.frame())
          assign('message', 'Sampling posterior', envir = parent.frame())
     }

     if(get('num.temps',parent.frame()) == 1){
          assign('temperatures', 1, envir = parent.frame()) #default to posterior sampling
     }

     if(is.list(get('data',parent.frame())[[1]])){
          stop('Expected data[[1]] to be a vector not a list',call. = FALSE)
     }

     if(is.null(get('n.chains',parent.frame()))){
          n.pars = get('n.pars',parent.frame())
          assign('n.chains', 3*n.pars,  envir =parent.frame())
     }

     if(get('n.samples',parent.frame()) <= 0){
          stop('n.samples must be > 0',call. = FALSE)
     }

     if(get('num.temps',parent.frame()) <= 0){
          stop('num.temps must be > 0',call. = FALSE)
     }

     if(get('burnin',parent.frame()) <= 0){
          stop('burnin must be > 0',call. = FALSE)
     }

     if(get('meltin',parent.frame()) <= 0 & !get('sample.posterior',parent.frame())){
          stop('meltin must be > 0',call. = FALSE)
     }

     if(get('n.sequences',parent.frame()) > get('num.temps',parent.frame())){
          stop('n.sequences must be <= num.temps',call. = FALSE)
     }

     if(is.null(get("current.sequence",parent.frame()))){
          stop('please choose current.sequence',call. = FALSE)
     }


}
