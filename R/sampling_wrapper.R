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
#' @param method A character vector that specifies the type of sampling to be performed and accepts one of the following:
#' * \code{standard} This options samples from each power posterior along the temperature schedule in sequence
#' * \code{parallel} This options samples from each power posterior in parallel, where the target density of each chain is
#' a power posterior at a given temperature. Using this method will cause n.chains to default to num.temps.
#' * \code{sample.posterior} This options samples from the posterior (i.e. the power posterior where temperature = 1).
#' Although this option is useful for parameter estimation, it is not possible to obtain marginal likelihoods via this option.
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
                   n.samples,n.chains,burnin,meltin,de_params=list(),method,
                   return.samples,verbose,update) UseMethod("powder")

#' @export
powder.Model.Individual = function(model,data,num.temps=NULL,alpha=.3,high.temps.first=FALSE,n.sequences=1,current.sequence=1,
                                   n.samples=1000,n.chains=NULL,burnin=500,meltin=250,de_params=list(),method='standard',
                                   return.samples=TRUE,verbose=TRUE,update=100){

     theta.names = model$theta.names
     n.pars = length(theta.names)

     #check and set default parameters
     check_pars()

     if (is.list(data[[1]])) {
          stop('Expected data[[1]] to be a vector not a list',call. = FALSE)
     }

     if (method == 'standard' | method == 'posterior') {
          samples = standard.sampling.individual(model,data,theta.names,n.pars,temperatures,de_params,
                                                 burnin,meltin,n.samples,n.chains,method,message,
                                                 return.samples,verbose,update)
     } else if (method == 'parallel') {
          samples = parallel.sampling.individual(model,data,theta.names,n.pars,temperatures,de_params,
                                                 burnin,n.samples,n.chains,method,message,
                                                 return.samples,verbose,update)
     }

     #options list for output
     opt = list(num.temps=num.temps,alpha=alpha,high.temps.first=high.temps.first,
                n.sequences=n.sequences,current.sequence=current.sequence,n.pars=n.pars,
                n.samples=n.samples,n.chains=n.chains,burnin=burnin,meltin=meltin,
                de_params=de_params,method=method,return.samples=return.samples,
                temperatures=temperatures)

     if(verbose){
          cat('\n')
          cat('Sampling Completed')
          cat('\n')
     }

     if (return.samples) {
          out = list(log.like.list=samples$log.like.list, theta=samples$theta, options=opt, model=model)
          class(out) = 'Powder.Individual'
          return(out)
     } else {
          out = list(log.like.list=samples$log.like.list, options=opt, model=model)
          class(out) = 'Powder.Individual'
          return()
     }


}

#' @export
powder.Model.Hierarchical = function(model,data,num.temps=NULL,alpha=.3,high.temps.first=FALSE,n.sequences=1,current.sequence=1,
                  n.samples=1000,n.chains=NULL,burnin=500,meltin=250,
                  de_params=list(),method='standard',return.samples=TRUE,verbose=TRUE,update=10){

     theta.names = model$theta.names
     n.pars = length(theta.names)

     check_pars()

     if (is.null(de_params$migration)) {
          de_params$migration = FALSE
          de_params$migration.freq = 0
          de_params$migration.start = -1
          de_params$migration.end = -1
     } else {
          if (is.null(de_params$migration.freq)) {
               warning('migration.freq not specified. Defaulting to 20',call. = FALSE,immediate. = TRUE)
               de_params$migration.freq = 20
          }

          if (is.null(de_params$migration.start)) {
               warning('migration.start not specified. Defaulting to burnin',call. = FALSE,immediate. = TRUE)
               de_params$migration.start = burnin
          }

          if (is.null(de_params$migration.end)) {
               warning('migration.end not specified. Defaulting to migration.start + round(n.samples/5)',
                       call. = FALSE,immediate. = TRUE)
               de_params$migration.end = de_params$migration.start + round(n.samples/5)
          }
     }

     if (method == 'standard' | method == 'posterior') {
          samples = standard.sampling.hierarchical(model,data,theta.names,n.pars,temperatures,de_params,
                                                   burnin,meltin,n.samples,n.chains,method,message,
                                                   return.samples,verbose,update)

     } else if (method == 'parallel') {
          samples = parallel.sampling.hierarchical(model,data,theta.names,n.pars,temperatures,de_params,
                                                    burnin,n.samples,n.chains,method,message,
                                                    return.samples,verbose,update)
     }


     opt = list(num.temps=num.temps, alpha=alpha, high.temps.first=high.temps.first,
                n.sequences=n.sequences, current.sequence=current.sequence,
                n.subj=length(data),n.pars=n.pars,n.hpars=length(model$phi.names),
                n.samples=n.samples, n.chains=n.chains,burnin=burnin, meltin=meltin,
                de_params=de_params, method=method, return.samples=return.samples,
                temperatures=temperatures)


     if(verbose){
         cat('\n')
          cat('Sampling Completed')
          cat('\n')
     }
     if (return.samples) {
          out = list(log.like.list=samples$log.like.list, theta=samples$theta, phi=samples$phi, options=opt, model=model)
          class(out) = 'Powder.Hierarchical'
          return(out)
     } else {
          out = list(log.like.list=samples$log.like.list, options=opt, model=model)
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
          assign('de_params', value = de_params, envir = parent.frame())
     }

     if (is.null(get('num.temps',parent.frame()))) {
          assign('num.temps', 30, envir = parent.frame())
     }

     if (get('method',parent.frame()) != 'posterior') {
          temperatures = get_temperatures(get('num.temps',parent.frame()),
                                          get('alpha',parent.frame()),
                                          get('high.temps.first',parent.frame()),
                                          get('n.sequences',parent.frame()),
                                          get('current.sequence',parent.frame()))

          assign('temperatures',temperatures,envir = parent.frame())
          if (get('method',parent.frame()) == 'standard') {
               assign('message', 'Sampling power posterior @ temperature', envir = parent.frame())
          } else {
               assign('message', paste('Sampling', length(temperatures),'power posteriors in parallel'), envir = parent.frame())
               assign('meltin', 0, envir = parent.frame())
          }
     } else {
          assign('temperatures', 1, envir = parent.frame())
          assign('meltin', 0, envir = parent.frame())
          assign('message', 'Sampling posterior', envir = parent.frame())
     }

     if (get('num.temps',parent.frame()) == 1) {
          assign('temperatures', 1, envir = parent.frame()) #default to posterior sampling
     }


     if (is.null(get('n.chains',parent.frame()))) {

          if(get('method',parent.frame()) == 'parallel'){
               assign('n.chains',length(temperatures), envir = parent.frame())
          }else{
               n.pars = get('n.pars',parent.frame())
               assign('n.chains', 3*n.pars,  envir =parent.frame())
          }

     } else {

          if(get('method',parent.frame()) == 'parallel' & get('n.chains',parent.frame()) != length(temperatures)){
               assign('n.chains',length(temperatures), envir = parent.frame())
               warning('n.chains != num.temps. Setting n.chains = num.temps.',
                       call. = FALSE, immediate. = TRUE)
          }

     }

     method = get('method',parent.frame())
     if( !(method == 'parallel' |  method == 'posterior' |  method == 'standard')){
          stop(paste('method', method, 'not found. Select \'parallel\', \'standard\', or \'posterior\''),
               call.=FALSE)
     }

     if (get('n.samples',parent.frame()) <= 0) {
          stop('n.samples must be > 0',call. = FALSE)
     }

     if (get('num.temps',parent.frame()) <= 0) {
          stop('num.temps must be > 0',call. = FALSE)
     }

     if (get('burnin',parent.frame()) < 0) {
          stop('burnin must be >= 0',call. = FALSE)
     }

     if (get('meltin',parent.frame()) < 0 & get('method',parent.frame()) == 'standard') {
          stop('meltin must be >= 0',call. = FALSE)
     }

     if (get('n.sequences',parent.frame()) > get('num.temps',parent.frame())) {
          stop('n.sequences must be <= num.temps',call. = FALSE)
     }

     if (is.null(get("current.sequence",parent.frame()))) {
          stop('please choose current.sequence',call. = FALSE)
     }


}
