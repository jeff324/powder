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
#' * \code{randomize_phi} Assume independence between subject and group-level parameters. Useful for improved sampling of group-level parameters.
#' * \code{zLag} When \code{randomize_phi = TRUE} and \code{method = 'parallel'}, zLag is the number of iterations that can be reached into from the past for the z-Update.
#' * \code{zStart} When \code{randomize_phi = TRUE} and \code{method = 'parallel'}, zStart is the iteration to begin z-Updating.
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
#' model = LBA.Individual$new(b=T)
#' data('individual',package='powder')
#' pow.out = powder(data=individual, model=model, num.temps=30)
#' est = marginal.likelihood(pow.out)
#' print(est)
#' }
#' @export
powder <- function(model,data,num.temps=NULL,alpha=.3,high.temps.first=FALSE,n.sequences=1,current.sequence=1,
                   n.samples=1000,n.chains=NULL,burnin=500,meltin=250,de_params=list(),method='standard',
                   return.samples=TRUE,verbose=TRUE,update=100) UseMethod("powder")

#' @export
powder.Model.Individual = function(model,data,num.temps=NULL,alpha=.3,high.temps.first=FALSE,n.sequences=1,current.sequence=1,
                                   n.samples=1000,n.chains=NULL,burnin=500,meltin=250,de_params=list(),method='standard',
                                   return.samples=TRUE,verbose=TRUE,update=100){

     theta.names = model$theta.names
     n.pars = length(theta.names)

     #check and set default parameters

     if (is.null(de_params$b)) {
          de_params$b = .001
     }

     if (is.null(num.temps)) {
          num.temps = 30
     }

     if (method != 'posterior') {
          temperatures = get_temperatures(num.temps,
                                          alpha,
                                          high.temps.first,
                                          n.sequences,
                                          current.sequence)

          if (method == 'standard') {
               message = 'Sampling power posterior @ temperature'
          } else {
               message = paste('Sampling', length(temperatures),'power posteriors in parallel')
               metlin = 0
          }
     } else {
          temperatures = 1
          meltin = 0
          message = 'Sampling posterior'
     }

     if (num.temps == 1) {
          temperatures = 1 #default to posterior sampling
     }


     if (is.null(n.chains)) {

          if(method == 'parallel'){
               n.chains = length(temperatures)
          }else{
               n.pars = n.pars
               n.chains = 3*n.pars
          }
     } else {

          if(method == 'parallel' & n.chains != length(temperatures)){
               n.chains = length(temperatures)
               warning('n.chains != num.temps. Setting n.chains = num.temps.',
                       call. = FALSE, immediate. = TRUE)
          }
     }

     if( !(method == 'parallel' |  method == 'posterior' |  method == 'standard')){
          stop(paste('method', method, 'not found. Select \'parallel\', \'standard\', or \'posterior\''),
               call.=FALSE)
     }

     if (n.samples <= 0) {
          stop('n.samples must be > 0',call. = FALSE)
     }

     if (num.temps <= 0) {
          stop('num.temps must be > 0',call. = FALSE)
     }

     if (burnin < 0) {
          stop('burnin must be >= 0',call. = FALSE)
     }

     if (meltin < 0 & method == 'standard') {
          stop('meltin must be >= 0',call. = FALSE)
     }

     if (n.sequences > num.temps) {
          stop('n.sequences must be <= num.temps',call. = FALSE)
     }

     if (is.null(current.sequence)) {
          stop('please choose current.sequence',call. = FALSE)
     }

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

     if (is.null(de_params$b)) {
          de_params$b = .001
     }

     if (is.null(num.temps)) {
          num.temps = 30
     }

     if (method != 'posterior') {
          temperatures = get_temperatures(num.temps,
                                          alpha,
                                          high.temps.first,
                                          n.sequences,
                                          current.sequence)

          if (method == 'standard') {
               message = 'Sampling power posterior @ temperature'
          } else {
               message = paste('Sampling', length(temperatures),'power posteriors in parallel')
               metlin = 0
          }
     } else {
          temperatures = 1
          meltin = 0
          message = 'Sampling posterior'
     }

     if (num.temps == 1) {
          temperatures = 1 #default to posterior sampling
     }


     if (is.null(n.chains)) {

          if(method == 'parallel'){
               n.chains = length(temperatures)
          }else{
               n.pars = n.pars
               n.chains = 3*n.pars
          }
     } else {

          if(method == 'parallel' & n.chains != length(temperatures)){
               n.chains = length(temperatures)
               warning('n.chains != num.temps. Setting n.chains = num.temps.',
                       call. = FALSE, immediate. = TRUE)
          }
     }

     if( !(method == 'parallel' |  method == 'posterior' |  method == 'standard')){
          stop(paste('method', method, 'not found. Select \'parallel\', \'standard\', or \'posterior\''),
               call.=FALSE)
     }

     if (n.samples <= 0) {
          stop('n.samples must be > 0',call. = FALSE)
     }

     if (num.temps <= 0) {
          stop('num.temps must be > 0',call. = FALSE)
     }

     if (burnin < 0) {
          stop('burnin must be >= 0',call. = FALSE)
     }

     if (meltin < 0 & method == 'standard') {
          stop('meltin must be >= 0',call. = FALSE)
     }

     if (n.sequences > num.temps) {
          stop('n.sequences must be <= num.temps',call. = FALSE)
     }

     if (is.null(current.sequence)) {
          stop('please choose current.sequence',call. = FALSE)
     }

     if (is.null(de_params$randomize_phi)){
          de_params$randomize_phi = TRUE
          warning('Assuming independence between group and subject-level parameters. Set de_params$randomize_phi = FALSE to model correlation between group and subject-level parameters.',call.=FALSE,immediate.=TRUE)
     }

     if (method == 'parallel' & de_params$randomize_phi) {
          if (is.null(de_params$zLag)) {
               de_params$zLag = 200
               warning('No zLag value supplied. Setting de_params$zLag = 200',call.=FALSE,immediate.=TRUE)
          }

          if(is.null(de_params$zStart)){
               de_params$zStart = 1500
               warning('No zStart value supplied. Setting de_params$zStart = 1500.',call.=FALSE,immediate.=TRUE)
               if(de_params$zStart > burnin){
                    warning('burnin must be greater than de_params$zStart. Setting burnin to 2*zStart.',call.=FALSE,immediate.=TRUE)
                    burnin = de_params$zStart * 2
               }
          }
     }

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

