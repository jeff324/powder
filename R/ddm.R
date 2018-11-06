#' Individual subject Drift Diffusion Model (DDM)
#'
#' @docType class
#' @importFrom R6 R6Class
#' @importFrom stats dnorm
#' @format \code{\link{R6Class}} object.
#' @param a threshold seperation
#' @param v drift rate
#' @param t0 non-decision time
#' @param z starting point
#' @param st0 inter-trial-variability of non-decision time
#' @param sv inter-trial-variability of drift rate
#' @param sz inter-trial-variability of starting point
#' @param prior list containing priors
#' @param conds numeric vector containing values corresponding to conditions
#' @param correct.response a character vector. Specifies the value corresponding to the correct response
#' @param error.response a character vector. Specifies the value corresponding to the error response
#' @param rt.column a character vector. Specifies the name of the column in the data for response times.
#' @param response.column a character vector. Specifies the name of the column in the data for respones.
#' @param condition.column a character vector. Specifies the name of the column in the data for the conditions.
#' @examples
#' \dontrun{
#' # DDM model that varies threshold and drift rate across 3 conditions
#' model = DDM.Individual$new(a=T,v=T,conds=1:3)
#'
#' # Note, inter-trial-variability components, st0, sv, and sz
#' # have 3 options: TRUE, FALSE, or 0. When TRUE,
#' # they will vary across conditions. When FALSE
#' # they are fit, but not allowed to vary across conditions.
#' By default, they are set to 0.
#'
#' # st0 is a parameter, other variability parameters are 0
#' model_st0 = DDM$new(st0=FALSE)
#'
#' # sv is allowed to vary over conditions
#' model_sv = DDM$new(sv=TRUE)
#'
#' # Setting priors:
#' prior = list(
#'     a = list(mu = c(mu = 1, sigma = 1),
#'              sigma = c(mu = 1, sigma = 1)),
#'     v = list(mu = c(mu = 1, sigma = 1),
#'              sigma = c(mu = 1, sigma = 1)),
#'     t0 = list(
#'          mu = c(mu = .3, sigma = .3),
#'          sigma = c(mu = .3, sigma = .3)),
#'     st0 = list(
#'          mu = c(mu = .1, sigma = .1),
#'          sigma = c(mu = .1, sigma = .1)))
#'
#'  # a and v vary across conditions
#'  # sz and sv are 0
#'  # t0 and st0 are parameters that do not vary across conditions
#'  model = DDM$new(a=TRUE,v=TRUE,prior=prior)
#'
#'  # Priors can also be changed after creating the model
#'  model$prior$a$mu = c(mu=2,sigma=2)
#'
#' }
#' @field theta.names a character vector containing the names of the subject-level parameters
#' @field theta.init a function that provides a random initial value for each subject-level parameter
#' @field theta.start.point a numeric vector containing means of start points used to initialize in theta.init
#' @field vary.parameter a logical vector containing parameters to vary
#' @field prior a list containing priors on all parameters
#' @section Methods:
#' \describe{
#' \item{\code{log.dens.prior(x,hyper)}}{likelihood of subject-level parameters given group-level parameters}
#' \item{\code{log.dens.like(x,data,par.names)}}{LBA likelihood function}}
#' @export
DDM = R6::R6Class(

     'Model.Hierarchical',

     inherit = DDM.Individual,

     public = list(
          prior = list(
               a = list(mu = c(mu = 1, sigma = 1),
                        sigma = c(mu = 1, sigma = 1)),

               v = list(mu = c(mu = 1, sigma = 1),
                        sigma = c(mu = 1, sigma = 1)),

               t0 = list(
                    mu = c(mu = .3, sigma = .3),
                    sigma = c(mu = .3, sigma = .3)
               ),

               z = list(
                    mu = c(mu = .5, sigma = .5),
                    sigma = c(mu = .5, sigma = .5)
               ),

               sz = list(
                    mu = c(mu = .2, sigma = .2),
                    sigma = c(mu = .2, sigma = .2)
               ),

               sv = list(
                    mu = c(mu = .1, sigma = .1),
                    sigma = c(mu = .1, sigma = .1)
               ),

               st0 = list(
                    mu = c(mu = .1, sigma = .1),
                    sigma = c(mu = .1, sigma = .1)
               )
          ),

          phi.names = NULL,

          phi.start.points = list(
               a = c(mu = 1, sd = 1),
               v = c(mu = 1, sd = 1),
               t0 = c(mu = .3, sd = .3),
               z = c(mu = 5, sd = .5),
               sz = c(mu = .2, sd = .2),
               sv = c(mu = .1, sd = .1),
               st0 = c(mu = .1, sd = .1)
          ),

          phi.init = function() {
               par_names = names(self$vary.parameter)
               phi = rep(NA, length(self$phi.names))
               for (i in 1:length(par_names)) {
                    tmp = grep(par_names[i], self$phi.names)
                    phi[tmp] = msm::rtnorm(length(tmp),
                                           self$phi.start.points[[par_names[i]]]['mu'],
                                           self$phi.start.points[[par_names[i]]]['sd'],0,1)
               }
               return(phi)

          },

          log.dens.prior = function(x, hyper) {
               out = 0
               v_idx = grep('^v', names(x))
               trunc_norm_pars = names(x)[-v_idx]

               for (p in trunc_norm_pars) {
                    out =  msm::dtnorm(x[p], hyper[paste(p, "mu", sep = ".")],
                                       hyper[paste(p, "sigma", sep = ".")],
                                       0, Inf, log = TRUE) + out
               }

               norm_pars = names(x)[v_idx]

               for (p in norm_pars) {
                    out =  dnorm(x[p], hyper[paste(p, "mu", sep = ".")],
                                 hyper[paste(p, "sigma", sep = ".")], log = TRUE) + out
               }

               return(out)
          },

          log.dens.hyper = function(theta, phi, prior) {
               if (length(grep('^v', names(phi))) > 0) {
                    sum((dnorm(theta, phi[1], phi[2], log = TRUE))) +
                         (dnorm(phi[1], prior$mu[1], prior$mu[2], log = TRUE)) +
                         (dnorm(phi[2], prior$mu[1], prior$sigma[2], log = TRUE))
               } else {
                    sum((msm::dtnorm(theta, phi[1], phi[2], 0, Inf, log = TRUE))) +
                         (msm::dtnorm(phi[1], prior$mu[1], prior$mu[2], 0, Inf, log = TRUE)) +
                         (msm::dtnorm(phi[2], prior$mu[1], prior$sigma[2], 0, Inf, log = TRUE))
               }
          },

          predict = function(pow.out,conds=NULL,thin=1,n=1,subjects=NULL){

               if(is.null(subjects)){
                    subjects = 1:length(pow.out$theta[1,1,,1])
               }

               if(is.null(conds)){
                    conds = self$conds
               }

               out = plyr::llply(subjects, function(x) private$predict_theta(theta=pow.out$theta[,,x,],
                                                                             conds=conds,
                                                                             thin=thin,
                                                                             n=n,
                                                                             verbose=FALSE),
                                 .progress='text'
               )
               out = Map(cbind,out,subject=subjects)
               out = do.call(rbind,out)

               return(out)

          },

          initialize = function(a = FALSE,
                                v = FALSE,
                                t0 = FALSE,
                                z = 0,
                                sz = 0,
                                sv = 0,
                                st0 = 0,
                                conds = NULL,
                                prior = NULL,
                                correct.response = NULL,
                                error.response = NULL,
                                rt.column = NULL,
                                response.column = NULL,
                                condition.column = NULL) {

               self$vary.parameter['a'] = a
               self$vary.parameter['v'] = v
               self$vary.parameter['t0'] = t0
               self$vary.parameter['z'] = z
               self$vary.parameter['sz'] = sz
               self$vary.parameter['sv'] = sv
               self$vary.parameter['st0'] = st0

               if (is.null(correct.response)) {
                    self$correct.response = 'upper'
                    warning('Value for correct response not specified, using \"upper\" as default. If this is not correct, please change in data.',call.=F)
               } else {
                    self$correct.response = correct.response
               }

               if (is.null(error.response)) {
                    self$error.response = 'lower'
                    warning('Value for error response not specified, using \"lower\" as default. If this is not correct, please change in data.', call. = F)
               } else {
                    self$error.response = error.response
               }

               if (is.null(rt.column)) {
                    self$rt.column = 'Time'
                    warning('Name of RT column not specified, using \"Time\" as default. If this is not correct, please change in data.', call. = F)
               } else {
                    self$rt.column = rt.column
               }

               if (is.null(response.column)) {
                    self$response.column = 'Correct'
                    warning('Name of response column not specified, using \"Correct\" as default. If this is not correct, please change in data.', call. = F)
               } else {
                    self$response.column = response.column
               }

               if (is.null(condition.column)) {
                    self$condition.column = 'Cond'
                    warning('Name of condition column not specified, using \"Cond\" as default. If this is not correct, please change in data.', call. = F)
               } else {
                    self$condition.column = condition.column
               }


               par_list = list('a'=a,'v'=v,'t0'=t0,'z'=z,'sz'=sz,'sv'=sv,'st0'=st0)

               for(i in 1:length(par_list)){
                    par_name = names(par_list)[i]
                    if (as.character(par_list[[i]]) == '0') {
                         if (is.null(prior[[par_name]])) {
                              self$prior[[par_name]] = 0
                              self$constant[[par_name]] = 0
                         } else{
                              if (is.numeric(prior[[par_name]])) {
                                   if(prior[[par_name]] == 0){
                                        self$prior[[par_name]] = 0
                                        self$constant[[par_name]] = 0
                                   }
                              }
                         }
                    }
               }


               if (!is.null(conds)) {
                    self$conds = conds
               }

               private$make.theta.names()
               private$make.start.points()
               private$make.phi.names()
               if (is.null(prior)) {
                    private$make.prior.default()
               } else {
                    private$make.prior(prior)
               }

          }
     ),

     private = list(
          make.phi.names = function() {
               self$phi.names = paste(rep(self$theta.names, each = 2), c("mu", "sigma"), sep = ".")
          },

          make.prior.default = function() {
               par_names = names(self$vary.parameter)
               prior = list()

               for (i in 1:length(par_names)) {
                    tmp = grep(par_names[i], self$theta.names, value = TRUE)
                    if(length(tmp) != 0){
                         for (n in 1:length(tmp)) {
                              tmp2 = tmp[n]
                              prior[[tmp2]] = list(
                                   mu = c(self$prior[[paste(par_names[i])]]$mu[['mu']],
                                          self$prior[[paste(par_names[i])]]$mu[['sigma']]),
                                   sigma = c(self$prior[[paste(par_names[i])]]$sigma[['mu']],
                                             self$prior[[paste(par_names[i])]]$sigma[['sigma']])
                              )
                         }
                    }
               }

               self$prior = prior
          },

          make.prior = function(prior) {
               par_names = names(self$vary.parameter)
               prior_tmp = list()

               for (i in 1:length(par_names)) {
                    tmp = grep(par_names[i], self$theta.names, value = TRUE)
                    if(length(tmp) != 0){
                         for (n in 1:length(tmp)) {
                              tmp2 = tmp[n]
                              prior_tmp[[tmp2]] = list(
                                   mu = c(prior[[paste(par_names[i])]]$mu[['mu']],
                                          prior[[paste(par_names[i])]]$mu[['sigma']]),
                                   sigma = c(prior[[paste(par_names[i])]]$sigma[['mu']],
                                             prior[[paste(par_names[i])]]$sigma[['sigma']])
                              )
                         }
                    }
               }

               self$prior = prior_tmp
          }


     )
)
