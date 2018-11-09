#' Hierarchical Linear Ballistic Accumulator model (Brown & Heathcote, 2008)
#'
#' @docType class
#' @importFrom R6 R6Class
#' @return Object of \code{\link{R6Class}} with methods for doing sampling with \code{\link{powder}}
#' @format \code{\link{R6Class}} object.
#' @examples
#' \dontrun{
#' #LBA model that varies threshold across 3 conditions
#' model_threshold = LBA$new(b=T,conds=1:3)
#' #LBA model that varies correct drift rate over 2 conditions
#' model_drift = LBA$new(vc=T,conds=1:2)
#' #LBA model that varies start point variability and non-decision over 2 conditions
#' model_sp_t0 = LBA$new(t0=T,A=T,conds=1:2)
#' }
#' @field theta.names a character vector containing the names of the subject-level parameters
#' @field theta.init a function that provides a random initial value for each subject-level parameter
#' @field theta.start.point a numeric vector containing means of start points used to initialize in theta.init
#' @field phi.names a character vector containing the names of the group-level parameters
#' @field phi.init a function that initializes the group-level parameters
#' @field vary.parameter a logical vector containing parameters to vary
#' @field prior a list containing priors on all parameters
#' @field contaminant a list specifying two values
#' \describe{
#' \item{\code{pct}}{the percentage of the LBA distribution assumed to be due to random contaminants.}
#' \item{\code{contaminant_bound}}{the upper bound of the contaminant distribution.}
#' The contaminant distribution is assumed to be a uniform spanning from 0 to \code{contaminant_bound}.
#' }
#' @section Methods:
#' \describe{
#' \item{\code{log.dens.prior(x,hyper)}}{likelihood of subject-level parameters given group-level parameters}
#' \item{\code{log.dens.hyper(theta,phi,prior)}}{summed log-likelihood for subject given group-level parameters and group given priors}
#' \item{\code{log.dens.like(x,data,par.names)}}{LBA likelihood function}}
#' @export
LBA = R6::R6Class('Model.Hierarchical',
    inherit = LBA.Individual,
    public = list(

        prior = list(
             A = list(mu=c(mu=1,sigma=3),
                      sigma=c(mu=1,sigma=3)),

             b = list(mu=c(mu=.4,sigma=3),
                      sigma=c(mu=.4,sigma=3)),

             t0 = list(mu=c(mu=.3,sigma=1),
                       sigma=c(mu=.3,sigma=1)),

             vc = list(mu=c(mu=3,sigma=5),
                       sigma=c(mu=3,sigma=5)),

             ve = list(mu=c(mu=2,sigma=3),
                       sigma=c(mu=3,sigma=3)),

             sve = list(mu=c(mu=1,sigma=3),
                        sigma=c(mu=1,sigma=3))
        ),

        phi.names = NULL,

        phi.start.points = list(A=c(mu=2,sd=1),
                                b=c(mu=1,sd=.5),
                                t0=c(mu=.1,sd=.05),
                                vc=c(mu=4,sd=1),
                                ve=c(mu=1,sd=.5),
                                sve=c(mu=.8,sd=.4)),

        phi.init = function(){

             par_names = names(self$vary.parameter)
             phi = rep(NA,length(self$phi.names))
             for(i in 1:length(par_names)){
                  tmp=grep(par_names[i],self$phi.names)
                  phi[tmp]=msm::rtnorm(length(tmp),self$phi.start.points[[par_names[i]]]['mu'],
                                       self$phi.start.points[[par_names[i]]]['sd'],0,1)
             }
             return(phi)

        },

        log.dens.prior=function(x,hyper){
             out=0
             for (p in names(x)) out =
                       out+msm::dtnorm(x[p],hyper[paste(p,"mu",sep=".")],hyper[paste(p,"sigma",sep=".")],0,Inf,log=TRUE)
             out
        },


        log.dens.hyper=function(theta,phi,prior){
             sum((msm::dtnorm(theta,phi[1],phi[2],0,Inf,log=TRUE))) +
                  (msm::dtnorm(phi[1],prior$mu[1],prior$mu[2],0,Inf,log=TRUE)) +
                  (msm::dtnorm(phi[2],prior$mu[1],prior$sigma[2],0,Inf,log=TRUE))
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

        initialize=function(A=F,b=F,vc=F,ve=F,t0=F,sve=F,conds=NULL,prior=NULL,contaminant=list()){
             self$vary.parameter['A'] = A
             self$vary.parameter['b'] = b
             self$vary.parameter['t0'] = t0
             self$vary.parameter['vc'] = vc
             self$vary.parameter['ve'] = ve
             self$vary.parameter['sve'] = sve

             if (length(contaminant)!=0) {
                  if(is.null(contaminant$pct)){
                       stop('contaminant list must contain an element named \"pct\". Example: list(pct=1,upper.bound=5)',call. = FALSE)
                  }
                  if(is.null(contaminant$upper.bound)){
                       stop('contaminant list must contain an element named \"upper.bound\". Example: list(pct=1,upper.bound=5)',call. = FALSE)
                  }
                  if (contaminant$pct >= 0 & contaminant$pct <= 100) {
                       self$contaminant$pct = contaminant$pct
                  } else {
                       stop('Contaminant percentatge must range from 0 to 100',call. = FALSE)
                  }
                  if (contaminant$upper.bound > 0){
                       self$contaminant$upper.bound = contaminant$upper.bound
                  } else {
                       stop('Contaminant upper bound must be > 0',call. = FALSE)
                  }
             }

             if(!is.null(conds)){
                  self$conds = conds
             }
             private$make.start.points()
             private$make.theta.names()
             if(is.null(prior)){
               private$make.prior()
             }

             private$make.phi.names()
        }
        ),

    private = list(

         make.phi.names = function(){
              self$phi.names = paste(rep(self$theta.names,each=2),c("mu","sigma"),sep=".")
         },

         make.prior = function(){

              par_names = names(self$vary.parameter)
              prior = list()

              for(i in 1:length(par_names)){
                  tmp=grep(par_names[i],self$theta.names,value=TRUE)
                  for (n in 1:length(tmp)) {
                       tmp2=tmp[n]
                       prior[[tmp2]]=list(mu=c(self$prior[[paste(par_names[i])]]$mu[['mu']],
                                               self$prior[[paste(par_names[i])]]$mu[['sigma']]),
                                          sigma=c(self$prior[[paste(par_names[i])]]$sigma[['mu']],
                                                  self$prior[[paste(par_names[i])]]$sigma[['sigma']]))
                  }
              }

             self$prior = prior
         }


     )
)






