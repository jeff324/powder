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
#' @examples
#' \dontrun{
#' # DDM model that varies threshold and drift rate across 3 conditions
#' model = DDM.Individual$new(a=T,v=T,conds=1:3)
#'
#' # Note, inter-trial-variability components, st0, sv, and sz
#' have 3 options: TRUE, FALSE, or 0. When TRUE,
#' # they will vary across conditions. When FALSE
#' they are fit, but not allowed to vary across conditions.
#' By default, they are set to 0.
#'
#' #st0 is a parameter, other variability parameters are 0
#' model_st0 = DDM.Individual$new(st0=FALSE)
#'
#' #sv is allowed to vary over conditions
#' model_sv = DDM.Individual$new(sv=TRUE)
#'
#'  # Setting priors:
#'  prior = list(
#'  a = c(mu=1,sigma=1),
#'  v = c(mu=1,sigma=1),
#'  t0 = c(mu=.3,sigma=.3),
#'  z = c(mu=.5,sigma=.5),
#'  sz = c(mu=.2,sigma=.2),
#'  sv = c(mu=.1,sigma=.1))
#'
#'  #a and v vary across conditions
#'  #sz and sv are parameters and do not vary across conditions
#'  #st0 is 0
#'  model = DDM.Individual$new(a=TRUE,v=TRUE)
#'
#'  # Priors can also be changed after creating the model
#'  model$prior$a = c(mu=2,sigma=2)
#'
#' }
#' @section Internal Fields:
#' @field theta.names a character vector containing the names of the subject-level parameters
#' @field theta.init a function that provides a random initial value for each subject-level parameter
#' @field theta.start.point a numeric vector containing means of start points used to initialize in theta.init
#' @field vary.parameter a logical vector containing parameters to vary
#' @field prior a list containing priors on all parameters
#' @section Internal Methods:
#' \describe{
#' \item{\code{log.dens.prior(x,hyper)}}{likelihood of subject-level parameters given group-level parameters}
#' \item{\code{log.dens.like(x,data,par.names)}}{LBA likelihood function}}
#' @export
DDM.Individual = R6::R6Class('Model.Individual',


    public = list(
         conds = c(1,2),

         prior = list(
              a = c(mu=1,sigma=1),
              v = c(mu=1,sigma=1),
              t0 = c(mu=.3,sigma=.3),
              z = c(mu=.5,sigma=.5),
              sz = c(mu=.2,sigma=.2),
              sv = c(mu=.1,sigma=.1),
              st0 = c(mu=.1,sigma=.1)
         ),

         constant = list(),

         vary.parameter = c(a=FALSE,
                            v=FALSE,
                            t0=FALSE,
                            z=FALSE,
                            sz=FALSE,
                            sv=FALSE,
                            st0=FALSE),

         theta.names = NULL,

         theta.start.points = c(a=1,v=1,t0=.3,z=.5,sz=.2,sv=1,st0=.1),

         theta.init = function(){
              msm::rtnorm(length(self$theta.start.points),
                          mean=self$theta.start.points,
                          sd=self$theta.start.points/5,0,Inf)
         },

         log.dens.prior=function(x,hyper){
              out=0
              trunc_norm_pars = names(x)[-grep('^v',names(x))]

              for (p in trunc_norm_pars) {
                   out = msm::dtnorm(x[p],hyper[[p]][1],hyper[[p]][2],0,Inf,log=TRUE) + out
              }

              out = dnorm(x['v'],hyper[['v']][1],hyper[['v']][2]) + out

              return(out)
         },

         log.dens.like=function(x,data,par.names){
              names(x)=par.names
              if (!self$vary.parameter['a']) {
                   a=x['a']
              }

              if (!self$vary.parameter['v']) {
                   v=x['v']
              }

              if (!self$vary.parameter['t0']) {
                   t0=x['t0']
              }

              if (!self$vary.parameter['z']) {
                   if(is.null(self$constant$z)){
                        z=x['z']
                   }else{
                        z=self$constant$z
                   }
              }

              if (!self$vary.parameter['sz']) {
                   if(is.null(self$constant$sz)){
                        sz=x['sz']
                   }else{
                        sz=self$constant$sz
                   }
              }

              if (!self$vary.parameter['sv']) {
                   if(is.null(self$constant$sv)){
                        sv=x['sv']
                   }else{
                        sv=self$constant$sv
                   }
              }

              if (!self$vary.parameter['st0']) {
                   if(is.null(self$constant$st0)){
                        st0=x['st0']
                   }else{
                        st0=self$constant$st0
                   }
              }

              out=0

              for (cond in self$conds) {

                   if (self$vary.parameter['a']) {
                        a=x[paste('a',cond,sep=".")]
                   }
                   if (self$vary.parameter['v']) {
                        v=x[paste('v',cond,sep=".")]
                   }
                   if (self$vary.parameter['t0']) {
                        t0=x[paste('t0',cond,sep=".")]
                   }
                   if (self$vary.parameter['z']) {
                        z=x[paste('z',cond,sep=".")]
                   }
                   if (self$vary.parameter['sz']) {
                        sz=x[paste('sz',cond,sep=".")]
                   }
                   if (self$vary.parameter['sv']) {
                        sv=x[paste('sv',cond,sep=".")]
                   }
                   if (self$vary.parameter['st0']) {
                        st0=x[paste('st0',cond,sep=".")]
                   }

                   tmp=data$Cond==cond
                   tmp=rtdists::ddiffusion(rt=data$Time[tmp],response=data$Correct[tmp],
                                           a=a,v=v,t0=t0,z=z+.5*a,sz=sz,sv=sv,st0=st0)

                   out=out+sum(log(pmax(tmp,1e-10)))
              }
              return(out)
         },

         initialize = function(a=FALSE,v=FALSE,t0=FALSE,z=0,sz=0,sv=0,st0=0,conds=NULL,prior=NULL){

              self$vary.parameter['a'] = a
              self$vary.parameter['v'] = v
              self$vary.parameter['t0'] = t0
              self$vary.parameter['z'] = z
              self$vary.parameter['sz'] = sz
              self$vary.parameter['sv'] = sv
              self$vary.parameter['st0'] = st0

              if (as.character(z) == '0') {
                   if(is.null(prior$z)){
                        self$prior[['z']] = 0
                        self$constant$z = 0
                   }else{
                        if (length(prior$z) == 1 & prior$z[1] == 0) {
                             self$prior[['z']] = 0
                             self$constant$z = 0
                        }
                   }
              }


              if (as.character(sz) == '0') {
                   if(is.null(prior$sz)){
                        self$prior[['sz']] = 0
                        self$constant$sz = 0
                   }else{
                        if (length(prior$sz) == 1 & prior$sz[1] == 0) {
                             self$prior[['sz']] = 0
                             self$constant$sz = 0
                        }
                   }
              }

              if (as.character(sv) == '0') {
                   if(is.null(prior$sv)){
                        self$prior[['sv']] = 0
                        self$constant$sv = 0
                   }else{
                        if (length(prior$sv) == 1 & prior$sv[1] == 0) {
                             self$prior[['sv']] = 0
                             self$constant$sv = 0
                        }
                   }
              }

              if (as.character(st0) == '0') {
                   if(is.null(prior$st0)){
                        self$prior[['st0']] = 0
                        self$constant$st0 = 0
                   }else{
                        if (length(prior$st0) == 1 & prior$st0[1] == 0) {
                             self$prior[['st0']] = 0
                             self$constant$st0 = 0
                        }
                   }
              }

              if (!is.null(conds)) {
                   self$conds = conds
              }

              private$make.theta.names()
              private$make.start.points()
              if (is.null(prior)) {
                   private$make.prior.default()
              } else {
                   private$make.prior(prior)
              }

         }
    ),

    private = list(

         make.start.points = function(){

              theta.start.points = rep(NA,length(self$theta.names))

              for(i in 1:length(self$theta.start.points)){
                   idx = grep(paste0('^',names(self$theta.start.points)[i]),self$theta.names)
                   theta.start.points[idx] = self$theta.start.points[i]
              }

              self$theta.start.points = theta.start.points
         },

         make.prior.default = function(){

              par_names = names(self$vary.parameter)
              prior = list()

              for(i in 1:length(par_names)){
                   tmp=grep(par_names[i],self$theta.names,value=TRUE)
                   if(length(tmp) != 0){
                        for (n in 1:length(tmp)) {
                             tmp2=tmp[n]
                             prior[[tmp2]]=c(mu=self$prior[[paste(par_names[i])]][['mu']],
                                             sigma=self$prior[[paste(par_names[i])]][['sigma']])
                        }
                   }
              }

              self$prior = prior
         },

         make.prior = function(prior){

              par_names = names(self$vary.parameter)
              prior_tmp = list()

              for(i in 1:length(par_names)){
                   tmp=grep(par_names[i],self$theta.names,value=TRUE)
                   if(length(tmp) != 0){
                        for (n in 1:length(tmp)) {
                             tmp2=tmp[n]
                             prior_tmp[[tmp2]]=c(mu=prior[[paste(par_names[i])]][['mu']],
                                             sigma=prior[[paste(par_names[i])]][['sigma']])
                        }
                   }
              }

              self$prior = prior_tmp
         },

         make.theta.names = function(){
              theta.names=NULL

              zeroed = sapply(1:length(self$prior),function(x) any(length(self$prior[[names(self$prior)[x]]]) == 1 &
                                   self$prior[[names(self$prior)[x]]] == 0))

              for (i in 1:length(self$vary.parameter)) {

                   if (self$vary.parameter[i] & !zeroed[i]) {
                        theta.names=c(theta.names,paste(names(self$vary.parameter)[i],self$conds,sep="."))
                   } else if (self$vary.parameter[i] & zeroed[i]) {
                        warning(paste(names(self$vary.parameter)[i],'\n'))
                        warning('Prior is a constant, but parameter is to be varied. Parameter will not be varied unless prior is changed.')
                   } else if (!self$vary.parameter[i] & !zeroed[i]) {
                        theta.names=c(theta.names,names(self$vary.parameter)[i])
                   }

              }

              self$theta.names = theta.names

         }
    )

)
