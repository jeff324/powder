#' Individual subject Linear Ballistic Accumulator model (Brown & Heathcote, 2008)
#'
#' @docType class
#' @importFrom R6 R6Class
#' @importFrom stats pnorm
#' @importFrom stats dnorm
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
#' @field vary.parameter a logical vector containing parameters to vary
#' @field prior a list containing priors on all parameters
#' @field contaminant a list specifying two values
#' \describe{
#' \item{\code{pct}} the percentage (ranging from 0 to 100) of the LBA distribution assumed to be due to random contaminants.
#' \item{\code{contaminant_bound}} the upper bound of the contaminant distribution.
#' The contaminant distribution is assumed to be a uniform spanning from 0 to \code{contaminant_bound}.
#' }
#' @section Methods:
#' \describe{
#' \item{\code{log.dens.prior(x,hyper)}}{likelihood of subject-level parameters given group-level parameters}
#' \item{\code{log.dens.like(x,data,par.names)}}{LBA likelihood function}}
#' @export
LBA.Individual = R6::R6Class('Model.Individual',

    public = list(

         contaminant = list('pct' = 0, 'upper.bound' = NULL),

         conds = c(1,2),

         prior = list(
              A = c(mu=1,sigma=1),
              b = c(mu=.4,sigma=.4),
              t0 = c(mu=.3,sigma=.3),
              vc = c(mu=3,sigma=3),
              ve = c(mu=1,sigma=1),
              sve = c(mu=1,sigma=1)
         ),

         vary.parameter = c(A=FALSE,
                            b=FALSE,
                            t0=FALSE,
                            vc=FALSE,
                            ve=FALSE,
                            sve=FALSE),

         theta.names = NULL,

         theta.start.points = c(A=1,b=1,t0=.2,vc=3,ve=2,sve=1),

         theta.init = function(){
              msm::rtnorm(length(self$theta.start.points),
                          mean=self$theta.start.points,
                          sd=self$theta.start.points/5,0,Inf)
         },

         log.dens.prior=function(x,hyper){
              out=0
              for (p in names(x)) out =
                        out+msm::dtnorm(x[p],hyper[[p]][1],hyper[[p]][2],0,Inf,log=TRUE)
              out
         },

         log.dens.like=function(x,data,par.names){
              names(x)=par.names
              if (!self$vary.parameter['A']) {
                   A=x["A"]
              }
              if (!self$vary.parameter['b'] & !self$vary.parameter['A']) {
                   b=x["b"] + x["A"]
              }
              if (!self$vary.parameter['t0']) {
                   t0=x["t0"]
              }
              if (!self$vary.parameter['vc'] & !self$vary.parameter['ve']) {
                   vs=c(x["vc"],x["ve"])
              }
              if (!self$vary.parameter['sve']) {
                   s=c(1,x["sve"])
              }

              out=0

              for (cond in self$conds) {

                   if (self$vary.parameter['A']) {
                        A=x[paste("A",cond,sep=".")]
                   }
                   if (self$vary.parameter['b'] & self$vary.parameter['A']) {
                        b=x[paste("b",cond,sep=".")] + x[paste("A",cond,sep=".")]
                   } else if (self$vary.parameter['b'] & !self$vary.parameter['A']) {
                        b=x[paste("b",cond,sep=".")] + x["A"]
                   } else if (!self$vary.parameter['b'] & self$vary.parameter['A']) {
                        b=x["b"] + x[paste("A",cond,sep=".")]
                   }
                   if (self$vary.parameter['t0']) {
                        t0=x[paste("t0",cond,sep=".")]
                   }
                   if (self$vary.parameter['vc'] & self$vary.parameter['ve']) {
                        vs=c(x[paste("vc",cond,sep=".")],x[paste("ve",cond,sep=".")])
                   } else if (!self$vary.parameter['vc'] & self$vary.parameter['ve']) {
                        vs=c(x["vc"],x[paste("ve",cond,sep=".")])
                   } else if (self$vary.parameter['vc'] & !self$vary.parameter['ve']) {
                        vs=c(x[paste("vc",cond,sep=".")],x["ve"])
                   }
                   if (self$vary.parameter['sve']) {
                        s=c(1,x[paste("sve",cond,sep=".")])
                   }

                   tmp=data$Cond==cond
                   tmp=private$get.dens.2choice(rt=data$Time[tmp],crct=data$Correct[tmp]==1,b=b,A=A,v=vs,s=s,t0=t0)

                   ## density contamination
                   if (self$contaminant$pct > 0) {
                        contaminant_pct = self$contaminant$pct * .01 #percent to decimal
                        tmp = (1-contaminant_pct) * tmp + ((1/self$contaminant$upper.bound) * contaminant_pct)/2
                        #tmp=0.98*tmp+0.01/5
                   }

                   out=out+sum(log(pmax(tmp,1e-10)))
              }
              return(out)
         },

        initialize = function(A=F,b=F,vc=F,ve=F,t0=F,sve=F,conds=NULL,prior=NULL,contaminant=list()){
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

             if (!is.null(conds)) {
                  self$conds = conds
             }
             private$make.start.points()
             private$make.theta.names()
             if (is.null(prior)) {
                  private$make.prior()
             }
        }
        ),

    private = list(

         make.start.points = function(){
              n.cond = length(self$conds)

              start.points = sapply(names(self$vary.parameter),
                                    function(x){
                                         rep(self$theta.start.points[[x]],ifelse(self$vary.parameter[[x]],n.cond,1))
                                    })
              self$theta.start.points = unlist(start.points)
         },

         make.prior = function(){

              par_names = names(self$vary.parameter)
              prior = list()

              for (i in 1:length(par_names)) {
                   tmp=grep(par_names[i],self$theta.names,value=TRUE)
                   for (n in 1:length(tmp)) {
                        tmp2=tmp[n]
                        prior[[tmp2]]=c(self$prior[[paste(par_names[i])]][['mu']],
                                        self$prior[[paste(par_names[i])]][['sigma']])
                   }
              }

              self$prior = prior
         },

         make.theta.names = function(){
              theta.names=NULL

              for (i in 1:length(self$vary.parameter)) {
                   if (self$vary.parameter[i]) {
                        theta.names=c(theta.names,paste(names(self$vary.parameter)[i],self$conds,sep="."))
                   } else {
                        theta.names=c(theta.names,names(self$vary.parameter)[i])
                   }

              }

              self$theta.names = theta.names

         },

         fptcdf=function(z,x0max,chi,driftrate,sddrift) {
              if (x0max<1e-10) return(pnorm(chi/z,mean=driftrate,sd=sddrift,lower.tail=F))
              zs=z*sddrift ; zu=z*driftrate ; chiminuszu=chi-zu ; xx=chiminuszu-x0max
              chizu=chiminuszu/zs ; chizumax=xx/zs
              tmp1=zs*(dnorm(chizumax)-dnorm(chizu))
              tmp2=xx*pnorm(chizumax)-chiminuszu*pnorm(chizu)
              1+(tmp1+tmp2)/x0max
         },

         fptpdf=function(z,x0max,chi,driftrate,sddrift) {
              if (x0max<1e-10) return( (chi/z^2)*dnorm(chi/z,mean=driftrate,sd=sddrift))
              zs=z*sddrift ; zu=z*driftrate ; chiminuszu=chi-zu
              chizu=chiminuszu/zs ; chizumax=(chiminuszu-x0max)/zs
              (driftrate*(pnorm(chizu)-pnorm(chizumax)) +
                        sddrift*(dnorm(chizumax)-dnorm(chizu)))/x0max
         },

         n1PDFfixedt0=function(t,x0max,chi,drift,sdI,truncdrifts=TRUE) {
              # Generates defective PDF for responses on node #1.
              # "truncdrifts" sets whether that part of the multi-variate
              # normal distribution on drift rates which would otherwise
              # lead to non-terminating trials is truncated.
              N=length(drift) # Number of responses.
              if (N>2) {
                   tmp=array(dim=c(length(t),N-1))
                   for (i in 2:N) tmp[,i-1]=fptcdf(z=t,x0max=x0max[i],chi=chi[i],
                                                   driftrate=drift[i],sddrift=sdI[i])
                   G=apply(1-tmp,1,prod)
              } else {
                   G=1-private$fptcdf(z=t,x0max=x0max[2],chi=chi[2],driftrate=drift[2],sddrift=sdI[2])
              }
              out=G*private$fptpdf(z=t,x0max=x0max[1],chi=chi[1],driftrate=drift[1],sddrift=sdI[1])
              if (truncdrifts) {
                   out=out/(1-prod(pnorm(-drift/sdI)))
                   out[t<=0]=0
                   return(out)
              } else {
                   return(out)
              }
         },

         get.dens.2choice=function(rt,crct,b,A,v,s,t0){
              out=numeric(length(rt))
              out[crct] =private$n1PDFfixedt0(rt[crct]-t0,x0max=c(A,A),chi=c(b,b),drift=v,sdI=s)
              out[!crct]=private$n1PDFfixedt0(rt[!crct]-t0,x0max=c(A,A),chi=c(b,b),drift=v[2:1],sdI=s[2:1])
              pmax(out,1e-10)
         }
    )

)
