#' Individual subject Linear Ballistic Accumulator model (Brown & Heathcote, 2008)
#'
#' @docType class
#' @importFrom R6 R6Class
#' @importFrom stats pnorm
#' @importFrom stats dnorm
#' @return Object of \code{\link{R6Class}} with methods for doing sampling with \code{\link{powder}}
#' @format \code{\link{R6Class}} object.
#' @param A threshold seperation
#' @param b threshold
#' @param vc correct drift rate
#' @param ve error drift rate
#' @param t0 non-decision time
#' @param sve inter-trial-variability of error drift rate. Note that correct drift rate variability is fixed at 1.
#' @param prior list containing priors
#' @param conds numeric vector containing values corresponding to conditions
#' @param correct.response a character vector. Specifies the value corresponding to the correct response
#' @param error.response a character vector. Specifies the value corresponding to the error response
#' @param rt.column a character vector. Specifies the name of the column in the data for response times.
#' @param response.column a character vector. Specifies the name of the column in the data for respones.
#' @param condition.column a character vector. Specifies the name of the column in the data for the conditions.
#' @param contaminant a list specifying two values
#' \describe{
#' \item{\code{pct}}{the percentage of the LBA distribution assumed to be due to random contaminants.}
#' \item{\code{contaminant_bound}}{the upper bound of the contaminant distribution.}
#' The contaminant distribution is assumed to be a uniform spanning from 0 to \code{contaminant_bound}.
#' }
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
#' @field correct.response is a character or numeric that corresponds to the correct response in data
#' \describe{
#' \item{\code{pct}}{the percentage (ranging from 0 to 100) of the LBA distribution assumed to be due to random contaminants.}
#' \item{\code{contaminant_bound}}{the upper bound of the contaminant distribution.}
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

         correct.response = NULL,

         error.response = NULL,

         rt.column = NULL,

         response.column = NULL,

         condition.column = NULL,

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
                   tmp=private$get.dens.2choice(rt=data$Time[tmp],crct=data$Correct[tmp]==self$correct.response,b=b,A=A,v=vs,s=s,t0=t0)

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

         predict = function(pow.out,conds=NULL,thin=1,n=1){

              if(is.null(conds)){
                   conds = self$conds
              }

              out = private$predict_theta(theta=pow.out$theta,conds=conds,thin=thin,n=n)

              return(out)

         },

        initialize = function(A=F,b=F,vc=F,ve=F,t0=F,sve=F,conds=NULL,prior=NULL,contaminant=list(),
                              correct.response = NULL,
                              error.response = NULL,
                              rt.column = NULL,
                              response.column = NULL,
                              condition.column = NULL){

             self$vary.parameter['A'] = A
             self$vary.parameter['b'] = b
             self$vary.parameter['t0'] = t0
             self$vary.parameter['vc'] = vc
             self$vary.parameter['ve'] = ve
             self$vary.parameter['sve'] = sve

             if (is.null(correct.response)) {
                  self$correct.response = 1
                  warning('Value for correct response not specified, using 1 as default. If this is not correct, please change in data.',call.=F)
             } else {
                  self$correct.response = correct.response
             }

             if (is.null(error.response)) {
                  self$error.response = 2
                  warning('Value for error response not specified, using 2 as default. If this is not correct, please change in data.', call. = F)
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

         predict_theta = function(theta,conds=NULL,thin=1,n=1,verbose=TRUE){

              theta = private$prettify_theta(theta,conds,thin)

              if(verbose){
                   progress_type = 'text'
              } else {
                   progress_type = 'none'
              }
              preds = lapply(theta,function(x) plyr::ldply(1:nrow(x), function(y)
                   rtdists::rLBA(n=n,
                                 A=x[y,'A'],
                                 b=x[y,'A']+x[y,'b'],
                                 t0=x[y,'t0'],
                                 mean_v=c(x[y,'vc'],x[y,'ve']),
                                 sd_v=c(1,x[y,'sve']),
                                 silent = T),
                   .progress = progress_type)
              )
              preds = Map(cbind,preds,condition=conds)
              preds = do.call(rbind,preds)
              return(preds)

         },

         prettify_theta = function(theta,conds=NULL,thin=1){
              theta = theta[,,seq(1,length(theta[1,1,]),by=thin)]
              pars = colnames(theta)
              n.iter = length(theta[1,1,])
              n.chains = length(theta[,1,1])
              theta.list = lapply(1:length(pars), function(x) cbind('iteration'=1:n.iter,'parameter'=pars[x],t(theta[,pars[x],])))
              theta = data.frame(do.call(rbind,theta.list),stringsAsFactors = FALSE)
              theta.chains = theta[,3:ncol(theta)]
              theta.chains[] = lapply(theta.chains,as.numeric)
              theta[,3:ncol(theta)] = theta.chains
              colnames(theta) = c('iteration','parameter',1:n.chains)
              theta$iteration = as.numeric(theta$iteration)
              theta = reshape2::melt(theta,id.vars = c('parameter','iteration'), value.name='value', variable.name = 'chain')
              theta = reshape2::dcast(theta,iteration + chain ~ parameter, value.var = 'value')
              constants = c(1:length(colnames(theta)))[-grep('[.]',colnames(theta))]
              if (length(constants) > 0) {
                   varied_pars = lapply(conds,function(x) grep(paste0('.',x),colnames(theta)))
                   theta = lapply(varied_pars,function(x) theta[,c(x,constants)])
                   cond_df_names = lapply(theta,function(x) colnames(x))
                   cond_df_names = lapply(1:length(conds),function(x) gsub(paste0('.',conds[x]),'',cond_df_names[[x]]))
                   for(i in 1:length(cond_df_names)){
                        colnames(theta[[i]]) = cond_df_names[[i]]
                   }
              } else {
                   theta = list(theta)
              }
              return(theta)
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
                   for (i in 2:N) tmp[,i-1]=private$fptcdf(z=t,x0max=x0max[i],chi=chi[i],
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
