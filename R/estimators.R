#' Estimates the marginal likelihood and summarizes the samples
#'
#' Summarizes the samples or estimates the marginal likelihood depending on the type of samples
#'
#' @importFrom stats var
#' @importFrom stats aggregate
#' @param  object a Powder.Hierarchical object from \code{\link{powder}}
#' @return A \code{\link{tibble}}. If the Powder.* object was generated via the \code{standard} or \code{parallel} method from \code{\link{powder}} then the returned value will be a \code{\link{tibble}} with the following rows:
#' \describe{
#' \item{\code{TI}}{log marginal likelihood estimate based on thermodynamic integration}
#' \item{\code{TI Corrected}}{log marginal likelihood estimate based on modified thermodynamic integration using variance of log likelihood}
#' \item{\code{Harmonic Mean}}{log marginal likelihood estimate based on harmonic mean}
#' \item{\code{Steppingstone}}{marginal likelihood esimate based on steppingstone estimator. Because of underflow, this will usually not be meaningful.}
#' \item{\code{Log Steppingstone}}{log marginal likelihood estimate based on steppingstone estimator.}
#' \item{\code{TI Variance}}{estimated variance of \code{TI}}
#' \item{\code{Steppingstone Variance}}{estimated variance of \code{Log Steppingstone}}
#' }
#' @return If the object was generated via the \code{posterior} method from \code{\link{powder}} then the returned value will be a \code{\link{tibble}} with the following columns:
#' \describe{
#' \item{\code{mean}}{the mean of the parameter}
#' \item{\code{95\% HDI UB}}{upper bound of 95\% highest density interval}
#' \item{\code{95\% HDI LB}}{lower bound of 95\% highest density interval}
#' }
#' @examples
#' \dontrun{
#' model = LBA$new()
#' data('null',package='powder')
#' #low samples for illustration purposes only
#' pow.out = powder(data=null, model=model, method='parallel', burnin=10, n.samples=20)
#' summary(pow.out,options=list(sig.digits=6))
#'
#' # # A tibble: 7 x 2
#' # Method                       Value
#' # <fct>                        <dbl>
#' # 1 TI                       1822.12
#' # 2 TI Corrected             1920.92
#' # 3 Harmonic Mean            2128.10
#' # 4 Steppingstone            NaN
#' # 5 Log Steppingstone        2086.49
#' # 6 TI Variance              4144.37
#' # 7 Steppingstone Variance   0.58262
#'
#' pow.out = powder(data=null, model=model, method='posterior',burnin=10,n.samples=20)
#' summary(pow.out,options=list(sig.digits=3,num.rows=20))
#'
#' # $group.level
#' # # A tibble: 12 x 4
#' # parameter  mean         95% HDI LB         95% HDI UB
#' # <fct>     <dbl>              <dbl>               <dbl>
#' # 1 A.mu      0.769            0.0803               1.33
#' # 2 A.sigma   0.781            0.275                1.25
#' # 3 b.mu      0.750            0.0459               1.29
#' # 4 b.sigma   0.724            0.177                1.46
#' # 5 sve.mu    0.924            0.115                1.41
#' # 6 sve.sigma 0.687            0.242                1.23
#' # 7 t0.mu     0.150            0.00687              0.227
#' # 8 t0.sigma  0.122            0.0488               0.240
#' # 9 vc.mu     3.22             0.427                4.61
#' # 10 vc.sigma 1.84             0.731                3.98
#' # 11 ve.mu    1.57             0.262                3.10
#' # 12 ve.sigma 1.49             0.278                2.97
#' #
#' # $subject.level
#' # # A tibble: 60 x 5
#' # parameter subject  mean   95% HDI LB    95% HDI UB
#' # <fct>     <fct>   <dbl>        <dbl>         <dbl>
#' # 1 A         1       0.850      0.194          1.61
#' # 2 b         1       1.14       0.125          1.89
#' # 3 sve       1       1.05       0.549          1.49
#' # 4 t0        1       0.141      0.00486        0.271
#' # 5 vc        1       3.61       2.42           4.59
#' # 6 ve        1       1.60       0.238          2.84
#' # 7 A         10      0.943      0.232          1.68
#' # 8 b         10      0.950      0.340          1.44
#' # 9 sve       10      1.27       0.866          1.70
#' # 10 t0       10      0.205      0.0877         0.321
#' # # ... with 50 more rows
#' }
#' @export
summary.Powder.Hierarchical = function(object,options=list(), ...){

     pow.out =  object

     if (!is.null(options$sig.digits)) {
          options(pillar.sigfig = options$sig.digits)
     }

     if (!is.null(options$num.rows)) {
          options(tibble.print_max = options$num.rows)
     }

     prob_list = pow.out$log.like.list
     meltin = pow.out$options$meltin
     burnin = pow.out$options$burnin
     threads = pow.out$options$n.sequences
     temperatures = pow.out$options$temperatures
     n.chains = pow.out$options$n.chains
     K = pow.out$options$num.temps

     if (pow.out$options$method == 'posterior') {
          theta = reshape.theta(pow.out)
          phi = reshape.phi(pow.out)
          theta.agg = do.call(data.frame,aggregate(value ~ parameter * subject, theta, function(x) c(mean(x),hdi(x))))
          colnames(theta.agg) = c('parameter','subject','mean','95% HDI LB','95% HDI UB')
          phi.agg = do.call(data.frame,aggregate(value ~ parameter, phi, function(x) c(mean(x),hdi(x))))
          colnames(phi.agg) = c('parameter','mean','95% HDI LB','95% HDI UB')
          theta.agg = tibble::as.tibble(theta.agg)
          phi.agg = tibble::as.tibble(phi.agg)
          out = list(group.level=phi.agg,subject.level=theta.agg)
          return(out)
     }

     if (pow.out$options$method == 'standard') {
          burnin = rep(c(burnin,rep(meltin,(K/threads)-1)),threads)
          prob_list = lapply(1:length(prob_list),function(x) prob_list[[x]][burnin[x]:dim(prob_list[[x]])[1],,])
          prob_list = lapply(1:length(prob_list),function(x) c(do.call(rbind,lapply(1:length(prob_list[[x]][1,,1]),function(c) apply(prob_list[[x]][,c,],1,sum)))))

          if (pow.out$options$high.temps.first) {
               prob_list = rev(prob_list)
               temperatures = rev(temperatures)
          }

          df = tibble::as.tibble(get_estimates(prob_list,temperatures))
     }

     if(pow.out$options$method == 'parallel'){

          prob_list = sapply(1:n.chains,function(c) apply(prob_list[,c,],1,sum))
          prob_list = lapply(1:n.chains,function(x) prob_list[burnin:nrow(prob_list),x])
          prob_list = lapply(1:n.chains,function(x)prob_list[[x]][!is.na(prob_list[[x]])])

          if (pow.out$options$high.temps.first) {
               prob_list = rev(prob_list)
               temperatures = rev(temperatures)
          }

          df = tibble::as.tibble(get_estimates(prob_list,temperatures))
     }

     if (pow.out$options$method == 'wbic') {
          prob_list = prob_list[[1]]
          prob_list = sapply(1:n.chains,function(c) apply(prob_list[,c,],1,sum))
          prob_list = prob_list[burnin:nrow(prob_list),]
          df = data.frame('wbic' = mean(prob_list), 'SE' = sd(prob_list) / sqrt(length(prob_list)))
     }

     return(df)

}

#' Estimates marginal likelihood from power posterior samples
#'
#' Estimates marginal likelihood using output from \code{powder}
#'
#' @importFrom stats var
#' @importFrom stats aggregate
#' @param  object A Powder.Individual or Powder.Hierarchical object from \code{\link{powder}}
#' @return A \code{\link{tibble}}. If the Powder.* object was generated via the \code{standard} or \code{parallel} method from \code{\link{powder}} then the returned value will be a \code{\link{tibble}} with the following rows:
#' \describe{
#' \item{\code{TI}}{log marginal likelihood estimate based on thermodynamic integration}
#' \item{\code{TI Corrected}}{log marginal likelihood estimate based on modified thermodynamic integration using variance of log likelihood}
#' \item{\code{Harmonic Mean}}{log marginal likelihood estimate based on harmonic mean}
#' \item{\code{Steppingstone}}{marginal likelihood esimate based on steppingstone estimator. Because of underflow, this will usually not be meaningful.}
#' \item{\code{Log Steppingstone}}{log marginal likelihood estimate based on steppingstone estimator.}
#' \item{\code{TI Variance}}{estimated variance of \code{TI}}
#' \item{\code{Steppingstone Variance}}{estimated variance of \code{Log Steppingstone}}
#' }
#' @return If the object was generated via the \code{posterior} method from \code{\link{powder}} then the returned value will be a \code{\link{tibble}} with the following columns:
#' \describe{
#' \item{\code{mean}}{the mean of the parameter}
#' \item{\code{95\% HDI UB}}{upper bound of 95\% highest density interval}
#' \item{\code{95\% HDI LB}}{lower bound of 95\% highest density interval}
#' }
#' @examples
#' \dontrun{
#' model = LBA.Individual$new()
#' data('individual',package='powder')
#' # low number of samples for illustration purposes only
#' pow.out = powder(data=individual, model=model, method='parallel',burnin=100,n.samples=500)
#' summary(pow.out,options=list(sig.digits=5))
#'
#' # A tibble: 7 x 2
#' # Method                       Value
#' # <fct>                        <dbl>
#' # 1 TI                     5.6759e+  2
#' # 2 TI Corrected           5.6799e+  2
#' # 3 Harmonic Mean          5.6184e+  2
#' # 4 Steppingstone          6.1592e+245
#' # 5 Log Steppingstone      5.6955e+  2
#' # 6 TI Variance            3.6981e-  2
#' # 7 Steppingstone Variance 9.8905e-  2
#'
#' # do posterior estimation
#' pow.out = powder(data=individual, model=model, method='posterior',burnin=100,n.samples=500)
#' summary(pow.out,options=list(sig.digits=3))
#'
#' #   parameter  mean   95% HDI LB   95% HDI UB
#' #   <fct>     <dbl>        <dbl>        <dbl>
#' # 1 A         0.827        0.553        1.02
#' # 2 b         0.575        0.261        1.01
#' # 3 sve       1.18         0.876        1.46
#' # 4 t0        0.269        0.178        0.320
#' # 5 vc        3.44         2.95         3.90
#' # 6 ve        1.12         0.655        1.60
#' }
#' @export
summary.Powder.Individual = function(object,options=list(), ...){

     pow.out = object

     if (!is.null(options$sig.digits)) {
          options(pillar.sigfig = options$sig.digits)
     }

     if (!is.null(options$num.rows)) {
          options(tibble.print_max = options$num.rows)
     }

     prob_list = pow.out$log.like.list
     meltin = pow.out$options$meltin
     burnin = pow.out$options$burnin
     threads = pow.out$options$n.sequences
     temperatures = pow.out$options$temperatures
     K = pow.out$options$num.temps

     if (pow.out$options$method == 'posterior') {
          theta = reshape.theta(pow.out)
          theta.agg = do.call(data.frame,aggregate(value ~ parameter, theta, function(x) c(mean(x),hdi(x))))
          colnames(theta.agg) = c('parameter','mean','95% HDI LB','95% HDI UB')
          theta.agg = tibble::as.tibble(theta.agg)
          return(theta.agg)
     }

     if (pow.out$options$method == 'standard') {
          burnin = rep(c(burnin,rep(meltin,(K/threads)-1)),threads)
          prob_list = lapply(1:length(prob_list),function(x) prob_list[[x]][burnin[x]:length(prob_list[[x]][,1]),])

          if (pow.out$options$high.temps.first) {
               prob_list = rev(prob_list)
               temperatures = rev(temperatures)
          }
          df = tibble::as.tibble(get_estimates(prob_list,temperatures))
     }

     if (pow.out$options$method == 'parallel') {

          prob_list = lapply(1:ncol(prob_list),function(x) prob_list[burnin:nrow(prob_list),x])
          if (pow.out$options$high.temps.first) {
               prob_list = rev(prob_list)
               temperatures = rev(temperatures)
          }

          df = tibble::as.tibble(get_estimates(prob_list,temperatures))
     }

     if (pow.out$options$method == 'wbic') {

          prob_list = prob_list[[1]]
          prob_list = prob_list[burnin:nrow(prob_list),]
          df = data.frame('wbic' = mean(prob_list), 'SE' = sd(prob_list) / sqrt(length(prob_list)))

     }

     return(df)

}



get_estimates = function(prob_list,t){
     mean_lp = sapply(1:length(t),function(x)mean(prob_list[[x]][!is.infinite(prob_list[[x]])]))
     mean_var = sapply(1:length(t),function(x)var(prob_list[[x]][!is.infinite(prob_list[[x]])]))
     ti = thermo_int(t,mean_lp,mean_var)
     mean_lp_list = lapply(1:length(prob_list[[1]]),function(y) sapply(1:length(t),function(x) prob_list[[x]][y]))
     ti_all = do.call(rbind,lapply(mean_lp_list,function(x) thermo_int(t,x,mean_var)$ti))
     ti_var = var(ti_all[!is.infinite(ti_all)])/(length(ti_all[!is.infinite(ti_all)]))
     ss_var = var_stepping_stone(prob_list,t)
     hm = harmonic_mean(prob_list[[which(t==1)]])
     step_stone = stepping_stone(prob_list,t)
     step_stone_stable = stepping_stone_stable(prob_list,t)
     df = data.frame('Method'=c('TI','TI Corrected','Harmonic Mean','Steppingstone','Log Steppingstone','TI Variance','Steppingstone Variance'),
                     'Value'=c(ti$ti,ti$ti_cor,hm,step_stone,step_stone_stable,ti_var,ss_var))
     return(df)
}

hdi = function( sampleVec , credMass=0.95 ) {
     # Computes highest density interval from a sample of representative values,
     #   estimated as shortest credible interval.
     # Arguments:
     #   sampleVec
     #     is a vector of representative values from a probability distribution.
     #   credMass
     #     is a scalar between 0 and 1, indicating the mass within the credible
     #     interval that is to be estimated.
     # Value:
     #   HDIlim is a vector containing the limits of the HDI
     sortedPts = sort( sampleVec )
     ciIdxInc = ceiling( credMass * length( sortedPts ) )
     nCIs = length( sortedPts ) - ciIdxInc
     ciWidth = rep( 0 , nCIs )
     for ( i in 1:nCIs ) {
          ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
     }
     HDImin = sortedPts[ which.min( ciWidth ) ]
     HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
     HDIlim = c( HDImin , HDImax )
     return( HDIlim )
}

thermo_int = function(t,lp,var){

     # Computes the log evidence using thermodynamic integration
     # Computes a corrected estimate using the variance (Friel, Hurn, Wyse, 2014)
     #
     # Args:
     #   t: vector of temperatures from 0 to 1
     #   lp: mean likelihood at each value of t
     #   var: variance of likelihod
     #
     # Returns:
     #   A list containing the estimated log evidence and a corrected estimate

     I = sum(sapply(1:(length(t)-1),function(x) (t[x+1] - t[x]) * ((lp[x] + lp[x+1])/2)))
     cor = sum(sapply(1:(length(t)-1),function(x) (((t[x+1] - t[x])^2)/12) * (var[x+1] - var[x])))

     return(list(ti=I,ti_cor=I-cor))
}


harmonic_mean = function(lp){
     # Computes the log evidence using the harmonic mean estimator
     #
     # Args:
     #   lp: the likelihood under each posterior sample
     #
     # Returns:
     #   The harmonic mean estimate of the marginal likelihood
     hm = 1/mean(1/lp)
     return(hm)
}


stepping_stone = function(prob_list,t){
     pD = prod(sapply(1:(length(prob_list)-1),function(k) mean(exp(prob_list[[k]])^(t[k+1]-t[k]),na.rm=T) ))
     return(pD)
}

stepping_stone_2 = function(prob_list,t){
     prob_list = lapply(prob_list,function(x) exp(x))
     rss = sapply(2:length(t),
                  function(k)(max(prob_list[[k-1]])^(t[k]-t[k-1])*sum((prob_list[[k-1]]/max(prob_list[[k-1]]))^(t[k]-t[k-1]))) / length(prob_list[[1]]))
     return(rss)
}

var_stepping_stone = function(prob_list,t){
     rss = stepping_stone_2(prob_list,t)
     prob_list = lapply(prob_list,function(x) exp(x))
     K = which(!is.nan(rss))
     if(length(K)==0){
          return(NA)
     }else{
          out=sum(sapply(K+1,function(k)sum((((prob_list[[k-1]]^(t[k]-t[k-1]))/rss[k-1])-1)^2))) / (length(prob_list[[1]])^2)
     }
     return(out)
}



stepping_stone_stable = function(prob_list,t){

     # Computes the log evidence using the stepping stone estimator
     #
     # Args:
     #   prob_list: a list of length t, each element containing a
     #              vector of log-likelihoods under the posterior at
     #              each temperature in t
     #
     # Returns:
     #   The stepping stone estimate of the marginal likelihood

     ss = sum(sapply(1:(length(prob_list)-1),
                     function(k) log(mean(exp((prob_list[[k]] -
                                                    max(prob_list[[k]])) *
                                                   (t[k+1]-t[k])))) +
                          (t[k+1]-t[k]) *
                          max(prob_list[[k]])))
     return(ss)
}



