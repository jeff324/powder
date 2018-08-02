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

#' Estimates marginal likelihood from power posterior samples
#'
#' \code{marginal.likelihood} Estimates marginal likelihood using output from \code{powder}
#'
#' @importFrom stats var
#' @param  pow.out The output from the \code{powder} function
#' @return A data frame with the following columns
#' @return \code{TI} Log marginal likelihood estimate based on thermodynamic integration
#' @return \code{TI Corrected} Log marginal likelihood estimate based on modified thermodynamic integration using variance of log likelihood
#' @return \code{Harmonic Mean} Log marginal likelihood estimate based on harmonic mean
#' @return \code{Steppingstone} Marginal likelihood esimate based on steppingstone estimator. Because of underflow, this will usually not be meaningful.
#' @return \code{Log Steppingstone} Log marginal likelihood estimate based on steppingstone estimator.
#' @return \code{TI Variance} Estimated variance of \code{TI}
#' @return \code{Steppingstone Variance} Estimated variance of \code{Log Steppingstone}
#'
#' @examples
#' \dontrun{
#' model = LBA$new(b=T)
#' data('null',package='powder')
#' num.temps = 30
#' pow_out = powder(data=null, model=model, num.temps=num.temps)
#' }
#' @export
marginal.likelihood = function(pow.out) UseMethod('marginal.likelihood')

#' @export
marginal.likelihood.Powder.Hierarchical = function(pow.out){

     prob_list = pow.out$log.like.list
     meltin = pow.out$options$meltin
     burnin = pow.out$options$burnin
     threads = pow.out$options$n.sequences
     temperatures = pow.out$options$temperatures
     n.chains = pow.out$options$n.chains
     K = pow.out$options$num.temps

     if (pow.out$options$method == 'standard') {
          burnin = rep(c(burnin,rep(meltin,(K/threads)-1)),threads)
          prob_list = lapply(1:length(prob_list),function(x) prob_list[[x]][burnin[x]:dim(prob_list[[x]])[1],,])
          prob_list = lapply(1:length(prob_list),function(x) c(do.call(rbind,lapply(1:length(prob_list[[x]][1,,1]),function(c) apply(prob_list[[x]][,c,],1,sum)))))

          if (pow.out$options$high.temps.first) {
               prob_list = rev(prob_list)
               temperatures = rev(temperatures)
          }

          df = get_estimates(prob_list,temperatures)
     }

     if(pow.out$options$method == 'parallel'){
          prob_list = sapply(1:n.chains,function(c) apply(prob_list[,c,],1,sum))
          prob_list = lapply(1:n.chains,function(x) prob_list[burnin:nrow(prob_list),x])
          prob_list = lapply(1:n.chains,function(x)prob_list[[x]][!is.na(prob_list[[x]])])

          if (pow.out$options$high.temps.first) {
               prob_list = rev(prob_list)
               temperatures = rev(temperatures)
          }

          df = get_estimates(prob_list,temperatures)
     }

     return(df)

}

#' @export
marginal.likelihood.Powder.Individual = function(pow.out){

     prob_list = pow.out$log.like.list
     meltin = pow.out$options$meltin
     burnin = pow.out$options$burnin
     threads = pow.out$options$n.sequences
     temperatures = pow.out$options$temperatures
     K = pow.out$options$num.temps

     if (pow.out$options$method == 'standard') {
          burnin = rep(c(burnin,rep(meltin,(K/threads)-1)),threads)
          prob_list = lapply(1:length(prob_list),function(x) prob_list[[x]][burnin[x]:length(prob_list[[x]][,1]),])

          if (pow.out$options$high.temps.first) {
               prob_list = rev(prob_list)
               temperatures = rev(temperatures)
          }
          df = get_estimates(prob_list,temperatures)
     }

     if (pow.out$options$method == 'parallel') {
          prob_list = lapply(1:ncol(prob_list),function(x) prob_list[burnin:nrow(prob_list),x])
          if (pow.out$options$high.temps.first) {
               prob_list = rev(prob_list)
               temperatures = rev(temperatures)
          }
          df = get_estimates(prob_list,temperatures)
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


