#' Differential Evolution functions
#'
#' @importFrom stats runif
#' @keywords internal
crossover=function(i,n.chains,b,pars,use.theta,use.like,data,hyper,par.names,temperature,model){
     #require(msm)
     use.weight=use.like[i]*temperature + model$log.dens.prior(use.theta[i,],hyper[i,])
     gamma = 2.38/sqrt(2*length(pars))
     index=sample(c(1:n.chains)[-i],2,replace=F)
     theta=use.theta[i,]
     theta[pars]=use.theta[i,pars] + gamma*(use.theta[index[1],pars]-use.theta[index[2],pars]) + runif(1,-b,b)
     #  theta=matrix(theta,1,length(theta))
     like=model$log.dens.like(theta,data,par.names=par.names)*temperature
     weight=like + model$log.dens.prior(theta,hyper[i,])
     if(is.na(weight))weight=-Inf
     if(runif(1) < exp(weight-use.weight)) {
          use.theta[i,]=theta
          use.like[i]=like
     }
     c(use.like[i],use.theta[i,])
}

crossover_hyper=function(i,n.chains,b,pars,use.theta,use.phi,prior,model){
     #require(msm)
     use.weight=model$log.dens.hyper(use.theta[i,],use.phi[i,pars],prior)
     gamma = 2.38/sqrt(2*length(pars))
     index=sample(c(1:n.chains)[-i],2,replace=F)
     phi=use.phi[i,]
     phi[pars]=use.phi[i,pars] + gamma*(use.phi[index[1],pars]-use.phi[index[2],pars]) + runif(1,-b,b)
     #  phi=matrix(phi,1,length(phi))
     weight=model$log.dens.hyper(use.theta[i,],phi[pars],prior)
     if(is.na(weight))weight=-Inf
     if(runif(1) < exp(weight-use.weight)) use.phi[i,]=phi
     use.phi[i,]
}

migration.crossover=function(pars,n.chains,use.theta,use.like,data,hyper,par.names,model){
     # migration step
     n.migration.chains=ceiling(runif(1,0,n.chains))
     use.chains=sample(1:n.chains,n.migration.chains)
     migration.use.weight=rep(NA,n.migration.chains)
     migration.weight=rep(NA,n.migration.chains)
     for (mi in 1:n.migration.chains) {
          migration.use.weight[mi]=use.like[use.chains[mi]] + model$log.dens.prior(use.theta[use.chains[mi],],hyper[use.chains[mi],])
          newChain = mi - 1
          if (newChain == 0) newChain = n.migration.chains
          migration.weight[mi]=use.like[use.chains[newChain]] + model$log.dens.prior(use.theta[use.chains[newChain],],hyper[use.chains[mi],])
          if(runif(1) < exp(migration.weight[mi]-migration.use.weight[mi])) {
               use.theta[use.chains[mi],]=use.theta[use.chains[newChain],]
               use.like[use.chains[mi]]=use.like[use.chains[newChain]]
          }
     }
     cbind(use.like,use.theta)
}

migration.crossover_hyper=function(pars,n.chains,use.theta,use.phi,prior,model){
     # migration step
     n.migration.chains=ceiling(runif(1,0,n.chains))
     use.chains=sample(1:n.chains,n.migration.chains)
     migration.use.weight=rep(NA,n.migration.chains)
     migration.weight=rep(NA,n.migration.chains)
     for (mi in 1:n.migration.chains) {
          migration.use.weight[mi]=model$log.dens.hyper(use.theta[use.chains[mi],],use.phi[use.chains[mi],pars],prior)
          newChain = mi - 1
          if (newChain == 0) newChain = n.migration.chains
          migration.weight[mi]=model$log.dens.hyper(use.theta[use.chains[mi],],use.phi[use.chains[newChain],pars],prior)
          if(runif(1) < exp(migration.weight[mi]-migration.use.weight[mi])) {
               use.phi[use.chains[mi],pars]=use.phi[use.chains[newChain],pars]
          }
     }
     use.phi
}
