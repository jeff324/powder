standard.sampling.individual = function(model,data,theta.names,n.pars,temperatures,de_params,
                                        burnin,meltin,n.samples,n.chains,method,message,
                                        return.samples,verbose,update){

     theta = array(NA,c(n.chains, n.pars, meltin*length(temperatures) + burnin + n.samples*length(temperatures)))
     weight = array(-Inf,c(meltin*length(temperatures) + burnin + n.samples*length(temperatures),n.chains))
     prior = model$prior
     log.like.list = list()
     colnames(theta) = theta.names
     b = de_params$b

     idx = 2
     for (t in 1:length(temperatures)) {
          if(verbose){
               cat('\n')
               if(method=='posterior'){
                    cat(paste(message))
               }else{
                    cat(paste(message, format(temperatures[t],scientific=T,digits=4),':',t,'/', length(temperatures)))
               }
          }
          if (t==1) {
               n.iter = burnin + n.samples
               for (i in 1:n.chains) {
                    while (weight[1, i]==-Inf) {
                         theta[i, ,1] = model$theta.init()
                         weight[1, i] = model$log.dens.like(theta[i, ,1],data=data,par.names=theta.names)
                    }
               }
          }else{
               n.iter = meltin + n.samples
               log.like[1, ] = sapply(1:n.chains, function(x) model$log.dens.like(theta[x, ,idx-1], data=data, par.names=theta.names))
               weight[idx-1, ] = log.like[1, ]
          }

          log.like = array(NA,c(n.iter, n.chains))
          for (i in 2:n.iter) {
               if(verbose){
                    if (i%%update==0)cat("\n ",i,'/',n.iter)
               }
               temp = t(sapply(1:n.chains, crossover_ind, n.chains=n.chains, b=b, pars=1:n.pars, use.theta=theta[ , ,idx-1], use.like=weight[idx-1, ],
                               data=data, hyper=prior, par.names=theta.names, temperature=temperatures[t], model=model))
               weight[idx, ] = temp[ ,1]
               theta[ , ,idx] = temp[ ,2:(n.pars+1)]
               log.like[i, ] = temp[ ,1]
               idx = idx + 1
          }
          log.like.list[[t]] = log.like
     }
     return(list(log.like.list=log.like.list,theta=theta))
}

parallel.sampling.individual = function(model,data,theta.names,n.pars,temperatures,de_params,
                                        burnin,n.samples,n.chains,method,message,
                                        return.samples,verbose,update){

     theta = array(NA,c(n.chains, n.pars, burnin + n.samples))
     weight = array(-Inf,c(burnin + n.samples,n.chains))
     prior = model$prior
     log.like.list = list()
     colnames(theta) = theta.names
     b = de_params$b

     idx = 2

     if(verbose){
          cat('\n')
          cat(message)
     }

     n.iter = burnin + n.samples
     for (i in 1:n.chains) {
          while (weight[1, i]==-Inf) {
               theta[i, ,1] = model$theta.init()
               weight[1, i] = model$log.dens.like(theta[i, ,1],data=data,par.names=theta.names)
          }
     }

     for (i in 2:n.iter) {
          if(verbose){
               if (i%%update==0)cat("\n ",i,'/',n.iter)
          }
          temp = t(sapply(1:n.chains, crossover_ind_parallel, n.chains=n.chains, b=b, pars=1:n.pars, use.theta=theta[ , ,idx-1], use.like=weight[idx-1, ],
                          data=data, hyper=prior, par.names=theta.names, temperatures=temperatures, model=model))
          weight[idx, ] = temp[ ,1]
          theta[ , ,idx] = temp[ ,2:(n.pars+1)]
          idx = idx + 1
     }

     return(list(log.like.list=weight,theta=theta))
}

standard.sampling.hierarchical = function(model,data,theta.names,n.pars,temperatures,de_params,
                                          burnin,meltin,n.samples,n.chains,method,message,
                                          return.samples,verbose,update){

     prior = model$prior
     n.subj = length(data)
     phi.names = model$phi.names
     n.hpars = length(phi.names)

     theta = array(NA, c(n.chains, n.pars, n.subj, meltin*length(temperatures) + burnin + n.samples*length(temperatures)))
     phi = array(NA, c(n.chains, n.hpars, meltin*length(temperatures) + burnin + n.samples*length(temperatures)))
     weight = array(-Inf, c(meltin*length(temperatures) + burnin + n.samples*length(temperatures),n.chains,n.subj))
     log.like.list = list()

     colnames(theta) = theta.names
     colnames(phi) = phi.names

     b = de_params$b
     migration.freq = de_params$migration.freq
     migration.start = de_params$migration.start
     migration.end = de_params$migration.end

     for (i in 1:n.chains) {
          phi[i, ,1] = model$phi.init()
     }

     idx = 2
     for (t in 1:length(temperatures)) {
          if(verbose){
               cat('\n')
               if(method=='posterior'){
                    cat(paste(message))
               }else{
                    cat(paste(message, format(temperatures[t],scientific=T,digits=4),':',t,'/', length(temperatures)))
               }
          }
          if (t==1) {
               n.iter = burnin + n.samples
          } else {
               n.iter = meltin + n.samples
          }
          log_like = array(NA,dim=c(n.iter, n.chains, n.subj))
          if (t==1) {
               for (i in 1:n.chains) {
                    for (j in 1:n.subj) {
                         while (weight[1,i,j]==-Inf) {
                              theta[i, ,j,1] = model$theta.init()
                              weight[1,i,j] = model$log.dens.like(theta[i, ,j,1],data=data[[j]],par.names=theta.names)
                         }
                    }
               }
          } else {
               for (j in 1:n.subj) {
                    log_like[1, ,j] = sapply(1:n.chains,function(x)model$log.dens.like(theta[x, ,j,idx-1],data=data[[j]],par.names=theta.names))
                    weight[idx-1, ,j] = log_like[1, ,j]
               }
          }

          for (i in 2:n.iter) {
               if(verbose){
                    if (i%%update == 0)cat("\n ", i, '/', n.iter, ' ')
               }
               phi[,,idx] = phi[,,idx-1]
               rand.samp = sample(1:n.chains,n.chains)
               for (p in theta.names) {
                    which.theta = match(x=p, table=theta.names)
                    which.phi = grep(paste0('^',p),phi.names)
                    if (idx %% migration.freq == 0 & idx > migration.start & idx < migration.end) {
                         phi[ , ,idx] = migration.crossover_hyper(pars=which.phi, n.chains=n.chains, use.theta=theta[rand.samp,which.theta,,idx-1],
                                                                  use.phi=phi[ , ,idx], prior=prior[[p]], model)
                    } else {
                         phi[ , ,idx] = t(sapply(1:n.chains, crossover_hyper, pars=which.phi, n.chains=n.chains, b=b,
                                                 use.theta=theta[rand.samp, which.theta, ,idx-1], use.phi=phi[ , ,idx], prior=prior[[p]], model))
                    }
               }
               rand.samp = sample(1:n.chains, n.chains)
               hyper = phi[rand.samp, , idx]

               for (j in 1:n.subj) {
                    #cat(j)
                    if (idx %% migration.freq == 0 & idx > migration.start & idx < migration.end) {
                         temp = migration.crossover(pars=1:n.pars, n.chains=n.chains, use.theta=theta[,,j,idx-1], use.like=weight[idx-1,,j],
                                                    data=data[[j]],hyper=hyper, par.names=theta.names, model)
                    } else {
                         temp = t(sapply(1:n.chains, crossover, pars=1:n.pars, n.chains=n.chains, b=b, use.theta=theta[ , ,j,idx-1],
                                         use.like=weight[idx-1, ,j], data=data[[j]], hyper=hyper,par.names=theta.names,temperatures[t],model))
                    }
                    weight[idx, ,j] = temp[ ,1]
                    theta[ , ,j,idx] = temp[ ,2:(n.pars+1)]
                    log_like[i, ,j] = temp[ ,1]
               }
               idx = idx + 1
          }
          log.like.list[[t]] = log_like
     }

     return(list(log.like.list=log.like.list,theta=theta,phi=phi))
}


parallel.sampling.hierarchical = function(model,data,theta.names,n.pars,temperatures,de_params,
                                          burnin,n.samples,n.chains,method,message,
                                          return.samples,verbose,update){

     prior = model$prior
     n.subj = length(data)
     phi.names = model$phi.names
     n.hpars = length(phi.names)

     theta = array(NA, c(n.chains, n.pars, n.subj, burnin + n.samples))
     phi = array(NA, c(n.chains, n.hpars, burnin + n.samples))
     weight = array(-Inf, c(burnin + n.samples,n.chains,n.subj))
     log.like.list = list()

     colnames(theta) = theta.names
     colnames(phi) = phi.names

     b = de_params$b
     migration.freq = de_params$migration.freq
     migration.start = de_params$migration.start
     migration.end = de_params$migration.end

     for (i in 1:n.chains) {
          phi[i, ,1] = model$phi.init()
     }

     idx = 2

     if(verbose){
          cat('\n')
          cat(paste(message))
     }

     n.iter = burnin + n.samples
     log_like = array(NA,dim=c(n.iter, n.chains, n.subj))
     for (i in 1:n.chains) {
          for (j in 1:n.subj) {
               while (weight[1,i,j]==-Inf) {
                    theta[i, ,j,1] = model$theta.init()
                    weight[1,i,j] = model$log.dens.like(theta[i, ,j,1],data=data[[j]],par.names=theta.names)
               }
          }
     }

     for (i in 2:n.iter) {
          if(verbose){
               if (i%%update == 0)cat("\n ", i, '/', n.iter, ' ')
          }
          phi[,,idx] = phi[,,idx-1]
          rand.samp = sample(1:n.chains,n.chains)
          for (p in theta.names) {
               which.theta = match(x=p, table=theta.names)
               which.phi = grep(paste0('^',p),phi.names)
               phi[ , ,idx] = t(sapply(1:n.chains, crossover_hyper, pars=which.phi, n.chains=n.chains, b=b,
                                       use.theta=theta[rand.samp, which.theta, ,idx-1], use.phi=phi[ , ,idx], prior=prior[[p]], model))
          }
          rand.samp = sample(1:n.chains, n.chains)
          hyper = phi[rand.samp, , idx]

          for (j in 1:n.subj) {
               temp = t(sapply(1:n.chains, crossover_parallel, pars=1:n.pars, n.chains=n.chains, b=b, use.theta=theta[ , ,j,idx-1],
                               use.like=weight[idx-1, ,j], data=data[[j]], hyper=hyper,par.names=theta.names,temperatures,model))
               weight[idx, ,j] = temp[ ,1]
               theta[ , ,j,idx] = temp[ ,2:(n.pars+1)]
          }
          idx = idx + 1
     }

     return(list(log.like.list=weight,theta=theta,phi=phi))
}
