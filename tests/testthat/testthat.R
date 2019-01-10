library(testthat)
library(powder)
context('powder')


test_that('powder test',{
     data('null',package='powder')
     model = LBA$new()
     out = powder(data=null, model = model, num.temps = 1, n.samples = 1, burnin=1, meltin=1,verbose=F)
     expect_equal(length(out),4)
     expect_is(out,'Powder.Hierarchical')
     expect_is(out$log.like.list,'list')
     expect_is(out$theta,'array')
     expect_is(out$phi,'array')
     expect_is(out$options,'list')

     expect_equal(length(out$log.like.list),1)
     opts = out$options
     theta_dim = c(opts$n.chains,
                   opts$n.pars,opts$n.subj,
                   opts$meltin*length(opts$temperatures) +
                        opts$burnin + opts$n.samples*length(opts$temperatures))
     expect_equal(dim(out$theta),theta_dim)

     phi_dim = c(opts$n.chains,opts$n.hpars,opts$meltin*length(opts$temperatures) + opts$burnin +
            opts$n.samples*length(opts$temperatures))
     expect_equal(dim(out$phi),phi_dim)

     out = powder(data=null, model = model, num.temps = 2, n.samples = 1, burnin=1, meltin=1,
                  n.sequences = 2, current.sequence = 1, verbose=F)
     expect_equal(length(out$log.like.list),1)
     expect_equal(out$options$temperatures,0)
     out = powder(data=null, model = model, num.temps = 2, n.samples = 1, burnin=1, meltin=1,
                  n.sequences = 2, current.sequence = 1,high.temps.first = TRUE, verbose=F)
     expect_equal(length(out$log.like.list),1)
     expect_equal(out$options$temperatures,1)

})

test_that('powder individual', {
     data('null',package='powder')
     dat = null
     model = LBA.Individual$new()
     out = powder(model=model, data=dat[[1]], num.temps=1, n.samples=1, burnin=1, meltin=1, verbose=F)
     expect_is(out,'Powder.Individual')
     expect_equal(length(out),3)
     expect_is(out$log.like.list,'list')
     expect_is(out$theta,'array')
     expect_is(out$options,'list')

     out = powder(model=model, data=dat[[1]], num.temps=2, n.samples=2, burnin=2, meltin=2, verbose=F)
     expect_equal(length(out),3)
     expect_is(out$log.like.list,'list')
     expect_is(out$theta,'array')
     expect_is(out$options,'list')
})


test_that('LBA class expectations',{
     model = LBA$new()
     expect_is(model$log.dens.hyper,'function')
     expect_is(model$log.dens.like,'function')
     expect_is(model$log.dens.prior,'function')
     expect_is(model$phi.init,'function')
     expect_is(model$phi.names,'character')
     expect_is(model$phi.start.points,'list')
     expect_is(model$prior,'list')
     expect_is(model$theta.init,'function')
     expect_is(model$theta.start.points,'numeric')
     expect_is(model$vary.parameter,'logical')
})

test_that('LBA parameter vector lengths',{
     model = LBA$new()
     expect_equal(length(model$phi.names),12)
     expect_equal(length(model$phi.start.points),6)
     expect_equal(length(model$theta.names),6)
     expect_equal(length(model$theta.start.points),6)
     expect_equal(names(model$phi.start.points),model$theta.names)
     model = LBA$new(b=T)
     expect_equal(length(model$phi.names),14)
     expect_equal(length(model$phi.start.points),6)
     expect_equal(length(model$theta.names),7)
     expect_equal(length(model$theta.start.points),7)
     model = LBA$new(A=T,b=T)
     expect_equal(length(model$phi.names),16)
     expect_equal(length(model$phi.start.points),6)
     expect_equal(length(model$theta.names),8)
     expect_equal(length(model$theta.start.points),8)
})

test_that('LBA name collision',{
     model = LBA$new()
     which.phi = lapply(model$theta.names,function(p)grep(paste0('^',p),model$phi.names))
     which.phi.len = sapply(which.phi,function(x)length(x))
     expect_false(any(which.phi.len != 2))
     model = LBA$new(b=T,sve=T,ve=T,vc=T,A=T)
     which.phi = lapply(model$theta.names,function(p)grep(paste0('^',p),model$phi.names))
     which.phi.len = sapply(which.phi,function(x)length(x))
     expect_false(any(which.phi.len != 2))
})

test_that('marginal likelihood', {
     model = LBA$new()
     data('null')
     out = powder(model=model,data=null,num.temps=3,burnin=5,meltin=5,n.samples=5,verbose=F)
     ml = summary(out)
     data('null',package='powder')
     model = LBA$new()
     out = powder(data=null, model = model, num.temps = 1, n.samples = 1, burnin=1, meltin=1,verbose=F)

     ml = summary(out)
     expect_is(ml,'data.frame')

     model = LBA.Individual$new()
     data('individual')
     out = powder(model=model,data=individual,num.temps=3,burnin=5,meltin=5,n.samples=5,verbose=F)
     ml = summary(out)
     expect_is(ml,'data.frame')
})

test_that('parallel',{
     model = LBA$new()
     data('null')
     out = powder(model=model,data=null,n.samples=5,burnin=1,verbose=F,method = 'parallel')
     ml = summary(out)
     expect_is(ml,'data.frame')
})
