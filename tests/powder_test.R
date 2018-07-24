rm(list=ls())
library(powder)
model = LBA$new(b=T)
data('null',package='powder')
pow.out.null = powder(data=null,model = model,num.temps = 10,
                        n.samples = 10,burnin=10,meltin=2)
ml = marginal_likelihood(pow.out.null)
