
<!-- README.md is generated from README.Rmd. Please edit that file -->
powder
------

**pow**er posteriors via **d**ifferential **e**volution in **r**

When presented with several competing formal models, one is required to select between these different explanations: a process commonly known as model selection. One of the more robust methods of performing model selection is through Bayes factors. However, Bayes factors require the computation of the marginal likelihood for each model, which, in turn, requires integrating over the entire parameter space. This makes the estimation of the Bayes factor for models with more than a few dimensions problematic for numerical integration techniques. Early Monte Carlo methods such as arithmetic/harmonic mean estimators or reversible-jump MCMC have been used with some success, but still face problems (Friel & Wyse, 2012). Other solutions exist, such as using GPUs (Evans & Brown, 2017), but are challenging to implement and are done on a model-by-model basis. In this vignette, we present two recent advancements, thermodynamic integration (Friel & Pettitt, 2008; Lartillot & Philippe, 2006) and stepping stone sampling (Xie, Lewis, Fan, Kuo, & Chen, 2011), and show how the `powder` package is used to estimate them. We use the Linear Ballistic Accumulator model (Brown & Heathcote, 2008) as an example.

Installing powder
-----------------

``` r
#install.packages('devtools')
devtools::install_github('jeff324/powder',build_vignettes = TRUE)
```

Vignette
--------

View the tutorial on thermodynamic integration and steppingstone sampling.

``` r
vignette("powder-lba")
```

Obtaining the power posteriors and estimating the marginal likelihood
---------------------------------------------------------------------

Now that we have an overview of the methods, let's use `powder` to define an LBA model, obtain the power posteriors, and compute the estimates for the marginal likelihood. To define an LBA model that varies threshold across conditions we can do the following:

``` r
library(powder)
model = LBA$new(b=T,conds = c(1,2))
```

There are many other LBA models that can be defined, but for now we will stick with a simple one. Next, we load our data. This dataset is called "null" and is a simulated dataset from the LBA with two conditions in which no parameters were varied.

``` r
data('null',package='powder')
```

Next, we will run the `powder` function to obtain the posterior samples:

``` r
pow.out.null = powder(data=null, model = model, num.temps = 30, n.samples = 500, burnin=1000, meltin=250)
```

This runs Differential Evolution MCMC (Turner et al. 2014) over `num.temps` power posteriors, each containing `n.samples` samples after burning in `burnin` samples. It runs each power posterior sequentially, using the last sample of each power posterior as the starting point for the new power posterior. When switching to the next power posterior, it takes time to adapt. We have set the adaptation time to `meltin` iterations. This will take a fairly long time to run.

Lastly, we can estimate the marginal likelihood with:

``` r
ml1 = marginal_likelihood(pow.out.null)
```

This will produce a data frame with the following format (note the values are for illustration only):

``` r
print(ml1)
```

    ##              Method     Value
    ## 1                TI  164.8057
    ## 2      TI Corrected 5242.6978
    ## 3     Harmonic Mean 1568.9184
    ## 4     Steppingstone       Inf
    ## 5 Log Steppingstone 1036.3150

The Steppingstone estimate will usually be `Inf` because of underflow. Use the log steppingstone estimate instead.
