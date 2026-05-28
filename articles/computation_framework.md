# Computation framework

## JAGS vs. NIMBLE vs. Stan

We use both JAGS (Plummer 2003) and NIMBLE (de Valpine et al. 2017).
NIMBLE uses a derivative of the [JAGS/BUGS
language](https://r-nimble.org/manual/cha-writing-models.html#sec:supp-feat-bugs),
so that the models can be used interchangeably with only minor
modifications. Still, JAGS and NIMBLE are conceptually quite different.
While it is theoretically possible to extend the functionality of JAGS
by custom functions, this process is not straightforward and requires
advanced programming skills. Therefore, the available functions and
samplers have to be taken as they are. In contrast, NIMBLE offers a
great suite of customization possibilities. This is especially true for
the samplers used even for individual nodes. We tried to tune the model
parameters, but with only limited success. We found that the default
samplers of NIMBLE exhibit poorer mixing than those of JAGS and thus
need longer chains to reach the same level in terms of quality measures.
Overall, we observed that the same model with the same data takes longer
to run in NIMBLE than in JAGS. Additionally, NIMBLE is less stable than
JAGS, especially in the multicore setup.

For the reasonse given above, we propose to use JAGS for the univariate
cases and cases with multiple normal likelihood, while for multivariate
applications with multinormal likelihood function, NIMBLE is currently
the only viable option.

Despite the fact that the models of JAGS and NIMBLE are almost
identical, the differences are subtle enough so that the same data and
parameters will produce slightly different results. This is because JAGS
and NIMBLE use different samplers and the way the ordinal thresholds are
ordered works differently.

We interface JAGS with a wrapper function which in turn relies on
`runjags` (Denwood 2016) to run JAGS. NIMBLE is a package in its own
right (de Valpine et al. 2017). The processing of the chains is done
with the `coda` package (Plummer et al. 2006). For further analysis,
`baytaAAR` provides custom functions but in the vignettes we also show
how to analyze the data with the packages `bayesplot` (Gabry et al.
2019) and `tidybayes` (Kay 2024).

We do not use [`Stan`](https://mc-stan.org) which is currently probably
the most popular choice when it comes to working with Bayesian models in
R. There are two main reasons for this: Firstly, `Stan` is much more
restrictive with respect to unobserved nodes. While JAGS and NIMBLE
simply treat those as parameters to be estimated, they have to be
explicitly declared in `Stan`. Secondly, and more importantly, `Stan`
does not provide the interval-distribution (`dinterval`) which is
instrumental to our application because only with the
`dinterval`-distribution is it possible to implement the multivariate
model.

## Parallel processing and number of chains

On modern hardware with usually more than one core, parallel processing
is very useful as it cuts down time significantly. In most cases, you
would use as many cores as chains. The default for JAGS and NIMBLE is
three chains. More chains can help in reaching quality criteria like
`PSRF` faster. However, it also increases the size of the MCMC samples
so that post-processing can become limited by available memory.
Therefore, and following J. Kruschke (2015), we control the number of
saved iterations by dividing `numSavedSteps` by the number of chains.
So, increasing the number of chains will automatically shorten the
number of saved iterations which will usually lead to poorer mixing. In
this respect, we observed that a higher number of chains does not
outweigh the disadvantage of shorter chains. We therefore do not
recommend increasing the number of chains.

Parallel processing is straightforward in
[`bay.ta()`](https://isaakiel.github.io/baytaAAR/reference/bay.ta.md).
You simply specify the mode of operation with the parameter `multicore`
which can be `FALSE` (single-core-processing) or `TRUE` (parallel
processing). This runs out-of-the-box with JAGS but parallel processing
in NIMBLE is much more complicated. The [NIMBLE
manual](https://r-nimble.org/examples/parallelizing_NIMBLE.html) gives a
general account of the necessary steps. For a package like ours, there
is the additional difficulty that we provide custom functions which have
to be supplied in the environment of each cluster. However, simply
copying the functions to the global environment would be unacceptable
behaviour of a package. We use the R-package `parallel`, which is part
of base R, for multi-core processing in NIMBLE.

In the end, we hope that we succeeded to make parallel processing with
NIMBLE as unobtrusive as possible by providing a wrapper function which
takes care of the different modes of operation. Please note, though,
that you will get no feedback while your model is running in parallel.
Therefore, it is good practice to start with a low number of steps and
see how long these take, ‘guesstimate’ the required number of steps from
the quality measures (especially `ESS`) and change the settings
accordingly.

Generally, JAGS is very stable in parallel processing so apart from your
fans possibly spinning up, you will notice very little difference
between single- and multicore-processing. For NIMBLE, this is
unfortunately not always the case. Apart from the lack of feedback,
there is also the possibility of a memory leak. This makes the cores
consume more and more memory until the machine becomes unusable and
might eventually crash. Therefore, we would strongly advise to monitor
the cores closely to see if their memory consumption remains stable or
steadily increases. If the latter is the case, the model should be
stopped and the R session should be restarted. Furthermore, it might be
necessary to terminate the processes manually using your system monitor.

------------------------------------------------------------------------

## References

de Valpine, Perry, Daniel Turek, Christopher Paciorek, Cliff
Anderson-Bergman, Duncan Temple Lang, and Ras Bodik. 2017. “Programming
with Models: Writing Statistical Algorithms for General Model Structures
with NIMBLE.” *Journal of Computational and Graphical Statistics* 26:
403–13. <https://doi.org/10.1080/10618600.2016.1172487>.

Denwood, Matthew J. 2016. “runjags: An R Package Providing Interface
Utilities, Model Templates, Parallel Computing Methods and Additional
Distributions for MCMC Models in JAGS.” *Journal of Statistical
Software* 71 (9): 1–25. <https://doi.org/10.18637/jss.v071.i09>.

Gabry, Jonah, Daniel Simpson, Aki Vehtari, Michael Betancourt, and
Andrew Gelman. 2019. “Visualization in Bayesian Workflow.” *Journal of
the Royal Statistical Society A* 182: 389–402.
<https://doi.org/10.1111/rssa.12378>.

Kay, Matthew. 2024. *tidybayes: Tidy Data and Geoms for Bayesian
Models*. <https://doi.org/10.5281/zenodo.1308151>.

Kruschke, John K. 2015. *Doing Bayesian data analysis: a tutorial with
R, JAGS, and Stan*. Academic Press.

Plummer, Martyn. 2003. “JAGS: A Program for Analysis of Bayesian
Graphical Models Using Gibbs Sampling.” In *Proceedings of the 3rd
International Workshop on Distributed Statistical Computing (DSC 2003),
Vienna, 20-22 March 2003*, edited by Kurt Hornik, Friedrich Leisch, and
Achim Zeileis. Technische Universität Wien.

Plummer, Martyn, Nicky Best, Kate Cowles, and Karen Vines. 2006. “CODA:
Convergence Diagnosis and Output Analysis for MCMC.” *R News* 6 (1):
7–11.
