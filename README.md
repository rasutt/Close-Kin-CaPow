Close Kin CaPow
================

This is a little web app in Shiny to illustrate close-kin
capture-recapture parameter estimation via simulation. It simulates
datasets from capture-recapture studies of a population of animals over
time. It takes the simulated population growth rate as input, and shows
the population size over time and the negative log-likelihood for the
first sample, and the MLEs over all samples.

### Things to do

-   Include outputs for the sim validation functions in the Check.R
    script.
-   Incorporate the number of animals for which we don’t know the
    parents into the calculations of the expected numbers of kin pairs.
-   Change proportions calculations to use numbers captured, not total
    numbers alive.
