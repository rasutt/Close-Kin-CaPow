---
editor_options: 
  markdown: 
    wrap: 72
---

# Close-kin CaPow - A web application for close-kin capture-recapture power analysis

This is a web app for close-kin capture-recapture study power analysis
by simulation. It simulates datasets from capture-recapture studies of
populations of animals over time, checks whether the simulations match
certain predictions, fits close-kin and/or conventional
capture-recapture models to the datasets, and analyses the resulting
estimates. It was developed to study the New Zealand Southern Right
Whale population, but may be useful for other similar populations.

Note - This application is no longer under active development and
further significant updates are unlikely, however requests for minor
assistance or adjustments may be considered and can be made in the
Issues section above.

## Using the app

The app starts up in the Simulate studies tab, where you can choose the
settings for the next simulation to run, see its features based on those
settings, and run it. Next is the Check simulation tab, where you can
check various features of the simulation, including focusing on the
first study in detail, or looking over all the studies together. Then
there is the Fit models tab, where you can choose one or more models and
fit them to the data from the simulated studies. Finally there is the
Save/load tab, where you can save the simulation and analysis, or load a
previous one.

Much of the application is self-explanatory given familiarity with the
area, and the best way to get to know it is to open it in a web browser
and interact with it, but for those with less familiarity some
explanation of each part of the app is provided below.

### Simulate studies

This tab contains the following inputs:

-   Survival rate - The probability that an animal survives from one
    year to the next

-   Per capita birth rate - The number of animals expected to be born
    each year as a proportion of the population size

-   Expected population size in base year

-   Base year for expected population size

-   Survey years - The years in which there are surveys, separated by
    commas, with eg. 2000:2005 indicating surveys in all years from 2000
    to 2005. App outputs may show errors while inputs are being changed,
    until inputs are correct.

-   Base level capture probability - The base-probability that an animal
    is captured in each survey, see additional capture probability below

-   Number of SNP loci - The number of independent binary loci in each
    genotype

-   Females breed in order of time since last breeding - Whether mature
    females are selected to breed randomly, or those that have waited
    longest breed first

-   Additional capture probability when calving - The extra probability
    that an animal is captured in each survey if it is a female with a
    newborn calf

-   Probability that males absent - The probability that a male animal
    is absent from the population each year \*Not accounted for in
    breeding\*

-   Age of sexual maturity - Minimum age at which an animal may have a
    calf

-   Length of population histories - Number of years for which the
    population should be simulated, ending with the final survey

-   Number of studies to simulate

It also contains the following outputs:

-   Expected population size over time - Plot over length of population
    histories, with base year and survey years indicated

-   Implied parameters

    -   Population growth rate - Ratio of expected population sizes of
        consecutive years
    -   Birthrate among mature females - The number of animals expected
        to be born each year as a proportion of the number of mature
        females in the population. If this number is greater than one
        the simulation will not be able to run because female Southern
        right whales are not able to have more than one calf per year.
    -   Initial population size - In the first year of the simulation
    -   Expected final population size - In the last year of the
        simulation (the final survey year)
    -   Expected superpopulation size - The expected number of animals
        that are alive in the population in at least one survey year

-   Predicted numbers of samples and kinpairs - Either expected numbers
    or an approximation of them

    -   Among sampled individuals - Including all pairs among samples
        from all surveys

    -   Among offset pairs - Including a randomly selected minimal set
        of pairs joining every sample within a survey, and every sample
        in each survey with one in every other survey

    -   In population - Including all pairs among all animals in the
        population in each survey year

        -   Number of samples - The expected number of animals sampled
            in each survey

        -   Total number of pairs - The predicted total number of pairs
            of animals sampled in each survey

        -   Self-pairs (all) - The predicted number of pairs of animals
            that are the same animal sampled in two different surveys

        -   Parent-offspring pairs - One animal the parent of the other

        -   Half-sibling pairs - Two animals sharing exactly one parent

        -   Population sizes

        -   All pairs - The predicted total number of pairs of animals
            in the population

        -   Same-mother pairs - Pairs sharing the same mother

        -   Same-father pairs - Pairs sharing the same father

        -   Full-sibling pairs - Same mother and same father

Select the input values corresponding to your population and study
design of interest. Outputs update whenever inputs are changed. Check
the outputs, and when you are happy with the settings click the Simulate
studies button to run the simulation.

### Check simulation

This tab contains the following sub-tabs:

-   Simulation features - Inputs selected, and parameters implied, for
    the last simulation that was run.

    -   Datasets removed and retained - Numbers of datasets removed for
        each reason, and number retained for analysis

        -   Number of extinct populations - Number of populations in
            which every animal died before the end of the simulation

        -   Number of studies with \< 2 samples - Fewer than two samples
            over all surveys, so that there are no pairs of sampled
            animals

        -   Number of datasets retained

-   First study - Analyses of first study retained from simulation

    -   First study samples - Example sample histories and sufficient
        statistics

        -   First sample-histories - Animal ID, mother and father IDs,
            sample history, and calving history, for first few animals
            sampled, by birthyear

        -   Sufficient statistics - Summary statistics of the
            sample-histories that are sufficient to fit the
            corresponding (Popan) model

    -   First study kin-pairs - Numbers of kin-pairs simulated in the
        first population and sampling study simulated

    -   First study genetics  - Genetic analysis of samples from first
        population and study simulated

        -   First sample genotypes - The maternally and paternally
            inherited binary genes at each locus for each of the first
            few animals sampled

        -   Allele frequencies - Relative frequencies, excluding samples
            from same animals in different surveys

        -   Genotype probabilities - Of possible genotypes at each locus

        -   Genopair probabilities at first locus - Probabilities of all
            possible pairs of genotypes given various family
            relationships between the animals

        -   Genopair log-probabilities at first locus

        -   First sample genopair log-probabilities

        -   All sample genopair log-probabilities - Log-probabilities
            corresponding to zero-probability genopairs are removed

        -   First sample genopair probabilities - For all pairs of
            samples, and just offset pairs

            -   Shown along with number of pairs from each pair of
                surveys

        -   Half-sibling vs unrelated pair PLODs - Pseudo log-likelihood
            ratios for each pair of samples being from half-siblings, or
            unrelated animals, given their observed genotypes

            -   Expected values are indicated for animals with various
                family relationships

    -   First study likelihoods - Negative log-likelihood surface over
        each parameter while others held at true values

        -   Shown for each model

    -   First study estimates

    -   

\

## Technical details

The app is written in R shiny, and uses C++ code with the TMB package
for automatic differentiation and fast model-fitting, and multicore
processing to speed up the most expensive computations (probabilities of
pairs of genotypes given various family relationships). 

The main parts of the app are:

-   global.R - global variables used throughout the app

-   ui.R - user interface and inputs

-   server.R - values to output

-   ckc_saved_objs.Rdata - dataset loaded when app started

-   Tabs - code specific to each tab

-   Functions - more general functions

-   TMB_files - likelihood function and files compiled using TMB package

There is also a Notes folder containing written material explaining some
of the maths and theory behind the app.

\
