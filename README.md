# Close-kin CaPow - A web application for close-kin capture-recapture power analysis

This is a web app for close-kin capture-recapture study power analysis by simulation. It simulates datasets from capture-recapture studies of populations of animals over time, checks whether the simulations match certain predictions, fits close-kin and/or conventional capture-recapture models to the datasets, and analyses the resulting estimates. It was developed to study the New Zealand Southern Right Whale population, but may be useful for other similar populations.

Note - This application is no longer under active development and further significant updates are unlikely, however requests for minor assistance or adjustments may be considered and can be made in the Issues section above.

## Docker branch

This branch of the repo is setup to run in a docker container from the [Rocker project](https://rocker-project.org/images/versioned/shiny.html). It contains the app files in a new folder called app, with a docker file outside it. It works locally with the instructions at <https://anderfernandez.com/en/blog/put-shiny-app-into-production/>, you may just have to remove any compiled files (.o and .dll) from the TMB folder before building it. Unfortunately on Google cloud it ran inconsistently, and shiny apps are not supported there.

## Using the app

The app starts up in the Simulate studies tab, where you can choose the settings for the next simulation to run, see its features based on those settings, and run it. Next is the Check simulation tab, where you can check various features of the simulation, including focusing on the first study in detail, or looking over all the studies together. Then there is the Fit models tab, where you can choose one or more models and fit them to the data from the simulated studies. Finally there is the Save/load tab, where you can save the simulation and analysis, or load a previous one.

Much of the application is self-explanatory given familiarity with the area, and the best way to get to know it is to open it in a web browser and interact with it, but for those with less familiarity some explanation of each part of the app is provided below.

## Technical details

The app is written in R shiny, and uses C++ code with the TMB package for automatic differentiation and fast model-fitting, and multicore processing to speed up the most expensive computations (probabilities of pairs of genotypes given various family relationships). 

The main parts of the app are:

-   global.R - global variables used throughout the app

-   ui.R - user interface and inputs

-   server.R - values to output

-   ckc_saved_objs.Rdata - dataset loaded when app started

-   Tabs - code specific to each tab

-   Functions - more general functions

-   TMB_files - likelihood function and files compiled using TMB package

There is also a Notes folder containing written material explaining some of the maths and theory behind the app.

\
