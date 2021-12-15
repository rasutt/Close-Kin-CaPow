### Things to do

-   Include expected super/final population size in inputs and update selected value so that it doesn't include survey in 2020.

-   Show parameters of last simulation to keep clear when the inputs been changed but not simulated from.

-   Split into separate tabs for simulation planning, checking datasets, and model fitting and power analysis.

-   Incorporate the number of animals for which we don't know the parents into the calculations of the expected numbers of kin pairs. And model-fitting? Warnings when too high?

    -   Check if bias in numbers of kinpairs due to finite sums vs closed form expressions. The finite sums going back exactly the length of the simulation should be unbiased. Using closed form expressions so that could be the reason. Should check how Mark handles the long-term implications of exponential growth. Maybe doesn't matter for short-lived animals?
    -   Could it also be because of the approximation of the expected total number of pairs? Check kinpair probabilities correct, as well as numbers. Could we get the analytic expression. Basically a binomial RV squared right?
    -   Check this stuff and then maybe try a proper logistic difference equation.

-   Find confidence interval coverage and show number of results discarded.

-   Try plotting relationships between captures?

-   Investigate when close kin does and does not do better than popan

    -   One interesting observation is that popan depends much more on having surveys over a good range of times, whereas close kin seems to do better when there's only a short range. Mainly for lambda (5-8x better), and phi (2-3x better), but also somewhat for N_final (15-20% better), when three consecutive surveys. Tried with two but strange results, need to actually check numbers converged now. Should check confidence interval coverage too.
