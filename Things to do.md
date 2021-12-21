### Things to do

-   Extend model analysis to popan models and all parameters.

-   Change to log-normal CI's for population parameters.

-   Include expected final population size in inputs and update selected value so that it doesn't include survey in 2020.

-   Allow for just one survey for close kin model

-   Show parameters simulated from in checking and analysis tabs.

-   Add output for differences between population size estimates and true values, rather than expected values.

-   Also average confidence interval width, since often models fail by returning large variances and CI's that still cover their true values. Also average bias, since biases can cancel out, although I guess that comes into CV. Yeah nah, CV covers all of that I think. You check the CI coverage and bias for problems, and then look at the CV as the actual performance metric, which will still be bad for failed models. Still, easy to report them all so you can see how they fail I guess.

-   Include the parameters controlling the biological scenarios Emma's interested in as inputs and investigate their effects.

-   Show parameters of last simulation to keep clear when the inputs been changed but not simulated from.

-   Include mathematical explanations. Can use includeMarkdown() or withMathJax().

-   Can I make it so when new model requested it doesn't fit old model again too?

-   Investigate when close kin does and does not do better than popan

    -   One interesting observation is that popan depends much more on having surveys over a good range of times, whereas close kin seems to do better when there's only a short range. Mainly for lambda (5-8x better), and phi (2-3x better), but also somewhat for N_final (15-20% better), when three consecutive surveys. Tried with two but strange results, need to actually check numbers converged now. Should check confidence interval coverage too.

-   Incorporate the number of animals for which we don't know the parents into the calculations of the expected numbers of kin pairs. And model-fitting? Warnings when too high?

    -   Check if bias in numbers of kinpairs due to finite sums vs closed form expressions. The finite sums going back exactly the length of the simulation should be unbiased. Using closed form expressions so that could be the reason. Should check how Mark handles the long-term implications of exponential growth. Maybe doesn't matter for short-lived animals?
    -   I think that's right. So I think it'll be equivalent and fast and easy to use the closed form expressions in cases where the constant dynamics apply for the whole period that matters, but we'll need a non-constant dynamics model, with no closed form solution for the kin pairs, in other cases, especially the case of whales. Probably wanna provide both options, and shouldn't be too difficult, since already derived most of the closed form expressions. Might be fast and convenient to do the rest of the testing of them with lower survival rate and age of maturity, and thus shorter necessary simulation lengths? Will wanna establish the expression for those lengths, not only for the simulations, but for the assumptions of the constant dynamics models, which people will wanna use when possible since they'll be simpler and allow greater statistical power. They're the same expression, the length of time the model has to have applied to account for all of the kin-ships in the data with high enough probability. What's high enough will depend on the sensitivity of the model estimates, which we can investigate with the app and maybe also mathematically.
    -   Could it also be because of the approximation of the expected total number of pairs? Check kinpair probabilities correct, as well as numbers. Could we get the analytic expression. Basically a binomial RV squared right? Can check our approximation against the sim, and any better attempt too. I think this is more likely to be the problem because very few animals are getting captured with unknown parents, which would distinguish the closed-form and finite sum expressions.
    -   Check this stuff and then maybe try a proper logistic difference equation.

<!-- -->

-   Try plotting relationships between captures?
