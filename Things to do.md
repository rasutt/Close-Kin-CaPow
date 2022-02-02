### Things to do

-   The expression for the expected numbers of half-sibling pairs within samples is wrong and the error depends on the age of maturity alpha. The one for same-father pairs is now correct, but the one for same-mother pairs is wrong when there are calving intervals, and the error depends on alpha. It would be good to derive the new expression for half-sibling pairs and check it. I should write up the derivation of the one for same-father pairs first though.

-   There are some outputs in the check tab to show the average numbers of kin-pairs that include an animal of ages up to the length of the population history in the final year of the simulation, and compare them to the expected values. They are just tables but it would be good to try plotting them.

-   When the age of maturity is set to zero none of the expressions for the expected numbers of kin-pairs match the simulation, and it would be good to find out why.

-   It would be good to add outputs to show the numbers of parent-offspring pairs within and between survey years in the whole population, show that they match the expressions for the expected numbers, and check the bias in the approximations for their probabilities. I think I changed those expressions as well as the simulation to allow animals to have calves as soon as they are mature, rather than requiring them to be mature when the calves are conceived, but I should check it and rewrite the derivations for them.

-   Add combined model, and more simulation checks.

-   Allow selection of kin-pair types to include in close-kin model. Add half-sibling pairs between surveys.

-   Allow for just one survey for close kin model.

-   Allow results to be downloaded, as both R objects and pdf report.

-   Include mathematical explanations. Can use includeMarkdown() or withMathJax().

-   Investigate when close kin does and does not do better than popan

    -   One interesting observation is that popan depends much more on having surveys over a good range of times, whereas close kin seems to do better when there's only a short range. Mainly for lambda (5-8x better), and phi (2-3x better), but also somewhat for N_final (15-20% better), when three consecutive surveys. Tried with two but strange results, need to actually check numbers converged now. Should check confidence interval coverage too.

-   Incorporate the number of animals for which we don't know the parents into the calculations of the expected numbers of kin pairs. And model-fitting? Warnings when too high?

    -   Check if bias in numbers of kinpairs due to finite sums vs closed form expressions. The finite sums going back exactly the length of the simulation should be unbiased. Using closed form expressions so that could be the reason. Should check how Mark handles the long-term implications of exponential growth. Maybe doesn't matter for short-lived animals?
    -   I think that's right. So I think it'll be equivalent and fast and easy to use the closed form expressions in cases where the constant dynamics apply for the whole period that matters, but we'll need a non-constant dynamics model, with no closed form solution for the kin pairs, in other cases, especially the case of whales. Probably wanna provide both options, and shouldn't be too difficult, since already derived most of the closed form expressions. Might be fast and convenient to do the rest of the testing of them with lower survival rate and age of maturity, and thus shorter necessary simulation lengths? Will wanna establish the expression for those lengths, not only for the simulations, but for the assumptions of the constant dynamics models, which people will wanna use when possible since they'll be simpler and allow greater statistical power. They're the same expression, the length of time the model has to have applied to account for all of the kin-ships in the data with high enough probability. What's high enough will depend on the sensitivity of the model estimates, which we can investigate with the app and maybe also mathematically.
    -   Could it also be because of the approximation of the expected total number of pairs? Check kinpair probabilities correct, as well as numbers. Could we get the analytic expression. Basically a binomial RV squared right? Can check our approximation against the sim, and any better attempt too. I think this is more likely to be the problem because very few animals are getting captured with unknown parents, which would distinguish the closed-form and finite sum expressions.
    -   If we get the distribution we might also be abel to put CI's on the population plots which would be cool.
    -   Check this stuff and then maybe try a proper logistic difference equation.

-   Try plotting relationships between captures?
