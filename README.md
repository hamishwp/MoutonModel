# Dissertation-Appendix
Supplementary materials for the dissertation for MSc Statistical Science, University of Oxford

The text can be found at https://www.overleaf.com/read/svskqkmrqpnd

Abstract:

We explore the Integral Projection Model (IPM), a common tool in statistical ecology for use with census data. IPMs attempt to draw inference on the demography of trait-structured populations by fitting functions which describe the relationship between the size (mass, length, or other) of an individual and certain vital rates (such as its probability of survival, reproduction or growth to a given size). Traditionally, these functions are fitted piecemeal using the same data, by maximum likelihood estimation (MLE). This approach is unsatisfactory, as it is unable to elucidate parameter correlation structures, which is highly important for effective inference in the ecological sphere, where many complicated biological behaviours are interdependent. We propose a method making use of Particle Markov Chain Monte Carlo to produce Bayesian analyses of data sets whose structure is traditionally associated with IPMs. We favour this approach due to Gibbs samplers' weakness to mix well when parameters are highly correlated, because parameters are only updated in small sets at a time \citep{Finke}. Further to this, the use of particle filtering in Integrated Population Models is a relatively recent addition to the toolbox of ecological statisticians, with a recent pioneering paper by \citet{Finke}, which inspires the model specification we use. We apply our methodology to the analysis of both simulated and real Soay Sheep (\textit{Ovis aries}) data from Scotland. We conclude that the methodology shows promise in obtaining estimates of vital rate parameters with smaller, easier to collect data sets, allowing a robust exploration of the posterior correlation, and that the use of prior information is extremely important, as certain parameters become unidentifiable in the new parameterisation.
