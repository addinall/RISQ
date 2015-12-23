# RISQ
## Representativity Indicators for Survey Quality
### Mark Addinall - January 2016

-----

Missing data due to non-response impose a serious threat to the quality of surveys and register-based statistics. Non-response is often not a random phenomenon . It usually depends on demographic and socio-economic characteristics of individuals or enterprises. Also the data collection process may have a substantial influence.

The response rate is often used as an indicator for survey quality. It has the advantage that is can be easily computed. However, low response rates will not necessarily cause estimates to be biased. There are ample examples in the literature where increased data collection efforts has lead to a higher response rate but also to a larger non-response bias.

To assess the effects of nonresponse on the quality of statistics, other quality indicators are needed. These indicators should measure the degree to which respondents and non-respondents differ from each other. In other words, such indicators should measure the degree to which the group of respondents in a survey or register resembles the population. The indicators are called Representativity Indicators or, for short, R-indicators.

It is the objective of the RISQ Project to develop R-indicators, to explore their characteristics and to show how to implement and use them in a practical data collection environment.

The project will demonstrate that R-indicators are not only very useful in the analysis of survey data, but also during fieldwork. They can be used to the monitor data collection processes, and therefore facilitate efficient allocation of interviewing resources.

## This Project

It is my intention to use work done by the RISQ team in my exploration of emerging statistical techniques in the field of Epidemiology.  This will be encompassed in more general R/Python tool formation as my code matures.

## Tools from RISQ

An extended version of the RISQ 2 code in SAS and R was released at September 14, 2015. RISQ 2.1 estimates partial CV at the variable-level and at the category level. Furthermore, it approximates standard errors for the overall CV and partial CV. The manual is extended and discusses the new indicators.

## Response rate and non-response bias

The response rate of a survey is frequently used as an indicator of the quality of the survey. However, it is well-known in the survey methodological literature that response rates by themselves can be poor indicators of non-response bias. Up until now no alternative indicators of non-response have been suggested that may be more useful to assess survey quality.

Non-response has two main consequences. First, it reduces the sample size. Therefore. it decreases the precision of the estimates. Second, it affects the sampling design. The probability to obtain an observation depends on both the selection probabilities specified in the sampling design and the unknown probabilities to respond (if selected). Just using selection probabilities in estimation procedures will lead to biased estimates.

The response rate is directly related to the precision. The higher the response rate, the higher the precision. The same does not hold for the bias. The bias does not only depend on the response rate, but also on the extent to which respondents on non-respondents differ. A higher response rate only minimises the maximal impact that non-response has under the worst-case scenario.

It can be shown that (under the Random Response Model) the bias of the response mean is equal to

```
     R(ρ, Y) x S(ρ) x S(Y) / A(ρ)
```

in which R(ρ,Y) is the correlation between response probabilities ρ and the values of the target variable (Y), S(ρ) is the standard error of the response probabilities, S(Y) is the standard error of the target variable (Y) and A(ρ) is the average of the response probabilities.

This implies that low response rates have a negative impact on bias if there is a linear relation between the survey variable and response behaviour. We will use response probabilities to define the concept of representativity.



## The concept of representativity

The concept of representativity is not one without debate. Kruskal and Mosteller (1979) give an extensive overview of the interpretations of this concept in the statistical and non-statistical literature. Here, we relate representativity to the missing-data-mechanism alone. We use the individual response probabilities ρ[k], for k = 1, 2, .., N. We define the response data set to be strongly representative if all response probabilities are identical.

This is an attractive definition, but not useful is practice. We can not test whether our survey data set satisfies this definition, simply because we can not estimate the response probabilities. Therefore, we introduce a weaker definition that relates representativity to the available auxiliary information. A response data set is said to be weakly representative with respect to the sample and an auxiliary variable X if the average response probability is the same in each category of X.

The weak definition of representativity implies that the auxiliary variable and response behaviour are independent. This assumption can be tested with a chi-squared test.

The two definitions given here relate to the concept Missing Completely At Random (MCAR). A missing-data-mechanism is MCAR in case the probability to respond does not depend on the value of the target variable. Hence, in case the response dataset is strongly representative, the missing-data-mechanism is MCAR for any target variable. In case the response subset is weakly representative for X, then the missing-data-mechanism is MCAR for X but not necessarily for other variables.

Response probabilities can be used for the construction of indicators that measure divergence from representativity.

## Measures of representativity

In practice response will hardly ever be fully representative with respect to all available auxiliary socio-demographic characteristics. In order to evaluate the quality of the response there is a need for indicators of how dissimilar the response data set is from the sample data set. One example of a Representativity Indicator (or: R-indicator) is given here.

This R-indicator is based on the standard deviation of estimated response probabilities. It is defined by

```
     M(ρ)= 1 - 2 S(ρ)
```

The response data set is representative if all response probabilities are equal. In this case the standard deviation is zero, and M(p) assumes the value 1.

The response data set is not representative if there is much variation in response probabilities. This is reflected by a large standard error. The maximum value the standard error can assume is 0.5. In this case the value of the R-indicator is equal to 0.

In practice, the values of the response probabilities are not known. Therefore, they are estimated using, for instance, a logistic or probit regression model.
