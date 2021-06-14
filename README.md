**Detrended Fluctuation Analysis (DFA)

Javascript port

available on https://github.com/deemeetree


based on the approach described on https://noduslabs.com/featured/fractal-variability-feedback-system/

based on the python script https://github.com/dokato/dfa

*Use
Include the dfa.js file into your browser or a node.js app. Then:

```javascript

let time_series = [8, 10, 6, 9, 7, 5, 5, 11, 11, 8, 6, 7, 9, 10, 7, 9]

let dfa = new DFA(time_series)

let alpha_component = dfa.compute()

console.log(alpha_component)

```


*Concept 

DFA is used to measure the behaviour of a time series. 

An obtained alpha component (closely related to the Hurst exponent) will indicate 
the presence of correlations in the time series. 

It is based on the relationship between the length of an observation and cumulated variability.


*Algorithm:

1) represent a time series as a one-dimensional vector

2) transform it into cumulative sum of variances from the mean

3) generate the scales of different window sizes for the time series — each containing the number of observations

4) split the cumulative variance vector into those chunks

5) for each chunk: generate polinomial in the range in a window (to detrend it)

6) calculate RMS of the fluctuations of the original from the fit

7) get the average RMS value for each window in the scale

8) take the average for each squared RMS for each scale

9) if they align along in a straight line in loglog plot, there is power law relation
i.e. as observations increase in length, the amplitude of average fluctuations increases as well 
that is, there is an exponential growth of fluctuations when there's an exponential growth of scale (2^n)
on the smaller scales the deviations are smaller. on the bigger scales they are much bigger.
