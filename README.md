## Detrended Fluctuation Analysis (DFA)

Javascript port with no external dependencies

Based on the approach described in the [Nodus Labs Fractal Variability Feedback System](https://noduslabs.com/featured/fractal-variability-feedback-system/) article. 

Currently used to analyze the fractal variability of movement in [EightOS](https://8os.io) practice and cognitive variability in [InfraNodus](https://infranodus.com) tool. 

If you have more use cases, please, let us know.

Based on the python script https://github.com/dokato/dfa and the polinomial script https://github.com/rfink/polyfit.js

### Use

Include the dfa.js file into your browser or a node.js app. Then:

```javascript

let time_series = [8, 10, 6, 9, 7, 5, 5, 11, 11, 8, 6, 7, 9, 10, 7, 9]

let dfa = new DFA(time_series)

let alpha_component = dfa.compute()

console.log(alpha_component)

alpha_component = 
  {
    scales: scales,
    fluctuations: fluctuations,
    scales_log: scales_log,
    fluctuations_log: flucts_log,
    coefficients: coefficients,
    alpha: alpha
  }

```


### Concept 

DFA is used to measure the behaviour of a time series. 

An obtained alpha component (closely related to the Hurst exponent) will indicate 
the presence of correlations in the time series. 

It is based on the relationship between the length of an observation and cumulated variability.


### Algorithm:

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


### Author

Created by [Dmitry Paranyushkin](https://deemeetree.com) in 2020


### GPL License

This open source, free software is available under the GNU Affero General Public License version 3 (AGPLv3) license.
You can make modifications to this code and binaries based on it, but only on the condition that you provide access to those modifications under the same license (including remotely  through a computer network).
It is provided as is, with no guarantees and no liabilities.
You can re-use it as long as you keep this notice inside the code.