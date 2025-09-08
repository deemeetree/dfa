## Detrended Fluctuation Analysis (DFA)

Javascript port with no external dependencies

Inspired by the approach described in the [Nodus Labs Fractal Variability Feedback System](https://noduslabs.com/featured/fractal-variability-feedback-system/) article.

Currently used to analyze the fractal variability of movement in [EightOS](https://8os.io) practice and cognitive variability in [InfraNodus](https://infranodus.com) tool.

If you have more use cases, please, let us know.

Based on the python script https://github.com/dokato/dfa

### Use

#### When using NPM

```
npm install dfa-variability
```

#### In the browser

Include the dfa.js file into your browser or a node.js app. Then use the Javascript below:

```javascript
// when running on the backend (Node.Js)
let DFA = require("dfa-variability");

// When running on the frontend (via browser) or backend:
let time_series = [8, 10, 6, 9, 7, 5, 5, 11, 11, 8, 6, 7, 9, 10, 7, 9];

let dfa = new DFA(time_series);

// Calculate DFA with default parameters
let alpha_component = dfa.compute();

// Or with custom parameters:
// dfa.compute(minWindow = 4, step = 2, expStep = 0.5, shortMax = 16, longMin = 17, longMaxFraction = 0.25)

console.log(alpha_component);
```

The object structure of the response is:

```javascript
alpha_component = {
    // Basic HRV statistics
    averageVariance: avgVar,      // Population variance of the time series
    meanValue: meanValue,          // Mean of the time series
    lengthOfData: N,              // Length of the time series
    SDNN: SDNN,                   // Standard deviation of NN intervals
    RMSSD: RMSSD,                 // Root mean square of successive differences
    lnRMSSD: lnRMSSD,            // Natural log of RMSSD
    PNN50: PNN50,                 // Percentage of successive differences > 50ms
    averageDifferences: avgDiff,  // Mean absolute successive differences
    
    // DFA results
    scales: scales,               // Array of window sizes used
    segments: segments,           // Array of forward segment counts per scale
    fluctuations: fluctuations,   // F(s) values for each scale
    scalesLog: scalesLog,        // Natural log of scales
    fluctuationsLog: fluctLog,    // Natural log of fluctuations (null for F=0)
    
    // Alpha components
    coefficients: coefficients,   // {slope, intercept} of global fit
    alpha: alpha,                 // Global scaling exponent (DFA slope)
    alpha1: alpha1,              // Short-term scaling (scales 4-16)
    alpha2: alpha2,              // Long-term scaling (scales 17+)
    alpha1Range: [min, max],     // Actual scale range used for α₁
    alpha2Range: [min, max],     // Actual scale range used for α₂
    
    // Scoring
    alphaScore: alphaScore,      // Categorical: "recovering"/"regular"/"resilient"/"tension"
    alphaScoreNumeric: score     // Numeric: 0-100 (distance from 1.0)
}
```

### Concept

DFA is used to measure the behaviour of a time series.

An obtained alpha component (closely related to the Hurst exponent) will indicate
the presence of correlations in the time series.

It is based on the relationship between the length of an observation and cumulated variability.

### Algorithm:

1. **Input Processing**: Represent the time series as a one-dimensional vector (RR intervals for HRV)

2. **Basic HRV Statistics**: Calculate population variance, SDNN, RMSSD, pNN50 for heart rate variability metrics

3. **Integrated Profile**: Transform the series into cumulative sum of deviations from the mean (detrending)

4. **Scale Generation**: 
   - Generate window sizes with α₁ anchors (4, 6, 8, ..., 16) using linear steps
   - For longer series, add geometric progression for α₂ scales (17+) using factor 2^0.5
   - Cap maximum scale at min(64, N/4) to ensure ≥4 forward segments

5. **Segmentation**: For each scale s:
   - Divide the profile into non-overlapping segments of size s (forward segments)
   - If remainder ≥ minWindow, also create backward segments from the end
   - This captures both forward and backward fluctuations

6. **Linear Detrending**: For each segment:
   - Fit a linear trend using least squares regression
   - Calculate residuals (deviations from the trend line)
   - Compute RMS of residuals for that segment

7. **Fluctuation Aggregation**: For each scale:
   - Calculate F(s) = √(mean(RMS²)) across all segments
   - This gives the characteristic fluctuation at that scale

8. **Log-Log Analysis**:
   - Take natural logarithm of scales and fluctuations (filtering F=0)
   - Perform linear regression on log(scales) vs log(fluctuations)
   - The slope gives the scaling exponent α

9. **Multi-Scale Analysis**:
   - **Global α**: Fit across all valid scales
   - **α₁**: Short-term correlations (scales 4-16)
   - **α₂**: Long-term correlations (scales 17+, with ≥4 segments)
   
10. **Interpretation**:
    - α ≈ 0.5: Uncorrelated (white noise)
    - α ≈ 1.0: Scale-invariant, fractal-like (1/f noise)
    - α ≈ 1.5: Brownian motion
    - α < 0.5: Anti-correlated
    - α > 1.0: Strong correlations

### Author

Created by [Dmitry Paranyushkin](https://deemeetree.com) in 2020-2021

### GPL License

This open source, free software is available under the GNU Affero General Public License version 3 (AGPLv3) license.
You can make modifications to this code and binaries based on it, but only on the condition that you provide access to those modifications under the same license (including remotely through a computer network).
It is provided as is, with no guarantees and no liabilities.
You can re-use it as long as you keep this notice inside the code.
