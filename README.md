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
// dfa.compute(minWindow, expStep, step, shortMax, longMin, longMaxFraction, level)

// where:
// minWindow is the minimum window size for scales (default: 4)
// expStep is the step for increasing scales for global alpha (default: 0.25)
// step is the step for increasing scales for short alpha1 and long alpha2 (default: 2)
// shortMax is the maximum window size for short alpha1 (default: 16)
// longMin is the minimum window size for long alpha2 (default: 16)
// longMaxFraction is the maximum fraction of the series length for long alpha2 (default: 0.25)
// level is the threshold sensitivity: "relaxed", "moderate", or "strict" (default: "moderate")

console.log(alpha_component);
```

The object structure of the response is:

```javascript
alpha_component = {
	// Basic HRV statistics
	averageVariance: avgVar, // Population variance of the time series
	meanValue: meanValue, // Mean of the time series
	lengthOfData: N, // Length of the time series
	SDNN: SDNN, // Standard deviation of NN intervals
	RMSSD: RMSSD, // Root mean square of successive differences
	lnRMSSD: lnRMSSD, // Natural log of RMSSD
	PNN50: PNN50, // Percentage of successive differences > 50ms
	averageDifferences: avgDiff, // Mean absolute successive differences

	// DFA results
	scales: scales, // Array of window sizes used
	segments: segments, // Array of forward segment counts per scale
	fluctuations: fluctuations, // F(s) values for each scale
	scalesLog: scalesLog, // Natural log of scales
	fluctuationsLog: fluctLog, // Natural log of fluctuations (null for F=0)

	// Alpha components
	coefficients: coefficients, // {slope, intercept} of global fit
	alpha: alpha, // Global scaling exponent (DFA slope)
	alpha1: alpha1, // Short-term scaling (scales 4–16), null if insufficient data
	alpha2: alpha2, // Long-term scaling (scales 16–N/4), null if insufficient data
	alpha1Range: [min, max], // Actual scale range used for α₁ (null when alpha1 is null)
	alpha2Range: [min, max], // Actual scale range used for α₂ (null when alpha2 is null)
	scalesAlpha1: [...], // Window sizes included in the α₁ fit
	scalesAlpha2: [...], // Window sizes included in the α₂ fit

	// Labels — simplified (4 categories)
	alphaLabel: "fractal", // "random" / "regular" / "fractal" / "complex"
	alpha1Label: "fractal", // Same for α₁ (null when alpha1 is null)
	alpha2Label: "regular", // Same for α₂ (null when alpha2 is null)

	// Labels — traditional DFA (6 categories)
	dfaLabel: "1/f noise", // Scientific DFA classification for global α
	dfa1Label: "1/f noise", // Same for α₁ (null when alpha1 is null)
	dfa2Label: "correlated", // Same for α₂ (null when alpha2 is null)

	// Labels — legacy HRV (same thresholds as alphaLabel, different names)
	alphaScore: "resilient", // "recovering" / "regular" / "resilient" / "tension"
	alpha1Score: "resilient", // Same for α₁ (null when alpha1 is null)
	alpha2Score: "regular", // Same for α₂ (null when alpha2 is null)

	// Scoring — numeric
	alphaScoreNumeric: 95, // 0–100 for global α (distance from 1.0)
	alpha1ScoreNumeric: 90, // 0–100 for α₁ (null when alpha1 is null)
	alpha2ScoreNumeric: 70, // 0–100 for α₂ (null when alpha2 is null)
};
```

### Alpha1 and Alpha2

The DFA algorithm produces three scaling exponents:

- **Global α** (`alpha`) — fitted across all scales, from `minWindow` up to the full series length.
- **α₁** (`alpha1`) — short-term correlations, fitted over scales 4–16 (controlled by `minWindow` and `shortMax`). Captures beat-to-beat dynamics in HRV or local patterns in other time series.
- **α₂** (`alpha2`) — long-term correlations, fitted over scales 16–N/4 (controlled by `longMin` and `longMaxFraction`). Captures slower trends and regulatory patterns.

The crossover point between α₁ and α₂ is at scale 16 by default (both ranges include it).

**When alpha1 or alpha2 is `null`:**

- `alpha1` requires at least 3 valid scales in the 4–16 range. For very short series (N < 8), there may not be enough.
- `alpha2` requires at least 3 valid scales in the 16–N/4 range, each with at least 4 forward segments. This means you need roughly **N ≥ 256** for a reliable alpha2 (since at scale 32, you need 4 segments → N ≥ 128, but the range 16–32 with the segment filter often leaves too few points for shorter series).

For short series, rely on the global `alpha` and `alpha1`. For long series (500+ points), all three exponents are typically available.

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
   - Generate window sizes with α₂ anchors (16, 18, 20, ..., N/4) using linear steps
   - For global α, add geometric progression from N/4 using factor 2^0.25 to the whole window length
   - Cap maximum scale at N/4 for α₂ to ensure ≥4 forward segments

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
   - **α₂**: Long-term correlations (scales 16+, with ≥4 segments)
10. **Interpretation** (see Scoring section below)

### Scoring and Interpretation

The output includes multiple label systems and numeric scores for each alpha (global, α₁, α₂). All categorical labels are affected by the `level` parameter passed to `compute()`.

#### Simplified labels (4 categories)

**`alphaLabel` / `alpha1Label` / `alpha2Label`** — general-purpose labels suitable for any time series.
**`alphaScore` / `alpha1Score` / `alpha2Score`** — legacy HRV-specific labels (same thresholds, different names), designed for biofeedback contexts where α ≈ 1.0 indicates a healthy, resilient heart.

| `alphaLabel` | `alphaScore` | `"moderate"` (default) | `"relaxed"` | `"strict"` |
|---|---|---|---|---|
| `"random"` | `"recovering"` | α ≤ 0.55 | α ≤ 0.65 | α ≤ 0.50 |
| `"regular"` | `"regular"` | 0.55 < α < 0.95 | 0.65 < α < 0.90 | 0.50 < α < 0.98 |
| `"fractal"` | `"resilient"` | 0.95 ≤ α ≤ 1.05 | 0.90 ≤ α ≤ 1.10 | 0.98 ≤ α ≤ 1.02 |
| `"complex"` | `"tension"` | α > 1.05 | α > 1.10 | α > 1.02 |

#### Traditional DFA labels (6 categories)

**`dfaLabel` / `dfa1Label` / `dfa2Label`** — scientific DFA classification with finer granularity, using zones around the three canonical alpha values (0.5, 1.0, 1.5):

| `dfaLabel` | `"moderate"` (default) | `"relaxed"` | `"strict"` | Signal characteristic |
|---|---|---|---|---|
| `"anti-correlated"` | α < 0.45 | α < 0.40 | α < 0.48 | Successive values tend to alternate |
| `"white noise"` | 0.45 ≤ α ≤ 0.55 | 0.40 ≤ α ≤ 0.60 | 0.48 ≤ α ≤ 0.52 | Uncorrelated, random |
| `"correlated"` | 0.55 < α < 0.95 | 0.60 < α < 0.90 | 0.52 < α < 0.98 | Long-range power-law correlations |
| `"1/f noise"` | 0.95 ≤ α ≤ 1.05 | 0.90 ≤ α ≤ 1.10 | 0.98 ≤ α ≤ 1.02 | Scale-invariant, fractal |
| `"strongly correlated"` | 1.05 < α < 1.45 | 1.10 < α < 1.40 | 1.02 < α < 1.48 | Strong persistence, trending |
| `"Brownian motion"` | α ≥ 1.45 | α ≥ 1.40 | α ≥ 1.48 | Random walk, integrated noise |

The simplified `alphaLabel` maps to `dfaLabel` as follows: `"random"` covers both `"anti-correlated"` and `"white noise"`, `"regular"` equals `"correlated"`, `"fractal"` equals `"1/f noise"`, and `"complex"` covers both `"strongly correlated"` and `"Brownian motion"`.

#### Threshold levels

- **`"relaxed"`** — wider zones around key values, useful for noisy real-world data or short recordings where alpha estimates have higher variance.
- **`"moderate"`** — balanced thresholds, suitable for most applications.
- **`"strict"`** — narrow zones, for research contexts requiring precise classification or long, clean recordings.

#### Numeric scores

**`alphaScoreNumeric` / `alpha1ScoreNumeric` / `alpha2ScoreNumeric`** — Numeric score from 0 to 100, where 100 means α = 1.0 (perfect fractal scaling) and 0 means α is ≥1.0 away from the ideal. Not affected by the `level` parameter.

### References

- Peng, C.-K., Havlin, S., Stanley, H.E., & Goldberger, A.L. (1995). Quantification of scaling exponents and crossover phenomena in nonstationary heartbeat time series. _Chaos_, 5(1), 82–87. [doi:10.1063/1.166141](https://doi.org/10.1063/1.166141) — The foundational DFA paper that named the method, applied it to HRV, and introduced the α1/α2 crossover.

- Shaffer, F. & Ginsberg, J.P. (2017). An overview of heart rate variability metrics and norms. _Frontiers in Public Health_, 5, 258. [doi:10.3389/fpubh.2017.00258](https://doi.org/10.3389/fpubh.2017.00258) — Comprehensive review of HRV metrics including DFA, with normative values for clinical and healthy populations.

- Kantelhardt, J.W., Zschiegner, S.A., Koscielny-Bunde, E., Havlin, S., Bunde, A., & Stanley, H.E. (2002). Multifractal detrended fluctuation analysis of nonstationary time series. _Physica A_, 316(1–4), 87–114. [doi:10.1016/S0378-4371(02)01383-3](<https://doi.org/10.1016/S0378-4371(02)01383-3>) — Extension to multifractal DFA (MF-DFA) with detailed discussion of scale selection and segment requirements.

- Hardstone, R., Poil, S.-S., Schiavone, G., Jansen, R., Nikulin, V.V., Mansvelder, H.D., & Linkenkaer-Hansen, K. (2012). Detrended fluctuation analysis: A scale-free view on neuronal oscillations. _Frontiers in Physiology_, 3, 450. [doi:10.3389/fphys.2012.00450](https://doi.org/10.3389/fphys.2012.00450) — Practical guide to DFA parameter choices including minimum segment counts and scale ranges.

- Paranyushkin D (2021). Fractal Variability Movement Feedback System. Nodus Labs. [https://noduslabs.com/featured/fractal-variability-feedback-system/](https://noduslabs.com/featured/fractal-variability-feedback-system/)

### Author

Created by [Dmitry Paranyushkin](https://deemeetree.com) in 2020-2026

### GPL License

This open source, free software is available under the GNU Affero General Public License version 3 (AGPLv3) license.
You can make modifications to this code and binaries based on it, but only on the condition that you provide access to those modifications under the same license (including remotely through a computer network).
It is provided as is, with no guarantees and no liabilities.
You can re-use it as long as you keep this notice inside the code.
