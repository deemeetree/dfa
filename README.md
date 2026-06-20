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

### Multifractal DFA (MFDFA)

Standard `compute()` characterises a signal with a single scaling exponent — it assumes the same fractal behaviour governs both small and large fluctuations (a _monofractal_). Many real signals are **multifractal**: their scaling depends on the magnitude of fluctuations. `computeMultifractal()` captures this by generalising the fluctuation function with a **q-order parameter**.

```js
let dfa = new DFA(time_series);

let mf = dfa.computeMultifractal({
	qMin: -5, // smallest q (emphasises small fluctuations)
	qMax: 5, // largest q (emphasises large fluctuations)
	qStep: 0.5, // step between q values
});
```

> **Note:** the base of the logarithm does _not_ control multifractality. α is the slope of `log F(s)` vs `log s`, and changing the log base scales both axes equally, leaving the slope unchanged. The q parameter — not the log base — is what generalises DFA to the multifractal case.

A single call sweeps the whole q range and returns the generalised Hurst exponent for each q. Mirroring the monofractal `compute()`, three h(q) curves are returned — a **global** one plus **α₁** (short-scale) and **α₂** (long-scale) curves:

```js
{
	q: [-5, -4.5, ..., 0, ..., 4.5, 5],       // the q values swept

	hq: [...],         // GLOBAL h(q) — fitted across all scales, one per q
	hq1: [...],        // α₁ h(q)    — fitted over short scales (4–16), one per q
	hq2: [...],        // α₂ h(q)    — fitted over long scales (16–N/4); null where insufficient

	tau: [...],        // mass exponent τ(q) = q·h(q) − 1 (from the global curve)
	alpha: [...],      // Hölder exponent α / singularity strength (from the global curve)
	falpha: [...],     // singularity spectrum f(α) (from the global curve)

	hMin: 0.55,        // global h(q) at qMax — scaling of large fluctuations
	hMax: 0.92,        // global h(q) at qMin — scaling of small fluctuations
	multifractalWidth: 0.37, // width of the GLOBAL curve (hMax − hMin)
	width1: 0.87,      // width of the α₁ curve
	width2: 0.26,      // width of the α₂ curve

	monofractal: {     // h(2) per curve — matches compute() alpha/alpha1/alpha2
		alpha: 0.71, alpha1: 0.68, alpha2: 0.73,
	},
	ranges: {          // scale ranges each curve was fitted over
		global: [4, N], alpha1: [4, 16], alpha2: [16, N/4],
	},

	scales: [...],            // window sizes used (same as compute())
	fluctuationsByQ: [[...]], // F_q(s) per q, for plotting/debugging
	lengthOfData: N,
}
```

**How to read it:**

- `q`, `hq`, `hq1`, `hq2`, `tau`, `alpha`, `falpha` are all **index-aligned** — `hq[i]` is the global Hurst exponent for `q[i]`, `hq1[i]` the short-scale one, `hq2[i]` the long-scale one.
- **Three curves, same meaning as monofractal α/α₁/α₂:** `hq` is the global multifractal scaling, `hq1` captures fast/short-scale dynamics, `hq2` slow/long-scale dynamics. Comparing their widths (`multifractalWidth`, `width1`, `width2`) shows whether the multifractality lives mainly in the fast or the slow part of the signal.
- Plot any curve against `q`: a roughly **flat** line means **monofractal** behaviour in that range; a **downward slope** means **multifractal**.
- Negative q emphasise small fluctuations, positive q emphasise large ones. `q = 2` reproduces standard DFA, so `monofractal.alpha`/`alpha1`/`alpha2` match `compute().alpha`/`alpha1`/`alpha2` exactly. `q = 0` uses a logarithmic-averaging special case internally (the `1/q` form is undefined there).
- `multifractalWidth` (= `hMax − hMin`) is the single-number summary most often reported: ≈ 0 indicates monofractal, larger values indicate a richer multifractal structure.
- `alpha`/`falpha` are the standard **multifractal spectrum** (derived from the global curve) — for a multifractal signal, plotting `falpha` against `alpha` gives an inverted-parabola shape.
- **`hq2`/`width2` are `null` for short series** (same N ≥ ~256 requirement as `compute().alpha2`). Fall back to `hq` and `hq1` in that case.

**Typical parameter values:** q in the range **−5 to +5** (sometimes −10 to +10), with a step of 0.25, 0.5, or 1. `computeMultifractal()` accepts the same scale-generation parameters as `compute()` (`minWindow`, `step`, `expStep`, `shortMax`, `longMin`, `longMaxFraction`) and reuses the same scales and the same α₁/α₂ fit ranges, so the q = 2 curves match the monofractal `alpha`/`alpha1`/`alpha2` exactly. `compute()` itself is unchanged and remains the q = 2 special case.

### Understanding h(q), multifractalWidth, and the multifractal curves

Standard DFA gives you **one** number. Multifractal DFA gives you a **spectrum** of numbers — one scaling exponent for each q.

- In standard DFA, `alpha` is a single scaling exponent.
- In MFDFA, `h(q)` is the scaling exponent _for each value of q_.
- `h(2)` is the closest equivalent to standard DFA alpha, because ordinary DFA is built on the second-order (squared) fluctuation function.
- **Negative q** values emphasize the **small** fluctuations in the signal.
- **Positive q** values emphasize the **large** fluctuations.
- If `h(q)` is almost **flat** across q, the signal is closer to **monofractal** — small and large fluctuations scale the same way.
- If `h(q)` **changes strongly** across q, the signal is **multifractal** — small and large fluctuations scale differently.

A small example of the `q → hq` curve:

```js
q:  [-5, -4.5, -4, ..., 0, ..., 4, 4.5, 5]
hq: [1.2, 1.18, 1.15, ..., 1.0, ..., 0.82, 0.79, 0.76]
```

Plotting `hq` against `q` draws this curve. The package gives you **three** such curves:

```js
q → hq   // global multifractal curve (all scales)
q → hq1  // short-scale multifractal curve
q → hq2  // long-scale multifractal curve
```

`hq2` (and `width2`) can be `null` if the input series is too short for reliable long-scale fitting — see the short-series note above.

#### multifractalWidth

`multifractalWidth` is a one-number summary of **how different** the scaling behaviour is between small and large fluctuations. It is calculated as:

```js
multifractalWidth = hMax - hMin
```

where, for the global curve:

```js
hMax = h(qMin)   // small fluctuations, e.g. q = -5
hMin = h(qMax)   // large fluctuations, e.g. q = 5
```

Note that `multifractalWidth` only compares the **extremes** of the h(q) curve. The intermediate `h(q)` values still matter — they reveal the **shape** of the multifractal spectrum, which the width alone cannot show:

```js
Signal A:
h(q) smoothly decreases from 1.2 to 0.8
width = 0.4

Signal B:
h(q) stays near 1.2, then suddenly drops near q = 4
width = 0.4
```

Both signals have the **same** `multifractalWidth`, but very different internal dynamics. So treat `multifractalWidth` as a useful summary, and read the full `hq` curve when you need the actual structure.

#### Curve-shape metrics

Because the width alone hides the shape (as shown above), the result also includes four numbers that summarise the **shape** of the global `hq` curve. They are derived from the global `q → hq` curve:

```js
result.multifractalWidth   // difference between the extremes (hMax − hMin)
result.hCurveSlope         // general direction of h(q) vs q
result.hCurveCurvature     // U-shape vs inverted-U shape
result.hCurveNonlinearity  // how much the curve deviates from a straight line
result.hMinLocation        // the q where persistence (h) is lowest
```

What each one means:

- **`hCurveSlope`** — the slope of the best-fit straight line through `h(q)`. It captures the **overall direction**. Almost always negative (h decreases as q increases); a steeper negative slope means a stronger overall difference between small- and large-fluctuation scaling. A slope near 0 means a flat, near-monofractal curve.
- **`hCurveCurvature`** — the second derivative (`2c`) of a quadratic fit to `h(q)`. It tells you the **bend** of the curve: **> 0** is a **U-shape** (convex), **< 0** is an **inverted-U** (concave), and **≈ 0** is roughly straight. This separates the two example signals above that share the same width.
- **`hCurveNonlinearity`** — the RMS deviation of `h(q)` from the straight-line fit. **0** means a perfectly straight curve; larger values mean the multifractal spectrum has real structure (bends, kinks, asymmetry) that the slope and width don't capture.
- **`hMinLocation`** — the q value at which `h(q)` is smallest, i.e. **where persistence is lowest**. For a monotonically decreasing curve this is simply `qMax`; an interior minimum signals a more complex, non-monotonic spectrum.

The four fields above describe the **global** curve. The same four metrics are also provided for the **α₁** (short-scale) and **α₂** (long-scale) curves, using the `1`/`2` suffix (the same convention as `width1`/`width2`):

```js
// α₁ (short-scale) curve shape
result.hCurveSlope1
result.hCurveCurvature1
result.hCurveNonlinearity1
result.hMinLocation1

// α₂ (long-scale) curve shape
result.hCurveSlope2
result.hCurveCurvature2
result.hCurveNonlinearity2
result.hMinLocation2
```

This lets you see **where** the multifractal structure lives: for example, a large `hCurveNonlinearity1` with a near-zero `hCurveNonlinearity2` means the rich multifractal behaviour is concentrated in the fast/short-scale dynamics, while the slow/long-scale part stays close to monofractal.

All shape metrics are `null` for curves too short to fit (fewer than 3 valid q points). The α₂ shape fields are also `null` for short series, the same way `hq2`/`width2` are.

### Interpreting MFDFA in practice

- **DFA alpha** tells you whether the signal has _one overall_ fractal scaling pattern.
- **MFDFA** tells you whether small, medium, and large fluctuations follow the _same_ or _different_ scaling patterns.
- A **low** `multifractalWidth` means the signal is closer to monofractal — one dominant scaling logic.
- A **higher** `multifractalWidth` means the signal is more heterogeneous — small and large fluctuations behave differently across scales.
- **Very high** width can indicate rich adaptive complexity, but it can also reflect noise, bursts, outliers, or nonstationarity, depending on the data.

> **Warning:** do not interpret a larger multifractal width as automatically "better." Its meaning depends on the domain, the preprocessing, the data length, and whether the signal contains artifacts or outliers.

### Example: financial time series

In finance, MFDFA is usually applied to **log returns**, **absolute returns**, or **squared returns** — not to raw prices:

```js
price → log returns
price → absolute returns / squared returns
```

Financial markets often contain volatility clustering, fat tails, regime changes, and different behaviour between calm periods and crisis-like periods. MFDFA is useful here because it can separate the scaling behaviour of small fluctuations from large fluctuations.

| Output                       | Possible financial interpretation                                |
| ---------------------------- | ---------------------------------------------------------------- |
| `h(2) ≈ 0.5`                 | returns are close to random-walk-like                            |
| `h(2) > 0.5`                 | persistence / trend-like memory                                  |
| `h(2) < 0.5`                 | anti-persistence / mean-reversion tendency                       |
| large `multifractalWidth`    | small and large moves scale differently                          |
| rising `multifractalWidth`   | market dynamics becoming more heterogeneous or regime-fragmented |
| strong positive-q distortion | large moves / jumps / crashes have different scaling             |
| strong negative-q distortion | calm-period microstructure has its own scaling                   |

This is **not** a trading signal by itself. It is a descriptive tool for market structure, volatility regimes, and risk dynamics.

A simple rolling-window idea:

```js
// 1. Convert prices to log returns
// 2. Run computeMultifractal() on a rolling window
// 3. Track h(2), multifractalWidth, hq, width1, and width2 over time
```

A rolling increase in `multifractalWidth` can indicate that the market is becoming less uniform and more regime-dependent. But interpretation requires comparison against the **same asset, same timeframe, and same preprocessing method**.

### Recommended outputs to inspect

```js
result.monofractal.alpha   // h(2), closest to standard DFA alpha
result.hq                  // full q → h(q) curve
result.multifractalWidth   // summary width between q extremes
result.hCurveSlope         // general direction of the h(q) curve
result.hCurveCurvature     // U-shape (>0) vs inverted-U (<0)
result.hCurveNonlinearity  // deviation from a straight line
result.hMinLocation        // q where persistence is lowest
result.hq1                 // short-scale h(q)
result.hq2                 // long-scale h(q), if available
result.width1              // short-scale multifractal width
result.width2              // long-scale multifractal width, if available
result.hCurveSlope1        // α₁ curve shape (also Curvature1/Nonlinearity1/hMinLocation1)
result.hCurveSlope2        // α₂ curve shape (also Curvature2/Nonlinearity2/hMinLocation2)
result.alpha               // singularity strengths for f(alpha)
result.falpha              // multifractal singularity spectrum
```

**Important naming clarification:** in the multifractal output, `alpha` is **not** the same as standard DFA `alpha`.

- `monofractal.alpha` is the DFA-like alpha, approximately `h(2)`.
- `alpha` together with `falpha` is the **singularity spectrum** used in formal multifractal analysis.

MFDFA is a **descriptive** analysis method — it characterises structure in data. It is not a prediction engine.

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
