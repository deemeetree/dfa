/* 
    Detrended Fluctuation Analysis (DFA)
    Javascript port
    available on https://github.com/deemeetree/dfa

    based on the approach described on https://noduslabs.com/featured/fractal-variability-feedback-system/
    based on the python script https://github.com/dokato/dfa


    ~ Concept 

    DFA is used to measure the behaviour of a time series. 
    
    An obtained alpha component (closely related to the Hurst exponent) will indicate 
    the presence of correlations in the time series. 

    It is based on the relationship between the length of an observation and cumulated variability.
    

    ~ Algorithm:

    1) represent a time series as a one-dimensional vector
    2) transform it into cumulative sum of variances from the mean
    3) generate the scales of different window sizes for the time series — each containing the number of observations
    4) split the cumulative variance vector into those chunks
    5) for each chunk: generate polynomial in the range in a window (to detrend it)
    6) calculate RMS of the fluctuations of the original from the fit
    7) get the average RMS value for each window in the scale
    8) take the average for each squared RMS for each scale
    9) if they align along in a straight line in loglog plot, there is power law relation
    i.e. as observations increase in length, the amplitude of average fluctuations increases as well 
    that is, there is an exponential growth of fluctuations when there's an exponential growth of scale (2^n)
    on the smaller scales the deviations are smaller. on the bigger scales they are much bigger.

    ~ Simple use

    let time_series = [8, 10, 6, 9, 7, 5, 5, 11, 11, 8, 6, 7, 9, 10, 7, 9]
    let dfa = new DFA(time_series)
    let alpha_component = dfa.compute()
    console.log(alpha_component)

    */

let DFA = (function () {
	// Helper functions
	function mean(x) {
		return x.reduce((a, b) => a + b, 0) / x.length;
	}

	function variancePopulation(x) {
		const m = mean(x);
		let s = 0;
		for (let i = 0; i < x.length; i++) {
			const d = x[i] - m;
			s += d * d;
		}
		return s / x.length; // population variance
	}

	function sdnn(x) {
		return Math.sqrt(variancePopulation(x));
	}

	function rmssd(x) {
		if (x.length < 2) return 0;
		let s = 0;
		for (let i = 1; i < x.length; i++) {
			const d = x[i] - x[i - 1];
			s += d * d;
		}
		return Math.sqrt(s / (x.length - 1));
	}

	function pnn50(x) {
		if (x.length < 2) return 0;
		let c = 0;
		for (let i = 1; i < x.length; i++) {
			if (Math.abs(x[i] - x[i - 1]) > 50) c++;
		}
		return (c / (x.length - 1)) * 100;
	}

	function integratedProfile(rr) {
		const m = mean(rr);
		const y = new Array(rr.length);
		let acc = 0;
		for (let i = 0; i < rr.length; i++) {
			acc += rr[i] - m;
			y[i] = acc;
		}
		return { profile: y, meanValue: m };
	}

	function generateScales(
		N,
		{
			minWindow = 4,
			step = 2, // linear step for small scales (and short series)
			expStep = 0.5, // log2 spacing for geometric progression (≈√2)
			shortMax = 16,
			longMin = 17,
			longMaxFraction = 0.25,
			hardLongMax = 128, // Changed from 64 to 128
		} = {}
	) {
		// Always use the full data length (N) up to hardLongMax, no fraction limitation
		const maxReli = Math.min(N, hardLongMax);
		const S = new Set();

		// α1 anchors: minWindow..shortMax inclusive, in +step increments
		for (let s = minWindow; s <= Math.min(shortMax, maxReli); s += step)
			S.add(s);

		// α2 candidates:
		const a2Lo = Math.max(longMin, minWindow + step); // ensure > α1
		// geometric progression over [a2Lo..maxReli]
		let g = [];
		if (a2Lo <= maxReli) {
			let s = a2Lo;
			const mult = Math.pow(2, expStep);
			while (s <= maxReli) {
				g.push(Math.round(s));
				s *= mult;
			}
		}

		// Short-series fallback: if too few geometric points, fill with linear increments
		if (g.length < 3) {
			for (let s = a2Lo; s <= maxReli; s += step) S.add(s);
		} else {
			for (const s of g) S.add(s);
		}

		// Ensure at least 3 scales overall (global fit needs ≥3 points)
		let scales = Array.from(S).filter((s) => s >= minWindow && s <= maxReli);
		scales.sort((a, b) => a - b);
		if (scales.length < 3) {
			// minimal fallback
			scales = [
				minWindow,
				Math.max(minWindow + step, Math.floor((minWindow + maxReli) / 2)),
				maxReli,
			]
				.filter((v, i, a) => a.indexOf(v) === i)
				.sort((a, b) => a - b);
		}
		return scales;
	}

	function linearRegression(x, y) {
		const n = x.length;
		if (n === 0) return { slope: 0, intercept: 0 };
		let sx = 0,
			sy = 0,
			sxx = 0,
			sxy = 0;
		for (let i = 0; i < n; i++) {
			sx += x[i];
			sy += y[i];
			sxx += x[i] * x[i];
			sxy += x[i] * y[i];
		}
		const denom = n * sxx - sx * sx;
		if (denom === 0) return { slope: 0, intercept: sy / n };
		const slope = (n * sxy - sx * sy) / denom;
		const intercept = (sy - slope * sx) / n;
		return { slope, intercept };
	}

	function detrendLinear(ySeg) {
		// returns fitted values for the segment
		const n = ySeg.length;
		const x = new Array(n);
		for (let i = 0; i < n; i++) x[i] = i;
		const { slope, intercept } = linearRegression(x, ySeg);
		const trend = new Array(n);
		for (let i = 0; i < n; i++) trend[i] = slope * i + intercept;
		return trend;
	}

	function fluctuationForScale(profile, s, minWindow = 4) {
		const n = profile.length;
		const numSegments = Math.floor(n / s);
		if (numSegments <= 0) return { F: 0, segments: 0 };

		const rms = [];

		// forward segments
		for (let seg = 0; seg < numSegments; seg++) {
			const start = seg * s,
				end = start + s;
			const segData = profile.slice(start, end);
			const trend = detrendLinear(segData);
			let ss = 0;
			for (let i = 0; i < s; i++) {
				const r = segData[i] - trend[i];
				ss += r * r;
			}
			rms.push(Math.sqrt(ss / s));
		}

		// backward segments if remainder is sufficiently long
		const remainder = n % s;
		if (remainder >= minWindow) {
			for (let seg = 0; seg < numSegments; seg++) {
				const end = n - seg * s;
				const start = end - s;
				if (start < 0) break;
				const segData = profile.slice(start, end);
				const trend = detrendLinear(segData);
				let ss = 0;
				for (let i = 0; i < s; i++) {
					const r = segData[i] - trend[i];
					ss += r * r;
				}
				rms.push(Math.sqrt(ss / s));
			}
		}

		// canonical DFA aggregation: F(s) = sqrt( mean( RMS^2 ) )
		if (rms.length === 0) return { F: 0, segments: numSegments };
		let meanSq = 0;
		for (let i = 0; i < rms.length; i++) meanSq += rms[i] * rms[i];
		meanSq /= rms.length;

		return { F: Math.sqrt(meanSq), segments: numSegments };
	}

	function fitAlphaInRange(
		scales,
		logs,
		range,
		segments = null,
		minSegments = 0
	) {
		const xs = [],
			ys = [],
			used = [];
		for (let i = 0; i < scales.length; i++) {
			const s = scales[i];
			const ylog = logs.fluct[i];
			if (s < range[0] || s > range[1]) continue;
			if (ylog === null) continue; // filtered out (e.g., F==0)
			if (segments && segments[i] < minSegments) continue;
			xs.push(logs.scale[i]);
			ys.push(ylog);
			used.push(s);
		}
		if (xs.length < 3) return { alpha: null, usedRange: null };
		const { slope } = linearRegression(xs, ys);
		return {
			alpha: slope,
			usedRange: [Math.min(...used), Math.max(...used)],
		};
	}

	function alphaScoreNumeric(alpha) {
		const dev = Math.abs(alpha - 1.0);
		return Math.max(0, (1.0 - dev) * 100);
	}

	// DFA class definition
	let DFA = function (x) {
		if (!(this instanceof DFA)) {
			return new DFA(x);
		}

		if (
			!(
				x instanceof Array ||
				x instanceof Float32Array ||
				x instanceof Float64Array
			)
		) {
			throw new Error("x must be an array");
		}

		this.x = x;
	};

	// Legacy methods kept for backward compatibility
	DFA.prototype.cumsumVariance = function (x, mean) {
		let cumulativeSumVariance = ((sum) => (value) => (sum += value - mean))(
			0
		);
		return x.map(cumulativeSumVariance);
	};

	DFA.prototype.averageVariance = function (x, mean) {
		// Note: This is actually mean absolute deviation, kept for compatibility
		let totalDeviation = 0;
		for (let i = 0; i < x.length; i++) {
			totalDeviation += Math.abs(x[i] - mean);
		}
		return totalDeviation / x.length;
	};

	DFA.prototype.deviationsFromMean = function (x, mean) {
		const deviationsVector = x.map((value) => value - mean);
		return deviationsVector;
	};

	DFA.prototype.squareVector = function (x) {
		return x.map((value) => Math.pow(value, 2));
	};

	DFA.prototype.SDNN = function (x) {
		return sdnn(x);
	};

	DFA.prototype.averageDifferences = function (x) {
		let sumOfDifferences = 0;
		for (let i = 0; i < x.length - 1; i++) {
			sumOfDifferences += Math.abs(x[i + 1] - x[i]);
		}
		return sumOfDifferences / (x.length - 1);
	};

	DFA.prototype.RMSSD = function (x) {
		return rmssd(x);
	};

	DFA.prototype.naturalLog = function (value) {
		return Math.log(value);
	};

	DFA.prototype.PNN50 = function (rrIntervals) {
		return pnn50(rrIntervals);
	};

	DFA.prototype.meanOfVector = function (x) {
		let y = x || this.x;
		return mean(y);
	};

	DFA.prototype.meanOfVectorNoBias = function (x) {
		let y = x || this.x;
		let meanNoBias = (arr) => arr.reduce((p, c) => p + c, 0) / (arr.length - 1);
		return meanNoBias(y);
	};

	DFA.prototype.alphaScore = function (alpha, level = "moderate") {
		// Legacy categorical score - kept for compatibility
		let thresholds = [0.55, 0.95, 1.05];
		if (level == "relaxed") {
			thresholds = [0.65, 0.9, 1.1];
		} else if (level == "strict") {
			thresholds = [0.5, 0.98, 1.02];
		}

		if (alpha <= thresholds[0]) {
			return "recovering";
		} else if (alpha > thresholds[0] && alpha < thresholds[1]) {
			return "regular";
		} else if (alpha >= thresholds[1] && alpha <= thresholds[2]) {
			return "resilient";
		} else if (alpha > thresholds[2]) {
			return "tension";
		}
	};

	DFA.prototype.compute = function (
		minWindow = 4,
		step = 2, // linear increment for α1 and for short series fallback
		expStep = 0.5, // geometric log2 step for longer series (≈√2)
		shortMax = 16,
		longMin = 17,
		longMaxFraction = 0.25
	) {
		const rr = this.x;
		const N = rr.length;
		if (N <= minWindow) {
			return {
				averageVariance: 0,
				meanValue: 0,
				lengthOfData: N,
				SDNN: 0,
				RMSSD: 0,
				lnRMSSD: 0,
				PNN50: 0,
				averageDifferences: 0,
				scales: [],
				segments: [],
				fluctuations: [],
				scalesLog: [],
				fluctuationsLog: [],
				scales_log: [], // Legacy alias
				fluctuations_log: [], // Legacy alias
				coefficients: { slope: 0, intercept: 0 },
				alpha: 1.0,
				alphaScore: 0,
				alpha1: null,
				alpha2: null,
				alpha1Range: null,
				alpha2Range: null,
			};
		}

		// Basic stats (parity with Swift)
		const meanValue = mean(rr);
		const avgVar = variancePopulation(rr);
		const SDNN = Math.sqrt(avgVar);
		const RMSSD = rmssd(rr);
		const lnRMSSD = RMSSD > 0 ? Math.log(RMSSD) : 0;
		const PNN50 = pnn50(rr);
		// Keep averageDifferences to match legacy JS result
		let sumAbs = 0;
		for (let i = 1; i < N; i++) sumAbs += Math.abs(rr[i] - rr[i - 1]);
		const averageDifferences = N > 1 ? sumAbs / (N - 1) : 0;

		// Profile
		const { profile } = integratedProfile(rr);

		// Scales
		const scales = generateScales(N, {
			minWindow,
			step,
			expStep,
			shortMax,
			longMin,
			longMaxFraction,
		});
		const segments = new Array(scales.length);
		const fluctuations = new Array(scales.length);

		// F(s) per scale (forward+backward); record forward segments only
		for (let i = 0; i < scales.length; i++) {
			const s = scales[i];
			const { F, segments: segs } = fluctuationForScale(profile, s, minWindow);
			fluctuations[i] = F;
			segments[i] = segs; // parity with Swift (forward only)
		}

		// Filter zeros before logs
		const scalesLog = new Array(scales.length);
		const fluctuationsLog = new Array(scales.length).fill(null);
		for (let i = 0; i < scales.length; i++) {
			scalesLog[i] = Math.log(scales[i]); // natural log
			if (fluctuations[i] > 0)
				fluctuationsLog[i] = Math.log(fluctuations[i]);
		}

		// Global α on valid (non-null) points
		const xAll = [],
			yAll = [];
		for (let i = 0; i < scales.length; i++) {
			if (fluctuationsLog[i] !== null) {
				xAll.push(scalesLog[i]);
				yAll.push(fluctuationsLog[i]);
			}
		}
		const coefficients = linearRegression(xAll, yAll);
		const alpha = coefficients.slope;
		const alphaScoreVal = alphaScoreNumeric(alpha);
		const alphaScoreCat = this.alphaScore(alpha); // Legacy categorical score

		// α1 (minWindow..shortMax)
		const { alpha: alpha1, usedRange: alpha1Range } = fitAlphaInRange(
			scales,
			{ scale: scalesLog, fluct: fluctuationsLog },
			[minWindow, Math.min(shortMax, N)],
			null,
			0
		);

		// α2 (longMin..maxReliable), require ≥4 forward segments
		// Cap at N/4 for alpha2 reliability, but still max at 128
		const maxReliableAlpha2 = Math.min(Math.floor(N * 0.25), 128);
		const { alpha: alpha2, usedRange: alpha2Range } = fitAlphaInRange(
			scales,
			{ scale: scalesLog, fluct: fluctuationsLog },
			[longMin, maxReliableAlpha2],
			segments,
			4
		);

		return {
			averageVariance: avgVar,
			meanValue,
			lengthOfData: N,
			SDNN,
			RMSSD,
			lnRMSSD,
			PNN50,
			averageDifferences,
			scales,
			segments,
			fluctuations,
			scalesLog,
			fluctuationsLog,
			scales_log: scalesLog, // Legacy alias
			fluctuations_log: fluctuationsLog, // Legacy alias
			coefficients,
			alpha,
			alphaScore: alphaScoreCat, // Legacy categorical
			alphaScoreNumeric: alphaScoreVal, // New numeric score
			alpha1,
			alpha2,
			alpha1Range,
			alpha2Range,
		};
	};

	// Deprecated/unused methods kept for backward compatibility
	DFA.prototype.generateRange = function (array, startAt = 0, step = 1) {
		let size = array.length;
		let result = [];
		for (let i = 0; i < size; i++) {
			result.push(i * step + startAt);
		}
		return result;
	};

	DFA.prototype.generateZeroArray = function (size) {
		let result = [];
		for (let i = 0; i < size; i++) {
			result.push(0);
		}
		return result;
	};

	DFA.prototype.generateScaleRange = function (array, startAt = 4, step = 0.5) {
		// Legacy method - deprecated but kept for compatibility
		let startPow = Math.sqrt(startAt);
		let size = Math.floor((Math.log2(array.length) - startPow) / step + 1);
		let result = [];
		for (let i = 0; i < size; i++) {
			result.push(Math.floor(Math.pow(2, i * step + startPow)));
		}
		return result;
	};

	DFA.prototype.splitArrayIntoChunks = function (arr, size, strict = false) {
		// Legacy method - deprecated but kept for compatibility
		let chunkedArray = [];
		for (let i = 0; i < arr.length; i += size) {
			chunkedArray.push(arr.slice(i, i + size));
		}
		if (strict) {
			if (chunkedArray[chunkedArray.length - 1].length < size) {
				chunkedArray.pop();
			}
		}
		return chunkedArray;
	};

	DFA.prototype.log2Vector = function (x) {
		// Legacy method - deprecated but kept for compatibility
		if (!(x instanceof Array)) {
			throw new Error("log2Vector input must be an array");
		}
		let y = [];
		for (let i in x) {
			y[i] = Math.log2(x[i]);
		}
		return y;
	};

	DFA.prototype.RMS = function (cut, fit) {
		// Legacy method - deprecated but kept for compatibility
		if (cut.length != fit.length) {
			throw new Error(
				"for calculating the RMS of the two vectors they should be equal"
			);
		}
		let variance_vector = [];
		for (let i = 0; i < cut.length; i++) {
			variance_vector[i] = Math.pow(cut[i] - fit[i], 2);
		}
		return Math.sqrt(mean(variance_vector));
	};

	return DFA;
})();

module.exports = DFA;