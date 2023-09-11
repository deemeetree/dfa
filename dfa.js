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
    5) for each chunk: generate polinomial in the range in a window (to detrend it)
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

var DFA = (function () {
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

	DFA.prototype.cumsumVariance = function (x, mean) {
		let cumulativeSumVariance = (
			(sum) => (value) =>
				(sum += value - mean)
		)(0);
		return x.map(cumulativeSumVariance);
	};

	DFA.prototype.averageVariance = function (x, mean) {
		// Sum the absolute deviations
		let totalDeviation = 0;
		for (let i = 0; i < x.length; i++) {
			totalDeviation += Math.abs(x[i] - mean);
		}

		// Return the average deviation
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
		const deviationsFromMean = this.deviationsFromMean(x, this.meanOfVector(x));
		const squaredVector = this.squareVector(deviationsFromMean);
		const meanSquaredVector = this.meanOfVectorNoBias(squaredVector);
		return Math.pow(meanSquaredVector, 0.5);
	};

	DFA.prototype.averageDifferences = function (x) {
		let sumOfDifferences = 0;

		// Loop through the intervals and calculate the sum of differences
		for (let i = 0; i < x.length - 1; i++) {
			sumOfDifferences += Math.abs(x[i + 1] - x[i]);
		}

		// Return the average difference
		return sumOfDifferences / (x.length - 1);
	};

	DFA.prototype.RMSSD = function (x) {
		let sumOfSquaredDifferences = 0;

		// Loop through the intervals and calculate the sum of squared differences
		for (let i = 0; i < x.length - 1; i++) {
			const difference = Math.abs(x[i + 1] - x[i]);
			sumOfSquaredDifferences += difference * difference;
		}

		// Calculate RMSSD
		const rmssd = Math.sqrt(sumOfSquaredDifferences / (x.length - 1));

		return rmssd;
	};

	DFA.prototype.naturalLog = function (value) {
		return Math.log(value);
	};

	DFA.prototype.PNN50 = function (rrIntervals) {
		let differences = [];
		let countDifferencesOver50ms = 0;

		// Calculate differences between successive RR intervals
		for (let i = 1; i < rrIntervals.length; i++) {
			differences.push(Math.abs(rrIntervals[i] - rrIntervals[i - 1])); // absolute difference
			if (differences[i - 1] > 50) {
				countDifferencesOver50ms++;
			}
		}

		// Calculate pNN50
		let pNN50 = (countDifferencesOver50ms / (rrIntervals.length - 1)) * 100;

		return pNN50;
	};

	DFA.prototype.generateRange = function (array, startAt = 0, step = 1) {
		let size = array.length;
		return [...Array(size).keys()].map((i) => i * step + startAt);
	};

	DFA.prototype.generateZeroArray = function (size) {
		return [...Array(size).keys()].map((i) => 0);
	};

	DFA.prototype.generateScaleRange = function (array, startAt = 4, step = 0.5) {
		// the size of the scale should be not longer than the array
		let startPow = Math.sqrt(startAt);
		let size = Math.floor((Math.log2(array.length) - startPow) / step + 1);

		return [...Array(size).keys()].map((i) =>
			Math.floor(Math.pow(2, i * step + startPow))
		); // Math.round could be used, but .floor makes it possible to have more windows
	};

	DFA.prototype.splitArrayIntoChunks = function (arr, size, strict = false) {
		let chunkedArray = [];
		for (var i = 0; i < arr.length; i += size) {
			chunkedArray.push(arr.slice(i, i + size));
		}
		if (strict) {
			if (chunkedArray[chunkedArray.length - 1].length < size) {
				chunkedArray.pop();
			}
		}
		return chunkedArray;
	};

	DFA.prototype.meanOfVector = function (x) {
		!x ? (y = this.x) : (y = x);
		let mean = (y) => y.reduce((p, c) => p + c, 0) / y.length;
		return mean(y);
	};

	DFA.prototype.meanOfVectorNoBias = function (x) {
		!x ? (y = this.x) : (y = x);
		let mean = (y) => y.reduce((p, c) => p + c, 0) / (y.length - 1);
		return mean(y);
	};

	DFA.prototype.log2Vector = function (x) {
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
		if (cut.length != fit.length) {
			throw new Error(
				"for calculating the RMS of the two vectors they should be equal"
			);
		}

		// let's build the variance vector for the 2 vectors
		let variance_vector = [];
		for (let i = 0; i < cut.length; i++) {
			variance_vector[i] = Math.pow(cut[i] - fit[i], 2);
		}

		// let's calculate the average for that vector and take a square root of it
		return Math.sqrt(this.meanOfVector(variance_vector));
	};

	DFA.prototype.alphaScore = function (alpha, level = "moderate") {
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

	DFA.prototype.compute = function (min_window = 4, step = 0.5) {
		/* 
        Transform the series into the cumulative sum of variances (deviations from the mean).
        We need this to only look at the fluctuations and not the absolute values.
        */
		let cumsum = this.cumsumVariance(this.x, this.meanOfVector());

		let averageVariance = this.averageVariance(
			this.x,
			this.meanOfVector(this.x)
		);

		let meanValue = this.meanOfVector(this.x);

		let SDNN = this.SDNN(this.x);

		let RMSSD = this.RMSSD(this.x);

		let lnRMSSD = this.naturalLog(RMSSD);

		let PNN50 = this.PNN50(this.x);

		let averageDifferences = this.averageDifferences(this.x);

		let lengthOfData = this.x.length;

		/* 
        Create observation windows of length, e.g scales = [4,5,8,11,16]
        Scale of the window size from min_window to max series length.
        Each increment is made at the rate of 2^step.
        This is done to have exponential increase of the scale windows.
        */
		let scales = this.generateScaleRange(cumsum, min_window, step);

		/* 
        Split the cumulative sum into the window chunks of scales.
        E.g. for a series of 16 values, the first window of size 4 will split it into 4 windows of 4 values in each.
        */
		let split_scales = {};

		let fluctuations = [];

		for (let s in scales) {
			// split cumsum into the window sizes of scales[s]
			let window_size = scales[s];

			split_scales[window_size] = this.splitArrayIntoChunks(
				cumsum,
				window_size,
				true
			);

			let num_windows = split_scales[window_size].length;

			/* 
            For each window let's build a polynomial in order to detrend it.
            Then we will calculate RMS (root mean square) for each window (standard deviation)
            */

			let rms = this.generateZeroArray(num_windows);

			for (let chunk in split_scales[window_size]) {
				// a chunk window of the scale w
				let xcut = split_scales[window_size][chunk];

				// values vector
				let vector = this.generateRange(xcut, 0, 1);

				// calculate polynomial fit of the cumsum in the split scale
				// e.g. make a straight line that fits values found in that window row
				let poly = new Polyfit(vector, xcut);

				// use polyfit coefficients to build a straight line xfit along the split_scales[w] data
				// using the coefficiens that can be retried via poly.computeCoefficients(1)
				// this is how we detrend the data
				let coeff = poly.getPolynomial(1); // 1 is the degree
				let xfit = [];
				for (let c in vector) {
					xfit[c] = coeff(c);
				}

				// let's now calculate the root mean square off the differences
				// between the values in this window chunk (xcut) and the fit (xfit)
				rms[chunk] = Math.pow(this.RMS(xcut, xfit), 2);
			}

			// let's get the average of the variance for each window in this scale
			// and calculate its deviation
			fluctuations[s] = Math.sqrt(this.meanOfVector(rms));
		}
		/* 
        Calculate best-fit for the Scales to Fluctuations (RMS)
        We now have the vector with scales that shows the number of values in each window.
        We also have the fluctuations vector, which shows the square root of the mean RMS.
        We plot them against each other on the log base 2 scales and find the coefficient of the best-fit.
        */

		let scales_log = this.log2Vector(scales);
		let flucts_log = this.log2Vector(fluctuations);
		let alpha_poly = new Polyfit(scales_log, flucts_log);

		let coefficients = alpha_poly.computeCoefficients(1);

		let alpha = coefficients[1];

		let arraySize = scales.length;

		let alphaScore = this.alphaScore(alpha);

		let result = {
			averageVariance: averageVariance,
			meanValue: meanValue,
			lengthOfData: lengthOfData,
			SDNN: SDNN,
			RMSSD: RMSSD,
			lnRMSSD: lnRMSSD,
			PNN50: PNN50,
			averageDifferences: averageDifferences,
			scales: scales,
			segments: arraySize,
			fluctuations: fluctuations,
			scales_log: scales_log,
			fluctuations_log: flucts_log,
			coefficients: coefficients,
			alpha: alpha,
			alphaScore: alphaScore,
		};

		return result;
	};

	/* 
    POLINOMIAL CALCULATION
    taken from https://github.com/rfink/polyfit.js
    
    ~ Concept:
    
    What is a function that best describes a relationship? 
    n is the degree, so if n = 1, what is the straight line that describes the relationship

    ~ Use: 
    
    let poly = new Polyfit([0,1,2,3],[0,0,-2,-3]) 
    poly.computeCoefficients(1) // 1 is the degree
    // result: [0.4, -1.1] - python DFA algorithm yields reverse order
    let func = poly.getPolynomial(1)
    func(0)
    func(1)
    etc. for the fit

    */

	let Polyfit = function (x, y) {
		// Make sure we return an instance
		if (!(this instanceof Polyfit)) {
			return new Polyfit(x, y);
		}
		// Check that x any y are both arrays of the same type
		if (
			!(
				(x instanceof Array && y instanceof Array) ||
				(x instanceof Float32Array && y instanceof Float32Array) ||
				(x instanceof Float64Array && y instanceof Float64Array)
			)
		) {
			throw new Error("x and y must be arrays");
		}
		if (x instanceof Float32Array) {
			this.FloatXArray = Float32Array;
		} else if (x instanceof Float64Array) {
			this.FloatXArray = Float64Array;
		}
		// Make sure we have equal lengths
		if (x.length !== y.length) {
			throw new Error("x and y must have the same length");
		}
		this.x = x;
		this.y = y;
	};
	/**
	 * Perform gauss-jordan division
	 *
	 * @param {number[][]|Float32Array[]|Float64Array[]} matrix - gets modified
	 * @param {number} row
	 * @param {number} col
	 * @param {number} numCols
	 * @returns void
	 */
	Polyfit.gaussJordanDivide = function (matrix, row, col, numCols) {
		for (var i = col + 1; i < numCols; i++) {
			matrix[row][i] /= matrix[row][col];
		}
		matrix[row][col] = 1;
	};
	/**
	 * Perform gauss-jordan elimination
	 *
	 * @param {number[][]|Float64Array[]} matrix - gets modified
	 * @param {number} row
	 * @param {number} col
	 * @param {number} numRows
	 * @param {number} numCols
	 * @returns void
	 */
	Polyfit.gaussJordanEliminate = function (matrix, row, col, numRows, numCols) {
		for (var i = 0; i < numRows; i++) {
			if (i !== row && matrix[i][col] !== 0) {
				for (var j = col + 1; j < numCols; j++) {
					matrix[i][j] -= matrix[i][col] * matrix[row][j];
				}
				matrix[i][col] = 0;
			}
		}
	};
	/**
	 * Perform gauss-jordan echelon method
	 *
	 * @param {number[][]|Float32Array[]|Float64Array[]} matrix - gets modified
	 * @returns {number[][]|Float32Array[]|Float64Array[]} matrix
	 */
	Polyfit.gaussJordanEchelonize = function (matrix) {
		var rows = matrix.length;
		var cols = matrix[0].length;
		var i = 0;
		var j = 0;
		var k;
		var swap;
		while (i < rows && j < cols) {
			k = i;
			// Look for non-zero entries in col j at or below row i
			while (k < rows && matrix[k][j] === 0) {
				k++;
			}
			// If an entry is found at row k
			if (k < rows) {
				// If k is not i, then swap row i with row k
				if (k !== i) {
					swap = matrix[i];
					matrix[i] = matrix[k];
					matrix[k] = swap;
				}
				// If matrix[i][j] is != 1, divide row i by matrix[i][j]
				if (matrix[i][j] !== 1) {
					Polyfit.gaussJordanDivide(matrix, i, j, cols);
				}
				// Eliminate all other non-zero entries
				Polyfit.gaussJordanEliminate(matrix, i, j, rows, cols);
				i++;
			}
			j++;
		}
		return matrix;
	};
	/**
	 * Perform regression
	 *
	 * @param {number} x
	 * @param {number[]|Float32Array[]|Float64Array[]} terms
	 * @returns {number}
	 */
	Polyfit.regress = function (x, terms) {
		var a = 0;
		var exp = 0;
		for (var i = 0, len = terms.length; i < len; i++) {
			a += terms[i] * Math.pow(x, exp++);
		}
		return a;
	};
	/**
	 * Compute correlation coefficient
	 *
	 * @param {number[]|Float32Array[]|Float64Array[]} terms
	 * @returns {number}
	 */
	Polyfit.prototype.correlationCoefficient = function (terms) {
		var r = 0;
		var n = this.x.length;
		var sx = 0;
		var sx2 = 0;
		var sy = 0;
		var sy2 = 0;
		var sxy = 0;
		var x;
		var y;
		for (var i = 0; i < n; i++) {
			x = Polyfit.regress(this.x[i], terms);
			y = this.y[i];
			sx += x;
			sy += y;
			sxy += x * y;
			sx2 += x * x;
			sy2 += y * y;
		}
		var div = Math.sqrt((sx2 - (sx * sx) / n) * (sy2 - (sy * sy) / n));
		if (div !== 0) {
			r = Math.pow((sxy - (sx * sy) / n) / div, 2);
		}
		return r;
	};
	/**
	 * Run standard error function
	 *
	 * @param {number[]|Float32Array[]|Float64Array[]} terms
	 * @returns number
	 */
	Polyfit.prototype.standardError = function (terms) {
		var r = 0;
		var n = this.x.length;
		if (n > 2) {
			var a = 0;
			for (var i = 0; i < n; i++) {
				a += Math.pow(Polyfit.regress(this.x[i], terms) - this.y[i], 2);
			}
			r = Math.sqrt(a / (n - 2));
		}
		return r;
	};
	/**
	 * Compute coefficients for given data matrix
	 *
	 * @param {number} p
	 * @returns {number[]|Float32Array|Float64Array}
	 */
	Polyfit.prototype.computeCoefficients = function (p) {
		var n = this.x.length;
		var r;
		var c;
		var rs = 2 * ++p - 1;
		var i;
		var m = [];
		// Initialize array with 0 values
		if (this.FloatXArray) {
			// fast FloatXArray-Matrix init
			var bytesPerRow = (p + 1) * this.FloatXArray.BYTES_PER_ELEMENT;
			var buffer = new ArrayBuffer(p * bytesPerRow);
			for (i = 0; i < p; i++) {
				m[i] = new this.FloatXArray(buffer, i * bytesPerRow, p + 1);
			}
		} else {
			var zeroRow = [];
			for (i = 0; i <= p; i++) {
				zeroRow[i] = 0;
			}
			m[0] = zeroRow;
			for (i = 1; i < p; i++) {
				// copy zeroRow
				m[i] = zeroRow.slice();
			}
		}
		var mpc = [n];
		for (i = 1; i < rs; i++) {
			mpc[i] = 0;
		}
		for (i = 0; i < n; i++) {
			var x = this.x[i];
			var y = this.y[i];
			// Process precalculation array
			for (r = 1; r < rs; r++) {
				mpc[r] += Math.pow(x, r);
			}
			// Process RH column cells
			m[0][p] += y;
			for (r = 1; r < p; r++) {
				m[r][p] += Math.pow(x, r) * y;
			}
		}
		// Populate square matrix section
		for (r = 0; r < p; r++) {
			for (c = 0; c < p; c++) {
				m[r][c] = mpc[r + c];
			}
		}
		Polyfit.gaussJordanEchelonize(m);
		var terms = (this.FloatXArray && new this.FloatXArray(m.length)) || [];
		for (i = m.length - 1; i >= 0; i--) {
			terms[i] = m[i][p];
		}
		return terms;
	};
	/**
	 * Using given degree of fitment, return a function that will calculate
	 * the y for a given x
	 *
	 * @param {number} degree  > 0
	 * @returns {Function}     f(x) =
	 */
	Polyfit.prototype.getPolynomial = function (degree) {
		if (isNaN(degree) || degree < 0) {
			throw new Error("Degree must be a positive integer");
		}
		var terms = this.computeCoefficients(degree);
		var eqParts = [];
		eqParts.push(terms[0].toPrecision());
		for (var i = 1, len = terms.length; i < len; i++) {
			eqParts.push(terms[i] + " * Math.pow(x, " + i + ")");
		}
		var expr = "return " + eqParts.join(" + ") + ";";
		/* jshint evil: true */
		return new Function("x", expr);
		/* jshint evil: false */
	};
	/**
	 * Convert the polynomial to a string expression, mostly useful for visual
	 * debugging
	 *
	 * @param {number} degree
	 * @returns {string}
	 */
	Polyfit.prototype.toExpression = function (degree) {
		if (isNaN(degree) || degree < 0) {
			throw new Error("Degree must be a positive integer");
		}
		var terms = this.computeCoefficients(degree);
		var eqParts = [];
		var len = terms.length;
		eqParts.push(terms[0].toPrecision());
		for (var i = 1; i < len; i++) {
			eqParts.push(terms[i] + "x^" + i);
		}
		return eqParts.join(" + ");
	};

	return DFA;
})();

module.exports = DFA;
