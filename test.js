let DFA = require("./dfa");

let GenerateFractal = require("generate-fractal-signal");

const signalGenerator = GenerateFractal;
const signalConfig = {
	signalLength: 1024, //4096, 1024, 256, 128
	minWindow: 4,
	scaleGrow: 0.5,
	signalRange: [1, 4],
	signalType: "Fractal",
};

let DFA_NEW = require("./dfa-new");

const testLengths = [128, 200, 256, 512, 65536]; // 65536, 262144

let randomArrays = [];

let fractalArrays = [];

console.log("generating signals");
testLengths.forEach((length) => {
	randomArrays.push(
		Array.from({ length: length }, () => Math.floor(Math.random() * 1000) + 1)
	);
	fractalArrays.push(
		signalGenerator.generateSignal({ ...signalConfig, signalLength: length })
			.timeSeries
	);
});

const before = Date.now();

randomArrays.forEach((data) => {
	let dfa = new DFA(data);
	let alpha_component = dfa.compute();
	console.log(alpha_component.alpha);
	console.log(alpha_component.RMSSD);
	console.log(alpha_component.lengthOfData);
});

fractalArrays.forEach((data) => {
	let dfa = new DFA(data);
	let alpha_component = dfa.compute();
	console.log(alpha_component.alpha);
	console.log(alpha_component.RMSSD);
	console.log(alpha_component.lengthOfData);
});

const after = Date.now();

console.log(`Time taken for the old algo: ${after - before}ms`);

const before_new = Date.now();

randomArrays.forEach((data) => {
	let dfa = new DFA_NEW(data);
	let alpha_component = dfa.compute();
	console.log(alpha_component.alpha);
	console.log(alpha_component.RMSSD);
	console.log(alpha_component.lengthOfData);
});

fractalArrays.forEach((data) => {
	let dfa = new DFA_NEW(data);
	let alpha_component = dfa.compute();
	console.log(alpha_component.alpha);
	console.log(alpha_component.RMSSD);
	console.log(alpha_component.lengthOfData);
});
const after_new = Date.now();

console.log(`Time taken for the old algo: ${after_new - before_new}ms`);
