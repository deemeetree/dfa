let DFA = require("./dfa");

const testLengths = [128, 200, 256, 512, 65536, 262144];

let randomArrays = [];

console.log("generating signals");

// testLengths.forEach((length) => {
// 	randomArrays.push(
// 		Array.from({ length: length }, () => Math.floor(Math.random() * 1000) + 1)
// 	);
// });

randomArrays.push([
	3, 3, 3, 3, 3, 3, 3, 3, 1, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 3, 2, 2, 2, 1, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 2, 3, 3, 3, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2,
	3, 3, 2, 2, 4, 4, 4, 3, 3, 3, 4, 3, 3, 4, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 2, 2,
	3, 2, 3, 3, 3, 3, 3, 3, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
	3, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
]);
const before = Date.now();

randomArrays.forEach((data) => {
	console.log(JSON.stringify(data));
	let dfa = new DFA(data);
	let alpha_component = dfa.compute();
	console.log(
		`N=${alpha_component.lengthOfData}, scales:`,
		alpha_component.scales
	);
	console.log(
		`Alpha: ${alpha_component.alpha}, Alpha1: ${alpha_component.alpha1}, Alpha2: ${alpha_component.alpha2}`
	);
	console.log(`scalesAlpha1:`, alpha_component.scalesAlpha1);
	console.log(`scalesAlpha2:`, alpha_component.scalesAlpha2);
	console.log("---");
});

const after = Date.now();

console.log(`Time taken: ${after - before}ms`);

// ── Multifractal DFA (MFDFA) ──
console.log("\n=== Multifractal DFA ===");

// Longer correlated signal so all three h(q) curves are available.
let acc = 0;
let seed = 42;
const mfData = Array.from({ length: 4000 }, () => {
	seed = (seed * 1103515245 + 12345) & 0x7fffffff;
	acc += seed / 0x7fffffff - 0.5;
	return acc;
});

const mfDfa = new DFA(mfData);
const mono = mfDfa.compute();
const mf = mfDfa.computeMultifractal({ qMin: -5, qMax: 5, qStep: 1 });

console.log(`N=${mf.lengthOfData}, q values:`, mf.q);
console.log(
	"global h(q):",
	mf.hq.map((v) => (v === null ? "null" : v.toFixed(3)))
);
console.log(
	"alpha1 h(q):",
	mf.hq1.map((v) => (v === null ? "null" : v.toFixed(3)))
);
console.log(
	"alpha2 h(q):",
	mf.hq2.map((v) => (v === null ? "null" : v.toFixed(3)))
);
console.log(
	`widths — global: ${mf.multifractalWidth.toFixed(
		3
	)}, alpha1: ${mf.width1.toFixed(3)}, alpha2: ${mf.width2.toFixed(3)}`
);

// q=2 must reproduce the monofractal compute() result exactly.
const close = (a, b) => Math.abs(a - b) < 1e-9;
const ok =
	close(mono.alpha, mf.monofractal.alpha) &&
	close(mono.alpha1, mf.monofractal.alpha1) &&
	close(mono.alpha2, mf.monofractal.alpha2);
console.log(
	`q=2 sanity check (matches compute alpha/alpha1/alpha2): ${
		ok ? "PASS" : "FAIL"
	}`
);
