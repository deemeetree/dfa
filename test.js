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

// ── Multifractal verdicts ───────────────────────────────────────────────────
// Three reference signals: white noise (monofractal), a binomial multiplicative
// cascade (multifractal through its ordering) and the same cascade shuffled
// (same values, ordering destroyed — only the distribution-driven width stays).
function seededRng(seed) {
	let a = seed >>> 0;
	return () => {
		a = (a + 0x6d2b79f5) >>> 0;
		let t = a;
		t = Math.imul(t ^ (t >>> 15), t | 1);
		t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
		return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
	};
}
const rng = seededRng(7);
const gauss = () => {
	const u = rng() || 1e-12;
	return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * rng());
};
const cascade = (levels, a = 0.75) => {
	let s = [1];
	for (let i = 0; i < levels; i++) {
		const next = [];
		for (const v of s) next.push(v * a, v * (1 - a));
		s = next;
	}
	return s;
};
const white = Array.from({ length: 1024 }, gauss);
const casc = cascade(10);
const shuffledCasc = casc.slice();
for (let i = shuffledCasc.length - 1; i > 0; i--) {
	const j = Math.floor(rng() * (i + 1));
	[shuffledCasc[i], shuffledCasc[j]] = [shuffledCasc[j], shuffledCasc[i]];
}

const verdictOpts = { qMin: -5, qMax: 5, qStep: 0.5 };
const vWhite = new DFA(white).computeMultifractal(verdictOpts).verdict.global;
const vCasc = new DFA(casc).computeMultifractal(verdictOpts).verdict.global;
const vShuf = new DFA(shuffledCasc).computeMultifractal(verdictOpts).verdict
	.global;
const show = (v) =>
	`${v.label} (source: ${v.source}, z: ${v.z.toFixed(2)}, p: ${v.p.toFixed(
		3
	)}, width: ${v.width.toFixed(3)}, surrogate: ${v.surrogateWidth.toFixed(3)})`;
console.log("white noise:      ", show(vWhite));
console.log("cascade:          ", show(vCasc));
console.log("shuffled cascade: ", show(vShuf));

const verdictsOk =
	vWhite.label === "monofractal" &&
	vWhite.significant === false &&
	vCasc.label === "multifractal" &&
	vCasc.significant === true &&
	vShuf.significant === false &&
	vShuf.source === "distribution";
console.log(`verdict check: ${verdictsOk ? "PASS" : "FAIL"}`);

// Same seed → same surrogates → same verdict.
const again = new DFA(casc).computeMultifractal(verdictOpts).verdict.global;
console.log(
	`verdict determinism: ${
		JSON.stringify(again) === JSON.stringify(vCasc) ? "PASS" : "FAIL"
	}`
);

// surrogateCount: 0 falls back to width thresholds only.
const noSurr = new DFA(casc).computeMultifractal({
	...verdictOpts,
	surrogateCount: 0,
}).verdict.global;
console.log(
	`no-surrogate fallback: ${
		noSurr.label === "multifractal" && noSurr.significant === null
			? "PASS"
			: "FAIL"
	}`
);

if (!ok || !verdictsOk) process.exitCode = 1;

// Global alpha cap (2.5.0): default fits up to N/4, globalMaxFraction = 1 up to N
const capSeries = Array.from({ length: 512 }, (_, i) => Math.sin(i * 0.37) + ((i * 7919) % 101) / 101);
const capDefault = new DFA(capSeries).compute();
const capFull = new DFA(capSeries).compute(4, 0.25, 2, 16, 16, 0.25, "moderate", 1);
console.log(
	`global cap: default range ${JSON.stringify(capDefault.globalRange)} ${
		capDefault.globalRange[1] <= 512 / 4 ? "PASS" : "FAIL"
	}, full range ${JSON.stringify(capFull.globalRange)} ${capFull.globalRange[1] === 512 ? "PASS" : "FAIL"}, alpha1 equal ${
		capDefault.alpha1 === capFull.alpha1 ? "PASS" : "FAIL"
	}, alpha2 equal ${capDefault.alpha2 === capFull.alpha2 ? "PASS" : "FAIL"}`
);
