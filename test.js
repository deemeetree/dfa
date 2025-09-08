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
