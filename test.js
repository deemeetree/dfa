let DFA = require("./dfa");

const testLengths = [128, 200, 256, 512, 65536, 262144];

let randomArrays = [];

console.log("generating signals");

testLengths.forEach((length) => {
	randomArrays.push(
		Array.from({ length: length }, () => Math.floor(Math.random() * 1000) + 1)
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

const after = Date.now();

console.log(`Time taken: ${after - before}ms`);
