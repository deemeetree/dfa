let DFA = require("./dfa");

const testLengths = [128, 200, 256, 512, 65536, 262144];

let randomArrays = [];

// console.log("generating signals");

// testLengths.forEach((length) => {
// 	randomArrays.push(
// 		Array.from({ length: length }, () => Math.floor(Math.random() * 1000) + 1)
// 	);
// });

randomArrays.push([
	959.0, 1156.2, 1084.0, 1041.0, 1095.7, 1093.8, 1113.3, 1091.8, 933.6, 1121.1,
	1039.1, 1056.6, 1103.5, 1080.1, 1058.6, 1072.3, 1031.2, 906.2, 1046.9, 1015.6,
	1037.1, 1078.1, 1054.7, 1095.7, 1064.5, 1029.3, 1009.8,
]);
const before = Date.now();

randomArrays.forEach((data) => {
	let dfa = new DFA(data);
	let alpha_component = dfa.compute();
	console.log(alpha_component);
	console.log(alpha_component.alpha);
	console.log(alpha_component.RMSSD);
	console.log(alpha_component.lengthOfData);
});

const after = Date.now();

console.log(`Time taken: ${after - before}ms`);
