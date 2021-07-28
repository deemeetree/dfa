let DFA = require('./dfa')

let time_series = [8, 10, 6, 9, 7, 5, 5, 11, 11, 8, 6, 7, 9, 10, 7, 9]
let dfa = new DFA(time_series)
let alpha_component = dfa.compute()
console.log(alpha_component)