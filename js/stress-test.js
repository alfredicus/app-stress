stress = globalThis['@alfredo-taboada/stress']

const f = new stress.Fault({strike: 45, dipDirection: 'SE', dip: 60})
f.setStriation({rake: 30, strikeDirection: 'NE', sensMouv: stress.SensOfMovement.I})

console.log(f)