stress = globalThis['@alfredo-taboada/stress']


function test() {
    const f = new stress.Fault({strike: 45, dipDirection: stress.Direction.SE, dip: 60})
    f.setStriation({rake: 30, strikeDirection: stress.Direction.NE, sensMouv: stress.SensOfMovement.I})

    console.log(f)
}
