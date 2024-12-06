stress = globalThis['@alfredo-taboada/stress']


function test() {
    const f = new stress.Fault({strike: 180, dipDirection: stress.Direction.S, dip: 90})
    f.setStriation({rake: 30, strikeDirection: stress.Direction.S, typeMov: stress.TypeOfMovement.RL})

    console.log(f)
}
