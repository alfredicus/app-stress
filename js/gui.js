
function connectGUI() {
    kepler.getColorMapNames().forEach(name => {
        kepler.ColorMap.addColorMap(name, kepler.getColorMap(name, 40, false).colors)
    })
    const colorTables = kepler.colorMapNames()

    let attributes = []
    if (lineDataframe.length !== 0) {
        const manager = new dataframe.Manager(lineDataframe[0], [
            new math.PositionDecomposer,       // x y z
            new math.ComponentDecomposer,      // Ux Uy Uz Sxx Sxy Sz Syy Syz Szz
            new math.VectorNormDecomposer,     // U
            new math.EigenValuesDecomposer,    // S1 S2 S3
            new math.EigenVectorsDecomposer,   // S1 S2 S3
        ])
        attributes = manager.names(1)
    }

    // =========================================================

    const upload = document.getElementById('upload')
    upload.onchange = () => {
        upload.files[0].arrayBuffer().then(arrayBuffer => {
            console.log(new TextDecoder().decode(arrayBuffer))
        })
    }

    const GUI = lil.GUI
    const parent = document.getElementById('menu')
    const gui = new GUI({ container: parent })

    const myObject = {
        screenshot: function () {
            takeScreenshot()
        },
        showBBox: false,
        showSphere: true,
        colorTable: 'Insar',

        showLines: true,
        loadData: function () {
            document.getElementById('upload').click()
        },
        clearData: function () {

        },
        attribute: 'x',

        sigma1: 1,
        sigma2: 0,
        sigma3: 0,

        showData: true,
        dataWidth: 0.005,

        run: function () {
            changeSigma(1, 2, 3)
        }
    }
    const colorFormats = {
        string: '#ffffff',
        int: 0xffffff,
        object: { r: 1, g: 1, b: 1 },
        array: [1, 1, 1]
    }

    // ---------------------------------

    const general = gui.addFolder('General');
    general.addColor(colorFormats, 'string').name('Background').onChange(value => {
        extra.changeBackground({ scene, color: value })
    })
    general.add(myObject, 'showBBox').name('Bounding box').onChange(value => {
        bbox.visible = value
    })
    general.add(myObject, 'screenshot').name('Take screenshot') // Button
    general.close()

    // ---------------------------------

    const sphere = gui.addFolder('Sphere');
    sphere.add(myObject, 'showSphere').name('Visible').onChange(value => {
        surfs.visible = value
    })
    sphere.addColor(colorFormats, 'string').name('Color').onChange(value => {
        surfs.traverse(node => {
            if (node.name === 'sphere') {
                node.material.color = new three.Color(value)
            }
        })
    })
    sphere.close()

    // ---------------------------------

    const stressTensor = new stress.StressTensor({
        trendS1 : 0, 
        plungeS1: 0, 
        trendS3 : 0, 
        plungeS3: 0, 
        masterStress: 'σ1'
    })

    const tensor = {
        master      : 'σ1',
        trendMaster : 0,
        plungeMaster: 0,
        trendSlave  : 0,
        plungeSlave : 0
    }
    const tensorFolder = gui.addFolder('Stress tensor orientation')

    function updateStressTensor() {
        tensor.plungeSlave = stressTensor.plunge
    }

    tensorFolder.add(tensor, 'master', ['σ1', 'σ3']).name('Master principal stress').onChange( value => {
        tm.name(`    Trend ${value}`)
        pm.name(`    Plunge ${value}`)
        if (value === 'σ1') {
            ts.name(`    Trend σ3`)
            ps.name(`    Plunge σ3`)
        }
        else {
            ts.name(`    Trend σ1`)
            ps.name(`    Plunge σ1`)
        }
        stressTensor.changeMasterStress(value)
        updateStressTensor()
    })
    const tm = tensorFolder.add(tensor, 'trendMaster', 0, 360, 1).name('    Trend σ1').onChange( v => {
        if (tensor.master === 'σ1') {
            stressTensor.trendS1 = v
        }
        else {
            stressTensor.trendS3 = v
        }
        updateStressTensor()
    })
    const pm = tensorFolder.add(tensor, 'plungeMaster', -90, 90, 1).name('    Plunge σ1').onChange( v => {
        if (tensor.master === 'σ1') {
            stressTensor.plungeS1 = v
        }
        else {
            stressTensor.plungeS3 = v
        }
        updateStressTensor()
    })
    const ts = tensorFolder.add(tensor, 'trendSlave', 0, 360, 1).name('    Trend σ3').onChange( v => {
        if (tensor.master === 'σ1') {
            stressTensor.trendS3 = v
        }
        else {
            stressTensor.trendS1 = v
        }
        updateStressTensor()
    })
    const ps = tensorFolder.add(tensor, 'plungeSlave').name('    Plunge σ3').disable().listen()

    // ---------------------------------

    let integralBuilder      = new stress.IntegralCurve     ([-myObject.sigma1, -myObject.sigma3, -myObject.sigma2], 1.001)
    let equipotentialBuilder = new stress.EquipotentialCurve([-myObject.sigma1, -myObject.sigma3, -myObject.sigma2], 1.001)
    let mohrCoulombBuilder   = new stress.MohrCoulombCurve  ([-myObject.sigma1, -myObject.sigma3, -myObject.sigma2], 1.001)

    const ic = {
        visible  : true,
        lineWidth: 0.005,
        theta    : 0,
        phi      : 0
    }

    const mc = {
        visible  : true,
        lineWidth: 0.005,
        cohesion : 0,
        frictionAngle : 0
    }

    gui.add(myObject, 'sigma2', 0, 1, 0.01).name('σ2').onChange( value => {
        integralBuilder      = new stress.IntegralCurve     ([-myObject.sigma1, -myObject.sigma3, -value], 1.001)
        equipotentialBuilder = new stress.EquipotentialCurve([-myObject.sigma1, -myObject.sigma3, -value], 1.001)
        mohrCoulombBuilder   = new stress.MohrCoulombCurve  ([-myObject.sigma1, -myObject.sigma3, -value], 1.001)
        rebuildIntegralCurves     ()
        rebuildEquipotentialCurves()
        rebuildMohrCoulombCurves  ()
        redrawMohrCircle          ()
    })

    function rebuildIntegralCurves() {
        integrals.clear()
        if (ic.visible) {
            const buffer = integralBuilder.generate(ic.theta, ic.phi)
            const lines = createLineFromPl(buffer, color = '#0000ff', ic.lineWidth)
            lines.forEach( line => integrals.add(line) )
        }
    }

    function rebuildEquipotentialCurves() {
        equipotentials.clear()
        if (ic.visible) {
            const buffer = equipotentialBuilder.generate(ic.theta, ic.phi)
            const lines = createLineFromPl(buffer, color = '#0000ff', ic.lineWidth)
            lines.forEach( line => equipotentials.add(line) )
        }
    }

    function rebuildMohrCoulombCurves() {
        mohrCoulombs.clear()
        if (mc.visible) {
            const buffer = mohrCoulombBuilder.generate(mc.frictionAngle, mc.cohesion)
            const lines = createLineFromPl(buffer, color = '#FF0000', mc.lineWidth)
            lines.forEach( line => mohrCoulombs.add(line) )
        }
    }

    function redrawMohrCircle() {
        mohrCircle({ element: "mohr", width: 150, height: 150, S1: myObject.sigma1, S2: myObject.sigma2, S3: myObject.sigma3, scale: 50 })
        console.log(myObject.sigma1, myObject.sigma2, myObject.sigma3)
    }

    // ---------------------------------
 
    const guiC = gui.addFolder('Integral curves')

    guiC.add(ic, "visible").onChange( value => {
        rebuildIntegralCurves     ()
        rebuildEquipotentialCurves()
    })
    guiC.add(ic, 'lineWidth', 0, 0.01, 0.0001).name('Width').onChange( value => {
        ic.lineWidth = value
        rebuildIntegralCurves     ()
        rebuildEquipotentialCurves()
    })

    guiC.add(ic, 'theta', 0, 90, 1).name('Theta').onChange( value => {
        ic.theta = value
        rebuildIntegralCurves     ()
        rebuildEquipotentialCurves()
    })

    guiC.add(ic, 'phi', 0, 90, 1).name('Phi').onChange( value => {
        ic.phi = value
        rebuildIntegralCurves     ()
        rebuildEquipotentialCurves()
    })

    // ---------------------------------

    const guiMC = gui.addFolder('Mohr Coulomb')
    
    guiMC.add(mc, "visible").onChange( value => {
        rebuildMohrCoulombCurves()
    })
    guiMC.add(mc, 'frictionAngle', 0, 45, 1).name('Friction angle').onChange( value => {
        mc.frictionAngle = value
        const maxCohe = mohrCoulombBuilder.maxCohesion(mc.frictionAngle)
        if (mc.cohesion > maxCohe) {
            mc.cohesion = maxCohe
        }
        rebuildMohrCoulombCurves()
    }).listen()

    guiMC.add(mc, 'cohesion', 0, 0.5, 0.01).name('Cohesion').onChange( value => {
        mc.cohesion = value
        const maxFric = mohrCoulombBuilder.maxFrictionAngle(mc.cohesion)
        if (mc.frictionAngle > maxFric) {
            mc.frictionAngle = maxFric
        }
        rebuildMohrCoulombCurves()
    }).listen()

    // ---------------------------------

    const data = gui.addFolder('Data')
    data.add(myObject, 'loadData').name('Load') // Button
    data.add(myObject, 'clearData').name('clear') // Button
    data.add(myObject, 'showData').name('Visible')
    data.addColor(colorFormats, 'string').name('Color')
    data.add(myObject, 'dataWidth', 0, 0.01, 0.0001).name('Width')
    data.close()

    // ---------------------------------

    const run = gui.addFolder('Simulation')
    run.add(myObject, 'run').name('run')
    run.close()
}
