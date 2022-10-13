function postInit() {
    // extra.changeBackground( {scene, color: '#888888'} )
    if (model.backgroundColor) {
        extra.changeBackground( {scene, color: model.backgroundColor} )
    } else {
        extra.changeBackground( {scene, color: '#fff'} )
    }
    
    addLights()

    const keyboard = new extra.Keyboard(document, 'keydown')
    // keyboard.setUpEvent(e => {controls.constraint = extra.CONSTRAINT.NONE})

    keyboard.addKey({key:'u', cb:e => extra.changeView('up', {scene, camera, controls, selection: sphereFake}) })
    keyboard.addKey({key:'d', cb:e => extra.changeView('down', {scene, camera, controls, selection: sphereFake}) })
    keyboard.addKey({key:'s', cb:e => extra.changeView('south', {scene, camera, controls, selection: sphereFake}) })
    keyboard.addKey({key:'n', cb:e => extra.changeView('north', {scene, camera, controls, selection: sphereFake}) })
    keyboard.addKey({key:'e', cb:e => extra.changeView('east', {scene, camera, controls, selection: sphereFake}) })
    keyboard.addKey({key:'w', cb:e => extra.changeView('west', {scene, camera, controls, selection: sphereFake}) })

    keyboard.addKey({key:' ', cb: e => {
        if (cube) cube.restoreView()
    } })
    
    keyboard.addKey({key:'f', cb:e => extra.zoomToModel({scene, camera, controls, duration:300}) })

    // keyboard.addKey({key:'x', cb:e => controls.constraint = extra.CONSTRAINT.X})
    // keyboard.addKey({key:'y', cb:e => controls.constraint = extra.CONSTRAINT.Y})
    // keyboard.addKey({key:'z', cb:e => controls.constraint = extra.CONSTRAINT.Z})

    renderer.domElement.addEventListener('dblclick', event => {
        extra.zoomToIntersection( {scene, event, camera, controls} )
    })

    const LOADER = document.getElementById('js-loader')
    if (LOADER) LOADER.remove()
}
