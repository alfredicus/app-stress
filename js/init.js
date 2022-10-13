function init() {
    const three = globalThis['THREE']
    const extra = globalThis['@youwol/three-extra']
    
    scene = new three.Scene

    sphereFake = undefined
    
    /* Organization of the scene
        scene
            group
                geodesics
                surfs (sphere)
                    lines
            lights
            bbox (model)
    */
    group = new three.Group

    surfDataframe = []
    surfs = new three.Group
    group.add(surfs)

    lineDataframe = []
    lines = new three.Group
    group.add(lines)

    scene.add(group)

    bbox   = undefined
    lights = undefined

    camera = new three.PerspectiveCamera( 30, window.innerWidth / window.innerHeight, 0.01, 100000 )
    camera.position.z = 100

    renderer = new three.WebGLRenderer({alpha: true})
    renderer.setPixelRatio( window.devicePixelRatio )
    renderer.setSize( window.innerWidth, window.innerHeight )
    renderer.shadowMap.enabled = true
    document.body.appendChild( renderer.domElement )

    renderFct = new extra.RenderFunctions({renderer, scene, camera})

    window.addEventListener( 'resize', onWindowResize )

    controls = new three.TrackballControls( camera, renderer.domElement )
    controls.rotateSpeed = 3.0
    // controls.zoomSpeed = 1.2
    // controls.panSpeed = 0.8
    renderFct.add( controls.update )

    cube = new extra.installNavigationCube(
        new extra.NavigationCubeParameters({
            scene, 
            camera, 
            renderer,
            controls, 
            renderFunctions: renderFct, // will also add the cube in renderFct
            labels: ['Right', 'Left', 'Up', 'Down', 'Front', 'Back'],
            // labels: ['East', 'West', 'Up', 'Down', 'South', 'North'],
            // labels: ['Y', '-Y', 'Z', '-Z', 'X', '-X'],

            domElement : document.getElementById('orientCubeWrapper'),
            domHome    : document.getElementById('goHome'), 
            domSaveHome: document.getElementById('saveHome')
        })
    )
}

