class LinesBuilder {
    constructor() {
        this.clear()
    }
    clear() {
        this.s       = ''
        this.id      = 0
        this.beginId = 0
        this.endId   = 0
    }
    startLine() {
        this.id = 0
        this.s += `GOCAD PLine 1.0
        HEADER {
            name: circle
        }
        `
    }
    add(x,y,z) {
        this.s += `VRTX ${this.id} ${x} ${y} ${z}\n`
        this.id++
    }
    endLine() {
        for (let i=0; i<this.id-1; ++i) {
            this.s += `SEG ${i} ${i+1}\n`
        }
        this.s += 'END\n'
    }
    get buffer() {
        return this.s
    }
}

function generateSphere() {
    const geometry = new three.SphereGeometry( 1, 64, 32 )
    const material = new THREE.MeshPhongMaterial( { color: 0xaaaaaa } )
    sphereFake = new THREE.Mesh( geometry, material )
    const sphere = new THREE.Mesh( geometry, material )
    sphere.name = "sphere"
    surfs.add(sphere)

    // --------------------------------------

    const r = 1.001
    const l = new LinesBuilder()

    l.startLine()
    for (let i=0; i<360; i++) {
        const a = i*Math.PI/180
        l.add(r*Math.cos(a), r*Math.sin(a), 0)
    }
    l.endLine()

    l.startLine()
    for (let i=0; i<360; i++) {
        const a = i*Math.PI/180
        l.add(r*Math.cos(a), 0, r*Math.sin(a))
    }
    l.endLine()

    l.startLine()
    for (let i=0; i<360; i++) {
        const a = i*Math.PI/180
        l.add(0, r*Math.cos(a), r*Math.sin(a))
    }
    l.endLine()

    createLineFromPl(l.buffer, '#000000', 0.002).forEach( skin => sphere.add(skin) )
    
    arrow([0,0,1], [0,0,1], [0,0,1.6], sphere, 0x0000ff, 'Up') // up
    arrow([1,0,0], [1,0,0], [1.6,0,0], sphere, 0xff0000, 'East') // east
    arrow([0,1,0], [0,1,0], [0,1.6,0], sphere, 0x00ff00, 'North') // north
}

function arrow(d, o, pos, parent, color, name) {
    const dir = new three.Vector3( d[0], d[1], d[2] )
    dir.normalize()
    const origin = new three.Vector3( o[0], o[1], o[2] );
    const length = 0.5
    const hex = color
    const arrowHelper = new three.ArrowHelper( dir, origin, length, hex, length*0.2, length*0.3 )
    parent.add( arrowHelper )

    addObjectLabel(parent, name, pos[0], pos[1], pos[2])
}
