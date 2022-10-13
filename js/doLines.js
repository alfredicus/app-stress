
function doLines(plines) {
    const promises = []

    const promise = doPLine(plines)
    if (promise) {
        if (Array.isArray(promise)) {
            promises.push(...promise)
        }
        else {
            promises.push(promise)
        }
    }

    return promises
}

function doPLine(plineInfo) {
    if (plineInfo.url === undefined) return []
    if (Array.isArray(plineInfo.url) && plineInfo.url.length===0) return []
    if (Array.isArray(plineInfo.url)) {
        const promises = []
        plineInfo.url.forEach( url => {
            const promise = doOnePline(url)
            if (promise) promises.push(promise)
        })
        return promises
    }
    else {
        return doOnePline(plineInfo.url)
    }

    function doOnePline(url) {
        const promise = fetch(url)
            .then( res => {
                if ( res.ok ) return res.text()
                return undefined
            })
            .then( buffer => {
                if (! buffer) return undefined

                const filter = io.IOFactory.getFilter(url)
                if (filter) {
                    const dfs = filter.decode(buffer, {shared: false, merge: true})
                    dfs.forEach( df => {
                        lineDataframe.push(df)
                        createGlLine(df, plineInfo)

                        console.log(math.minMax(df.series.positions))
                    })
                }
            })
        return promise
    }
}

function updateLines() {
    lines.clear()
    lineDataframe.forEach( df => createGlLine(df, plines) )
}

function createGlLine(df, plineInfo) {
    const manager = new dataframe.Manager(df, [
        new math.PositionDecomposer,       // x y z
        new math.ComponentDecomposer,      // Ux Uy Uz Sxx Sxy Sz Syy Syz Szz
        new math.VectorNormDecomposer,     // U
        new math.EigenValuesDecomposer,    // S1 S2 S3
        new math.EigenVectorsDecomposer,   // S1 S2 S3
    ])

    let skin = kepler.createLineset2({
        position: df.series.positions,
        parameters: {
            width  : plineInfo.width,
            color  : plineInfo.color,
            opacity: plineInfo.opacity
        }
    })
    lines.add( skin )

    if (plineInfo.attr) {
        const attrName = plineInfo.attr
        const attr = manager.serie(1, attrName)
        if (attr) {
            kepler.paintAttribute(skin, attr, new kepler.PaintParameters({
                atVertex  : true,
                lut       : plineInfo.lut!==undefined?plineInfo.lut: 'insar',
                reverseLut: plineInfo.reverseLut!==undefined?plineInfo.reverseLut: false
            }))
        }
    }
}