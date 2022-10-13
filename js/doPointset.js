
function doPointsets(pointsets) {
    const promises = []

    if (pointsets) {
        pointsets.forEach( pointset => {
            if (pointset.show) {
                const promise = doPointset(pointset)
                if (promise) {
                    if (Array.isArray(promise)) {
                        promises.push(...promise)
                    }
                    else {
                        promises.push(promise)
                    }
                }
            }
        })
    }

    return promises
}

// ------------------------------------------

function doPointset(psetInfo) {
    if (psetInfo.url === undefined || psetInfo.show === false) return []
    if (Array.isArray(psetInfo.url)) {
        const promises = []
        psetInfo.url.forEach( url => {
            const promise = doOnePset(url)
            if (promise) promises.push(promise)
        })
        return promises
    }
    else {
        return doOnePset(psetInfo.url)
    }

    function doOnePset(url) {
        const promise = fetch(url)
            .then( res => {
                if ( res.ok ) return res.text()
                return undefined
            })
            .then( buffer => {
                if (! buffer) return undefined

                // const scolor = randColor()

                const filter = io.IOFactory.getFilter(url)
                if (filter === undefined) {
                    return undefined
                }

                const dfs = filter.decode(buffer, { shared: false, merge: true })
                
                dfs.forEach( df => {

                    let position = df.series['positions']
                    console.log('min-max position pointset:', math.minMax(position) )

                    if (psetInfo.translation) {
                        const x = psetInfo.translation[0]
                        const y = psetInfo.translation[1]
                        const z = psetInfo.translation[2]
                        position = position.map(v => [v[0] + x, v[1] + y, v[2] + z])
                    }

                    const manager = new dataframe.Manager(df, [
                        new math.PositionDecomposer,
                        new math.ComponentDecomposer,
                        new math.VectorNormDecomposer,
                        new math.EigenValuesDecomposer
                    ])

                    let attr = manager.serie(1, psetInfo.attr)

                    const SKIN = kepler.createPointset({
                        position: position,
                        parameters: new kepler.PointsetParameters({
                            size: psetInfo.size,
                            color: psetInfo.color!==undefined?psetInfo.color:undefined,
                            sizeAttenuation: psetInfo.sizeAttenuation!==undefined?psetInfo.sizeAttenuation:true
                        })
                    })

                    if (psetInfo.show) {
                        group.add( SKIN )
                    }

                    if (psetInfo.attr !== undefined && attr && psetInfo.show) {
                        console.log('min-max attr '+psetInfo.attr+' pointset:', math.minMax(attr) )

                        kepler.paintAttribute( SKIN, attr, new kepler.PaintParameters({
                            atVertex: true,
                            lut: psetInfo.lut,
                            duplicateLut: psetInfo.duplicateLut,
                            reverseLut: psetInfo.reverseLut,
                            min: psetInfo.useMinMax ? psetInfo.min : 0,
                            max: psetInfo.useMinMax ? psetInfo.max : 1
                        }) )
                    }
                    
                    if (psetInfo.vectors !== undefined && psetInfo.vectors.show === true) {
                        const vattr = manager.serie(3, psetInfo.vectors.attr)
                        if (vattr) {
                            if (psetInfo.vectors.useTube) {
                                group.add( kepler.createTubeVectors({
                                    geometry: SKIN.geometry,
                                    vectorField: vattr,
                                    // attribute: attr,
                                    parameters: new kepler.TubeVectorsParameters({
                                        scale: psetInfo.vectors.scale,
                                        color: psetInfo.vectors.color,
                                        radius: psetInfo.vectors.radius
                                    })
                                }) )
                            }
                            else {
                                group.add( kepler.createVectors({
                                    geometry: SKIN.geometry,
                                    vectorField: vattr,
                                    parameters: new kepler.TubeVectorsParameters({
                                        scale: psetInfo.vectors.scale,
                                        color: psetInfo.vectors.color
                                    })
                                }) )
                            }
                        }
                    }

                    if (psetInfo.failure !== undefined && psetInfo.failure.show === true) {
                        if (df.series[psetInfo.failure.stress] !== undefined) {
                            let skin = kepler.createFailurePlanes({
                                geometry: SKIN.geometry,
                                dataframe: df,
                                parameters: new kepler.FailurePlanesParameters({
                                    stress: psetInfo.failure.stress,
                                    size: psetInfo.failure.size,
                                    sizeAttribute: psetInfo.failure.sizeAttribute?psetInfo.failure.sizeAttribute:'',
                                    paintAttribute: psetInfo.failure.paintAttribute?psetInfo.failure.paintAttribute:'',
                                    color: psetInfo.failure.color,
                                    circle: psetInfo.failure.circle,
                                    borders: psetInfo.failure.borders,
                                    type: psetInfo.failure.type
                                })
                            })
                            group.add(skin)

                            if (psetInfo.failure.borders) {
                                skin = kepler.createEdges(skin.geometry, new kepler.EdgesParameters({
                                    thresholdAngle: 10,
                                    color: psetInfo.failure.borderColor?psetInfo.failure.borderColor:"#000000"
                                }))
                                group.add(skin)
                            }
                        }
                    }

                })
            })
        return promise
    }
}
