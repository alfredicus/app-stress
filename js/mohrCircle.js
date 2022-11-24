const PIO2 = Math.PI / 2
const O2 = 1 / 2

function mohrCircle({ element, width, height, S1, S2, S3, scale = 100 }) {
    const transformX = x => x * scale + tx
    const transformY = y => height + y * scale - ty
    const makeCircle = (x, y, r, fill = 'none', color = 'black') => {
        svg.append("path")
            .attr("transform", `translate(${transformX(x)}, ${transformY(y)})`)
            .attr("d", d3.arc()
                .innerRadius(0)
                .outerRadius(r * scale)
                .startAngle(-PIO2)
                .endAngle(PIO2)
            )
            .attr('fill', fill) // or black
            .attr('stroke-width', '1')
            .attr('stroke', color)
    }
    const drawSigma = (x, id) => svg.append("text")
        .attr('x', transformX(x / m - 0.1))
        .attr('y', transformY(0.15))
        .attr("dy", ".35em")
        .text(`σ${id}`)

    // let scale = 120
    let tx = width / 2
    let ty = 30

    let c1 = O2 * (S2 + S3)
    let c2 = O2 * (S1 + S3)
    let c3 = O2 * (S1 + S2)
    let r1 = O2 * (S2 - S3)
    let r2 = O2 * (S1 - S3)
    let r3 = O2 * (S1 - S2)

    // Normalizing everything between 0 and 1
    const m = Math.max(r1, r2, r3)
    if (m === 0) m === 1
    r1 /= m
    r2 /= m
    r3 /= m
    c1 /= m
    c2 /= m
    c3 /= m

    tx = width / 2 - c2 * scale - 10

    const svg = d3.select(`#${element}`)
        .html(null) // clear the element
        .append("svg")
        .attr("width", width)
        .attr("height", height)

    makeCircle(c2, 0, r2, 'white')
    makeCircle(c1, 0, r1, 'gray')
    makeCircle(c3, 0, r3, 'gray')

    drawSigma(S1, 1)
    drawSigma(S2, 2)
    drawSigma(S3, 3)

    // Axis X
    // svg.append("line")
    //     .attr("x1", transformX(0))
    //     .attr("x2", transformX(1.75))
    //     .attr("y1", transformY(0))
    //     .attr("y2", transformY(0))
    //     .attr("stroke", "black")
    // svg.append("text")
    //     .attr('x', transformX(1.75))
    //     .attr('y', transformY(0))
    //     .attr("dy", ".35em")
    //     .text('σ')

    // Axis Y
    svg.append("line")
        .attr("x1", transformX(0))
        .attr("x2", transformX(0))
        .attr("y1", transformY(0))
        .attr("y2", transformY(-1.2))
        .attr("stroke", "black")
    svg.append("text")
        .attr('x', transformX(- 0.1))
        .attr('y', transformY(-1.2))
        .attr("dy", ".35em")
        .text('τ')
}

