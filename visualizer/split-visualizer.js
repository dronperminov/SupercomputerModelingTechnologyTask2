function Visualizer() {
    this.plot = document.getElementById('plot')
    this.width = this.plot.clientWidth
    this.height = this.plot.clientHeight

    this.gridSizeBox = document.getElementById('gridSizeBox')
    this.gridSizeBox.addEventListener('change', () => this.Draw())

    this.processorsBox = document.getElementById('processorsBox')
    this.processorsBox.addEventListener('change', () => this.Draw())

    this.info = document.getElementById('info')
    this.Draw()
}

Visualizer.prototype.MakeParallelepiped = function(xmin, xmax, ymin, ymax, zmin, zmax) {
    return {
        type: "mesh3d",
        x: [xmin, xmin, xmax, xmax, xmin, xmin, xmax, xmax],
        y: [ymin, ymax, ymax, ymin, ymin, ymax, ymax, ymin],
        z: [zmin, zmin, zmin, zmin, zmax, zmax, zmax, zmax],
        i: [7, 0, 0, 0, 4, 4, 6, 1, 4, 0, 3, 6],
        j: [3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3],
        k: [0, 7, 2, 3, 6, 7, 1, 6, 5, 5, 7, 2],
        text:  `[${xmin}, ${xmax}] x [${ymin}, ${ymax}] x [${zmin}, ${zmax}] = ${(xmax - xmin) * (ymax - ymin) * (zmax - zmin)}`,
        opacity: 1,
        flatshading: true
    }
}

Visualizer.prototype.Split = function(xmin, xmax, ymin, ymax, zmin, zmax, size, axis) {
    if (size == 1)
        return [this.MakeParallelepiped(xmin, xmax, ymin, ymax, zmin, zmax)]

    let splited = []

    if (size % 2 == 1) {
        if (axis == 'x') {
            let x = xmin + Math.floor((xmax - xmin) / size)
            splited = [this.MakeParallelepiped(xmin, x, ymin, ymax, zmin, zmax)]
            xmin = x + 1
            axis = 'y'
        }
        else if (axis == 'y') {
            let y = ymin + Math.floor((ymax - ymin) / size)
            splited = [this.MakeParallelepiped(xmin, xmax, ymin, y, zmin, zmax)]
            ymin = y + 1
            axis = 'z'
        }
        else {
            let z = zmin + Math.floor((zmax - zmin) / size)
            splited = [this.MakeParallelepiped(xmin, xmax, ymin, ymax, zmin, z)]
            zmin = z + 1
            axis = 'x'
        }

        size--
    }

    if (axis == 'x') {
        let x = Math.floor((xmin + xmax) / 2)
        splited = splited.concat(this.Split(xmin, x, ymin, ymax, zmin, zmax, size / 2, 'y'))
        splited = splited.concat(this.Split(x + 1, xmax, ymin, ymax, zmin, zmax, size / 2, 'y'))
    }
    else if (axis == 'y') {
        let y = Math.floor((ymin + ymax) / 2)
        splited = splited.concat(this.Split(xmin, xmax, ymin, y, zmin, zmax, size / 2, 'z'))
        splited = splited.concat(this.Split(xmin, xmax, y + 1, ymax, zmin, zmax, size / 2, 'z'))
    }
    else {
        let z = Math.floor((zmin + zmax) / 2)
        splited = splited.concat(this.Split(xmin, xmax, ymin, ymax, zmin, z, size / 2, 'x'))
        splited = splited.concat(this.Split(xmin, xmax, ymin, ymax, z + 1, zmax, size / 2, 'x'))
    }

    return splited
}

Visualizer.prototype.Draw = function() {
    let N = +this.gridSizeBox.value
    let P = +this.processorsBox.value

    this.info.innerHTML = '<b>Параметры:</b><br>'
    this.info.innerHTML += '<b>N</b>: ' + N + '<br>'
    this.info.innerHTML += '<b>P</b>: ' + P + '<br>'

    let meshes = this.Split(0, N, 0, N, 0, N, P, 'x')

    for (mesh of meshes) {
        let xmin = mesh.x[0]
        let xmax = mesh.x[2]

        let ymin = mesh.y[0]
        let ymax = mesh.y[1]

        let zmin = mesh.z[0]
        let zmax = mesh.z[4]

        this.info.innerHTML += `[${xmin}, ${xmax}] x [${ymin}, ${ymax}] x [${zmin}, ${zmax}]<br>`
    }

    let layout = {
        width: this.width,
        height: this.height,
        autosize: true,
        margin: { l: 0, r: 0, b: 0, t: 0 }
    }

    Plotly.newPlot(this.plot, meshes, layout);
}
