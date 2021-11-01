function Visualizer() {
    this.plot = document.getElementById('plot')
    this.width = this.plot.clientWidth
    this.height = this.plot.clientHeight

    this.gridSizeBox = document.getElementById('gridSizeBox')
    this.gridSizeBox.addEventListener('change', () => this.Draw())

    this.splitBox = document.getElementById('splitBox')
    this.splitBox.addEventListener('change', () => this.Draw())

    this.processorsBox = document.getElementById('processorsBox')
    this.processorsBox.addEventListener('change', () => this.Draw())

    this.neighbourBox = document.getElementById('neighbourBox')
    this.neighbourBox.addEventListener('change', () => this.Draw())

    this.drawNeighbourBox = document.getElementById('drawNeighbourBox')
    this.drawNeighbourBox.addEventListener('change', () => this.Draw())

    this.info = document.getElementById('info')
    this.Draw()
}

Visualizer.prototype.MakeParallelepipedMesh = function(p, isNeighbour = false) {
    if (isNeighbour) {
        if (p.xmin == p.xmax) {
            p.xmax += 0.1
            p.xmin -= 0.1
        }

        if (p.ymin == p.ymax) {
            p.ymax += 0.1
            p.ymin -= 0.1
        }

        if (p.zmin == p.zmax) {
            p.zmax += 0.1
            p.zmin -= 0.1
        }
    }

    let mesh = {
        type: "mesh3d",
        x: [p.xmin, p.xmin, p.xmax, p.xmax, p.xmin, p.xmin, p.xmax, p.xmax],
        y: [p.ymin, p.ymax, p.ymax, p.ymin, p.ymin, p.ymax, p.ymax, p.ymin],
        z: [p.zmin, p.zmin, p.zmin, p.zmin, p.zmax, p.zmax, p.zmax, p.zmax],
        i: [7, 0, 0, 0, 4, 4, 6, 1, 4, 0, 3, 6],
        j: [3, 4, 1, 2, 5, 6, 5, 2, 0, 1, 6, 3],
        k: [0, 7, 2, 3, 6, 7, 1, 6, 5, 5, 7, 2],
        text: p.text,
        opacity: 1,
        flatshading: true
    }

    if (isNeighbour)
        mesh.color = '#000'

    return mesh
}

Visualizer.prototype.MakeParallelepiped = function(xmin, xmax, ymin, ymax, zmin, zmax) {
    let volume = (xmax - xmin + 1) * (ymax - ymin + 1) * (zmax - zmin + 1)

    return {
        xmin: xmin, xmax: xmax,
        ymin: ymin, ymax: ymax,
        zmin: zmin, zmax: zmax,
        volume: volume,
        text: `[${xmin}, ${xmax}] x [${ymin}, ${ymax}] x [${zmin}, ${zmax}] = ${volume}`
    }
}

Visualizer.prototype.SplitTape = function(xmin, xmax, ymin, ymax, zmin, zmax, size, axis) {
    let splited = []

    let prev

    if (axis == 'x') {
        prev = xmin
    }
    else if (axis == 'y') {
        prev = ymin
    }
    else if (axis == 'z') {
        prev = zmin
    }

    for (let i = 0; i < size; i++) {
        if (axis == 'x') {
            let x = Math.min(xmax, xmin + Math.floor((xmax - xmin) * (i + 1) / size))
            splited.push(this.MakeParallelepiped(prev, x, ymin, ymax, zmin, zmax))
            prev = x + 1
        }
        else if (axis == 'y') {
            let y = Math.min(ymax, ymin + Math.floor((ymax - ymin) * (i + 1) / size))
            splited.push(this.MakeParallelepiped(xmin, xmax, prev, y, zmin, zmax))
            prev = y + 1
        }
        else if (axis == 'z') {
            let z = Math.min(zmax, zmin + Math.floor((zmax - zmin) * (i + 1) / size))
            splited.push(this.MakeParallelepiped(xmin, xmax, ymin, ymax, prev, z))
            prev = z + 1
        }
    }

    return splited
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

Visualizer.prototype.Decompose = function(n) {
    let dividers = []

    for (let i = 2; i <= n; i++) {
        while (n % i == 0) {
            dividers.push(i)
            n /= i
        }
    }

    while (dividers.length < 3)
        dividers.push(1)

    if (dividers.length > 3) {
        for (let i = 0; i < dividers.length && dividers.length > 3; i++) {
            let p = dividers[i] * dividers[i + 1]
            dividers.splice(i, 2, p)
        }

        dividers.sort()
    }

    return dividers
}

Visualizer.prototype.SplitProduct = function(xmin, xmax, ymin, ymax, zmin, zmax, size) {
    let dividers = this.Decompose(size)
    let px = dividers[0]
    let py = dividers[1]
    let pz = dividers[2]

    let dx = Math.max(1, Math.floor((xmax - xmin + 1) / px))
    let dy = Math.max(1, Math.floor((ymax - ymin + 1) / py))
    let dz = Math.max(1, Math.floor((zmax - zmin + 1) / pz))

    let volumes = []

    for (let i = 0; i < px; i++) {
        for (let j = 0; j < py; j++) {
            for (let k = 0; k < pz; k++) {
                let x1 = xmin + i * dx
                let x2 = i == px - 1 ? xmax : x1 + dx - 1

                let y1 = ymin + j * dy
                let y2 = j == py - 1 ? ymax : y1 + dy - 1

                let z1 = zmin + k * dz
                let z2 = k == pz - 1 ? zmax : z1 + dz - 1

                volumes.push(this.MakeParallelepiped(x1, x2, y1, y2, z1, z2))
            }
        }
    }

    return volumes
}

// 1 in 2
Visualizer.prototype.IsInside = function(xmin1, xmax1, ymin1, ymax1, xmin2, xmax2, ymin2, ymax2) {
    return xmin2 <= xmin1 && xmax1 <= xmax2 && ymin2 <= ymin1 && ymax1 <= ymax2
}

Visualizer.prototype.GetNeighbour = function(v1, v2) {
    if (v1.xmin == v2.xmax + 1 || v2.xmin == v1.xmax + 1) {
        let x = v1.xmin == v2.xmax + 1 ? v1.xmin : v1.xmax

        if (this.IsInside(v1.ymin, v1.ymax, v1.zmin, v1.zmax, v2.ymin, v2.ymax, v2.zmin, v2.zmax))
            return this.MakeParallelepiped(x, x, v1.ymin, v1.ymax, v1.zmin, v1.zmax)

        if (this.IsInside(v2.ymin, v2.ymax, v2.zmin, v2.zmax, v1.ymin, v1.ymax, v1.zmin, v1.zmax))
            return this.MakeParallelepiped(x, x, v2.ymin, v2.ymax, v2.zmin, v2.zmax)

        return null
    }

    if (v1.ymin == v2.ymax + 1 || v2.ymin == v1.ymax + 1) {
        let y = v1.ymin == v2.ymax + 1 ? v1.ymin : v1.ymax

        if (this.IsInside(v1.xmin, v1.xmax, v1.zmin, v1.zmax, v2.xmin, v2.xmax, v2.zmin, v2.zmax))
            return this.MakeParallelepiped(v1.xmin, v1.xmax, y, y, v1.zmin, v1.zmax)

        if (this.IsInside(v2.xmin, v2.xmax, v2.zmin, v2.zmax, v1.xmin, v1.xmax, v1.zmin, v1.zmax))
            return this.MakeParallelepiped(v2.xmin, v2.xmax, y, y, v2.zmin, v2.zmax)

        return null
    }

    if (v1.zmin == v2.zmax + 1 || v2.zmin == v1.zmax + 1) {
        let z = v1.zmin == v2.zmax + 1 ? v1.zmin : v1.zmax

        if (this.IsInside(v1.xmin, v1.xmax, v1.ymin, v1.ymax, v2.xmin, v2.xmax, v2.ymin, v2.ymax))
            return this.MakeParallelepiped(v1.xmin, v1.xmax, v1.ymin, v1.ymax, z, z)

        if (this.IsInside(v2.xmin, v2.xmax, v2.ymin, v2.ymax, v1.xmin, v1.xmax, v1.ymin, v1.ymax))
            return this.MakeParallelepiped(v2.xmin, v2.xmax, v2.ymin, v2.ymax, z, z)

        return null
    }

    return null
}

Visualizer.prototype.FindNeighbours = function(volumes) {
    let neighbours = []

    for (let i = 0; i < volumes.length; i++) {
        neighbours[i] = []

        for (let j = 0; j < volumes.length; j++) {
            if (i == j)
                continue

            let neighbour = this.GetNeighbour(volumes[i], volumes[j])

            if (neighbour == null)
                continue

            neighbours[i].push(neighbour)
        }
    }

    return neighbours
}

Visualizer.prototype.Draw = function() {
    let N = +this.gridSizeBox.value
    let P = +this.processorsBox.value
    let split = this.splitBox.value

    this.neighbourBox.max = P - 1

    if (+this.neighbourBox.value >= P)
        this.neighbourBox.value = P - 1

    let drawNeighbour = this.drawNeighbourBox.checked
    let neighbour = +this.neighbourBox.value

    this.info.innerHTML = '<b>Параметры:</b><br>'
    this.info.innerHTML += '<b>N</b>: ' + N + '<br>'
    this.info.innerHTML += '<b>P</b>: ' + P + '<br>'

    let volumes = []

    if (split == 'blocks') {
        volumes = this.Split(0, N, 0, N, 0, N, P, 'x')
    }
    else if (split == 'tapes') {
        volumes = this.SplitTape(0, N, 0, N, 0, N, P, 'x')
    }
    else if (split == 'product') {
        volumes = this.SplitProduct(0, N, 0, N, 0, N, P)
    }

    let neighbours = this.FindNeighbours(volumes)

    for (let i = 0; i < P; i++) {
        this.info.innerHTML += `${i}: ${volumes[i].text}<br>`
    }

    let layout = {
        width: this.width,
        height: this.height,
        autosize: true,
        margin: { l: 0, r: 0, b: 0, t: 0 }
    }

    let meshes = volumes.map((volume) => this.MakeParallelepipedMesh(volume))

    if (drawNeighbour) {
        this.neighbourBox.parentNode.style.display = ''
        meshes = meshes.concat(neighbours[neighbour].map((volume) => this.MakeParallelepipedMesh(volume, true)))
    }
    else {
        this.neighbourBox.parentNode.style.display = 'none'
    }

    Plotly.newPlot(this.plot, meshes, layout);
}
