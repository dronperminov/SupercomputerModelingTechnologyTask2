function Visualizer() {
    this.plot = document.getElementById('plot')
    this.width = this.plot.clientWidth
    this.height = this.plot.clientHeight

    this.input = document.getElementById('input')
    this.inputSpan = document.getElementById('input-span')
    this.input.addEventListener('change', () => this.ChangeFile())

    this.step = document.getElementById('step')
    this.step.addEventListener('change', () => this.Draw())

    this.info = document.getElementById('info')
}

Visualizer.prototype.InitWebGL = function() {
    let gl = null

    try {
        gl = this.plot.getContext("webgl") || this.plot.getContext("experimental-webgl")
    }
    catch (e) {

    }

    if (!gl) {
        alert("WebGL не поддерживается Вашим браузером...")
        gl = null
    }

    return gl
}

Visualizer.prototype.ChangeFile = function() {
    this.inputSpan.innerHTML = this.input.files[0].name

    let fr = new FileReader()
    fr.onload = (e) => {
        this.data = JSON.parse(e.target.result)
        this.step.parentNode.style.display = ''
        this.Draw()
    }

    fr.readAsText(this.input.files[0])
    this.input.files = null
}

Visualizer.prototype.Plot = function(x, y, z, u, u_min, u_max, size) {
    let trace = {
        x: x,
        y: y,
        z: z,
        text: u.map((v) => 'u: ' + v),
        mode: 'markers',
        marker: {
            color: u.map((ui) => `hsl(${240 - (ui - u_min) / (u_max - u_min) * 240}, 50, 70)`),
            size: size,
            line: {
            color: 'rgba(127, 127, 127, 0.5)',
            width: 0.5},
            opacity: 1
        },
        type: 'scatter3d'
    }

    let layout = {
        width: this.width,
        height: this.height,
        autosize: true,
        margin: { l: 0, r: 0, b: 0, t: 0 }
    }

    Plotly.newPlot(this.plot, [trace], layout)
}

Visualizer.prototype.Draw = function() {
    let N = this.data['N']
    let points = this.data['points']

    let x = []
    let y = []
    let z = []
    let u = []

    let u_min = Infinity
    let u_max = -Infinity

    let index = 0;
    let step = +this.step.value

    for (let i = 0; i <= N; i++) {
        for (let j = 0; j <= N; j++) {
            for (let k = 0; k <= N; k++) {
                let point = points[index++]

                if ((i % step == 0 || i == N) && (j % step == 0 || j == N) && (k % step == 0 || k == N)) {
                    x.push(point[0])
                    y.push(point[1])
                    z.push(point[2])
                    u.push(point[3])
                }

                u_min = Math.min(u_min, point[3])
                u_max = Math.max(u_max, point[3])
            }
        }
    }

    this.info.innerHTML = '<b>Параметры решения:</b><br>'
    this.info.innerHTML += '<b>L<sub>x</sub></b>: ' + this.data['Lx'] + '<br>'
    this.info.innerHTML += '<b>L<sub>y</sub></b>: ' + this.data['Ly'] + '<br>'
    this.info.innerHTML += '<b>L<sub>z</sub></b>: ' + this.data['Lz'] + '<br>'
    this.info.innerHTML += '<b>N</b>: ' + N + '<br>'
    this.info.innerHTML += '<b>t</b>: ' + this.data['t'] + '<br>'
    this.info.innerHTML += '<b>Число точек</b>: ' + u.length + '<br>'

    this.Plot(x, y, z, u, u_min, u_max, 3*step)
}