function Visualizer() {
    this.plot = document.getElementById('plot')
    this.width = this.plot.clientWidth
    this.height = this.plot.clientHeight

    this.input = document.getElementById('input')
    this.inputSpan = document.getElementById('input-span')
    this.input.addEventListener('change', () => this.ChangeFile())

    this.stepBox = document.getElementById('step')
    this.stepBox.addEventListener('change', () => this.Draw())

    this.opacityBox = document.getElementById('opacity')
    this.opacityBox.addEventListener('change', () => this.Draw())

    this.info = document.getElementById('info')
}

Visualizer.prototype.ChangeFile = function() {
    this.inputSpan.innerHTML = this.input.files[0].name

    let fr = new FileReader()
    fr.onload = (e) => {
        this.data = JSON.parse(e.target.result)
        this.stepBox.parentNode.style.display = ''
        this.opacityBox.parentNode.style.display = ''
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
        text: u.map((ui) => 'u: ' + ui),
        mode: 'markers',
        marker: {
            color: u.map((ui) => `hsl(${((u_max - ui) / (u_max - u_min)) * 240}, 50, 70)`),
            size: size,
            opacity: +this.opacityBox.value
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
    let points = this.data['u']

    let x = []
    let y = []
    let z = []
    let u = []

    let u_min = Infinity
    let u_max = -Infinity

    let index = 0;
    let step = +this.stepBox.value

    let hx = this.data['Lx'] / N
    let hy = this.data['Ly'] / N
    let hz = this.data['Lz'] / N

    for (let i = 0; i <= N; i++) {
        for (let j = 0; j <= N; j++) {
            for (let k = 0; k <= N; k++) {
                let ui = points[index++]

                if ((i % step == 0 || i == N) && (j % step == 0 || j == N) && (k % step == 0 || k == N)) {
                    x.push(i * hx)
                    y.push(j * hy)
                    z.push(k * hz)
                    u.push(ui)
                }

                u_min = Math.min(u_min, ui)
                u_max = Math.max(u_max, ui)
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

    this.Plot(x, y, z, u, u_min, u_max, 4*step)
}