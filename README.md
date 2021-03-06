# SupercomputerModelingTechnologyTask2
Реализация второго практического задания в рамках курса "Суперкомпьютерное моделирование и технологии"

## Программная реализация

Реализовано три программы: последовательная и две гибридных параллельных (MPI + OpenMP и MPI + CUDA). Каждая из программ является консольным приложением и принимает входные данные в виде аргументов командной строки. Используются следующие аргументы (доступно при вызове с единственным аргументом `--help`):

* `-d` – отладочный режим (по умолчанию не используется)
* `-Lx` – длина параллелепипеда вдоль оси X (по умолчанию 1)
* `-Ly` – длина параллелепипеда вдоль оси Y (по умолчанию 1)
* `-Lz` – длина параллелепипеда вдоль оси Z (по умолчанию 1)
* `-T` – конечное время сетки (по умолчанию 1)
* `-N` – количество точек пространственной сетки (по умолчанию 40)
* `-K` – количество точек временной сетки (по умолчанию 100)
* `-steps` – количество шагов для решения (по умолчанию 20)
* `-btx` – тип граничного условия вдоль оси X (по умолчанию однородные первого рода)
* `-bty` – тип граничного условия вдоль оси Y (по умолчанию периодические-численные)
* `-btz` – тип граничного условия вдоль оси Z (по умолчанию однородные первого рода)
* `-s` – стратегия разбиения на блоки (по умолчанию блочное, доступно только для MPI версии)
* `-on` – путь к json файлу для сохранения численного решения (по умолчанию не используется)
* `-oa` – путь к json файлу для сохранения аналитического решения (по умолчанию не используется)
* `-od` – путь к json файлу для сохранения погрешности (по умолчанию не используется)
* `-o` – путь к текстовому файлу для вывода результатов (по умолчанию output.txt)

### Типы граничных условий:
* `first-kind (f)` – однородные граничные условия первого рода
* `periodic-analytical (pa)` – аналитические периодические граничные условия
* `periodic-numerical (pn)` – численные периодические граничные условия

### Типы стратегий разбиения
* `blocks (b)` – блочное разбиение
* `tapes (t)` – ленточное разбиение
* `products (p)` – разбиение путём разложения на множители

Для визуализации получаемых решений написан визуализатор (принимает json файлы, генерируемые программой): https://programforyou.ru/tests/supercomputer-modeling-technology/solve-visualizer

<img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/solve.gif">

## Сборка

Сборки осуществляется через makefile. Для сборки последовательной версии выполните `make main`. Для сборки параллельной версии выполните `make main-mpi`.

### Пример запуска
```bash
>make main
>./main -d -N 64 -K 150 -btx f -bty pa -btz pn -Lx 0.5 -Ly 12.5 -Lz 100
```


## Особенности параллельной реализации

Для распараллеливания вся сетка разбивается на области (также прямоугольные параллелепипеды) в количестве используемых процессов. Операции, производимые с получаемыми параллелепипедами, обрабатываются с помощью OpenMP в первой гибридной версии и с помощью CUDA ядер во второй гибридной версиях.

### Блочное разбиение:
* начнём разбиение с параллелепипеда [0, N] x [0, N] x [0, N], выберем начальную ось (X) и запустим рекурсивный процесс
* если текущее количество областей (size) равно 1, вернём обрабатываемый параллелепипед
* если размер нечётный, то по текущей оси выберем область 1 / size и сделаем из неё параллелепипед, и продолжим разбивать область 1 - 1 / size
* по выбранной оси делим область пополам и рекурсивно запускаем для этих подобластей

### Ленточное разбиение:
* выбираем ось, вдоль которой будет производиться разбиение
* начинаем делить на параллелепипеды по max⁡(1,(x_max-x_min)/size)  элементов

### Разбиение разложением на множители
* Раскладываем число на простые сомножители
* Если получилось менее трёх сомножителей, добавляем единицы до трёх
* В противном случае, пока сомножителей больше 3, заменяем пары чисел их произведением, а затем сортируем последовательность.

Для большей наглядности реализован визуализатор разбиения с возможностью установки размера сетки и числа процессов: https://programforyou.ru/tests/supercomputer-modeling-technology/split-visualizer

<img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/split.gif">

## Примеры разбиения

<details>
  <summary>Показать/скрыть</summary>

### Блочное разбиение
<table>
    <tr>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/blocks_split2.png"><br>P=2</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/blocks_split3.png"><br>P=3</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/blocks_split4.png"><br>P=4</td>
    </tr>
    <tr>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/blocks_split5.png"><br>P=5</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/blocks_split10.png"><br>P=10</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/blocks_split14.png"><br>P=14</td>
    </tr>
    <tr>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/blocks_split25.png"><br>P=25</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/blocks_split30.png"><br>P=30</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/blocks_split60.png"><br>P=60</td>
    </tr>
</table>

### Ленточное разбиение

<table>
    <tr>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/tapes_split2.png"><br>P=2</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/tapes_split3.png"><br>P=3</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/tapes_split4.png"><br>P=4</td>
    </tr>
    <tr>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/tapes_split8.png"><br>P=8</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/tapes_split14.png"><br>P=14</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/tapes_split17.png"><br>P=17</td>
    </tr>
</table>

### Разбиение разложением на множители

<table>
    <tr>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/product_split6.png"><br>P=6</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/product_split11.png"><br>P=11</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/product_split12.png"><br>P=12</td>
    </tr>
    <tr>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/product_split18.png"><br>P=18</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/product_split40.png"><br>P=40</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/product_split54.png"><br>P=54</td>
    </tr>
</table>

</details>

## Графики аналитического и полученного решений

<details>
  <summary>Показать/скрыть</summary>

### Lx = Ly = Lz = 1
<table>
    <tr>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/analytical_1_n128.png"><br>аналитическое</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/numerical_1_n128.png"><br>полученное</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/difference_1_n128.png"><br>погрешность</td>
    </tr>
</table>

### Lx = Ly = Lz = π
<table>
    <tr>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/analytical_pi_n128.png"><br>аналитическое</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/numerical_pi_n128.png"><br>полученное</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/difference_pi_n128.png"><br>погрешность</td>
    </tr>
</table>

</details>

## Polus (MPI версия)

<details>
  <summary>Показать/скрыть таблицы</summary>

### Lx = Ly = Lz = 1, N = 128, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th><th colspan="3">Разбиение разложением на множители</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>3.37891</td><td>1.000</td><td>0.00541829</td><td>3.38371</td><td>1.000</td><td>0.00541829</td><td>3.54794</td><td>1.000</td><td>0.00541829</td></tr>
    <tr align="center"><td>2</td><td>1.97175</td><td>1.714</td><td>0.00541829</td><td>1.98044</td><td>1.709</td><td>0.00541829</td><td>2.06306</td><td>1.720</td><td>0.00541829</td></tr>
    <tr align="center"><td>4</td><td>1.13682</td><td>2.972</td><td>0.00541829</td><td>1.13015</td><td>2.994</td><td>0.00541829</td><td>1.13537</td><td>3.125</td><td>0.00541829</td></tr>
    <tr align="center"><td>8</td><td>0.761273</td><td>4.438</td><td>0.00541829</td><td>0.715661</td><td>4.728</td><td>0.00541829</td><td>0.698494</td><td>5.079</td><td>0.00541829</td></tr>
    <tr align="center"><td>10</td><td>0.697569</td><td>4.844</td><td>0.00541829</td><td>0.573838</td><td>5.897</td><td>0.00541829</td><td>0.551368</td><td>6.435</td><td>0.00541829</td></tr>
    <tr align="center"><td>16</td><td>0.430247</td><td>7.853</td><td>0.00541829</td><td>0.409231</td><td>8.268</td><td>0.00541829</td><td>0.468169</td><td>7.578</td><td>0.00541829</td></tr>
    <tr align="center"><td>20</td><td>0.389413</td><td>8.677</td><td>0.00541829</td><td>0.36037</td><td>9.390</td><td>0.00541829</td><td>0.373345</td><td>9.503</td><td>0.00541829</td></tr>
    <tr align="center"><td>32</td><td>0.302416</td><td>11.173</td><td>0.00541829</td><td>0.259181</td><td>13.055</td><td>0.00541829</td><td>0.320764</td><td>11.061</td><td>0.00541829</td></tr>
    <tr align="center"><td>40</td><td>0.24562</td><td>13.757</td><td>0.00541829</td><td>0.250944</td><td>13.484</td><td>0.00541829</td><td>0.257525</td><td>13.777</td><td>0.00541829</td></tr>
    <tr align="center"><td>64</td><td>0.175378</td><td>19.266</td><td>0.00541829</td><td>0.182393</td><td>18.552</td><td>0.00541829</td><td>0.171637</td><td>20.671</td><td>0.00541829</td></tr>
</table>

### Lx = Ly = Lz = 1, N = 256, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th><th colspan="3">Разбиение разложением на множители</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>26.8418</td><td>1.000</td><td>0.00542011</td><td>26.8246</td><td>1.000</td><td>0.00542011</td><td>28.2069</td><td>1.000</td><td>0.00542011</td></tr>
    <tr align="center"><td>2</td><td>15.5534</td><td>1.726</td><td>0.00542011</td><td>15.6095</td><td>1.718</td><td>0.00542011</td><td>16.7504</td><td>1.684</td><td>0.00542011</td></tr>
    <tr align="center"><td>4</td><td>8.5587</td><td>3.136</td><td>0.00542011</td><td>8.58564</td><td>3.124</td><td>0.00542011</td><td>8.9247</td><td>3.161</td><td>0.00542011</td></tr>
    <tr align="center"><td>8</td><td>5.74982</td><td>4.668</td><td>0.00542011</td><td>5.16543</td><td>5.193</td><td>0.00542011</td><td>5.3995</td><td>5.224</td><td>0.00542011</td></tr>
    <tr align="center"><td>10</td><td>5.12554</td><td>5.237</td><td>0.00542011</td><td>4.26427</td><td>6.291</td><td>0.00542011</td><td>3.84234</td><td>7.341</td><td>0.00542011</td></tr>
    <tr align="center"><td>16</td><td>3.6738</td><td>7.306</td><td>0.00542011</td><td>2.85054</td><td>9.410</td><td>0.00542011</td><td>3.41574</td><td>8.258</td><td>0.00542011</td></tr>
    <tr align="center"><td>20</td><td>3.17586</td><td>8.452</td><td>0.00542011</td><td>2.39508</td><td>11.200</td><td>0.00542011</td><td>2.91717</td><td>9.669</td><td>0.00542011</td></tr>
    <tr align="center"><td>32</td><td>2.13905</td><td>12.548</td><td>0.00542011</td><td>1.74538</td><td>15.369</td><td>0.00542011</td><td>2.34652</td><td>12.021</td><td>0.00542011</td></tr>
    <tr align="center"><td>40</td><td>1.84101</td><td>14.580</td><td>0.00542011</td><td>1.55653</td><td>17.234</td><td>0.00542011</td><td>1.75374</td><td>16.084</td><td>0.00542011</td></tr>
    <tr align="center"><td>64</td><td>1.23209</td><td>21.786</td><td>0.00542011</td><td>1.12446</td><td>23.856</td><td>0.00542011</td><td>1.27153</td><td>22.183</td><td>0.00542011</td></tr>
</table>


### Lx = Ly = Lz = 1, N = 512, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th><th colspan="3">Разбиение разложением на множители</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>212.339</td><td>1.000</td><td>0.00542056</td><td>212.174</td><td>1.000</td><td>0.00542056</td><td>225.694</td><td>1.000</td><td>0.00542056</td></tr>
    <tr align="center"><td>2</td><td>119.572</td><td>1.776</td><td>0.00542056</td><td>120.541</td><td>1.760</td><td>0.00542056</td><td>130.791</td><td>1.726</td><td>0.00542056</td></tr>
    <tr align="center"><td>4</td><td>66.7043</td><td>3.183</td><td>0.00542056</td><td>66.9953</td><td>3.167</td><td>0.00542056</td><td>70.579</td><td>3.198</td><td>0.00542056</td></tr>
    <tr align="center"><td>8</td><td>38.7001</td><td>5.487</td><td>0.00542056</td><td>35.9597</td><td>5.900</td><td>0.00542056</td><td>41.4528</td><td>5.445</td><td>0.00542056</td></tr>
    <tr align="center"><td>10</td><td>38.6523</td><td>5.494</td><td>0.00542056</td><td>31.6109</td><td>6.712</td><td>0.00542056</td><td>30.7711</td><td>7.335</td><td>0.00542056</td></tr>
    <tr align="center"><td>16</td><td>22.3261</td><td>9.511</td><td>0.00542056</td><td>20.8502</td><td>10.176</td><td>0.00542056</td><td>25.9581</td><td>8.695</td><td>0.00542056</td></tr>
    <tr align="center"><td>20</td><td>23.1022</td><td>9.191</td><td>0.00542056</td><td>19.5098</td><td>10.875</td><td>0.00542056</td><td>21.8811</td><td>10.315</td><td>0.00542056</td></tr>
    <tr align="center"><td>32</td><td>14.4698</td><td>14.675</td><td>0.00542056</td><td>12.6782</td><td>16.735</td><td>0.00542056</td><td>17.111</td><td>13.190</td><td>0.00542056</td></tr>
    <tr align="center"><td>40</td><td>13.1529</td><td>16.144</td><td>0.00542056</td><td>10.8413</td><td>19.571</td><td>0.00542056</td><td>11.5291</td><td>19.576</td><td>0.00542056</td></tr>
    <tr align="center"><td>64</td><td>8.52932</td><td>24.895</td><td>0.00542056</td><td>7.82263</td><td>27.123</td><td>0.00542056</td><td>9.88421</td><td>22.834</td><td>0.00542056</td></tr>
</table>


### Lx = Ly = Lz = π, N = 128, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th><th colspan="3">Разбиение разложением на множители</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>3.37747</td><td>1.000</td><td>0.000549676</td><td>3.38126</td><td>1.000</td><td>0.000549676</td><td>3.68838</td><td>1.000</td><td>0.000549676</td></tr>
    <tr align="center"><td>2</td><td>1.99226</td><td>1.695</td><td>0.000549676</td><td>1.98962</td><td>1.699</td><td>0.000549676</td><td>2.06472</td><td>1.786</td><td>0.000549676</td></tr>
    <tr align="center"><td>4</td><td>1.12936</td><td>2.991</td><td>0.000549676</td><td>1.13578</td><td>2.977</td><td>0.000549676</td><td>1.26399</td><td>2.918</td><td>0.000549676</td></tr>
    <tr align="center"><td>8</td><td>0.847514</td><td>3.985</td><td>0.000549676</td><td>0.579493</td><td>5.835</td><td>0.000549676</td><td>0.920891</td><td>4.005</td><td>0.000549676</td></tr>
    <tr align="center"><td>10</td><td>0.673376</td><td>5.016</td><td>0.000549676</td><td>0.552498</td><td>6.120</td><td>0.000549676</td><td>0.541634</td><td>6.810</td><td>0.000549676</td></tr>
    <tr align="center"><td>16</td><td>0.46852</td><td>7.209</td><td>0.000549676</td><td>0.402093</td><td>8.409</td><td>0.000549676</td><td>0.453177</td><td>8.139</td><td>0.000549676</td></tr>
    <tr align="center"><td>20</td><td>0.421829</td><td>8.007</td><td>0.000549676</td><td>0.378469</td><td>8.934</td><td>0.000549676</td><td>0.356322</td><td>10.351</td><td>0.000549676</td></tr>
    <tr align="center"><td>32</td><td>0.29494</td><td>11.451</td><td>0.000549676</td><td>0.257307</td><td>13.141</td><td>0.000549676</td><td>0.309398</td><td>11.921</td><td>0.000549676</td></tr>
    <tr align="center"><td>40</td><td>0.278028</td><td>12.148</td><td>0.000549676</td><td>0.227738</td><td>14.847</td><td>0.000549676</td><td>0.222918</td><td>16.546</td><td>0.000549676</td></tr>
    <tr align="center"><td>64</td><td>0.172591</td><td>19.569</td><td>0.000549676</td><td>0.195346</td><td>17.309</td><td>0.000549676</td><td>0.175563</td><td>21.009</td><td>0.000549676</td></tr>
</table>


### Lx = Ly = Lz = π, N = 256, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th><th colspan="3">Разбиение разложением на множители</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>26.8546</td><td>1.000</td><td>0.000549861</td><td>26.8581</td><td>1.000</td><td>0.000549861</td><td>28.103</td><td>1.000</td><td>0.000549861</td></tr>
    <tr align="center"><td>2</td><td>15.2069</td><td>1.766</td><td>0.000549861</td><td>14.7859</td><td>1.816</td><td>0.000549861</td><td>15.4402</td><td>1.820</td><td>0.000549861</td></tr>
    <tr align="center"><td>4</td><td>8.37291</td><td>3.207</td><td>0.000549861</td><td>9.6132</td><td>2.794</td><td>0.000549861</td><td>8.70444</td><td>3.229</td><td>0.000549861</td></tr>
    <tr align="center"><td>8</td><td>4.72224</td><td>5.687</td><td>0.000549861</td><td>4.50557</td><td>5.961</td><td>0.000549861</td><td>4.87293</td><td>5.767</td><td>0.000549861</td></tr>
    <tr align="center"><td>10</td><td>4.64202</td><td>5.785</td><td>0.000549861</td><td>4.07747</td><td>6.587</td><td>0.000549861</td><td>4.57985</td><td>6.136</td><td>0.000549861</td></tr>
    <tr align="center"><td>16</td><td>3.2318</td><td>8.309</td><td>0.000549861</td><td>2.55629</td><td>10.507</td><td>0.000549861</td><td>3.93083</td><td>7.149</td><td>0.000549861</td></tr>
    <tr align="center"><td>20</td><td>2.85524</td><td>9.405</td><td>0.000549861</td><td>2.35552</td><td>11.402</td><td>0.000549861</td><td>2.80902</td><td>10.005</td><td>0.000549861</td></tr>
    <tr align="center"><td>32</td><td>2.20191</td><td>12.196</td><td>0.000549861</td><td>1.76651</td><td>15.204</td><td>0.000549861</td><td>1.95025</td><td>14.410</td><td>0.000549861</td></tr>
    <tr align="center"><td>40</td><td>2.04728</td><td>13.117</td><td>0.000549861</td><td>1.57969</td><td>17.002</td><td>0.000549861</td><td>1.86461</td><td>15.072</td><td>0.000549861</td></tr>
    <tr align="center"><td>64</td><td>1.2482</td><td>21.515</td><td>0.000549861</td><td>1.11757</td><td>24.033</td><td>0.000549861</td><td>1.24787</td><td>22.521</td><td>0.000549861</td></tr>
</table>


### Lx = Ly = Lz = π, N = 512, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th><th colspan="3">Разбиение разложением на множители</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>212.134</td><td>1.000</td><td>0.000549907</td><td>212.326</td><td>1.000</td><td>0.000549907</td><td>224.574</td><td>1.000</td><td>0.000549907</td></tr>
    <tr align="center"><td>2</td><td>119.237</td><td>1.779</td><td>0.000549907</td><td>119.366</td><td>1.779</td><td>0.000549907</td><td>127.467</td><td>1.762</td><td>0.000549907</td></tr>
    <tr align="center"><td>4</td><td>66.5768</td><td>3.186</td><td>0.000549907</td><td>67.0074</td><td>3.169</td><td>0.000549907</td><td>69.6428</td><td>3.225</td><td>0.000549907</td></tr>
    <tr align="center"><td>8</td><td>40.2159</td><td>5.275</td><td>0.000549907</td><td>34.0599</td><td>6.234</td><td>0.000549907</td><td>40.21</td><td>5.585</td><td>0.000549907</td></tr>
    <tr align="center"><td>10</td><td>36.5986</td><td>5.796</td><td>0.000549907</td><td>27.4475</td><td>7.736</td><td>0.000549907</td><td>30.8729</td><td>7.274</td><td>0.000549907</td></tr>
    <tr align="center"><td>16</td><td>23.8981</td><td>8.877</td><td>0.000549907</td><td>21.9308</td><td>9.682</td><td>0.000549907</td><td>26.087</td><td>8.609</td><td>0.000549907</td></tr>
    <tr align="center"><td>20</td><td>23.9314</td><td>8.864</td><td>0.000549907</td><td>17.8668</td><td>11.884</td><td>0.000549907</td><td>21.9503</td><td>10.231</td><td>0.000549907</td></tr>
    <tr align="center"><td>32</td><td>14.2082</td><td>14.930</td><td>0.000549907</td><td>13.0011</td><td>16.331</td><td>0.000549907</td><td>15.1717</td><td>14.802</td><td>0.000549907</td></tr>
    <tr align="center"><td>40</td><td>13.68</td><td>15.507</td><td>0.000549907</td><td>10.1392</td><td>20.941</td><td>0.000549907</td><td>12.3448</td><td>18.192</td><td>0.000549907</td></tr>
    <tr align="center"><td>64</td><td>8.81245</td><td>24.072</td><td>0.000549907</td><td>7.60838</td><td>27.907</td><td>0.000549907</td><td>8.89526</td><td>25.246</td><td>0.000549907</td></tr>
</table>
</details>  

## Polus (MPI + CUDA версия)

<details>
  <summary>Показать/скрыть таблицы</summary>

### Lx = Ly = Lz = 1, N = 128, K = 2000 (без использования CUDA ядер)

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th><th colspan="3">Разбиение разложением на множители</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>10.8421</td><td>1.000</td><td>0.00541829</td><td>10.828</td><td>1.000</td><td>0.00541829</td><td>10.8417</td><td>1.000</td><td>0.00541829</td></tr>
    <tr align="center"><td>2</td><td>5.56431</td><td>1.949</td><td>0.00541829</td><td>5.54688</td><td>1.952</td><td>0.00541829</td><td>5.61019</td><td>1.933</td><td>0.00541829</td></tr>
    <tr align="center"><td>4</td><td>2.80274</td><td>3.868</td><td>0.00541829</td><td>2.87798</td><td>3.762</td><td>0.00541829</td><td>2.96835</td><td>3.652</td><td>0.00541829</td></tr>
    <tr align="center"><td>6</td><td>1.95402</td><td>5.549</td><td>0.00541829</td><td>2.03811</td><td>5.313</td><td>0.00541829</td><td>1.97781</td><td>5.482</td><td>0.00541829</td></tr>
</table>

### Lx = Ly = Lz = 1, N = 128, K = 2000 (с использования CUDA ядер)

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th><th colspan="3">Разбиение разложением на множители</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>0.625809</td><td>1.000</td><td>0.00541829</td><td>0.417379</td><td>1.000</td><td>0.00541829</td><td>0.41405</td><td>1.000</td><td>0.00541829</td></tr>
    <tr align="center"><td>2</td><td>0.491011</td><td>1.275</td><td>0.00541829</td><td>0.395455</td><td>1.055</td><td>0.00541829</td><td>0.398651</td><td>1.039</td><td>0.00541829</td></tr>
    <tr align="center"><td>4</td><td>0.235194</td><td>2.661</td><td>0.00541829</td><td>0.235013</td><td>1.776</td><td>0.00541829</td><td>0.234564</td><td>1.765</td><td>0.00541829</td></tr>
    <tr align="center"><td>6</td><td>0.23938</td><td>2.614</td><td>0.00541829</td><td>0.254023</td><td>1.643</td><td>0.00541829</td><td>0.193185</td><td>2.143</td><td>0.00541829</td></tr>
</table>

### Lx = Ly = Lz = 1, N = 256, K = 2000 (без использования CUDA ядер)

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th><th colspan="3">Разбиение разложением на множители</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>83.0349</td><td>1.000</td><td>0.00542011</td><td>83.024</td><td>1.000</td><td>0.00542011</td><td>83.0116</td><td>1.000</td><td>0.00542011</td></tr>
    <tr align="center"><td>2</td><td>42.5433</td><td>1.952</td><td>0.00542011</td><td>42.4604</td><td>1.955</td><td>0.00542011</td><td>42.8079</td><td>1.939</td><td>0.00542011</td></tr>
    <tr align="center"><td>4</td><td>21.3244</td><td>3.894</td><td>0.00542011</td><td>21.7445</td><td>3.818</td><td>0.00542011</td><td>22.1411</td><td>3.749</td><td>0.00542011</td></tr>
    <tr align="center"><td>6</td><td>14.7643</td><td>5.624</td><td>0.00542011</td><td>14.8906</td><td>5.576</td><td>0.00542011</td><td>14.9234</td><td>5.563</td><td>0.00542011</td></tr>
</table>


### Lx = Ly = Lz = 1, N = 256, K = 2000 (с использования CUDA ядер)

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th><th colspan="3">Разбиение разложением на множители</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>3.19642</td><td>1.000</td><td>0.00542011</td><td>3.18848</td><td>1.000</td><td>0.00542011</td><td>3.20082</td><td>1.000</td><td>0.00542011</td></tr>
    <tr align="center"><td>2</td><td>2.81246</td><td>1.137</td><td>0.00542011</td><td>2.81574</td><td>1.132</td><td>0.00542011</td><td>2.79625</td><td>1.145</td><td>0.00542011</td></tr>
    <tr align="center"><td>4</td><td>1.39997</td><td>2.283</td><td>0.00542011</td><td>1.39985</td><td>2.278</td><td>0.00542011</td><td>1.39953</td><td>2.287</td><td>0.00542011</td></tr>
    <tr align="center"><td>6</td><td>0.986884</td><td>3.239</td><td>0.00542011</td><td>0.968354</td><td>3.293</td><td>0.00542011</td><td>0.969411</td><td>3.302</td><td>0.00542011</td></tr>
</table>

### Lx = Ly = Lz = 1, N = 512, K = 2000 (без использования CUDA ядер)

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th><th colspan="3">Разбиение разложением на множители</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>644.438</td><td>1.000</td><td>0.00542056</td><td>663.803</td><td>1.000</td><td>0.00542056</td><td>625.334</td><td>1.000</td><td>0.00542056</td></tr>
    <tr align="center"><td>2</td><td>328.795</td><td>1.960</td><td>0.00542056</td><td>335.254</td><td>1.980</td><td>0.00542056</td><td>334.403</td><td>1.870</td><td>0.00542056</td></tr>
    <tr align="center"><td>4</td><td>170.1</td><td>3.789</td><td>0.00542056</td><td>169.7</td><td>3.912</td><td>0.00542056</td><td>201.045</td><td>3.110</td><td>0.00542056</td></tr>
    <tr align="center"><td>6</td><td>114.436</td><td>5.631</td><td>0.00542056</td><td>117.82</td><td>5.634</td><td>0.00542056</td><td>112.172</td><td>5.575</td><td>0.00542056</td></tr>
</table>


### Lx = Ly = Lz = 1, N = 512, K = 2000 (с использования CUDA ядер)

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th><th colspan="3">Разбиение разложением на множители</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>25.6305</td><td>1.000</td><td>0.00542056</td><td>25.4344</td><td>1.000</td><td>0.00542056</td><td>25.4935</td><td>1.000</td><td>0.00542056</td></tr>
    <tr align="center"><td>2</td><td>21.9411</td><td>1.168</td><td>0.00542056</td><td>22.2357</td><td>1.144</td><td>0.00542056</td><td>22.2011</td><td>1.148</td><td>0.00542056</td></tr>
    <tr align="center"><td>4</td><td>10.5393</td><td>2.432</td><td>0.00542056</td><td>10.5945</td><td>2.401</td><td>0.00542056</td><td>10.5416</td><td>2.418</td><td>0.00542056</td></tr>
    <tr align="center"><td>6</td><td>7.13817</td><td>3.591</td><td>0.00542056</td><td>7.512</td><td>3.386</td><td>0.00542056</td><td>7.47077</td><td>3.412</td><td>0.00542056</td></tr>
</table>

### Lx = Ly = Lz = 1, N = 128, K = 2000, блочное разбиение

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">MPI</th><th colspan="3">MPI + CUDA</th><th rowspan="2">Ускорение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>10.8421</td><td>1.000</td><td>0.00541829</td><td>0.625809</td><td>1.000</td><td>0.00541829</td><td>17.325</td></tr>
    <tr align="center"><td>2</td><td>5.56431</td><td>1.949</td><td>0.00541829</td><td>0.491011</td><td>1.275</td><td>0.00541829</td><td>11.332</td></tr>
    <tr align="center"><td>4</td><td>2.80274</td><td>3.868</td><td>0.00541829</td><td>0.235194</td><td>2.661</td><td>0.00541829</td><td>11.917</td></tr>
    <tr align="center"><td>6</td><td>1.95402</td><td>5.549</td><td>0.00541829</td><td>0.23938</td><td>2.614</td><td>0.00541829</td><td>8.163</td></tr>
</table>

### Lx = Ly = Lz = 1, N = 128, K = 2000, ленточное разбиение

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">MPI</th><th colspan="3">MPI + CUDA</th><th rowspan="2">Ускорение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>10.828</td><td>1.000</td><td>0.00541829</td><td>0.417379</td><td>1.000</td><td>0.00541829</td><td>25.943</td></tr>
    <tr align="center"><td>2</td><td>5.54688</td><td>1.952</td><td>0.00541829</td><td>0.395455</td><td>1.055</td><td>0.00541829</td><td>14.027</td></tr>
    <tr align="center"><td>4</td><td>2.87798</td><td>3.762</td><td>0.00541829</td><td>0.235013</td><td>1.776</td><td>0.00541829</td><td>12.246</td></tr>
    <tr align="center"><td>6</td><td>2.03811</td><td>5.313</td><td>0.00541829</td><td>0.254023</td><td>1.643</td><td>0.00541829</td><td>8.023</td></tr>
</table>


### Lx = Ly = Lz = 1, N = 128, K = 2000, разбиение разложением на множители

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">MPI</th><th colspan="3">MPI + CUDA</th><th rowspan="2">Ускорение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>10.8417</td><td>1.000</td><td>0.00541829</td><td>0.41405</td><td>1.000</td><td>0.00541829</td><td>26.185</td></tr>
    <tr align="center"><td>2</td><td>5.61019</td><td>1.933</td><td>0.00541829</td><td>0.398651</td><td>1.039</td><td>0.00541829</td><td>14.073</td></tr>
    <tr align="center"><td>4</td><td>2.96835</td><td>3.652</td><td>0.00541829</td><td>0.234564</td><td>1.765</td><td>0.00541829</td><td>12.655</td></tr>
    <tr align="center"><td>6</td><td>1.97781</td><td>5.482</td><td>0.00541829</td><td>0.193185</td><td>2.143</td><td>0.00541829</td><td>10.238</td></tr>
</table>


### Lx = Ly = Lz = 1, N = 256, K = 2000, блочное разбиение

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">MPI</th><th colspan="3">MPI + CUDA</th><th rowspan="2">Ускорение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>83.0349</td><td>1.000</td><td>0.00542011</td><td>3.19642</td><td>1.000</td><td>0.00542011</td><td>25.977</td></tr>
    <tr align="center"><td>2</td><td>42.5433</td><td>1.952</td><td>0.00542011</td><td>2.81246</td><td>1.137</td><td>0.00542011</td><td>15.127</td></tr>
    <tr align="center"><td>4</td><td>21.3244</td><td>3.894</td><td>0.00542011</td><td>1.39997</td><td>2.283</td><td>0.00542011</td><td>15.232</td></tr>
    <tr align="center"><td>6</td><td>14.7643</td><td>5.624</td><td>0.00542011</td><td>0.986884</td><td>3.239</td><td>0.00542011</td><td>14.961</td></tr>
</table>

### Lx = Ly = Lz = 1, N = 256, K = 2000, ленточное разбиение

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">MPI</th><th colspan="3">MPI + CUDA</th><th rowspan="2">Ускорение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>83.024</td><td>1.000</td><td>0.00542011</td><td>3.18848</td><td>1.000</td><td>0.00542011</td><td>26.039</td></tr>
    <tr align="center"><td>2</td><td>42.4604</td><td>1.955</td><td>0.00542011</td><td>2.81574</td><td>1.132</td><td>0.00542011</td><td>15.080</td></tr>
    <tr align="center"><td>4</td><td>21.7445</td><td>3.818</td><td>0.00542011</td><td>1.39985</td><td>2.278</td><td>0.00542011</td><td>15.533</td></tr>
    <tr align="center"><td>6</td><td>14.8906</td><td>5.576</td><td>0.00542011</td><td>0.968354</td><td>3.293</td><td>0.00542011</td><td>15.377</td></tr>
</table>


### Lx = Ly = Lz = 1, N = 256, K = 2000, разбиение разложением на множители

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">MPI</th><th colspan="3">MPI + CUDA</th><th rowspan="2">Ускорение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>83.0116</td><td>1.000</td><td>0.00542011</td><td>3.20082</td><td>1.000</td><td>0.00542011</td><td>25.934</td></tr>
    <tr align="center"><td>2</td><td>42.8079</td><td>1.939</td><td>0.00542011</td><td>2.79625</td><td>1.145</td><td>0.00542011</td><td>15.309</td></tr>
    <tr align="center"><td>4</td><td>22.1411</td><td>3.749</td><td>0.00542011</td><td>1.39953</td><td>2.287</td><td>0.00542011</td><td>15.820</td></tr>
    <tr align="center"><td>6</td><td>14.9234</td><td>5.563</td><td>0.00542011</td><td>0.969411</td><td>3.302</td><td>0.00542011</td><td>15.394</td></tr>
</table>


### Lx = Ly = Lz = 1, N = 512, K = 2000, блочное разбиение

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">MPI</th><th colspan="3">MPI + CUDA</th><th rowspan="2">Ускорение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>644.438</td><td>1.000</td><td>0.00542056</td><td>25.6305</td><td>1.000</td><td>0.00542056</td><td>25.143</td></tr>
    <tr align="center"><td>2</td><td>328.795</td><td>1.960</td><td>0.00542056</td><td>21.9411</td><td>1.168</td><td>0.00542056</td><td>14.985</td></tr>
    <tr align="center"><td>4</td><td>170.1</td><td>3.789</td><td>0.00542056</td><td>10.5393</td><td>2.432</td><td>0.00542056</td><td>16.140</td></tr>
    <tr align="center"><td>6</td><td>114.436</td><td>5.631</td><td>0.00542056</td><td>7.13817</td><td>3.591</td><td>0.00542056</td><td>16.032</td></tr>
</table>

### Lx = Ly = Lz = 1, N = 512, K = 2000, ленточное разбиение

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">MPI</th><th colspan="3">MPI + CUDA</th><th rowspan="2">Ускорение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>663.803</td><td>1.000</td><td>0.00542056</td><td>25.4344</td><td>1.000</td><td>0.00542056</td><td>26.099</td></tr>
    <tr align="center"><td>2</td><td>335.254</td><td>1.980</td><td>0.00542056</td><td>22.2357</td><td>1.144</td><td>0.00542056</td><td>15.077</td></tr>
    <tr align="center"><td>4</td><td>169.7</td><td>3.912</td><td>0.00542056</td><td>10.5945</td><td>2.401</td><td>0.00542056</td><td>16.018</td></tr>
    <tr align="center"><td>6</td><td>117.82</td><td>5.634</td><td>0.00542056</td><td>7.512</td><td>3.386</td><td>0.00542056</td><td>15.684</td></tr>
</table>

### Lx = Ly = Lz = 1, N = 512, K = 2000, разбиение разложением на множители

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">MPI</th><th colspan="3">MPI + CUDA</th><th rowspan="2">Ускорение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>625.334</td><td>1.000</td><td>0.00542056</td><td>25.4935</td><td>1.000</td><td>0.00542056</td><td>24.529</td></tr>
    <tr align="center"><td>2</td><td>334.403</td><td>1.870</td><td>0.00542056</td><td>22.2011</td><td>1.148</td><td>0.00542056</td><td>15.062</td></tr>
    <tr align="center"><td>4</td><td>201.045</td><td>3.110</td><td>0.00542056</td><td>10.5416</td><td>2.418</td><td>0.00542056</td><td>19.072</td></tr>
    <tr align="center"><td>6</td><td>112.172</td><td>5.575</td><td>0.00542056</td><td>7.47077</td><td>3.412</td><td>0.00542056</td><td>15.015</td></tr>
</table>

</details>


## Blue Gene (MPI версия)

<details>
  <summary>Показать/скрыть таблицы</summary>

### Lx = Ly = Lz = 1, N = 128, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>141.493</td><td>1.000</td><td>0.00541829</td><td>141.494</td><td>1.000</td><td>0.00541829</td></tr>
    <tr align="center"><td>2</td><td>77.0068</td><td>1.837</td><td>0.00541829</td><td>77.0072</td><td>1.837</td><td>0.00541829</td></tr>
    <tr align="center"><td>4</td><td>40.8615</td><td>3.463</td><td>0.00541829</td><td>40.0407</td><td>3.534</td><td>0.00541829</td></tr>
    <tr align="center"><td>8</td><td>21.9844</td><td>6.436</td><td>0.00541829</td><td>20.883</td><td>6.776</td><td>0.00541829</td></tr>
    <tr align="center"><td>16</td><td>11.5286</td><td>12.273</td><td>0.00541829</td><td>11.2159</td><td>12.615</td><td>0.00541829</td></tr>
    <tr align="center"><td>32</td><td>6.04969</td><td>23.388</td><td>0.00541829</td><td>7.17778</td><td>19.713</td><td>0.00541829</td></tr>
    <tr align="center"><td>64</td><td>3.23408</td><td>43.751</td><td>0.00541829</td><td>6.74018</td><td>20.993</td><td>0.00541829</td></tr>
    <tr align="center"><td>128</td><td>1.69645</td><td>83.405</td><td>0.00541829</td><td>9.62115</td><td>14.707</td><td>0.00541829</td></tr>
    <tr align="center"><td>256</td><td>0.901327</td><td>156.983</td><td>0.00541829</td><td>-</td><td>-</td><td>-</td></tr>
</table>

### Lx = Ly = Lz = 1, N = 256, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>1110.44</td><td>1.000</td><td>0.00542011</td><td>1110.44</td><td>1.000</td><td>0.00542011</td></tr>
    <tr align="center"><td>2</td><td>605.897</td><td>1.833</td><td>0.00542011</td><td>605.875</td><td>1.833</td><td>0.00542011</td></tr>
    <tr align="center"><td>4</td><td>322.288</td><td>3.445</td><td>0.00542011</td><td>315.013</td><td>3.525</td><td>0.00542011</td></tr>
    <tr align="center"><td>8</td><td>174.383</td><td>6.368</td><td>0.00542011</td><td>162.972</td><td>6.814</td><td>0.00542011</td></tr>
    <tr align="center"><td>16</td><td>90.9261</td><td>12.213</td><td>0.00542011</td><td>84.5998</td><td>13.126</td><td>0.00542011</td></tr>
    <tr align="center"><td>32</td><td>47.7614</td><td>23.250</td><td>0.00542011</td><td>48.3104</td><td>22.986</td><td>0.00542011</td></tr>
    <tr align="center"><td>64</td><td>25.336</td><td>43.829</td><td>0.00542011</td><td>36.0454</td><td>30.807</td><td>0.00542011</td></tr>
    <tr align="center"><td>128</td><td>13.4745</td><td>82.410</td><td>0.00542011</td><td>42.6008</td><td>26.066</td><td>0.00542011</td></tr>
    <tr align="center"><td>256</td><td>6.70683</td><td>165.569</td><td>0.00542011</td><td>70.5973</td><td>15.729</td><td>0.00542011</td></tr>
</table>

### Lx = Ly = Lz = 1, N = 512, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>8600.92</td><td>1.000</td><td>0.00542056</td><td>9423.54</td><td>1.000</td><td>0.00542056</td></tr>
    <tr align="center"><td>2</td><td>4751.89</td><td>1.810</td><td>0.00542056</td><td>4759.36</td><td>1.980</td><td>0.00542056</td></tr>
    <tr align="center"><td>4</td><td>2568.59</td><td>3.348</td><td>0.00542056</td><td>2491.81</td><td>3.782</td><td>0.00542056</td></tr>
    <tr align="center"><td>8</td><td>1380.79</td><td>6.229</td><td>0.00542056</td><td>1284.44</td><td>7.337</td><td>0.00542056</td></tr>
    <tr align="center"><td>16</td><td>722.163</td><td>11.910</td><td>0.00542056</td><td>656.146</td><td>14.362</td><td>0.00542056</td></tr>
    <tr align="center"><td>32</td><td>380.695</td><td>22.593</td><td>0.00542056</td><td>352.157</td><td>26.759</td><td>0.00542056</td></tr>
    <tr align="center"><td>64</td><td>202.1</td><td>42.558</td><td>0.00542056</td><td>223.976</td><td>42.074</td><td>0.00542056</td></tr>
    <tr align="center"><td>128</td><td>104.098</td><td>82.623</td><td>0.00542056</td><td>210.467</td><td>44.774</td><td>0.00542056</td></tr>
    <tr align="center"><td>256</td><td>52.0685</td><td>165.185</td><td>0.00542056</td><td>312.152</td><td>30.188</td><td>0.00542056</td></tr>
</table>


### Lx = Ly = Lz = π, N = 128, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>142.496</td><td>1.000</td><td>0.000549676</td><td>142.496</td><td>1.000</td><td>0.000549676</td></tr>
    <tr align="center"><td>2</td><td>77.4895</td><td>1.839</td><td>0.000549676</td><td>77.4894</td><td>1.839</td><td>0.000549676</td></tr>
    <tr align="center"><td>4</td><td>41.1135</td><td>3.466</td><td>0.000549676</td><td>40.3135</td><td>3.535</td><td>0.000549676</td></tr>
    <tr align="center"><td>8</td><td>22.1162</td><td>6.443</td><td>0.000549676</td><td>20.9752</td><td>6.794</td><td>0.000549676</td></tr>
    <tr align="center"><td>16</td><td>11.6038</td><td>12.280</td><td>0.000549676</td><td>11.2888</td><td>12.623</td><td>0.000549676</td></tr>
    <tr align="center"><td>32</td><td>6.08839</td><td>23.405</td><td>0.000549676</td><td>7.16593</td><td>19.885</td><td>0.000549676</td></tr>
    <tr align="center"><td>64</td><td>3.25425</td><td>43.788</td><td>0.000549676</td><td>6.74699</td><td>21.120</td><td>0.000549676</td></tr>
    <tr align="center"><td>128</td><td>1.70468</td><td>83.591</td><td>0.000549676</td><td>9.62002</td><td>14.812</td><td>0.000549676</td></tr>
    <tr align="center"><td>256</td><td>1.50589</td><td>94.626</td><td>0.000549676</td><td>-</td><td>-</td><td>-</td></tr>
</table>

### Lx = Ly = Lz = π, N = 256, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>1120.28</td><td>1.000</td><td>0.000549861</td><td>1120.28</td><td>1.000</td><td>0.000549861</td></tr>
    <tr align="center"><td>2</td><td>610.702</td><td>1.834</td><td>0.000549861</td><td>610.702</td><td>1.834</td><td>0.000549861</td></tr>
    <tr align="center"><td>4</td><td>324.514</td><td>3.452</td><td>0.000549861</td><td>317.597</td><td>3.527</td><td>0.000549861</td></tr>
    <tr align="center"><td>8</td><td>175.456</td><td>6.385</td><td>0.000549861</td><td>163.753</td><td>6.841</td><td>0.000549861</td></tr>
    <tr align="center"><td>16</td><td>91.5347</td><td>12.239</td><td>0.000549861</td><td>85.2042</td><td>13.148</td><td>0.000549861</td></tr>
    <tr align="center"><td>32</td><td>48.1049</td><td>23.288</td><td>0.000549861</td><td>48.4248</td><td>23.134</td><td>0.000549861</td></tr>
    <tr align="center"><td>64</td><td>25.5064</td><td>43.922</td><td>0.000549861</td><td>36.1936</td><td>30.952</td><td>0.000549861</td></tr>
    <tr align="center"><td>128</td><td>13.1408</td><td>85.252</td><td>0.000549861</td><td>42.6262</td><td>26.281</td><td>0.000549861</td></tr>
    <tr align="center"><td>256</td><td>7.2389</td><td>154.758</td><td>0.000549861</td><td>70.9196</td><td>15.796</td><td>0.000549861</td></tr>
</table>

### Lx = Ly = Lz = π, N = 512, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>8605.17</td><td>1.000</td><td>0.000549907</td><td>9431.81</td><td>1.000</td><td>0.000549907</td></tr>
    <tr align="center"><td>2</td><td>4754.13</td><td>1.810</td><td>0.000549907</td><td>4765.72</td><td>1.979</td><td>0.000549907</td></tr>
    <tr align="center"><td>4</td><td>2574.72</td><td>3.342</td><td>0.000549907</td><td>2497.13</td><td>3.777</td><td>0.000549907</td></tr>
    <tr align="center"><td>8</td><td>1392.04</td><td>6.182</td><td>0.000549907</td><td>1293.15</td><td>7.294</td><td>0.000549907</td></tr>
    <tr align="center"><td>16</td><td>728.384</td><td>11.814</td><td>0.000549907</td><td>661.571</td><td>14.257</td><td>0.000549907</td></tr>
    <tr align="center"><td>32</td><td>383.046</td><td>22.465</td><td>0.000549907</td><td>353.86</td><td>26.654</td><td>0.000549907</td></tr>
    <tr align="center"><td>64</td><td>202.853</td><td>42.421</td><td>0.000549907</td><td>225.365</td><td>41.851</td><td>0.000549907</td></tr>
    <tr align="center"><td>128</td><td>103.621</td><td>83.045</td><td>0.000549907</td><td>211.249</td><td>44.648</td><td>0.000549907</td></tr>
    <tr align="center"><td>256</td><td>52.5081</td><td>163.883</td><td>0.000549907</td><td>316.831</td><td>29.769</td><td>0.000549907</td></tr>
</table>

</details>

## Blue Gene (гибридная версия, блочное разбиение)

<details>
  <summary>Показать/скрыть таблицы</summary>

### Lx = Ly = Lz = 1, N = 128, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">MPI</th><th colspan="3">MPI+OpenMP</th><th rowspan="2">Ускорение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>141.493</td><td>1.000</td><td>0.00541829</td><td>87.9429</td><td>1.000</td><td>0.00541829</td><td>1.609</td></tr>
    <tr align="center"><td>2</td><td>77.0068</td><td>1.837</td><td>0.00541829</td><td>47.2583</td><td>1.861</td><td>0.00541829</td><td>1.629</td></tr>
    <tr align="center"><td>4</td><td>40.8615</td><td>3.463</td><td>0.00541829</td><td>24.6753</td><td>3.564</td><td>0.00541829</td><td>1.656</td></tr>
    <tr align="center"><td>8</td><td>21.9844</td><td>6.436</td><td>0.00541829</td><td>12.9544</td><td>6.789</td><td>0.00541829</td><td>1.697</td></tr>
    <tr align="center"><td>16</td><td>11.5286</td><td>12.273</td><td>0.00541829</td><td>6.78328</td><td>12.965</td><td>0.00541829</td><td>1.700</td></tr>
    <tr align="center"><td>32</td><td>6.04969</td><td>23.388</td><td>0.00541829</td><td>3.49826</td><td>25.139</td><td>0.00541829</td><td>1.729</td></tr>
    <tr align="center"><td>64</td><td>3.23408</td><td>43.751</td><td>0.00541829</td><td>1.83381</td><td>47.956</td><td>0.00541829</td><td>1.764</td></tr>
    <tr align="center"><td>128</td><td>1.69645</td><td>83.405</td><td>0.00541829</td><td>1.02677</td><td>85.650</td><td>0.00541829</td><td>1.652</td></tr>
    <tr align="center"><td>256</td><td>0.901327</td><td>156.983</td><td>0.00541829</td><td>0.603308</td><td>145.768</td><td>0.00541829</td><td>1.494</td></tr>
</table>

### Lx = Ly = Lz = 1, N = 256, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">MPI</th><th colspan="3">MPI+OpenMP</th><th rowspan="2">Ускорение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>1110.44</td><td>1.000</td><td>0.00542011</td><td>682.316</td><td>1.000</td><td>0.00542011</td><td>1.627</td></tr>
    <tr align="center"><td>2</td><td>605.897</td><td>1.833</td><td>0.00542011</td><td>366.57</td><td>1.861</td><td>0.00542011</td><td>1.653</td></tr>
    <tr align="center"><td>4</td><td>322.288</td><td>3.445</td><td>0.00542011</td><td>192.21</td><td>3.550</td><td>0.00542011</td><td>1.677</td></tr>
    <tr align="center"><td>8</td><td>174.383</td><td>6.368</td><td>0.00542011</td><td>100.854</td><td>6.765</td><td>0.00542011</td><td>1.729</td></tr>
    <tr align="center"><td>16</td><td>90.9261</td><td>12.213</td><td>0.00542011</td><td>52.1929</td><td>13.073</td><td>0.00542011</td><td>1.742</td></tr>
    <tr align="center"><td>32</td><td>47.7614</td><td>23.250</td><td>0.00542011</td><td>26.8175</td><td>25.443</td><td>0.00542011</td><td>1.781</td></tr>
    <tr align="center"><td>64</td><td>25.336</td><td>43.829</td><td>0.00542011</td><td>13.8613</td><td>49.225</td><td>0.00542011</td><td>1.828</td></tr>
    <tr align="center"><td>128</td><td>13.4745</td><td>82.410</td><td>0.00542011</td><td>7.21106</td><td>94.621</td><td>0.00542011</td><td>1.869</td></tr>
    <tr align="center"><td>256</td><td>6.70683</td><td>165.569</td><td>0.00542011</td><td>4.45263</td><td>153.239</td><td>0.00542011</td><td>1.506</td></tr>
</table>

### Lx = Ly = Lz = 1, N = 512, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">MPI</th><th colspan="3">MPI+OpenMP</th><th rowspan="2">Ускорение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>8600.92</td><td>1.000</td><td>0.00542056</td><td>5531.18</td><td>1.000</td><td>0.00542056</td><td>1.555</td></tr>
    <tr align="center"><td>2</td><td>4751.89</td><td>1.810</td><td>0.00542056</td><td>2895.61</td><td>1.910</td><td>0.00542056</td><td>1.641</td></tr>
    <tr align="center"><td>4</td><td>2568.59</td><td>3.348</td><td>0.00542056</td><td>1517.01</td><td>3.646</td><td>0.00542056</td><td>1.693</td></tr>
    <tr align="center"><td>8</td><td>1380.79</td><td>6.229</td><td>0.00542056</td><td>794.502</td><td>6.962</td><td>0.00542056</td><td>1.738</td></tr>
    <tr align="center"><td>16</td><td>722.163</td><td>11.910</td><td>0.00542056</td><td>407.951</td><td>13.558</td><td>0.00542056</td><td>1.770</td></tr>
    <tr align="center"><td>32</td><td>380.695</td><td>22.593</td><td>0.00542056</td><td>210.244</td><td>26.308</td><td>0.00542056</td><td>1.811</td></tr>
    <tr align="center"><td>64</td><td>202.1</td><td>42.558</td><td>0.00542056</td><td>108.913</td><td>50.785</td><td>0.00542056</td><td>1.856</td></tr>
    <tr align="center"><td>128</td><td>104.098</td><td>82.623</td><td>0.00542056</td><td>55.9623</td><td>98.838</td><td>0.00542056</td><td>1.860</td></tr>
    <tr align="center"><td>256</td><td>52.0685</td><td>165.185</td><td>0.00542056</td><td>29.8807</td><td>185.109</td><td>0.00542056</td><td>1.743</td></tr>
</table>


### Lx = Ly = Lz = π, N = 128, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">MPI</th><th colspan="3">MPI+OpenMP</th><th rowspan="2">Ускорение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>142.496</td><td>1.000</td><td>0.000549676</td><td>88.0216</td><td>1.000</td><td>0.000549676</td><td>1.619</td></tr>
    <tr align="center"><td>2</td><td>77.4895</td><td>1.839</td><td>0.000549676</td><td>47.2734</td><td>1.862</td><td>0.000549676</td><td>1.639</td></tr>
    <tr align="center"><td>4</td><td>41.1135</td><td>3.466</td><td>0.000549676</td><td>24.6981</td><td>3.564</td><td>0.000549676</td><td>1.665</td></tr>
    <tr align="center"><td>8</td><td>22.1162</td><td>6.443</td><td>0.000549676</td><td>12.9627</td><td>6.790</td><td>0.000549676</td><td>1.706</td></tr>
    <tr align="center"><td>16</td><td>11.6038</td><td>12.280</td><td>0.000549676</td><td>6.7755</td><td>12.991</td><td>0.000549676</td><td>1.713</td></tr>
    <tr align="center"><td>32</td><td>6.08839</td><td>23.405</td><td>0.000549676</td><td>3.49033</td><td>25.219</td><td>0.000549676</td><td>1.744</td></tr>
    <tr align="center"><td>64</td><td>3.25425</td><td>43.788</td><td>0.000549676</td><td>1.8336</td><td>48.005</td><td>0.000549676</td><td>1.775</td></tr>
    <tr align="center"><td>128</td><td>1.70468</td><td>83.591</td><td>0.000549676</td><td>1.02571</td><td>85.815</td><td>0.000549676</td><td>1.662</td></tr>
    <tr align="center"><td>256</td><td>1.50589</td><td>94.626</td><td>0.000549676</td><td>1.15484</td><td>76.220</td><td>0.000549676</td><td>1.304</td></tr>
</table>

### Lx = Ly = Lz = π, N = 256, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">MPI</th><th colspan="3">MPI+OpenMP</th><th rowspan="2">Ускорение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>1120.28</td><td>1.000</td><td>0.000549861</td><td>684.851</td><td>1.000</td><td>0.000549861</td><td>1.636</td></tr>
    <tr align="center"><td>2</td><td>610.702</td><td>1.834</td><td>0.000549861</td><td>367.672</td><td>1.863</td><td>0.000549861</td><td>1.661</td></tr>
    <tr align="center"><td>4</td><td>324.514</td><td>3.452</td><td>0.000549861</td><td>192.997</td><td>3.549</td><td>0.000549861</td><td>1.681</td></tr>
    <tr align="center"><td>8</td><td>175.456</td><td>6.385</td><td>0.000549861</td><td>101.139</td><td>6.771</td><td>0.000549861</td><td>1.735</td></tr>
    <tr align="center"><td>16</td><td>91.5347</td><td>12.239</td><td>0.000549861</td><td>52.277</td><td>13.100</td><td>0.000549861</td><td>1.751</td></tr>
    <tr align="center"><td>32</td><td>48.1049</td><td>23.288</td><td>0.000549861</td><td>26.8123</td><td>25.542</td><td>0.000549861</td><td>1.794</td></tr>
    <tr align="center"><td>64</td><td>25.5064</td><td>43.922</td><td>0.000549861</td><td>13.8648</td><td>49.395</td><td>0.000549861</td><td>1.840</td></tr>
    <tr align="center"><td>128</td><td>13.1408</td><td>85.252</td><td>0.000549861</td><td>7.21582</td><td>94.910</td><td>0.000549861</td><td>1.821</td></tr>
    <tr align="center"><td>256</td><td>7.2389</td><td>154.758</td><td>0.000549861</td><td>4.64308</td><td>147.499</td><td>0.000549861</td><td>1.559</td></tr>
</table>

### Lx = Ly = Lz = π, N = 512, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">MPI</th><th colspan="3">MPI+OpenMP</th><th rowspan="2">Ускорение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>8605.17</td><td>1.000</td><td>0.000549907</td><td>5541.26</td><td>1.000</td><td>0.000549907</td><td>1.553</td></tr>
    <tr align="center"><td>2</td><td>4754.13</td><td>1.810</td><td>0.000549907</td><td>2894.92</td><td>1.914</td><td>0.000549907</td><td>1.642</td></tr>
    <tr align="center"><td>4</td><td>2574.72</td><td>3.342</td><td>0.000549907</td><td>1519.12</td><td>3.648</td><td>0.000549907</td><td>1.695</td></tr>
    <tr align="center"><td>8</td><td>1392.04</td><td>6.182</td><td>0.000549907</td><td>798.428</td><td>6.940</td><td>0.000549907</td><td>1.743</td></tr>
    <tr align="center"><td>16</td><td>728.384</td><td>11.814</td><td>0.000549907</td><td>408.788</td><td>13.555</td><td>0.000549907</td><td>1.782</td></tr>
    <tr align="center"><td>32</td><td>383.046</td><td>22.465</td><td>0.000549907</td><td>210.271</td><td>26.353</td><td>0.000549907</td><td>1.822</td></tr>
    <tr align="center"><td>64</td><td>202.853</td><td>42.421</td><td>0.000549907</td><td>108.978</td><td>50.848</td><td>0.000549907</td><td>1.861</td></tr>
    <tr align="center"><td>128</td><td>103.621</td><td>83.045</td><td>0.000549907</td><td>56.0311</td><td>98.896</td><td>0.000549907</td><td>1.849</td></tr>
    <tr align="center"><td>256</td><td>52.5081</td><td>163.883</td><td>0.000549907</td><td>29.3193</td><td>188.997</td><td>0.000549907</td><td>1.791</td></tr>
</table>

</details>
