# SupercomputerModelingTechnologyTask2
Реализация второго практического задания в рамках курса "Суперкомпьютерное моделирование и технологии"

## Программная реализация

Реализовано две программы: последовательная и гибридная параллельная (MPI + OpenMP). Каждая из программ является консольным приложением и принимает входные данные в виде аргументов командной строки. Используются следующие аргументы (доступно при вызове с единственным аргументом `--help`):

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

Для визуализации получаемых решений написан визуализатор (принимает json файлы, генерируемые программой): https://programforyou.ru/tests/supercomputer-modeling-technology/solve-visualizer


## Сборка

Сборки осуществляется через makefile. Для сборки последовательной версии выполните `make main`. Для сборки параллельной версии выполните `make main-mpi`.

### Пример запуска
```bash
>make main
>./main -d -N 64 -K 150 -btx f -bty pa -btz pn -Lx 0.5 -Ly 12.5 -Lz 100
```


### Особенности параллельной реализации

Для распараллеливания вся сетка разбивается на области (также прямоугольные параллелепипеды) в количестве используемых процессов по следующему алгоритму:
* начнём разбиение с параллелепипеда [0, N] x [0, N] x [0, N], выберем начальную ось (X) и запустим рекурсивный процесс
* если текущее количество областей (size) равно 1, вернём обрабатываемый параллелепипед
* если размер нечётный, то по текущей оси выберем область 1 / size и сделаем из неё параллелепипед, и продолжим разбивать область 1 - 1 / size
* по выбранной оси делим область пополам и рекурсивно запускаем для этих подобластей

Для большей наглядности реализован визуализатор разбиения с возможностью установки размера сетки и числа процессов: https://programforyou.ru/tests/supercomputer-modeling-technology/split-visualizer

### Получающиеся разбиения

<table>
    <tr>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/split2.png"><br>P=2</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/split3.png"><br>P=3</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/split4.png"><br>P=4</td>
    </tr>
    <tr>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/split5.png"><br>P=5</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/split10.png"><br>P=10</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/split14.png"><br>P=14</td>
    </tr>
    <tr>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/split25.png"><br>P=25</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/split30.png"><br>P=30</td>
        <td align="center"><img src="https://github.com/dronperminov/SupercomputerModelingTechnologyTask2/blob/master/examples/split60.png"><br>P=60</td>
    </tr>
</table>

## Графики аналитического и полученного решений

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


## Polus (без использования OpenMP)

### Lx = Ly = Lz = 1, N = 128, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>3.37891</td><td>1.000</td><td>0.00541829</td><td>3.38371</td><td>1.000</td><td>0.00541829</td></tr>
    <tr align="center"><td>2</td><td>1.97175</td><td>1.714</td><td>0.00541829</td><td>1.98044</td><td>1.709</td><td>0.00541829</td></tr>
    <tr align="center"><td>4</td><td>1.13682</td><td>2.972</td><td>0.00541829</td><td>1.13015</td><td>2.994</td><td>0.00541829</td></tr>
    <tr align="center"><td>8</td><td>0.761273</td><td>4.438</td><td>0.00541829</td><td>0.715661</td><td>4.728</td><td>0.00541829</td></tr>
    <tr align="center"><td>10</td><td>0.697569</td><td>4.844</td><td>0.00541829</td><td>0.573838</td><td>5.897</td><td>0.00541829</td></tr>
    <tr align="center"><td>16</td><td>0.430247</td><td>7.853</td><td>0.00541829</td><td>0.409231</td><td>8.268</td><td>0.00541829</td></tr>
    <tr align="center"><td>20</td><td>0.389413</td><td>8.677</td><td>0.00541829</td><td>0.36037</td><td>9.390</td><td>0.00541829</td></tr>
    <tr align="center"><td>32</td><td>0.302416</td><td>11.173</td><td>0.00541829</td><td>0.259181</td><td>13.055</td><td>0.00541829</td></tr>
    <tr align="center"><td>40</td><td>0.24562</td><td>13.757</td><td>0.00541829</td><td>0.250944</td><td>13.484</td><td>0.00541829</td></tr>
    <tr align="center"><td>64</td><td>0.175378</td><td>19.266</td><td>0.00541829</td><td>0.182393</td><td>18.552</td><td>0.00541829</td></tr>
</table>

### Lx = Ly = Lz = 1, N = 256, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>26.8418</td><td>1.000</td><td>0.00542011</td><td>26.8246</td><td>1.000</td><td>0.00542011</td></tr>
    <tr align="center"><td>2</td><td>15.5534</td><td>1.726</td><td>0.00542011</td><td>15.6095</td><td>1.718</td><td>0.00542011</td></tr>
    <tr align="center"><td>4</td><td>8.5587</td><td>3.136</td><td>0.00542011</td><td>8.58564</td><td>3.124</td><td>0.00542011</td></tr>
    <tr align="center"><td>8</td><td>5.74982</td><td>4.668</td><td>0.00542011</td><td>5.16543</td><td>5.193</td><td>0.00542011</td></tr>
    <tr align="center"><td>10</td><td>5.12554</td><td>5.237</td><td>0.00542011</td><td>4.26427</td><td>6.291</td><td>0.00542011</td></tr>
    <tr align="center"><td>16</td><td>3.6738</td><td>7.306</td><td>0.00542011</td><td>2.85054</td><td>9.410</td><td>0.00542011</td></tr>
    <tr align="center"><td>20</td><td>3.17586</td><td>8.452</td><td>0.00542011</td><td>2.39508</td><td>11.200</td><td>0.00542011</td></tr>
    <tr align="center"><td>32</td><td>2.13905</td><td>12.548</td><td>0.00542011</td><td>1.74538</td><td>15.369</td><td>0.00542011</td></tr>
    <tr align="center"><td>40</td><td>1.84101</td><td>14.580</td><td>0.00542011</td><td>1.55653</td><td>17.234</td><td>0.00542011</td></tr>
    <tr align="center"><td>64</td><td>1.23209</td><td>21.786</td><td>0.00542011</td><td>1.12446</td><td>23.856</td><td>0.00542011</td></tr>
</table>

### Lx = Ly = Lz = 1, N = 512, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>212.339</td><td>1.000</td><td>0.00542056</td><td>212.174</td><td>1.000</td><td>0.00542056</td></tr>
    <tr align="center"><td>2</td><td>119.572</td><td>1.776</td><td>0.00542056</td><td>120.541</td><td>1.760</td><td>0.00542056</td></tr>
    <tr align="center"><td>4</td><td>66.7043</td><td>3.183</td><td>0.00542056</td><td>66.9953</td><td>3.167</td><td>0.00542056</td></tr>
    <tr align="center"><td>8</td><td>38.7001</td><td>5.487</td><td>0.00542056</td><td>35.9597</td><td>5.900</td><td>0.00542056</td></tr>
    <tr align="center"><td>10</td><td>38.6523</td><td>5.494</td><td>0.00542056</td><td>31.6109</td><td>6.712</td><td>0.00542056</td></tr>
    <tr align="center"><td>16</td><td>22.3261</td><td>9.511</td><td>0.00542056</td><td>20.8502</td><td>10.176</td><td>0.00542056</td></tr>
    <tr align="center"><td>20</td><td>23.1022</td><td>9.191</td><td>0.00542056</td><td>19.5098</td><td>10.875</td><td>0.00542056</td></tr>
    <tr align="center"><td>32</td><td>14.4698</td><td>14.675</td><td>0.00542056</td><td>12.6782</td><td>16.735</td><td>0.00542056</td></tr>
    <tr align="center"><td>40</td><td>13.1529</td><td>16.144</td><td>0.00542056</td><td>10.8413</td><td>19.571</td><td>0.00542056</td></tr>
    <tr align="center"><td>64</td><td>8.52932</td><td>24.895</td><td>0.00542056</td><td>7.82263</td><td>27.123</td><td>0.00542056</td></tr>
</table>

### Lx = Ly = Lz = π, N = 128, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>3.37747</td><td>1.000</td><td>0.000549676</td><td>3.38126</td><td>1.000</td><td>0.000549676</td></tr>
    <tr align="center"><td>2</td><td>1.99226</td><td>1.695</td><td>0.000549676</td><td>1.98962</td><td>1.699</td><td>0.000549676</td></tr>
    <tr align="center"><td>4</td><td>1.12936</td><td>2.991</td><td>0.000549676</td><td>1.13578</td><td>2.977</td><td>0.000549676</td></tr>
    <tr align="center"><td>8</td><td>0.847514</td><td>3.985</td><td>0.000549676</td><td>0.579493</td><td>5.835</td><td>0.000549676</td></tr>
    <tr align="center"><td>10</td><td>0.673376</td><td>5.016</td><td>0.000549676</td><td>0.552498</td><td>6.120</td><td>0.000549676</td></tr>
    <tr align="center"><td>16</td><td>0.46852</td><td>7.209</td><td>0.000549676</td><td>0.402093</td><td>8.409</td><td>0.000549676</td></tr>
    <tr align="center"><td>20</td><td>0.421829</td><td>8.007</td><td>0.000549676</td><td>0.378469</td><td>8.934</td><td>0.000549676</td></tr>
    <tr align="center"><td>32</td><td>0.29494</td><td>11.451</td><td>0.000549676</td><td>0.257307</td><td>13.141</td><td>0.000549676</td></tr>
    <tr align="center"><td>40</td><td>0.278028</td><td>12.148</td><td>0.000549676</td><td>0.227738</td><td>14.847</td><td>0.000549676</td></tr>
    <tr align="center"><td>64</td><td>0.172591</td><td>19.569</td><td>0.000549676</td><td>0.195346</td><td>17.309</td><td>0.000549676</td></tr>
</table>

### Lx = Ly = Lz = π, N = 256, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>26.8546</td><td>1.000</td><td>0.000549861</td><td>26.8581</td><td>1.000</td><td>0.000549861</td></tr>
    <tr align="center"><td>2</td><td>15.2069</td><td>1.766</td><td>0.000549861</td><td>14.7859</td><td>1.816</td><td>0.000549861</td></tr>
    <tr align="center"><td>4</td><td>8.37291</td><td>3.207</td><td>0.000549861</td><td>9.6132</td><td>2.794</td><td>0.000549861</td></tr>
    <tr align="center"><td>8</td><td>4.72224</td><td>5.687</td><td>0.000549861</td><td>4.50557</td><td>5.961</td><td>0.000549861</td></tr>
    <tr align="center"><td>10</td><td>4.64202</td><td>5.785</td><td>0.000549861</td><td>4.07747</td><td>6.587</td><td>0.000549861</td></tr>
    <tr align="center"><td>16</td><td>3.2318</td><td>8.309</td><td>0.000549861</td><td>2.55629</td><td>10.507</td><td>0.000549861</td></tr>
    <tr align="center"><td>20</td><td>2.85524</td><td>9.405</td><td>0.000549861</td><td>2.35552</td><td>11.402</td><td>0.000549861</td></tr>
    <tr align="center"><td>32</td><td>2.20191</td><td>12.196</td><td>0.000549861</td><td>1.76651</td><td>15.204</td><td>0.000549861</td></tr>
    <tr align="center"><td>40</td><td>2.04728</td><td>13.117</td><td>0.000549861</td><td>1.57969</td><td>17.002</td><td>0.000549861</td></tr>
    <tr align="center"><td>64</td><td>1.2482</td><td>21.515</td><td>0.000549861</td><td>1.11757</td><td>24.033</td><td>0.000549861</td></tr>
</table>

### Lx = Ly = Lz = π, N = 512, K = 2000

<table>
    <tr><th rowspan="2">Число MPI процессов (P)</th><th colspan="3">Блочное разбиение</th><th colspan="3">Ленточное разбиение</th></tr>
    <tr><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th><th>Время решения (с)</th><th>Ускорение</th><th>Погрешность</th></tr>
    <tr align="center"><td>1</td><td>212.134</td><td>1.000</td><td>0.000549907</td><td>212.326</td><td>1.000</td><td>0.000549907</td></tr>
    <tr align="center"><td>2</td><td>119.237</td><td>1.779</td><td>0.000549907</td><td>119.366</td><td>1.779</td><td>0.000549907</td></tr>
    <tr align="center"><td>4</td><td>66.5768</td><td>3.186</td><td>0.000549907</td><td>67.0074</td><td>3.169</td><td>0.000549907</td></tr>
    <tr align="center"><td>8</td><td>40.2159</td><td>5.275</td><td>0.000549907</td><td>34.0599</td><td>6.234</td><td>0.000549907</td></tr>
    <tr align="center"><td>10</td><td>36.5986</td><td>5.796</td><td>0.000549907</td><td>27.4475</td><td>7.736</td><td>0.000549907</td></tr>
    <tr align="center"><td>16</td><td>23.8981</td><td>8.877</td><td>0.000549907</td><td>21.9308</td><td>9.682</td><td>0.000549907</td></tr>
    <tr align="center"><td>20</td><td>23.9314</td><td>8.864</td><td>0.000549907</td><td>17.8668</td><td>11.884</td><td>0.000549907</td></tr>
    <tr align="center"><td>32</td><td>14.2082</td><td>14.930</td><td>0.000549907</td><td>13.0011</td><td>16.331</td><td>0.000549907</td></tr>
    <tr align="center"><td>40</td><td>13.68</td><td>15.507</td><td>0.000549907</td><td>10.1392</td><td>20.941</td><td>0.000549907</td></tr>
    <tr align="center"><td>64</td><td>8.81245</td><td>24.072</td><td>0.000549907</td><td>7.60838</td><td>27.907</td><td>0.000549907</td></tr>
</table>


## Bluegene

### Lx = Ly = Lz = 1, N = 128, K = 2000

| Число MPI процессов (P) | Время решения (с) | Ускорение | Погрешность |
|                     :-: |               :-: |       :-: |         :-: |
|                       1 |           107.683 |           |  0.00541829 |
|                       2 |           313.807 |           |  0.00541829 |
|                       4 |           438.996 |           |  0.00541829 |
|                       8 |           217.415 |           |  0.00541829 |
|                      16 |           230.718 |           |  0.00541829 |
|                      32 |           125.958 |           |  0.00541829 |
|                      64 |           56.3501 |           |  0.00541829 |
|                     128 |            24.533 |           |  0.00541829 |
|                     256 |           9.87815 |           |  0.00541829 |
