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
* `-on` – путь к json файлу для сохранения численного решения (по умолчанию не используется)
* `-oa` – путь к json файлу для сохранения аналитического решения (по умолчанию не используется)
* `-o` – путь к текстовому файлу для вывода результатов (по умолчанию output.txt)

### Типы граничных условий:
* `first-kind (f)` – однородные граничные условия первого рода
* `periodic-analytical (pa)` – аналитические периодические граничные условия
* `periodic-numerical (pn)` – численные периодические граничные условия

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

| Число MPI процессов (P) | Время решения (с) | Ускорение | Погрешность |
|                     :-: |               :-: |       :-: |         :-: |
|                       1 |           3.34033 |         1 |  0.00541829 |
|                       2 |           4.20456 |     0.794 |  0.00541829 |
|                       4 |           2.84689 |     1.173 |  0.00541829 |
|                       8 |           1.74623 |     1.913 |  0.00541829 |
|                      10 |           1.82461 |     1.831 |  0.00541829 |
|                      16 |           1.18504 |     2.819 |  0.00541829 |
|                      20 |           1.22239 |     2.733 |  0.00541829 |
|                      32 |          0.786723 |     4.246 |  0.00541829 |
|                      40 |          0.792015 |     4.218 |  0.00541829 |
|                      64 |          0.551487 |     6.057 |  0.00541829 |


### Lx = Ly = Lz = 1, N = 256, K = 2000

| Число MPI процессов (P) | Время решения (с) | Ускорение | Погрешность |
|                     :-: |               :-: |       :-: |         :-: |
|                       1 |           27.6079 |         1 |  0.00542011 |
|                       2 |           31.7398 |     0.870 |  0.00542011 |
|                       4 |           19.0134 |     1.452 |  0.00542011 |
|                       8 |           13.2252 |     2.088 |  0.00542011 |
|                      10 |           10.9807 |     2.514 |  0.00542011 |
|                      16 |           7.04918 |     3.916 |  0.00542011 |
|                      20 |           7.21181 |     3.828 |  0.00542011 |
|                      32 |           5.40517 |     5.108 |  0.00542011 |
|                      40 |           5.18453 |     5.325 |  0.00542011 |
|                      64 |           2.96466 |     9.312 |  0.00542011 |


### Lx = Ly = Lz = 1, N = 512, K = 2000

| Число MPI процессов (P) | Время решения (с) | Ускорение | Погрешность |
|                     :-: |               :-: |       :-: |         :-: |
|                       1 |           219.528 |         1 |  0.00542056 |
|                       2 |           247.193 |     0.888 |  0.00542056 |
|                       4 |           145.986 |     1.504 |  0.00542056 |
|                       8 |           90.5833 |     2.423 |  0.00542056 |
|                      10 |           87.2427 |     2.516 |  0.00542056 |
|                      16 |           55.9972 |     3.920 |  0.00542056 |
|                      20 |           50.7379 |     4.327 |  0.00542056 |
|                      32 |           32.6001 |     6.734 |  0.00542056 |
|                      40 |           33.5516 |     6.543 |  0.00542056 |
|                      64 |           19.1128 |    11.486 |  0.00542056 |


### Lx = Ly = Lz = π, N = 128, K = 2000

| Число MPI процессов (P) | Время решения (с) | Ускорение | Погрешность |
|                     :-: |               :-: |       :-: |         :-: |
|                       1 |           3.34223 |         1 | 0.000549676 |
|                       2 |            4.2048 |     0.795 | 0.000549676 |
|                       4 |           2.66891 |     1.252 | 0.000549676 |
|                       8 |            1.7306 |     1.931 | 0.000549676 |
|                      10 |           1.77632 |     1.882 | 0.000549676 |
|                      16 |           1.17823 |     2.837 | 0.000549676 |
|                      20 |           1.20883 |     2.765 | 0.000549676 |
|                      32 |          0.907298 |     3.684 | 0.000549676 |
|                      40 |          0.792064 |     4.220 | 0.000549676 |
|                      64 |          0.504406 |     6.626 | 0.000549676 |


### Lx = Ly = Lz = π, N = 256, K = 2000

| Число MPI процессов (P) | Время решения (с) | Ускорение | Погрешность |
|                     :-: |               :-: |       :-: |         :-: |
|                       1 |           27.7564 |         1 | 0.000549861 |
|                       2 |           31.7372 |     0.875 | 0.000549861 |
|                       4 |           19.1404 |     1.450 | 0.000549861 |
|                       8 |           11.3117 |     2.454 | 0.000549861 |
|                      10 |           11.1413 |     2.491 | 0.000549861 |
|                      16 |           7.48548 |     3.708 | 0.000549861 |
|                      20 |           8.70398 |     3.189 | 0.000549861 |
|                      32 |           4.93602 |     5.623 | 0.000549861 |
|                      40 |           5.06855 |     5.476 | 0.000549861 |
|                      64 |           2.97572 |     9.328 | 0.000549861 |


### Lx = Ly = Lz = π, N = 512, K = 2000

| Число MPI процессов (P) | Время решения (с) | Ускорение | Погрешность |
|                     :-: |               :-: |       :-: |         :-: |
|                       1 |           223.427 |         1 | 0.000549907 |
|                       2 |           247.067 |     0.904 | 0.000549907 |
|                       4 |           145.405 |     1.537 | 0.000549907 |
|                       8 |           84.3112 |     2.650 | 0.000549907 |
|                      10 |           86.7045 |     2.577 | 0.000549907 |
|                      16 |            55.649 |     4.015 | 0.000549907 |
|                      20 |           56.9751 |     3.921 | 0.000549907 |
|                      32 |            32.388 |     6.898 | 0.000549907 |
|                      40 |           30.8168 |     7.250 | 0.000549907 |
|                      64 |           19.6482 |    11.371 | 0.000549907 |
