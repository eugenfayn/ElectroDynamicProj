Задача 5: Электродинамика
==================================================

Группа: Файн Е.Е., Баринов Д.С., Мещеряков А.О.

TODO + предварительное распределение
------------

- [ ] Код-ревью веток
  - @alexmeshr
- [ ] Эскиз проекта на Python
  - @alexmeshr
- [x] Стандартизация проекта
  - @dsbarinov1
- [x] CmakeLists + README.md
  - @dsbarinov1
- [ ] Шаблон функции универсального интегрирования
  - @eugenfayn
- [x] Переход к барицентрическим координатам
  - @dsbarinov1
- [x] Квадратуры Гаусса-3/Гаусса-4 (Не забыть умножить веса w_i на 2)
  - @dsbarinov1
- [ ] Вычисление правой части задачи (Dep: 2,4)
  - @eugenfayn
- [ ] Вычисление матрицы A (Dep:2,4)
  - @eugenfayn
- [ ] Перегрузки операций на классы Vertex, Face, Edge (Базовые операции с векторами, получение объектов-точек)
  - @alexmeshr
- [ ] Визуализация в сферических координатах
  - @eugenfayn
- [ ] Визуализация в ParaView (?)
  - @alexmeshr
- [ ] Презентация []
  - @dsbarinov1, @eugenfayn, @alexmeshr

Структура проекта
------------

```bash
ElectroDynamicProj/
├── include/         # Заголовочные файлы
│   ├── geometry/    # Геометрия
│   ├── quadrature/  # Квадратуры Гаусса
│   └── config.h.in  # Конфигурационный файл
├── src/             # Исходный код
├── data/            # Данные
│   ├── input/       # Входные
│   └── output/      # Выходные
├── scripts/         # Скрипты для визуализации
└── validation       # Данные для проверки результата
```

Формат входных данных
------------

```bash
Points <vertex_count>
x1 y1 z1
x2 y2 z2
...
Frames <face_count>
v1 v2 v3
v4 v5 v6
...
```

Зависимости
------------

- CMake 3.10+
- C++17
- Python 3.x (для визуализации)
- ParaView (для визуализации, опционально)

Сборка
------------

По умолчанию (-O0):

```bash
mkdir build
cd build
cmake ..
make
cd ..
```

One-liner:

```bash
mkdir build; cd build; cmake ..; make; cd ..
```

С оптимизациями (-On, n=0,1,2,3):

```bash
mkdir build
cd build
cmake .. -DOPTIMIZE_O3=ON
make
cd ..
```

One-liner:

```bash
mkdir build; cd build; cmake .. -DOPTIMIZE_O3=ON; make; cd ..
```

Запуск
-----

```bash
./build/bin/integral_solver
```
