#include <iostream>
#include <cmath>

// Функция, которую мы хотим интегрировать
double f(double x) {
    return x * x;
}

// Функция, которую мы хотим интегрировать
double f(double x, double y, double z) {
    return x * x + y * y + z * z;
}

// Метод трапеций для численного интегрирования
double integrate(double a, double b, int n) {
    double h = (b - a) / n; // Ширина трапеции
    double sum = 0.5 * (f(a) + f(b)); // Начальная сумма

    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        sum += f(x);
    }

    return sum * h;
}

// Метод трапеций для численного интегрирования по объему
double integrateVolume(double ax, double bx, double ay, double by, double az, double bz, int nx, int ny, int nz) {
    double hx = (bx - ax) / nx; // Ширина трапеции по x
    double hy = (by - ay) / ny; // Ширина трапеции по y
    double hz = (bz - az) / nz; // Ширина трапеции по z

    double sum = 0.0;

    for (int i = 0; i <= nx; ++i) {
        for (int j = 0; j <= ny; ++j) {
            for (int k = 0; k <= nz; ++k) {
                double x = ax + i * hx;
                double y = ay + j * hy;
                double z = az + k * hz;

                double weight = 1.0;
                if (i == 0 || i == nx) weight *= 0.5;
                if (j == 0 || j == ny) weight *= 0.5;
                if (k == 0 || k == nz) weight *= 0.5;

                sum += weight * f(x, y, z);
            }
        }
    }

    return sum * hx * hy * hz;
}

int main() {
    double a = 0; // Начало интервала
    double b = 1; // Конец интервала
    int n = 1000; // Количество трапеций

    double result = integrate(a, b, n);
    std::cout << "Результат интегрирования: " << result << std::endl;

    double ax = 0, bx = 1; // Границы интервала по x
    double ay = 0, by = 1; // Границы интервала по y
    double az = 0, bz = 1; // Границы интервала по z
    int nx = 10, ny = 10, nz = 10; // Количество трапеций по каждой оси

    result = integrateVolume(ax, bx, ay, by, az, bz, nx, ny, nz);
    std::cout << "Результат интегрирования по объему: " << result << std::endl;

    return 0;
}