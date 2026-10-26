#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>

using namespace std;
using namespace chrono;

// Функция для проверки симметричности матрицы
bool is_symmetric(const vector<vector<double>>& matrix) {
    // Проверяем размерность матрицы A
    if (matrix.size() != matrix[0].size()) {
        return false;
    }
    // Проверяем симметричность матрицы
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size(); j++) {
            if (matrix[i][j] != matrix[j][i]) {
                return false;
            }
        }
    }
    return true;
}

// Функция реализации метода вращений (Якоби)
vector<double> jacobi_method(vector<vector<double>> matrix, double tol = 1e-6, int max_iterations = 1000) {
    // Проверяем симметричность матрицы A перед применением метода Якоби
    if (!is_symmetric(matrix)) {
        throw runtime_error("Матрица не симметрична.");
    }
    int n = matrix.size();
    vector<double> eigenvalues(n, 0.0);
    for (int iter = 0; iter < max_iterations; iter++) {
        // Находим индекс наибольшего по модулю элемента вне главной диагонали матрицы A
        int max_idx = -1;
        double max_value = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double abs_value = abs(matrix[i][j]);
                if (abs_value > max_value) {
                    max_value = abs_value;
                    max_idx = i * n + j;
                }
            }
        }
        if (max_idx == -1) {
            // Если все внедиагональные элементы уже близки к нулю, выходим из цикла
            break;
        }
        int p = max_idx / n;  // Индекс строки
        int q = max_idx % n;  // Индекс столбца
        double theta;
        if (matrix[p][p] == matrix[q][q]) {
            // Если элементы на диагонали p и q равны, устанавливаем угол theta равным pi/4
            theta = M_PI / 4;
        }
        else {
            // Иначе вычисляем угол theta по формуле
            theta = 0.5 * atan(2 * matrix[p][q] / (matrix[p][p] - matrix[q][q]));
        }
        // Вычисляем значения cos(theta) и sin(theta)
        double c = cos(theta);
        double s = sin(theta);
        // Создаем матрицу поворота J
        vector<vector<double>> J(n, vector<double>(n, 0.0));
        for (int i = 0; i < n; i++) {
            J[i][i] = 1.0;
        }
        J[p][p] = c;
        J[q][q] = c;
        J[p][q] = -s;
        J[q][p] = s;
        // Применяем матрицу поворота J к матрице A
        vector<vector<double>> tempMatrix(n, vector<double>(n, 0.0));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    tempMatrix[i][j] += J[k][i] * matrix[k][j];
                }
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                matrix[i][j] = 0.0;
                for (int k = 0; k < n; k++) {
                    matrix[i][j] += tempMatrix[i][k] * J[k][j];
                }
            }
        }
        // Вычисляем величину внедиагональных элементов матрицы A
        double max_off_diag = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double off_diag = abs(matrix[i][j]);
                if (off_diag > max_off_diag) {
                    max_off_diag = off_diag;
                }
            }
        }
        if (max_off_diag < tol) {
            // Если величина внедиагональных элементов меньше заданной погрешности, выходим из цикла
            break;
        }
    }
    // Извлекаем собственные значения из диагонали матрицы A
    for (int i = 0; i < n; i++) {
        eigenvalues[i] = matrix[i][i];
    }
    return eigenvalues;
}

// Функция для генерации произвольной симметричной матрицы 
vector<vector<double>> generateRandomMatrix(int size, double minVal, double maxVal) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(minVal, maxVal);
    vector<vector<double>> matrix(size, vector<double>(size));
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j <= i; ++j) {
            double value = dis(gen);
            matrix[i][j] = value;
            matrix[j][i] = value;
        }
    }
    return matrix;
}

// Функция для исследования (измерения) скорости сходимости в зависимости от заданной точности
vector<double> measureConvergenceSpeedOfJacobiMethod(vector<vector<double>>& matrix, const vector<double>& tolerances) {
    vector<double> eigenvalues;
    for (double tolerance : tolerances) {
        // Запускаем таймер
        auto start = high_resolution_clock::now();
        // Выполняем метод Якоби
        eigenvalues = jacobi_method(matrix, tolerance);
        // Завершаем таймер и выводим время выполнения
        auto end = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(end - start);
        cout << "Точность: " << tolerance << ", Время выполнения: " << duration.count() << " мс" << endl;
    }
    return eigenvalues;
}

int main() {
    // Разные значения точности для измерения 
    vector<double> tolerances = {1e-2, 1e-4, 1e-6, 1e-8};
    // Определяем матрицу A
    vector<vector<double>> A = {{1, 1.5, 3.5, 4.5}, {1.5, 0.5, 2, 1.6}, {3.5, 2, 2, 1.7}, {4.5, 1.6, 1.7, 2}};
    cout << "Матрица A:" << endl;
    for (const auto& row : A) {
        for (const auto& element : row) {
            cout << element << " ";
        }
        cout << endl;
    }
    vector<double> eigenvalues = measureConvergenceSpeedOfJacobiMethod(A, tolerances);
    for (double eigenvalue : eigenvalues) {
        // Выводим каждое собственное значение
        cout << "Собственное значение: " << eigenvalue << "\n\n";
    }
    vector<vector<double>> randomMatrix = generateRandomMatrix(4, 0.0, 10.0);
    cout << "Произвольная матрица:" << endl;
    for (const auto& row : randomMatrix) {
        for (const auto& element : row) {
            cout << element << " ";
        }
        cout << endl;
    }
    eigenvalues = measureConvergenceSpeedOfJacobiMethod(randomMatrix, tolerances);
    for (double eigenvalue : eigenvalues) {
        cout << "Собственное значение: " << eigenvalue << "\n\n";
    }
    vector<vector<double>> randomMatrixSize11 = generateRandomMatrix(11, 0.0, 10.0);
    cout << "Матрица размерности 11:" << endl;
    for (const auto& row : randomMatrixSize11) {
        for (const auto& element : row) {
            cout << element << " ";
        }
        cout << endl;
    }
    eigenvalues = measureConvergenceSpeedOfJacobiMethod(randomMatrixSize11, tolerances);
    for (double eigenvalue : eigenvalues) {
        cout << "Собственное значение: " << eigenvalue << "\n\n";
    }
    return 0;
}
