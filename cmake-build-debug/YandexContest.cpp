//Bulat Akhmatov
#include <iostream>
#include <vector>
#include <iomanip>
#include <valarray>

using namespace std;

template<typename T>
class Matrix {

public:
    Matrix() = default;

    Matrix(int r, int c) {
        rows = r;
        cols = c;
        data.resize(r, vector<T>(c, 0));
    }

    Matrix(const Matrix &other) {
        rows = other.rows;
        cols = other.cols;
        data = other.data;
    }

    friend ostream &operator<<(ostream &out, const Matrix &mat) {
        for (int i = 0; i < mat.rows; i++) {
            for (int j = 0; j < mat.cols; j++) {
                if (abs(mat.data[i][j]) < 1e-10) {
                    out << "0.0000";
                } else {
                    out << fixed << setprecision(4) << mat.data[i][j];
                }
                if (j != mat.cols - 1) {
                    out << " ";
                }
            }
            out << endl;
        }
        return out;
    }

    friend istream &operator>>(istream &in, Matrix &mat) {
        for (int i = 0; i < mat.rows; i++) {
            for (int j = 0; j < mat.cols; j++) {
                in >> mat.data[i][j];
            }
        }
        return in;
    }


    Matrix &operator=(const Matrix &other) = default;

    Matrix operator+(const Matrix &other) {
        if (rows != other.rows || cols != other.cols) {
            cout << "Error: the dimensional problem occurred" << endl;
            return {0, 0};
        }
        Matrix result(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return result;
    }

    Matrix operator-(const Matrix &other) {
        if (rows != other.rows || cols != other.cols) {
            cout << "Error: the dimensional problem occurred" << endl;
            return {0, 0};
        }
        Matrix result(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.data[i][j] = data[i][j] - other.data[i][j];
            }
        }
        return result;
    }

    Matrix operator*(const Matrix &other) {
        if (cols != other.rows) {
            cout << "Error: the dimensional problem occurred" << endl;
            return {0, 0};
        }
        Matrix result(rows, other.cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < other.cols; j++) {
                for (int k = 0; k < cols; k++) {
                    result.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        return result;
    }

    Matrix transpose() {
        Matrix result(cols, rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result.data[j][i] = data[i][j];
            }
        }
        return result;
    }

    vector<T> &operator[](int i) {
        return data[i];
    }

    int rows{};
    int cols{};
    vector<vector<T>> data;
};

template<typename T>
class PermutationMatrix;

template<typename T>
class EliminationMatrix;

template<typename T>
class SquareMatrix : public Matrix<T> {
public:
    SquareMatrix() : Matrix<T>() {}

    explicit SquareMatrix(int n) : Matrix<T>(n, n) {}

    using Matrix<T>::Matrix;
    using Matrix<T>::operator=;

    SquareMatrix(SquareMatrix &mat) : Matrix<T>(mat) {
        (*this) = mat;
    }

    explicit SquareMatrix(Matrix<T> mat) : Matrix<T>(0, 0) {
        (*this) = mat;
    }

    int index_of_max(SquareMatrix &matrix, int i) {
        int max = i;
        for (int j = i + 1; j < matrix.rows; j++) {
            if (abs(matrix.data[j][i]) > abs(matrix.data[max][i])) {
                max = j;
            }
        }
        return max;
    }

    void determinant() {
        int stepNumber = 1;
        int permutationNumber = 0;
        SquareMatrix mat = *this;
        for (int i = 0; i < mat.rows; i++) {
            int pivot = index_of_max(mat, i);
            if (pivot != i) {
                mat = PermutationMatrix<T>(mat, pivot, i) * mat;
                cout << "step #" + to_string(stepNumber) + ": permutation" << endl;
                stepNumber++;
                permutationNumber++;
                cout << mat;
            }
            for (int j = i + 1; j < mat.rows; j++) {
                if (abs(mat.data[j][i]) < 1e-10) {
                    continue;
                }
                mat = EliminationMatrix<T>(mat, j, i) * mat;
                cout << "step #" + to_string(stepNumber) + ": elimination" << endl;
                stepNumber++;
                cout << mat;
            }
        }
        double det = 1.00;
        for (int i = 0; i < mat.rows; i++) {
            det *= mat.data[i][i];
        }
        if (permutationNumber % 2 == 1) {
            det *= -1;
        }
        cout << "result:" << endl;
        if (abs(det) < 1e-10) {
            cout << "0.00" << endl;
        } else {
            cout << fixed << setprecision(2) << det << endl;
        }
    }
};

template<typename T>
class IdentityMatrix : public SquareMatrix<T> {
public:
    explicit IdentityMatrix(int n) : SquareMatrix<T>(n) {
        for (int i = 0; i < n; i++) {
            this->data[i][i] = 1;
        }
    }
};

template<typename T>
class EliminationMatrix : public IdentityMatrix<T> {
public:
    EliminationMatrix(Matrix<T> &A, int n, int m) : IdentityMatrix<T>(A.rows) {
        (*this).data[n][m] = -A[n][m] / A[m][m];
    }
};

template<typename T>
class PermutationMatrix : public IdentityMatrix<T> {
public:
    PermutationMatrix(Matrix<T> &A, int n, int m) : IdentityMatrix<T>(A.rows) {
        swap(this->data[n], this->data[m]);
    }
};

template<typename T>
class ColumnVector : public Matrix<T> {
public:
    explicit ColumnVector(int n) : Matrix<T>(n, 1) {}

    explicit ColumnVector(const Matrix<T> &mat) : Matrix<T>(mat.rows, 1) {
        for (int i = 0; i < mat.rows; i++) {
            this->data[i][0] = mat.data[i][0];
        }
    }

    using Matrix<T>::operator=;
};

template<typename T>
class AugmentedMatrix : public Matrix<T> {
public:
    AugmentedMatrix(Matrix<T> A, Matrix<T> B) : Matrix<T>(A.rows, A.cols + B.cols) {
        for (int i = 0; i < A.rows; i++) {
            for (int j = 0; j < A.cols; j++) {
                this->data[i][j] = A[i][j];
            }
        }
        for (int i = 0; i < B.rows; i++) {
            for (int j = 0; j < B.cols; j++) {
                this->data[i][j + A.cols] = B[i][j];
            }
        }
    }

    int index_of_max(SquareMatrix<T> &matrix, int i) {
        int max = i;
        for (int j = i + 1; j < matrix.rows; j++) {
            if (abs(matrix.data[j][i]) > abs(matrix.data[max][i])) {
                max = j;
            }
        }
        return max;
    }

    Matrix<T> solveInverse() {
        int stepNumber = 1;
        AugmentedMatrix aug = *this;
        SquareMatrix<T> A(aug.rows);
        SquareMatrix<T> b(aug.rows);
        for (int i = 0; i < aug.rows; i++) {
            for (int j = 0; j < aug.cols; j++) {
                if (j < aug.cols / 2) {
                    A[i][j] = aug[i][j];
                } else {
                    b[i][j - aug.cols / 2] = aug[i][j];
                }
            }
        }
        for (int i = 0; i < A.rows; i++) {
            int pivot = index_of_max(A, i);
            if (pivot != i) {
                PermutationMatrix<T> P = PermutationMatrix<T>(A, pivot, i);
                A = P * A;
                b = P * b;
                stepNumber++;
            }
            for (int j = i + 1; j < A.rows; j++) {
                if (A.data[i][i] == 0 || A.data[j][i] == 0) {
                    continue;
                }
                EliminationMatrix<T> E = EliminationMatrix<T>(A, j, i);
                A = E * A;
                b = E * b;
                stepNumber++;
            }
        }
        for (int i = A.rows - 1; i >= 0; i--) {
            for (int j = i - 1; j >= 0; j--) {
                if (A.data[i][i] == 0) {
                    continue;
                }
                EliminationMatrix<T> E = EliminationMatrix<T>(A, j, i);
                A = E * A;
                b = E * b;
                stepNumber++;
            }
        }
        for (int i = 0; i < A.rows; i++) {
            double temp = A.data[i][i];
            for (int j = 0; j < A.cols; j++) {
                A.data[i][j] /= temp;
            }
            for (int j = 0; j < b.cols; j++) {
                b.data[i][j] /= temp;
            }
        }
        return b;
    }

    void solveSystem() {
        int stepNumber = 1;
        AugmentedMatrix aug = *this;
        SquareMatrix<T> A(aug.rows);
        ColumnVector<T> b(aug.rows);
        for (int i = 0; i < aug.rows; i++) {
            for (int j = 0; j < aug.cols; j++) {
                if (j < aug.cols - 1) {
                    A[i][j] = aug[i][j];
                } else {
                    b[i][j - aug.cols + 1] = aug[i][j];
                }
            }
        }
        for (int i = 0; i < A.rows; i++) {
            int pivot = index_of_max(A, i);
            if (pivot != i) {
                PermutationMatrix<T> P = PermutationMatrix<T>(A, pivot, i);
                A = P * A;
                b = P * b;
                stepNumber++;
            }
            for (int j = i + 1; j < A.rows; j++) {
                if (A.data[i][i] == 0 || A.data[j][i] == 0) {
                    continue;
                }
                EliminationMatrix<T> E = EliminationMatrix<T>(A, j, i);
                A = E * A;
                b = E * b;
                stepNumber++;
            }
        }
        for (int i = A.rows - 1; i >= 0; i--) {
            for (int j = i - 1; j >= 0; j--) {
                if (A.data[i][i] == 0) {
                    continue;
                }
                EliminationMatrix<T> E = EliminationMatrix<T>(A, j, i);
                A = E * A;
                b = E * b;
                stepNumber++;
            }
        }
        for (int i = 0; i < A.rows; i++) {
            double temp = A.data[i][i];
            for (int j = 0; j < A.cols; j++) {
                A.data[i][j] /= temp;
            }
            for (int j = 0; j < b.cols; j++) {
                b.data[i][j] /= temp;
            }
        }
        for (int i = 0; i < b.rows; i++) {
            for (int j = b.cols / 2; j < b.cols; j++) {
                if (abs(b.data[i][j]) < 1e-10) {
                    b.data[i][j] = 0;
                }
                printf("%.4f ", b.data[i][j]);
            }
            cout << endl;
        }
    }
};


int main() {
    int m;
    cin >> m;
    double t[m];
    double b[m];
    for (int i = 0; i < m; i++) {
        cin >> t[i] >> b[i];
    }
    int n;
    cin >> n;
    SquareMatrix<double> A(m, n + 1);
    ColumnVector<double> B(m);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n + 1; j++) {
            A[i][j] = pow(t[i], j);
        }
        B[i][0] = b[i];
    }
    cout << "A:" << endl;
    cout << A;
    cout << "A_T*A:" << endl;
    Matrix<double> AtA = A.transpose() * A;
    cout << AtA;
    cout << "(A_T*A)^-1:" << endl;
    Matrix<double> AtA_1 = AugmentedMatrix<double>(AtA, IdentityMatrix<double>(n + 1)).solveInverse();
    cout << AtA_1;
    cout << "A_T*b:" << endl;
    Matrix<double> AtB = A.transpose() * B;
    cout << AtB;
    cout << "x~:" << endl;
    cout << AtA_1 * AtB;
    return 0;
}