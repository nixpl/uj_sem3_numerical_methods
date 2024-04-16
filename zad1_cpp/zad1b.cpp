#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

vector<vector<double>> choleskyDiag3(const vector<vector<double>> &A) // zwraca macierz C po faktoryzacji
{
    const int diag_count = 3;
    size_t size = A.size();
    vector<double> row(size, 0.0);
    vector<vector<double>> C(size, row);

    C[0][0] = sqrt(A[0][0]);

    for (int i = 0; i < (int)size; i++)
        for (int j = i; j < i + (diag_count + 1) / 2 && j < (int)size; j++)
        {
            if (j == i)
            {
                double buff = A[j][i];
                for (int k = 0; k < i; k++)
                    buff -= C[j][k] * C[j][k];
                C[j][i] = sqrt(buff);
            }
            else
                C[j][i] = A[j][i] / C[j - 1][i];
        }

    return C;
}

vector<double> solveLowerTriangularWith2Diag(const vector<vector<double>> &C, const vector<double> &b)
{
    size_t size = C.size();
    vector<double> x;
    x.push_back(b[0] / C[0][0]);
    for (int j = 1; j < (int)size; j++)
        x.push_back((b[j] - C[j][j - 1] * x[j - 1]) / C[j][j]);
    return x;
}

vector<double> solveLowerTriangularWith2DiagAsTransposed(const vector<vector<double>> &C, const vector<double> &b) // zamiast transponować macierz a dopiero potem liczyc to od razu liczy na macierzy nietransponowanej
{
    size_t size = C.size();
    vector<double> x(size);
    x[size - 1] = b[size - 1] / C[size - 1][size - 1];
    for (int i = size - 2; i >= 0; i--)
        x[i] = (b[i] - C[i + 1][i] * x[i + 1]) / C[i][i];
    return x;
}

vector<double> solveForCholesky(const vector<vector<double>> &C, const vector<double> b)
{
    vector<double> y = solveLowerTriangularWith2Diag(C, b);
    return solveLowerTriangularWith2DiagAsTransposed(C, y);
}

double vectorMultiplicationToNumber(vector<double> &vec1, vector<double> &vec2)
{
    size_t size = vec2.size();
    double num = 0;
    for (int i = 0; i < (int)size; i++)
        num += vec1[i] * vec2[i];
    return num;
}

vector<double> multiplyVectorByNumber(vector<double> &vec, const double num)
{
    size_t size = vec.size();
    vector<double> multiplied(size);
    for (int i = 0; i < (int)size; i++)
        multiplied[i] = vec[i] * num;
    return multiplied;
}

vector<double> vectorSubtraction(vector<double> &vec1, vector<double> &vec2)
{
    size_t size = vec1.size();
    vector<double> subtracted(size);
    for (int i = 0; i < (int)size; i++)
        subtracted[i] = vec1[i] - vec2[i];
    return subtracted;
}

vector<double> shermanMorrison(vector<double> &u, vector<double> &z, vector<double> &q)
{
    double buff = vectorMultiplicationToNumber(u, z) / (1 + vectorMultiplicationToNumber(u, q));
    vector<double> qMultipliedByBuff = multiplyVectorByNumber(q, buff);
    vector<double> w = vectorSubtraction(z, qMultipliedByBuff);
    return w;
}

void printMatrix(const vector<vector<double>> &matrix, int columnWidth)
{
    for (const vector<double> &row : matrix)
    {
        for (double value : row)
            cout << setw(columnWidth) << value;
        cout << endl;
    }
}

void printVector(const vector<double> &vec)
{
    for (const double &value : vec)
        cout << value << endl;
}

int main()
{
    vector<vector<double>> A{
        {3, 1, 0, 0, 0, 0, 0},
        {1, 4, 1, 0, 0, 0, 0},
        {0, 1, 4, 1, 0, 0, 0},
        {0, 0, 1, 4, 1, 0, 0},
        {0, 0, 0, 1, 4, 1, 0},
        {0, 0, 0, 0, 1, 4, 1},
        {0, 0, 0, 0, 0, 1, 3},
    };

    vector<double> b{1, 2, 3, 4, 5, 6, 7};

    vector<double> u{1, 0, 0, 0, 0, 0, 1};

    vector<vector<double>> C = choleskyDiag3(A);

    cout << "Macierz trujkatna dolna C po faktoryzacji Cholesky'ego macierzy trojdiagonalnej A:" << endl;
    printMatrix(C, 12);
    cout << endl;

    vector<double> z = solveForCholesky(C, b);

    vector<double> q = solveForCholesky(C, u);

    vector<double> w = shermanMorrison(u, z, q);

    cout << "Wektor bedacy rozwiazaniem ukladu rownan dla macierzy z przykladu b): " << endl;
    printVector(w);

    return 0;
}