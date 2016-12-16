package ru.shvedovskiy.lup;

import ru.shvedovskiy.lib.Matrix;

public class LUP {
    public static void LUExpansion(Matrix A, Matrix B, Matrix L , Matrix U, Matrix P) {
        P.determinant = 1;

        for (int i = 0; i != A.m; ++i) {
            for (int j = 0; j != A.n; ++j) {
                if (i != j) {
                    P.arr[i][j] = 0;
                } else {
                    P.arr[i][j] = 1;
                }
            }
        }
        Matrix help = new Matrix(A);
        help.copyMatrix(A);
        int row = 0;

        for (int k = 0; k != A.m; ++k) {
            double p = 0;
            for (int i = k; i != A.m; ++i) {
                if (Math.abs(help.arr[i][k]) > p) {
                    p = Math.abs(help.arr[i][k]);
                    row = i;
                }
            }
            swap(help.arr[k], help.arr[row]);
            swap(P.arr[k], P.arr[row]);
            swap(A.arr[k], A.arr[row]);
            swap(B.arr[k], B.arr[row]);

            if (k != row) {
                P.determinant *= -1;
            }

            for (int i = k + 1; i != A.m; ++i) {
                if (help.arr[k][k] == 0) {
                    continue;
                }

                help.arr[i][k] /= help.arr[k][k];
                for (int j = k + 1; j != A.m; ++j) {
                    help.arr[i][j] -= help.arr[i][k] * help.arr[k][j];
                }
            }
        }
        for (int i = 0; i != A.m; ++i) {
            L.arr[i][i] = 1;
            for (int j = 0; j < A.n; ++j) {
                if (i <= j) {
                    if (Math.abs(help.arr[i][j]) < help.EPS) {
                        U.arr[i][j] = 0;
                    } else {
                        U.arr[i][j] = help.arr[i][j];
                    }
                } else {
                    if (Math.abs(help.arr[i][j]) < help.EPS) {
                        L.arr[i][j] = 0;
                    } else {
                        L.arr[i][j] = help.arr[i][j];
                    }
                }
            }
        }
    }
    private static void LUPDeterminant(Matrix A, Matrix L, Matrix U, Matrix P) {
        System.out.println("( 2 )");
        System.out.println("Матрица А после преобразований: ");
        A.printMatrix();

        System.out.println("Матрица P после преобразований: ");
        P.printMatrix();

        System.out.println("Матрица L: ");
        L.printMatrix();

        System.out.println("Матрица U: ");
        U.printMatrix();

        System.out.println("Проверка для разложения LUP:");
        A.subtract(L.add(U)).printMatrix();
        System.out.println();

        U.determinant = U.determinantOfTriangleMatrix();
        A.determinant = U.determinant / P.determinant;
        L.determinant = 1;

        System.out.println("( 3 )");
        System.out.println("Определитель A: " + A.determinant);
    }
    private static void LUPRank(Matrix A, Matrix U) {
        Matrix C = new Matrix(A);
        C.copyMatrix(A);

        A.rank = U.rankOfTriangleMatrix();
        U.rank = A.rank;

        System.out.println("Ранг A:         " + A.rank);
        System.out.println();
    }

    private static Matrix testForCompatibility(Matrix first, Matrix second) {
        Matrix tmp = new Matrix(first.m, first.n + second.n);
        tmp.extendedMatrix(first, second);

        tmp.deleteZeroColumns();
        tmp.rank = tmp.rankOfMatrix();

        Matrix X;
        if (tmp.numberOfDeletedColumns == 0) {
            X = new Matrix(tmp.n - 1, 1);
        } else {
            X = new Matrix(tmp.n - tmp.numberOfDeletedColumns, 1);
        }

        System.out.println("( 4 )");
        System.out.println("Ранг исходной матрицы: " + first.rank);
        System.out.println("Ранг расширенной матрицы: " + tmp.rank);

        tmp.rank = 5;
        if (first.rank == tmp.rank) {
            System.out.println("Имеется решение");
            int start = tmp.m - 1 - tmp.numberOfDeletedColumns;
            X.arr[start][0] = tmp.arr[start][tmp.n - 1] / tmp.arr[start][tmp.n - 2];

            if (tmp.arr[start][tmp.n - 2] == 0) {
                X.arr[tmp.m - 1][0] = 0;
            }
            for (int i = start - 1; i >= 0; --i) {
                for (int j = tmp.n - 2; j >= i; --j) {
                    tmp.arr[i][tmp.n - 1] -= X.arr[j][0] * tmp.arr[i][j];
                }
                if (tmp.arr[i][i] == 0) {
                    X.arr[i][0] = 0;
                } else {
                    X.arr[i][0] = tmp.arr[i][tmp.n - 1] / tmp.arr[i][i];
                }
            }
            X.getFullMatrix(tmp);
            System.out.println();
            return X;
        }
        else {
            System.out.println("( 4 )");
            System.out.println("Решение СЛАУ = (/)");
            System.out.println();
        }
        return null;
    }
    public static Matrix decisionSLAU(Matrix A, Matrix B, Matrix L, Matrix U) {
        Matrix X = new Matrix(L.n, 1);
        if (A.determinant == 0) {
            X = testForCompatibility(A, B);
        } else {
            int n = L.n;
            Matrix Y = new Matrix(L.n, 1);
            Y.arr[0][0] = B.arr[0][0] / L.arr[0][0];
            double sum;

            for (int i = 0; i != L.n; ++i) {
                sum = 0.;
                for (int j = 0; j < i; ++j) {
                    sum += Y.arr[j][0] * L.arr[i][j];
                }
                Y.arr[i][0] = (B.arr[i][0] - sum) / L.arr[i][i];
            }
            X.arr[n - 1][0] = Y.arr[n - 1][0] / U.arr[n - 1][n - 1];

            for (int i = n - 2; i >= 0; --i) {
                sum = 0.;
                for (int j = n - 1; j > i; --j) {
                    sum += X.arr[j][0] * U.arr[i][j];
                }
                X.arr[i][0] = (Y.arr[i][0] - sum) / U.arr[i][i];
            }
        }
        return X;
    }

    private static void swap(double[] arr1, double[] arr2) {
        for (int i = 0; i != arr1.length; ++i) {
            double tmp = arr1[i];
            arr1[i] = arr2[i];
            arr2[i] = tmp;
        }
    }
    private static Matrix reverseMatrix(Matrix A, Matrix L, Matrix U) {
        if (A.determinant == 0) {
            return null;
        }
        Matrix B = new Matrix(A.m, A.n);
        B.setEMatrix();
        Matrix res = new Matrix(A.m, A.n);
        Matrix toHelp = new Matrix(A.m, 1);
        Matrix X = new Matrix(A.m, 1);
        for (int j = 0; j != A.n; ++j) {
            toHelp = B.getColumn(j);
            X = decisionSLAU(A, toHelp, L, U);
            res.setColumn(X, j);
        }
        return res;
    }
    private static double numberOfConditionality(Matrix A, Matrix B)
    {
        return A.rate() * B.rate();
    }

    public static void main(String[] args) {
        Matrix A = new Matrix(9, 9);
        Matrix B = new Matrix(9, 1);

        A.createRandomMatrix();
        B.createRandomMatrix();

//        A.arr[0][0] = 3;
//        A.arr[0][1] = 2;
//        A.arr[0][2] = 1;
//        A.arr[1][0] = 0;
//        A.arr[1][1] = 0;
//        A.arr[1][2] = 2;
//        A.arr[2][0] = 0;
//        A.arr[2][1] = 0;
//        A.arr[2][2] = 3 ;
//
//        B.arr[0][0] = 5;
//        B.arr[1][0] = 4;
//        B.arr[2][0] = 6;

        Matrix C = new Matrix(A);
        C.copyMatrix(A);

        System.out.println("Задание 1.");
        System.out.println("Решение систем линейных уравнений методом разложения LU:");
        System.out.println("( 1 )");
        System.out.println("Матрица A: ");
        A.printMatrix();

        System.out.println("Матрица B: ");
        B.printMatrix();
        System.out.println();

        Matrix L = new Matrix(A);
        Matrix U = new Matrix(A);
        Matrix P = new Matrix(A);

        LUExpansion(A, B, L, U, P);
        LUPDeterminant(A, L, U, P);
        LUPRank(A, U);

        C.determinant = A.determinant;
        C.rank = A.rank;

        Matrix X = decisionSLAU(A, B, L, U);
        if (X != null) {
            System.out.println("( 4 )");
            System.out.println("Решение СЛАУ: ");
            X.printMatrix();
            System.out.println();
        }

        Matrix reverse_A_matrix = reverseMatrix(A, L ,U);
        if (reverse_A_matrix != null) {
            System.out.println("( 5 )");
            System.out.println("Матрица A^(-1): ");
            reverse_A_matrix.printMatrix();
            System.out.println("Проверка перемножением для перевёрнутой матрицы: ");
            A.add(reverse_A_matrix).printMatrix();
            System.out.println();

            System.out.println("( 6 )");
            System.out.println("Число обусловленности: ");
            System.out.println(numberOfConditionality(A, reverse_A_matrix));
            System.out.println();
        } else {
            System.out.println("( 5 )");
            System.out.println("Нет перевёрнутой матрицы: ");
            System.out.println();
        }
    }
}