package ru.shvedovskiy.snau;

import ru.shvedovskiy.lib.Matrix;
import ru.shvedovskiy.lup.LUP;

public class Newtons {
    private static double EPS = 10e-9;
    private static Matrix J;
    private static Matrix J2;
    private static Matrix functions;
    private static Matrix functions2;
    private static Matrix X;
    private static Matrix X2;

    private static void setJakobiMatrix(double x1, double x2, double x3, double x4, double x5,
                                double x6, double x7, double x8, double x9, double x10) {
        J.arr[0][0] = -Math.sin(x1 * x2) * x2;
        J.arr[0][1] = -Math.sin(x1 * x2) * x1;
        J.arr[0][2] = 3 * Math.exp(-3 * x3);
        J.arr[0][3] = x5 * x5;
        J.arr[0][4] = 2 * x4 * x5;
        J.arr[0][5] = -1;
        J.arr[0][6] = 0;
        J.arr[0][7] = -2 * Math.cosh(2 * x8) * x9;
        J.arr[0][8] = -Math.sinh(2 * x8);
        J.arr[0][9] = 2;
        J.arr[1][0] = Math.cos(x1 * x2) * x2;
        J.arr[1][1] = Math.cos(x1 * x2) * x1;
        J.arr[1][2] = x9 * x7;
        J.arr[1][3] = 0;
        J.arr[1][4] = 6 * x5;
        J.arr[1][5] = -Math.exp(-x10 + x6) -  x8 - 1;
        J.arr[1][6] = x3 * x9;
        J.arr[1][7] = -x6;
        J.arr[1][8] = x3 * x7;
        J.arr[1][9] = Math.exp(-x10 + x6);
        J.arr[2][0] = 1;
        J.arr[2][1] = -1;
        J.arr[2][2] = 1;
        J.arr[2][3] = -1;
        J.arr[2][4] = 1;
        J.arr[2][5] = -1;
        J.arr[2][6] = 1;
        J.arr[2][7] = -1;
        J.arr[2][8] = 1;
        J.arr[2][9] = -1;
        J.arr[3][0] = -x5 * Math.pow(x3 + x1, -2);
        J.arr[3][1] = -2 * Math.cos(x2 * x2) * x2;
        J.arr[3][2] = -x5 * Math.pow( x3 + x1, -2);
        J.arr[3][3] = -2 * Math.sin(-x9 + x4);
        J.arr[3][4] = Math.pow((x3 + x1), -1);
        J.arr[3][5] = 0;
        J.arr[3][6] = -2* Math.cos(x7 * x10) * Math.sin(x7 * x10) * x10;
        J.arr[3][7] = -1;
        J.arr[3][8] = 2 * Math.sin(-x9 + x4);
        J.arr[3][9] = -2 * Math.cos(x7 * x10) * Math.sin(x7 * x10) * x7;
        J.arr[4][0] = 2 * x8;
        J.arr[4][1] = -2 * Math.sin(x2);
        J.arr[4][2] = 2 * x8;
        J.arr[4][3] = Math.pow(-x9 + x4, -2);
        J.arr[4][4] = Math.cos(x5);
        J.arr[4][5] = x7 * Math.exp(-x7 * (-x10 + x6));
        J.arr[4][6] = -(x10 - x6) * Math.exp(-x7 * (-x10 + x6));
        J.arr[4][7] = (2 * x3) + 2 * x1;
        J.arr[4][8] = -Math.pow(-x9 + x4, -2);
        J.arr[4][9] = -x7 * Math.exp(-x7 * (-x10 + x6));
        J.arr[5][0] = Math.exp(x1 -  x4 - x9);
        J.arr[5][1] = -1.5 * Math.sin(3 * x10 * x2) * x10;
        J.arr[5][2] = -x6;
        J.arr[5][3] = -Math.exp(x1 - x4 - x9);
        J.arr[5][4] = 2 * x5 / x8;
        J.arr[5][5] = -x3;
        J.arr[5][6] = 0;
        J.arr[5][7] = -x5 * x5 * Math.pow( x8, -2);
        J.arr[5][8] = -Math.exp(x1 - x4 - x9);
        J.arr[5][9] = -1.5 * Math.sin(3 * x10 * x2) * x2;
        J.arr[6][0] = Math.cos(x4);
        J.arr[6][1] = 3 * x2 * x2 * x7;
        J.arr[6][2] = 1;
        J.arr[6][3] = -(x1 - x6) * Math.sin(x4);
        J.arr[6][4] = Math.cos(x10 / x5 + x8) * x10 * Math.pow(x5, -2);
        J.arr[6][5] = -Math.cos( x4);
        J.arr[6][6] = Math.pow(x2, 3);
        J.arr[6][7] = -Math.cos(x10 / x5 + x8);
        J.arr[6][8] = 0;
        J.arr[6][9] = -Math.cos(x10 / x5 + x8) / x5;
        J.arr[7][0] = 2 * x5 * (x1 - 2 * x6);
        J.arr[7][1] = -x7 * Math.exp(x2 * x7 + x10);
        J.arr[7][2] = -2 * Math.cos(-x9 + x3);
        J.arr[7][3] = 1.5;
        J.arr[7][4] = Math.pow(x1 - 2 * x6, 2);
        J.arr[7][5] = -4 * x5 * (x1 - 2 * x6);
        J.arr[7][6] = -x2 * Math.exp(x2 * x7 + x10);
        J.arr[7][7] = 0;
        J.arr[7][8] = 2 * Math.cos(-x9 + x3);
        J.arr[7][9] = -Math.exp(x2 * x7 + x10);
        J.arr[8][0] = -3;
        J.arr[8][1] = -2 * x8 * x10 * x7;
        J.arr[8][2] = 0;
        J.arr[8][3] = Math.exp(x5 + x4);
        J.arr[8][4] = Math.exp(x5 + x4);
        J.arr[8][5] = -7 * Math.pow(x6, -2);
        J.arr[8][6] = -2 * x2 * x8 * x10;
        J.arr[8][7] = -2 * x2 * x10 * x7;
        J.arr[8][8] = 3;
        J.arr[8][9] = -2 * x2 * x8 * x7;
        J.arr[9][0] = x10;
        J.arr[9][1] = x9;
        J.arr[9][2] = -x8;
        J.arr[9][3] = Math.cos(x4 + x5 + x6) * x7;
        J.arr[9][4] = Math.cos(x4 + x5 + x6) * x7;
        J.arr[9][5] = Math.cos(x4 + x5 + x6) * x7;
        J.arr[9][6] = Math.sin(x4 + x5 + x6);
        J.arr[9][7] = -x3;
        J.arr[9][8] = x2;
        J.arr[9][9] = x1;
    }
    private static void setFunctionsMatrix(int n) {
        functions = new Matrix(n, 1);
        functions.setZeroMatrix();

        J = new Matrix(n, n);
        X = new Matrix(n, 1);
    }
    private static void setStartedMeanings(double x1, double x2, double x3, double x4, double x5,
                                   double x6,double x7,double x8,double x9, double x10) {
        X.arr[0][0] = x1;
        X.arr[1][0] = x2;
        X.arr[2][0] = x3;
        X.arr[3][0] = x4;
        X.arr[4][0] = x5;
        X.arr[5][0] = x6;
        X.arr[6][0] = x7;
        X.arr[7][0] = x8;
        X.arr[8][0] = x9;
        X.arr[9][0] = x10;
    }
    private static void setJakobiMatrix2(double x, double y) {
        J2.arr[0][0] = Math.cos(x) + 1;
        J2.arr[0][1] = 2;
        J2.arr[1][0] = -1;
        J2.arr[1][1] = -(Math.sin(y));
    }
    private static void setFunctionsMatrix2(int n) {
        functions2 = new Matrix(n, 1);
        functions2.setZeroMatrix();

        J2 = new Matrix(n, n);
        X2 = new Matrix(n, 1);
    }
    private static void setStartedMeanings2(double x0, double y0) {
        X2.arr[0][0] = x0;
        X2.arr[1][0] = y0;
    }

    private static void findMeaning(double x1, double x2, double x3, double x4, double x5,
                            double x6,double x7,double x8,double x9, double x10) {
        functions.arr[0][0] = -(Math.cos(x1 * x2) - Math.exp(-3 * x3)
                + x4 * x5 * x5
                - x6 - Math.sinh(2 * x8) * x9 + 2 * x10 + 2.0004339741653854440 );

        functions.arr[1][0] = -(Math.sin(x1*x2) + x3 * x9 * x7
                - Math.exp(-x10 + x6) + 3 * x5
                * x5 - x6 * (x8 + 1) + 10.886272036407019994);

        functions.arr[2][0] = -(x1 - x2 + x3
                - x4 + x5 - x6
                + x7 - x8 + x9 - x10 - 3.1361904761904761904);

        functions.arr[3][0] = -(2 * Math.cos(-x9 + x4) + x5 / (x3 + x1)
                - Math.sin(x2 * x2) + Math.pow(Math.cos(x7 * x10), 2)
                - x8 - 0.1707472705022304757 );

        functions.arr[4][0] = -(Math.sin(x5) + 2 * x8 * (x3 + x1)
                - Math.exp(-x7 * (-x10 + x6)) + 2 * Math.cos(x2)
                - 1 / (x4 - x9) - 0.3685896273101277862);

        functions.arr[5][0] = -(Math.exp(x1 - x4 - x9) +
                x5 * x5 / x8 + 0.5 * Math.cos(3 * x10 * x2)
                - x6 * x3 + 2.0491086016771875115 );

        functions.arr[6][0] = -(x2 * x2 * x2 * x7
                - Math.sin(x10 / x5 + x8)
                + (x1 - x6) * Math.cos(x4) + x3 - 0.7380430076202798014);

        functions.arr[7][0] = -(x5 * Math.pow((x1 - 2 * x6), 2)
                - 2 * Math.sin(-x9 + x3)
                + 1.5 * x4 - Math.exp(x2 * x7 + x10) + 3.5668321989693809040);

        functions.arr[8][0] = -(7 / x6 + Math.exp(x5 + x4)
                - 2 * x2 * x8 * x7 * x10
                + 3 * x9 - 3 * x1 - 8.4394734508383257499);

        functions.arr[9][0] = -(x10 * x1 + x9 * x2
                - x8 * x3
                + Math.sin(x4+x5+x6) * x7 - 0.78238095238095238096);
    }
    private static void findMeaning(Matrix X) {
        findMeaning(X.arr[0][0], X.arr[1][0], X.arr[2][0], X.arr[3][0],
                X.arr[4][0], X.arr[5][0], X.arr[6][0],
                X.arr[7][0], X.arr[8][0], X.arr[9][0]);
    }
    private static void findJakobi(Matrix X) {
        setJakobiMatrix(X.arr[0][0], X.arr[1][0], X.arr[2][0], X.arr[3][0],
                X.arr[4][0], X.arr[5][0], X.arr[6][0],
                X.arr[7][0], X.arr[8][0], X.arr[9][0]);
    }
    private static void findMeaning2(double x, double y) {
        functions2.arr[0][0] = -(Math.sin(X2.arr[0][0] + 1) - X2.arr[1][0] - 1);
        functions2.arr[1][0] = -((2 * X2.arr[0][0]) + Math.cos(X2.arr[1][0]) - 2);
    }
    private static void findMeaning2(Matrix X2) {
        findMeaning2(X2.arr[0][0], X2.arr[1][0]);
    }
    private static void findJakobi2(Matrix X2) {
        setJakobiMatrix2(X2.arr[0][0], X2.arr[1][0]);
    }

    private static void SNAUDecision(boolean isModified, int index) {
        int numberOfSteps = 0;
        int k = 0;
        long time = System.currentTimeMillis();
        Matrix L = null;
        Matrix U = null;
        Matrix P = null;
        Matrix prev_X_matrix = new Matrix(X);
        Matrix res = new Matrix(X);

        prev_X_matrix.copyMatrix(X);
        res.copyMatrix(X);

        double maxDim = X.maxInColumn(0);

        while(Math.abs(maxDim) > EPS) {
            findMeaning(X);
            if (k % index == 0) {
                findJakobi(X);
                U = new Matrix(J);
                L = new Matrix(J);
                P = new Matrix(J);
                LUP.LUExpansion(J, functions, L, U, P);
                J.determinant = U.determinantOfTriangleMatrix();
            } else {
                functions = P.add(functions);
            }

            res = LUP.decisionSLAU(J, functions, L, U);
            prev_X_matrix.copyMatrix(X);
            X.add(res);
            maxDim = (X.maxInColumn(0) - prev_X_matrix.maxInColumn(0)) / prev_X_matrix.maxInColumn(0);
            k++;
            numberOfSteps++;
            if (numberOfSteps >= 1000) {
                break;
            }
        }
        System.out.println("Задание 2.");
        System.out.println("Решение систем нелинейных алгебраических уравнений: ");
        System.out.println("( 1 )");
        X.printMatrix();
        System.out.println("Число шагов: " + numberOfSteps);
        System.out.println("Время в мс:  " + (System.currentTimeMillis() - time));
        System.out.println();
    }
    private static void SNAUDecision2(boolean isModified, int index) {
        int numberOfSteps = 0;
        int k = 0;
        long time = System.currentTimeMillis();
        Matrix L = null;
        Matrix U = null;
        Matrix P = null;
        Matrix prev_X_matrix = new Matrix(X2);
        Matrix res = new Matrix(X2);

        prev_X_matrix.copyMatrix(X2);
        res.copyMatrix(X2);

        double maxDim = X2.maxInColumn(0);

        while(Math.abs(maxDim) > EPS) {
            findMeaning2(X2);
            if (k % index == 0) {
                findJakobi2(X2);
                U = new Matrix(J2);
                L = new Matrix(J2);
                P = new Matrix(J2);
                LUP.LUExpansion(J2 ,functions2, L, U, P);
                J2.determinant = U.determinantOfTriangleMatrix();
            } else {
                functions2 = P.add(functions2);
            }

            res = LUP.decisionSLAU(J2, functions2, L, U);
            prev_X_matrix.copyMatrix(X2);
            X2.add(res);
            maxDim = (X2.maxInColumn(0) - prev_X_matrix.maxInColumn(0)) / prev_X_matrix.maxInColumn(0);
            k++;
            numberOfSteps++;

            if (numberOfSteps >= 1000) {
                break;
            }
        }
        System.out.println("( 2 )");
        X2.printMatrix();
        System.out.println("Число шагов: " + numberOfSteps);
        System.out.println("Время в мс:  " + (System.currentTimeMillis() - time));
        System.out.println();
    }

    public static void main(String args[]) {
        setFunctionsMatrix(10);
        setStartedMeanings(0.5, 0.5, 1.5, -1.0, -0.2, 1.5, 0.5, -0.5, 1.5, -1.5);
        SNAUDecision(true, 3);

        setFunctionsMatrix2(2);
        setStartedMeanings2(0.5, 0.0);
        SNAUDecision2(true, 3);
    }
}