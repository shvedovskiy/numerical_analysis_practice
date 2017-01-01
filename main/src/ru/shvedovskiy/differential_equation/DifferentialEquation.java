package ru.shvedovskiy.differential_equation;

public class DifferentialEquation {
    public static void main(String[] args) {
        double eps = 1e-4;
        double A = 5. / 14., B = 7. / 15., c2 = 9. / 14.;
        int step = 5;
        System.out.println("Задание 4. Решение обыкновенных дифференциальных уравнений методом Рунге-Кутты\n");
        double h0 = new StartStep(A, B, eps).h;
        System.out.println("Значение h0: " + h0 + "\n");
        new AutomaticODU(c2, A, B, h0, eps);
        new Runge(c2, A, B, step, eps, false);
        new ODU(c2, A, B, Math.PI / step, step, true).print();
    }
}

class StartStep {
    double h;
    private double A, B;

    private double f1(double x, double[] y) {
        return A * y[1];
    }

    private double f2 (double x, double[] y) {
        return -B * y[0];
    }

    StartStep(double A, double B, double eps) {
        double s = 2;
        this.A = A;
        this.B = B;

        double[] y = new double[2];
        double[] f = new double[2];

        y[0] = A * Math.PI;
        y[1] = B * Math.PI;

        f[0] = f1(0, y);
        f[1] = f2(0, y);

        double norm = Math.abs(f[0] - f[1]);

        double delta = Math.pow(1. / Math.PI, s + 1) + Math.pow(norm, s + 1);
        h = Math.pow((eps / delta), 1. / (s + 1));
    }
}