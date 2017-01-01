package ru.shvedovskiy.differential_equation;

public class ODU {
    private double A, B, C;
    private double y1, y2;
    private double x0;

    private double f1(double x, double[] y) {
        return A * y[1];
    }

    private double f2 (double x, double[] y) {
        return -B * y[0];
    }

    public ODU(double C, double A, double B, double h, int n, boolean isOpponent) {
        this.C = C;
        this.A = A;
        this.B = B;
        if (!isOpponent) {
            double pi = Math.PI;
            double c2 = C;
            double a21, b1, b2;

            a21 = c2;
            b2 = 1. / ( 2 * c2 );
            b1 = 1 - b2;

            y1 = A * pi;
            y2 = B * pi;
            for (int i = 1; i != n; ++i) {
                double [] Y = new double[2];
                Y[0] = y1;
                Y[1] = y2;

                double y1k1 = h * f1(x0, Y);
                double y2k1 = h * f2(x0, Y);

                double[] toHelp = new double[2];
                toHelp[0] = Y[0] + a21 * y1k1;
                toHelp[1] = Y[1] + a21 * y2k1;

                double y2k2 = h * f2(x0 + c2 * h, toHelp);
                double y1k2 = h * f1(x0 + c2 * h, toHelp);

                y1 = Y[0] + b1 * y1k1 + b2 * y1k2;
                y2 = Y[1] + b1 * y2k1 + b2 * y2k2;
            }
        } else {
            double pi = Math.PI;
            y1 = A * pi;
            y2 = B * pi;
            for (int i = 1; i != n; ++i) {
                double [] Y = new double[2];
                Y[0] = y1;
                Y[1] = y2;

                double y1k1 = h * f1(x0, Y);
                double y2k1 = h * f2(x0, Y);

                double[] toHelp = new double[2];
                toHelp[0] = Y[0] + y1k1;
                toHelp[1] = Y[1] + y2k1;

                double y2k2 = h * f2(x0 + h, toHelp);
                double y1k2 = h * f1(x0 + h, toHelp);

                y1 = Y[0] + 1. / 2. * (y1k1 + y1k2);
                y2 = Y[1] + 1. / 2. * (y2k1 + y2k2);
            }
        }
    }

    public void print() {
        System.out.println("Y1 в пи: " + y1);
        System.out.println("Y2 в пи: " + y2);
    }

    public double lastY1() {
        return y1;
    }

    public double  lastY2 () {
        return y2;
    }
}
