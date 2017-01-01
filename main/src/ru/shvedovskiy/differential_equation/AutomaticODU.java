package ru.shvedovskiy.differential_equation;

public class AutomaticODU {
    private double y1, y2, Y1, Y2;
    private double A, B, C;
    private double[] y;

    private double f1(double x, double[] y) {
        return A * y[1];
    }
    private double f2 (double x, double[] y) {
        return -B * y[0];
    }

    public AutomaticODU(double C, double A, double B, double h0, double eps) {
        this.C = C;
        this.A = A;
        this.B = B;

        y = new double[2];

        double pi = Math.PI;
        double c2 = C;
        double a21, b1, b2;
        int n = 1;
        a21 = c2;
        b2 = 1. / (2 * c2);
        b1 = 1 - b2;

        y[0] = A * pi;
        y[1] = B * pi;

        double x0 = 0;
        double h = h0;

        while (x0 < pi) {
            double y1k1 = h * f1(x0, y);
            double y2k1 = h * f2(x0, y);

            double[] toHelp = new double [2];
            toHelp[0] = y[0] + a21 * y1k1;
            toHelp[1] = y[1] + a21 * y2k1;

            double y2k2 = h * f2(x0 + c2 * h, toHelp);
            double y1k2 = h * f1(x0 + c2 * h, toHelp);

            y1 = y[0] + b1 * y1k1 + b2 * y1k2;
            y2 = y[1] + b1 * y2k1 + b2 * y2k2;

            y1k1 = h / 2 * f1(x0, y);
            y2k1 = h / 2 * f1(x0, y);

            toHelp = new double [2];
            toHelp[0] = y[0] + a21 * y1k1;
            toHelp[1] = y[1] + a21 * y2k1;

            y2k2 = h / 2 * f2(x0 + c2 * h / 2, toHelp);
            y1k2 = h / 2 * f1(x0 + c2 * h / 2, toHelp);

            Y1 = y[0] + b1 * y1k1 + b2 * y1k2;
            Y2 = y[1] + b1 * y2k1 + b2 * y2k2;

            double s = 2;
            double localError = Math.abs(Math.sqrt(Math.pow(y1 - Y1, 2) + Math.pow(y2 - Y2, 2))) / 3;
            double e2s = eps * Math.pow(2, s);
            double eDiv2s = eps / (Math.pow(2, s + 1));
            double prevh = h;
            if (localError > e2s) {
                h /= 2.;
                System.out.println("На шаге " + n + " h уменьшилось");
            } else if (localError > eps && localError <= e2s) {
                h /= 2.;
                x0 += h;
                y[0] = Y1;
                y[1] = Y2;
                System.out.println("На шаге " + n + " h уменьшилось");
            } else if (localError < eps && localError >= eDiv2s) {
                y[0] = y1;
                y[1] = y2;
                x0 += h;
                //System.out.println("Шаг: " + h);
            } else if (localError < eDiv2s) {
                x0 += h;
                h *= 2.;
                y[0] = y1;
                y[1] = y2;
                System.out.println("На шаге " + n + " h увеличилось");
            }
			/*
			System.out.println("На шаге " + n + ":");
			System.out.println("Y1: " + y[0]);
			System.out.println("Y3: " + y[1]);
			System.out.println("Локальная ошибка: " + localError);
			System.out.println("H: " + prevh);
			System.out.println();
			*/
            ++n;
        }
        System.out.println();
    }
}
