package ru.shvedovskiy.differential_equation;

public class Runge {
    public Runge(double C, double A, double B, int n, double eps, boolean isOptimal) {
        double s = 2;
        double h = Math.PI / (double) n;
        double errorForY1 = 1000, errorForY2 = 1000;

        ODU UnionODU = new ODU(C, A, B, h, n, false);

        double y1H1 = UnionODU.lastY1();
        double y1H2;

        double y2H1 = UnionODU.lastY2();
        double y2H2;
        if (!isOptimal) {
            do {
                n *= 2;

                UnionODU = new ODU(C, A, B, Math.PI / n, n, false);

                y1H2 = y1H1;
                y1H1 = UnionODU.lastY1();

                errorForY1 = (y1H1 - y1H2) / 3;

                y2H2 = y2H1;
                y2H1 = UnionODU.lastY2();
                errorForY2 = (y2H1 - y2H2) / 3;

                System.out.println("Значение Y1: " + y1H1);
                System.out.println("Значение Y2: " + y2H1);
                System.out.println();
            } while (Math.abs(errorForY1) > eps || Math.abs(errorForY2) > eps);

            System.out.println("Ошибка для Y1: " + errorForY1);
            System.out.println("Ошибка для Y2: " + errorForY2);

            //System.out.println("Y1: " + y1H1);
            //System.out.println("Y2: " + y2H1);
            System.out.println("Число шагов: "+ n);
        } else {
            UnionODU = new ODU(C, A, B, Math.PI/ 2 * n, 2 * n, false);

            y1H2 = UnionODU.lastY1();
            y2H2 = UnionODU.lastY2();

            double norm = Math.abs(Math.sqrt(Math.pow(y1H1 - y1H2, 2) + Math.pow(y2H1 - y2H2, 2)));
            double optH = (Math.PI / n / 2) * Math.pow((Math.pow(2, s) - 1) * eps / norm, 1. / s);

            UnionODU = new ODU(C, A, B, optH, (int) Math.round(Math.PI / optH), false);

            ODU opponent = new ODU(C, A, B, optH, (int) Math.round(Math.PI / optH), true);


            System.out.println("Y1(не оппонент) для оптимального H: " + UnionODU.lastY1());
            System.out.println("Y2(not opponent) для оптимального H: " + UnionODU.lastY2());
            System.out.println();


            y1H2 = UnionODU.lastY1();
            y2H2 = UnionODU.lastY2();

            UnionODU = new ODU(C, A, B, optH / 2, (int) Math.round(2 * Math.PI / optH), false);

            y2H1 = UnionODU.lastY2();
            y1H1 = UnionODU.lastY1();

            norm = Math.abs(Math.sqrt(Math.pow(y1H1 - y1H2, 2) + Math.pow(y2H1 - y2H2, 2)));

            System.out.println("Ошибка(не оппонент) для оптимального H: " + norm / 3);
            System.out.println();

            System.out.println("Y1(оппонент) для оптимального H: " + opponent.lastY1());
            System.out.println("Y2(оппонент) для оптимального H: " + opponent.lastY2());
            System.out.println();

            y1H2 = opponent.lastY1();
            y2H2 = opponent.lastY2();

            opponent = new ODU(C, A, B, optH / 2, (int) Math.round(2 * Math.PI / optH), true);

            y2H1 = opponent.lastY2();
            y1H1 = opponent.lastY1();

            norm = Math.abs(Math.sqrt(Math.pow(y1H1 - y1H2, 2) + Math.pow(y2H1 - y2H2, 2)));

            System.out.println("Ошибка(оппонент) для оптимального H: " + norm / 3);
            System.out.println();
        }
    }
}
