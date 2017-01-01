package ru.shvedovskiy.integral;

public class Runge {
	private double anIntegral = 78.395972369586;

	public Runge(double left, double right, boolean isNewton, double eps, boolean withEitken) {
		if (!withEitken) {
			int n = 2;
			double p;
			if (isNewton) {
				p = 1. / 7.;
			} else {
				p = 1. / 63.;
			}
			double prevInt, currInt, err;

			currInt = new UnionIntegral(left, right, n, isNewton).getResult();
			do {
				n *= 2;
				prevInt = currInt;
				currInt = new UnionIntegral(left, right, n, isNewton).getResult();
				err = Math.abs(prevInt - currInt) * p;
			} while (err > eps);

			System.out.println("( 2 )");
            System.out.println("Оценка погрешности: ");
			System.out.println("Результат: " + currInt);
			System.out.println("интеграл в Wolfram - полученный интеграл: " + (anIntegral - currInt));
			System.out.println("Погрешность: " + err);
			System.out.println("Число шагов: " + n + "\n");
		} else {
			int H1 = 8; // шаг разбиения
			double SH1, SH2, SH3, err;
			SH1 = new UnionIntegral(left, right, H1 / 4, isNewton).getResult();
			SH2 = new UnionIntegral(left, right, H1 / 2, isNewton).getResult();
			SH3 = new UnionIntegral(left, right, H1, isNewton).getResult();

			do {
				double m = new Eitken(SH1, SH2, SH3).getM();
				System.out.println("Погрешность по Эйткену M:" + m);
				
				if (Math.abs(m - 3) > 1) {
                    m = 4;
                }
				
				err = (SH3 - SH2) / (Math.pow(2, m) - 1);
				H1 = (int) (H1 * Math.pow(eps * Math.abs(1 - Math.pow(2, -m)) / (SH3-SH2) , -1. / m));
				H1 += 4 - H1 % 4;
				System.out.println("H1: " + H1);
				SH1 = new UnionIntegral(left, right, H1 / 4, isNewton).getResult();
				SH2 = new UnionIntegral(left, right, H1 / 2, isNewton).getResult();
				SH3 = new UnionIntegral(left, right, H1, isNewton).getResult();
			} while (err > eps);
                        
            System.out.println("( 2 )");
            System.out.println("Оценка погрешности: ");
			System.out.println("Результат: " + SH3);
			System.out.println("интеграл в Wolfram - полученный интеграл: " + (anIntegral - SH3));
			System.out.println("Погрешность: " + err);
			System.out.println("Число шагов: " + H1 + "\n");
		}
	}
}