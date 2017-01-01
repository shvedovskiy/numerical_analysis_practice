package ru.shvedovskiy.integral;

import ru.shvedovskiy.lib.Matrix;
import ru.shvedovskiy.lup.LUP;
import static java.lang.Math.*;

public class Gauss  {
	private double[] moments;
	private double res = 0;
	private Matrix a, b, x; // узлы
	private Matrix A; // коэффициенты
        
	private void setMomentsAndMatrix(double left, double right) {
	    moments = new double[6];
        /*
		moments[0] = -5. / 3. * pow((2.9  - right), 3. / 5.) + 5. / 3. * pow((2.9  - left), 3. / 5.);
		moments[1] = -5. / 48. * pow((29. / 10. - right), (3. / 5.)) * (6 * right + 29) +
		    5. / 48. * pow((29. / 10. - left), (3. / 5.)) * (6 * left + 29);
		moments[2] = -(pow(5, 2. / 5.) * pow(29. / 2. - 5 * right, 3. / 5.) * (48 * right * right + 174 * right + 841)) / 624. +
		    (pow(5, 2. / 5.) * pow(29. / 2. - 5 * left, 3. / 5.) * (48 * left * left + 174 * left +841)) / 624.;
		moments[3] = -(pow(5, 2. / 5.) * pow(29. / 2. - 5 * right, 3. / 5.) *
			(416 * pow(right, 3) + 1392 * pow(right, 2) +  5046 * right + 24389)) / 7488. +
			(pow(5, 2. / 5.) * pow(29. / 2. - 5 * left, 3. / 5.) *
			(416 * pow(left, 3) + 1392 * pow(left, 2) +  5046 * left + 24389)) / 7488.;
		moments[4] = -(pow(5, 2. / 5.) * pow(29. / 2. - 5 * right, 3. / 5.) * (3744 * pow(right, 4) +
			12064 * pow(right, 3) + 40368 * pow(right, 2) + 146334 * right + 707281)) / 86112. +
			(pow(5, 2. / 5.) * pow(29. / 2. - 5 * left, 3. / 5.) * (3744 * pow(left, 4) +
			12064 * pow(left, 3) + 40368 * pow(left, 2) + 146334 * left + 707281)) / 86112.;
		moments[5] = -(pow(5, 2. / 5.) * pow(29. / 2. - 5 * right, 3. / 5.) *
			(172224 * pow(right, 5) + 542880 * pow(right, 4) + 1749280 * pow(right, 3) + 5853360 * pow(right, 2) + 21218430 * right + 102555745)) / 4822272. +
			(pow(5, 2. / 5.) * pow(29. / 2. - 5 * left, 3. / 5.) * (172224 * pow(left, 5) +
			542880 * pow(left, 4) + 1749280 * pow(left, 3) + 5853360 * pow(left, 2) + 21218430 * left + 102555745)) / 4822272.;
		*/
        moments[0] = -4. * exp(log(16. / 5. - right) / 4.) + 4. * exp(log(16. / 5. - left) / 4.);
        moments[1] = -(4. / 25.) * exp(log(16. / 5. - right) / 4.) * ((5. * right) + 64.)
            + (4. / 25.) * exp(log(16. / 5. - left) / 4.) * ((5. * left) + 64.);
        moments[2] = -(4. / 1125.) * exp(log(16. / 5. - right) / 4.) * ((125. * pow(right, 2)) + (640 * right) + 8192)
            + (4. / 1125.) * exp(log(16. / 5. - left) / 4.) * ((125. * pow(left, 2)) + (640 * left) + 8192);
        moments[3] = -(4. / 24375.) * exp(log(16. / 5. - right) / 4.) * ((1875. * pow(right, 3)) + (8000. * pow(right, 2)) + (40960. * right) + 524288.)
            + (4. / 24375.) * exp(log(16. / 5. - right) / 4.) * ((1875. * pow(right, 3)) + (8000. * pow(right, 2)) + (40960. * right) + 524288.);
        moments[4] = -(4. / 2071875.) * exp(log(16. / 5. - right)/4.) * ((121875. * pow(right, 4)) + (480000. * pow(right, 3)) + (2048000. * pow(right, 2)) + (10485760.  * right) + 134217728.)
            + (4. / 2071875.) * exp(log(16. / 5. - left) / 4.) * ((121875. * pow(left, 4)) + (480000. * pow(left, 3)) + (2048000. * pow(left, 2)) + (10485760. * left) + 134217728.);
        moments[5] = -(4. / 43509375.) * exp(log(16. / 5. - right) / 4.) * ((2071875. * pow(right, 5)) + (7800000. * pow(right, 4)) + (30720000. * pow(right, 3)) + (131072000. * pow(right, 2)) + (671088640. * right) + 8589934592.)
            + (4. / 43509375.) * exp(log(16. / 5. - left) / 4.) * ((2071875. * pow(left, 5)) + (7800000. * pow(left, 4)) + (30720000. * pow(left, 3)) + (131072000. * pow(left, 2)) + (671088640. * left) + 8589934592.);
        
		A = new Matrix(3,3);
		A.arr[0][0] = moments[0];
		A.arr[0][1] = moments[1];
		A.arr[0][2] = moments[2];
		A.arr[1][0] = moments[1];
		A.arr[1][1] = moments[2];
		A.arr[1][2] = moments[3];
		A.arr[2][0] = moments[2];
		A.arr[2][1] = moments[3];
		A.arr[2][2] = moments[4];
		
		b = new Matrix(3, 1);
		b.arr[0][0]= -moments[3];
		b.arr[1][0]= -moments[4];
		b.arr[2][0]= -moments[5];
		
		a = new Matrix(3, 1);
		x = new Matrix(3, 1);
	}
	
	private void step2() { //сумма (a_i * moment_i+s) = -moment_n+s
		Matrix L, U, P;
		L = new Matrix(A);
		U = new Matrix(A);
		P = new Matrix(A);

        LUP.LUExpansion(A, b, L, U, P);
		A.determinant = U.determinantOfTriangleMatrix();
        a = LUP.decisionSLAU(A, b, L, U);
	}

	private double function(double mean) { // нахождение корней узлового многочлена (для 3 шага)
		return pow(mean, 3) + a.arr[2][0] * pow(mean, 2) + a.arr[1][0] * mean + a.arr[0][0];
	}

	private void findOtherSolutions(double solution) {
		a.arr[0][0] = a.arr[2][0] * solution + pow(solution, 2) + a.arr[1][0];
		a.arr[1][0] = a.arr[2][0] + solution;
		a.arr[2][0] = 1;
		
		double disk = pow(a.arr[1][0], 2) - 4 * a.arr[0][0];
		
		if (disk < 0) {
			System.out.println("Комплексное решение");
			System.exit(0);
		}
		double solution1 = (-a.arr[1][0] + pow(disk, 1. / 2.)) / 2;
		double solution2 = (-a.arr[1][0] - pow(disk, 1. / 2.)) / 2;
		x.arr[1][0] = solution1;
		x.arr[2][0] = solution2;
	}

	private void step3(boolean isNewton, double leftMeaning, double rightMeaning) { // нахождение корней узлового многочлена
		if (isNewton) {
			double eps = 1e-14;
			double right = leftMeaning;
			double left = rightMeaning;
			double solution = (right + left) / 2;
			
			while (function(left) * function(right) > 0) {
				left = solution;
				solution = (left + right)/2;
			}
			double d = function(solution);
			while (abs(d) > eps) {
				double rightF = function(right);
				if (d * rightF > 0) {
                    right = solution;
                } else {
                    left = solution;
                }
				solution = (right + left) / 2;
				d = function(solution);
			}
			x.arr[0][0] = solution;
			findOtherSolutions(solution);
			if (x.arr[0][0] > rightMeaning || x.arr[0][0] < leftMeaning
                    || x.arr[1][0] > rightMeaning || x.arr[1][0] < leftMeaning
				    || x.arr[2][0] > rightMeaning || x.arr[2][0] < leftMeaning) {
				System.out.println("Нарушение диапазона!");
				System.exit(0);
			}
		} else {
			double Q = (a.arr[2][0] * a.arr[2][0] - 3. * a.arr[1][0]) / 9.;
			double R = (2. * pow(a.arr[2][0], 3) - 9. * a.arr[2][0] * a.arr[1][0] + 27. * a.arr[0][0]) / 54.;
			if (pow(R, 2) < pow(Q, 3)) {
				double t = (1. / 3.) * acos(R / pow(Q, 3. / 2.));
				x.arr[0][0] = -2. * pow(Q, 1. / 2.) * cos(t) - a.arr[2][0]/3;
				x.arr[1][0] = -2. * pow(Q, 1. / 2.) * cos(t + (2. * PI / 3.)) - a.arr[2][0] / 3.;
				x.arr[2][0] = -2. * pow(Q, 1. / 2.) * cos(t - (2. * PI / 3.)) - a.arr[2][0] / 3.;
			}
		}
		/*
		double tmp = x._arr[2][0];
		x._arr[2][0] = x._arr[1][0];
		x._arr[1][0] = tmp;
		*/
		//System.out.println("X is: ");
		//x.print();
		if (x.arr[0][0] > rightMeaning || x.arr[0][0] < leftMeaning
                || x.arr[1][0] > rightMeaning || x.arr[1][0] < leftMeaning
			    || x.arr[2][0] > rightMeaning || x.arr[2][0] < leftMeaning) {
            System.out.println("Нарушение диапазона!");
            System.exit(0);
        }
	}

	private void step4() { // сумма (A_i * (x_i)^s = moment_s) (s=0..n-1)
		for (int s = 0; s != 3; ++s) {
			for (int i = 0; i != 3; ++i) {
				A.arr[s][i] = pow(x.arr[i][0], s);
			}
		}
		b.arr[0][0] = moments[0];
		b.arr[1][0] = moments[1];
		b.arr[2][0] = moments[2];
		
		Matrix L, U, P;
		L = new Matrix(A);
		U = new Matrix(A);
		P = new Matrix(A);

        LUP.LUExpansion(A, b, L, U, P);
		A.determinant = U.determinantOfTriangleMatrix();
		A.rank = A.rankOfMatrix();
        a = LUP.decisionSLAU(A, b, L, U);
    }
        
	private double functionMeaning(double x) {
        return 4 * cos(2.5 * x) * exp((5. * x) / 4.) + 2.5 * sin(1.5 * x) * exp((-2. * x)/ 7.) + 5 * x;
	}
        
	private void finalStep() {
		for (int i = 0; i != 3; ++i) {
			res += a.arr[i][0] * functionMeaning(x.arr[i][0]);
		}
		System.out.println("Ответ: " + res);
	}
        
	public double getResult()
	{
		return res;
	}
        
	public Gauss(double left, double right) {
		setMomentsAndMatrix(left, right);
		step2();
		step3(true, left, right);
		step4();
		finalStep();
	}
}

