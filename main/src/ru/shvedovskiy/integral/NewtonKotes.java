package ru.shvedovskiy.integral;

import ru.shvedovskiy.lib.Matrix;
import ru.shvedovskiy.lup.LUP;
import static java.lang.Math.*;

public class NewtonKotes {
	private Matrix a, b;
	private Matrix A; // коэффициенты
	private Matrix x; // узлы
	private double[] moments;
	private double res = 0;
        
	private void setMomentsAndMatrix(double left, double right) {
		x = new Matrix(3, 1); // трёхточечная функция
		for (int i = 0; i != 3; ++i) {
			x.arr[i][0] = left + i*(right - left)/2.0;
		}
		moments = new double[3];
		/*
		moments[0] = -5. / 3. * pow((2.9  - right), 3. / 5.) + 5. / 3. * pow((2.9  - left), 3. / 5.);
		moments[1] = -5. / 48. * pow((29. / 10. - right), (3. / 5.)) * (6 * right + 29) +
		        5. / 48. * pow((29. / 10. - left), (3. / 5.)) * (6 * left + 29);
		moments[2] = -(pow(5, 2. / 5.) * pow(29. / 2. - 5 * right, 3. / 5.) * (48 * right * right + 174 * right + 841)) / 624. +
			    (pow(5, 2. / 5.) * pow(29. / 2. - 5 * left, 3. / 5.) * (48 * left * left + 174 * left + 841)) / 624.;
	    */
        moments[0] = -4. * exp(log(16. / 5. - right) / 4.) + 4. * exp(log(16. /5. - left) / 4.);
        moments[1] = -(4. / 25.) * exp(log(16. / 5. - right) / 4.) * ((5. * right) + 64.)
                + (4. / 25.) * exp(log(16. / 5. - left) / 4.) * ((5. * left) + 64.);
        moments[2] = -(4. / 1125.) * exp(log(16. / 5. - right) / 4.) * ((125. * pow(right, 2)) + (640 * right) + 8192)
                + (4. / 1125.) * exp(log(16. / 5. - left) / 4.) * ((125. * pow(left, 2)) + (640 * left) + 8192);
                
		b = new Matrix(3, 1);
		b.arr[0][0] = moments[0];
		b.arr[1][0] = moments[1];
		b.arr[2][0] = moments[2];

		System.out.println("Матрица X: ");
		x.printMatrix();
		
        System.out.println("Моменты весовой функции: ");
		b.printMatrix();
		
		A = new Matrix(3, 3);
		A.arr[0][0] = 1;
		A.arr[1][0] = x.arr[0][0];
		A.arr[2][0] = x.arr[0][0] * x.arr[0][0];
		A.arr[0][1] = 1;
		A.arr[1][1] = x.arr[1][0];
		A.arr[2][1] = x.arr[1][0] * x.arr[1][0];
		A.arr[0][2] = 1;
		A.arr[1][2] = x.arr[2][0];
		A.arr[2][2] = x.arr[2][0] * x.arr[2][0];
		
		
		Matrix L, U, P;
		L = new Matrix(A);
		U = new Matrix(A);
		P = new Matrix(A);

        LUP.LUExpansion(A, b, L, U, P);
		System.out.println("Матрица U: ");
		U.printMatrix();
		A.determinant = U.determinantOfTriangleMatrix();
        a = LUP.decisionSLAU(A, b, L, U);
		System.out.println("Решение СЛАУ (матрица a): ");
		a.printMatrix();
		
		for (int i = 0; i != 3; ++i) {
			res += a.arr[i][0] * functionMeaning(x.arr[i][0]);
		}
		System.out.println("Конец шага, ответ: " + res + "\n");
	}
        
	private double functionMeaning(double x) {
        return 4 * cos(2.5 * x) * exp((5. * x) / 4.) + 2.5 * sin(1.5 * x) * exp((-2. * x) /  7.) + 5 * x;
	}
        
	public double getResult()
	{
		return res;
	}
        
	public NewtonKotes(double a, double b)
	{
		setMomentsAndMatrix(a, b);
	}
}