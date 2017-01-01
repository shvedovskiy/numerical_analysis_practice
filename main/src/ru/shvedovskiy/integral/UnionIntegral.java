package ru.shvedovskiy.integral;

public class UnionIntegral {
	private double res = 0;
        
	public UnionIntegral(double left, double right, int n, boolean isNewton) {
		double[] x = new double[n + 1];
		for (int i = 0; i != n + 1; ++i) {
			x[i] = left + (right - left) * i / (n);
		}
		for (int i = 1; i < n + 1; ++i) {
			if (isNewton) {
				res += new NewtonKotes(x[i - 1], x[i]).getResult();
			} else {
				res += new Gauss(x[i - 1], x[i]).getResult();
			}
			System.out.println("Шаг "+ i + " завершён\n");
		}
	}
        
	public double getResult()
	{
		return res;
	}
}