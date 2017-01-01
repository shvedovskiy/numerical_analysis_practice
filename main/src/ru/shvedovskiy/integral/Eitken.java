package ru.shvedovskiy.integral;

public class Eitken {
	private double m;
	/*
	public Eitken(double left, double right, boolean isNewton, int n) {
		int H1 = n / 2;
		int H2 = H1 * 2;
		int H3 = H2 * 2;
		
		double SH1 = new UnionIntegral(left, right, H1, isNewton).getResult();
		double SH2 = new UnionIntegral(left, right, H2, isNewton).getResult();
		double SH3 = new UnionIntegral(left, right, H3, isNewton).getResult();
		double Cm = Math.abs(SH1 - SH2) / Math.abs(SH2 - SH3);
		System.out.println("SH1 is: "+ SH1 + "\nSH2 is: " + SH2 + "\nSH3 is: " + SH3);
		System.out.println("Cm is: " + Cm);
		m = Math.log(Cm) / Math.log(2.);
	}
	*/
	public Eitken(double SH1, double SH2, double SH3) {
		double Cm = Math.abs(SH1 - SH2) / Math.abs(SH2 - SH3);
		m = Math.log(Cm) / Math.log(2.);
	}
	public double getM()
	{
		return m;
	}
}