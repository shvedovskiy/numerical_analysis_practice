package ru.shvedovskiy.integral;

public class Integral {
	public static void main(String[] args) {
        System.out.println("Задание 3. Решение интегралов в квадратурных формулах");
		long time = System.currentTimeMillis();
		new Runge(2.7, 3.2, true, 1e-8, true);
        System.out.println("( 3 )");
        System.out.println("Время расчётов (мс): " + (System.currentTimeMillis() - time) + "\n");
	}
}