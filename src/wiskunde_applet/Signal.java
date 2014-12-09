package wiskunde_applet;

public class Signal {

	public static final double INTERVAL = 0.05;
	
	/**
	 * Returns f(x) where f is the signal function
	 * 
	 * @param x
	 * @return
	 */
	public static double compute(double x) {
		// Dirac
//		if(x == 0) {
//			return Double.MAX_VALUE;
//		} else {
//			return 0;
//		}
		
		// finite-pulse function
//		double t = 1;
//		if((-t < x) && (x < t)) {
//			return 1/(2 * t);
//		} else {
//			return 0;
//		}
		
		// Cosinus
		//return Math.cos(x);
		
		// sinc
		return Math.sin(x) / x;
	}
	
	
	/**
	 * Returns a range for this function
	 * 
	 * @param from
	 * @param to
	 * @return array with the values
	 */
	public static double[] complexRangeY(int from, int to) {
		int size = (int) ((to - from) / INTERVAL) * 2;
		double[] array = new double[size];
		for (int i = 0; i < size; i++) {
			if((i % 2) == 0) {
				array[i] = compute(from + ((i / 2) *INTERVAL));				
			} else {
				array[i] = 0;
			}
		}
		return array;
	}
	
	public static double[] realRangeY(int from, int to) {
		int size = (int) ((to - from) / INTERVAL);
		double[] array = new double[size];
		for (int i = 0; i < size; i++) {
				array[i] = compute(from + (i *INTERVAL));				
		}
		return array;
	}
	
	public static double[] rangeX(int from, int to) {
		int size = (int) ((to - from) / INTERVAL);
		double[] array = new double[size];
		for (int i = 0; i < size; i++) {
			array[i] = from + (i*INTERVAL);
		}
		return array;
	}
	
	public static double[] getRealValues(double[] array) {
		double[] returnArray = new double[array.length / 2];
		for(int i = 0; i < returnArray.length; i++) {
			double real = array[2*i];
			double imaginary = array[(2*i) + 1];
			// waarom?
			returnArray[i] = Math.sqrt((real*real) + (imaginary*imaginary));
		}
		return returnArray;
	}
}
