package wiskunde_applet;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;

import javax.swing.JApplet;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import java.lang.Math.*;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.ListIterator;

import static org.math.array.StatisticSample.*;

import org.math.array.util.Function;
import org.math.plot.*;
import org.math.plot.plotObjects.*;
import org.math.plot.utils.Array;
import org.jtransforms.fft.*;

public class BlackmanApplet extends JApplet {

	private Double[] frequencies;
	private double fs = 8000;
	private double fc_min = 1; // Cannot be zero, calculating the blackman
									// filter will cause NaN values
	private double fc_max = 1000;
	private double delta_f = 500;

	private double[] windowArray;
	private double[] filterArray;

	private static final double ALPHA0 = 0.42;
	private static final double ALPHA1 = 0.5;
	private static final double ALPHA2 = 0.08;

	// Called when this applet is loaded into the browser.
	public void init() {
		// Execute a job on the event-dispatching thread; creating this applet's
		// GUI.
		try {
			SwingUtilities.invokeAndWait(new Runnable() {
				public void run() {
					
					// Run tests
					try {
						tests();
					} catch(Exception e) {
						System.out.println(e.getMessage());
						return;
					}
					
					double Ts = 1 / fs;
					double signal_frequency = 250;
					double[] y = new double[(int) (Math.ceil((1 - Ts) / Ts))];
					for (int i = 0; i < y.length; i++) {
						double t = i * Ts;
						y[i] = 2 * Math.sin(2 * Math.PI * signal_frequency * t);
					}

					// Make one main panel
					JPanel mainPanel = new JPanel();

					Plot2DPanel plot1 = new Plot2DPanel();
					plot1.addLinePlot("signal", Color.blue, y);
					Plot2DPanel plot2 = new Plot2DPanel();
					Plot2DPanel plot3 = new Plot2DPanel();
					double[] fourrierValues = centeredFFT(y, fs);
					plot2.addLinePlot("fourier", freqDouble(), fourrierValues);
					double[] windowdata = generateWindowData(false);
					double[] windowed = convolute(y, windowdata);
					double[] windowedFourrierValues = centeredFFT(windowed, fs);
					plot2.addLinePlot("filter", Color.red, freqDouble(),
							windowedFourrierValues);
					//plot3.addLinePlot("test", centeredFFT(windowArray,fs));
					plot3.addLinePlot("moretes", filterArray);

					plot1.setPreferredSize(new Dimension(350, 500));
					plot2.setPreferredSize(new Dimension(350, 500));
					plot3.setPreferredSize(new Dimension(350, 500));
					mainPanel.setLayout(new BorderLayout());
					mainPanel.add(plot1, BorderLayout.WEST);
					mainPanel.add(plot2);
					mainPanel.add(plot3, BorderLayout.EAST);
					setContentPane(mainPanel);

				}
			});
		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("createGUI didn't complete successfully");
		}
	}

	public double[] centeredFFT(double[] data, double sampling_frequency) {

		double N = data.length * 2;
		double[] returnData = new double[data.length * 2];
		ArrayList<Double> k = new ArrayList<Double>();
		if ((N % 2) == 0) {
			for (double i = -N / 2; i < (N - 1) / 2; i = i + 1) {
				k.add(i);
			}
		} else {
			for (double i = -(N - 1) / 2; i < (N - 1) / 2; i = i + 1) {
				k.add(i);
			}
		}

		double T = N / sampling_frequency;
		for (ListIterator<Double> iterator = k.listIterator(); iterator
				.hasNext();) {
			Double double1 = (Double) iterator.next();
			iterator.set(Double.valueOf(double1.doubleValue() / T));
		}
		this.frequencies = k.toArray(new Double[k.size()]);

		for (int i = 0; i < data.length; ++i) {
			returnData[i] = data[i];
		}

		DoubleFFT_1D fft = new DoubleFFT_1D(data.length);
		fft.realForwardFull(returnData);
		for (int i = 0; i < returnData.length; i++) {
			returnData[i] = Math.abs(returnData[i]) / N; // Normalize data
		}

		// Shift
		double[] temp = new double[returnData.length];
		for (int i = 0; i < returnData.length; i++) { // make temp array with
														// same contents as x
			temp[i] = returnData[i];
		}

		for (int i = 0; i < returnData.length / 2; i++) {
			returnData[i] = temp[returnData.length / 2 + i];
			returnData[returnData.length / 2 + i] = temp[i];
		}
		return returnData;
	}

	private double[] freqDouble() {
		double[] r = new double[frequencies.length];
		for (int i = 0; i < frequencies.length; i++) {
			r[i] = frequencies[i].doubleValue();
		}
		return r;
	}

	// Generates window data for blackman window
	private double[] generateWindowData(boolean test) {
		double d_f = delta_f / fs;
		double N = test ? 10 : Math.ceil(5.5 / d_f);// See wpo 6-8 document page 18 for
		double[] ret = new double[(int) N];
		windowArray = new double[(int) N];
		filterArray = new double[(int) N];
		double window, filter;
		int index = 0;
		for (double i = -(N - 1) / 2; i < (N - 1) / 2; i = i + 1) {
			window = ALPHA0 - (ALPHA1 * Math.cos(2 * Math.PI * index / (N - 1)))
					+ (ALPHA2 * Math.cos(4 * Math.PI * index / (N - 1)));
			// window = 1;
			filter = this.filter(i, d_f);
			//filter = filterLowPass(i);
			windowArray[index] = window;
			filterArray[index] = filter;
			ret[index] = window * filter;
			index++;
		}
		return ret;
	}

	// Zie wpo 6-8 document
	private double filter(double n, double d_f) {
		double pi = Math.PI;
		double f1 = (fc_min / fs) + (d_f/2);
		double f2 = (fc_max / fs) + (d_f/2);
		if (n == 0) {
			return 2 * (f2 - f1);
		} else {
			return (2 * f2 * (Math.sin(n * 2 * pi * f2) / (n * 2 * pi * f2)))
					- (2 * f1 * (Math.sin(n * 2 * pi * f1) / (n * 2
							* pi * f1)));

		}
	}

	// Zie wpo 6-8 document
	private double filterLowPass(double n) {
		double pi = Math.PI;
		if (n == 0) {
			return 2 * (fc_max);
		} else {
			return (2 * fc_max * (Math.sin(n * 2 * pi * fc_max) / (n * 2 * pi * fc_max)));

		}
	}

	// Check http://www.songho.ca/dsp/convolution/convolution.html#definition
	// Kernel = window
	private double[] convolute(double[] input, double[] kernel) {
		double[] ret = new double[input.length];
		int kernelSize = kernel.length;
		int dataSize = input.length;
		// start convolution from out[kernelSize-1] to out[dataSize-1] (last)
		for (int i = kernelSize - 1; i < dataSize; ++i) {
			ret[i] = 0;

			for (int j = i, k = 0; k < kernelSize; --j, ++k) {
				ret[i] += input[j] * kernel[k];
			}
		}
		// convolution from out[0] to out[kernelSize-2]
		for (int i = 0; i < kernelSize - 1; ++i) {
			ret[i] = 0; // init to 0 before sum

			for (int j = i, k = 0; j >= 0; --j, ++k)
				ret[i] += input[j] * kernel[k];
		}
		return ret;
	}

	private double[] convoluteNew(double[] signal, double[] filter) {
		double[] ret = new double[signal.length];

		for (int i = 0; i < signal.length; ++i) {

			double sum_hk_and_signal = 0;

			double M = Math.min(filter.length, i + 1);

			for (int k = 0; k < M; ++k) {

				double h_of_k = filter[k];

				sum_hk_and_signal += h_of_k * filter[i - k];

			}

			ret[i] = sum_hk_and_signal;
		}

		return ret;
	}

	private void tests() throws Exception {
		
		//Arrayaprox
		double[] i1 = {4,5.000000000001, 6.999999999999999};
		double[] i2 = {4,5,7};
		if(!arrayAprox(i1, i2)) {
			throw new Exception("arrayAprox failed");
		}
		
		// Convolute
		double[] input1 = { 1, 2, 3, 4, 5 };
		double[] input2 = { 6, 5, 7 };
		double[] output = { 6, 17, 35, 53, 71};
		if(!arrayAprox(convolute(input1, input2), output)) {
			throw new Exception("Convolution failed");
		}
		
		// window
		double[] output2 = {  -0.00000,
				  0.05087,
				  0.25800,
				  0.63000,
				  0.95113,
				  0.95113,
				  0.63000,
				  0.25800,
				  0.05087,
				  -0.00000};
		double[] tst = generateWindowData(true);
		if(!arrayAprox(windowArray, output2)) {
			throw new Exception("Window calculation failed");
		}

	}
	
	private boolean arrayAprox(double[] ar1,double[] ar2) {
		if(ar1.length != ar2.length) return false;
		for(int i=0;i<ar1.length;i++) {
			if(((double)Math.round(ar1[i] * 100000) / 100000) == ((double)Math.round(ar2[i] * 100000) / 100000)) {
				continue;
			}else {
				System.out.println(((double)Math.round(ar1[i] * 10000) / 10000) + " where " + ar2[i]);
				return false;
			}
		}
		return true;
	}

}