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
	private double fc_min = 250; // Cannot be zero, calculating the blackman
									// filter will cause NaN values
	private double fc_max = 750;
	private double delta_f = 250; // transition width
	
	private double signal_frequency = 500;

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

					
					// Generate signal
					double Ts = 1 / fs;
					double[] y = new double[(int) (Math.ceil((1 - Ts) / Ts))];
					for (int i = 0; i < y.length; i++) {
						double t = i * Ts;
						y[i] = 2 * Math.sin(2 * Math.PI * signal_frequency * t);
					}

					// Make one main panel and plots
					JPanel mainPanel = new JPanel();
					Plot2DPanel plot1 = new Plot2DPanel();
					Plot2DPanel plot2 = new Plot2DPanel();
					Plot2DPanel plot3 = new Plot2DPanel();

					// Plot signal
					plot1.addLinePlot("Signal", Color.blue, y);
					
					// Plot fourrier and filtered fourrier
					double[] fourrierValues = centeredFFT(y, fs);
					plot2.addLinePlot("Transformed signal", freqDouble(), fourrierValues);
					double[] filterData = generateFilterData();
					double[] filterecSignal = convolute(y, filterData);
					double[] filterdSignalFourrierValues = centeredFFT(filterecSignal, fs);
					plot2.addLinePlot("Filtered transformed signal", Color.red, freqDouble(),
							filterdSignalFourrierValues);
					
					
					double[] windowTransformed = centeredFFT(filterData, fs);
					plot3.addLinePlot("test", freqDouble(), windowTransformed);
					//plot3.addLinePlot("Filter", windowdata);

					// GUI specific
					plot1.setPreferredSize(new Dimension(350, 500));
					plot2.setPreferredSize(new Dimension(350, 500));
					plot3.setPreferredSize(new Dimension(350, 500));
					plot1.setAxisLabel(0, "Time");
					plot1.setAxisLabel(1, "Signal");
					plot2.setAxisLabel(0, "Frequency (Hz)");
					plot2.setAxisLabel(1, "Amplitude");
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

	// Centered FFT as described in the exercices
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

	// Convert frequencies array from Double[] to double[]
	private double[] freqDouble() {
		double[] r = new double[frequencies.length];
		for (int i = 0; i < frequencies.length; i++) {
			r[i] = frequencies[i].doubleValue();
		}
		return r;
	}

	// Generates filter data, this is a bandpass filter using the blackman window 
	private double[] generateFilterData() {
		// Normalize transition width
		double d_f = delta_f / fs;
		double N = Math.ceil(5.5 / d_f); // Window specification
		double[] ret = new double[(int) N];
		windowArray = new double[(int) N];
		filterArray = new double[(int) N];
		double window, filter;
		int index = 0;
		for (double i = -(N - 1) / 2; i < (N - 1) / 2; i = i + 1) {
			// Calculate window
			window = ALPHA0 - (ALPHA1 * Math.cos(2 * Math.PI * index / (N - 1)))
					+ (ALPHA2 * Math.cos(4 * Math.PI * index / (N - 1))); 
			// Calculate filter
			filter = this.filter(i, d_f);
			
			// Save them in seperate arrays for debugging purposes
			windowArray[index] = window;
			filterArray[index] = filter;
			
			// Multiply filter and window
			ret[index] = window * filter;
			index++;
		}
		return ret;
	}

	// Bandpass filter
	private double filter(double n, double d_f) {
		double pi = Math.PI;
		
		//
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


	// Convolution function
	// Original C version http://www.songho.ca/dsp/convolution/convolution.html#definition
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

}