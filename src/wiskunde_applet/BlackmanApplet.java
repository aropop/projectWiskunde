package wiskunde_applet;

import java.awt.Color;

import javax.swing.JApplet;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.SwingUtilities;

import java.lang.Math.*;

import static org.math.array.StatisticSample.*;

import org.math.array.util.Function;
import org.math.plot.*;
import org.math.plot.plotObjects.*;
import org.math.plot.utils.Array;
import org.jtransforms.fft.*;

public class BlackmanApplet extends JApplet {

	// Called when this applet is loaded into the browser.
	public void init() {
		// Execute a job on the event-dispatching thread; creating this applet's
		// GUI.
		try {
			SwingUtilities.invokeAndWait(new Runnable() {
				public void run() {
					int min = -100;
					int max = 100;
					int size = (int) ((Math.abs(min) + Math.abs(max)) / Signal.INTERVAL);
					Plot2DPanel plot = new Plot2DPanel();

					// Complex
					double[] x = Signal.rangeX(min, max);
					double[] y = Signal.complexRangeY(min, max);

					plot.addLinePlot("signal", Color.blue, x,
							Signal.getRealValues(y));
					DoubleFFT_1D fft = new DoubleFFT_1D(size);
					double[] y2 = Array.copy(y);
					fft.complexForward(y2);
					double[] y3 = Array.copy(y2);
					fft.complexInverse(y3, true);
					plot.addLinePlot("fourier", x, Signal.getRealValues(y2));
					plot.addLinePlot("original", x, Signal.getRealValues(y3));

					// Real
//					double[] x = Signal.rangeX(min, max);
//					double[] y = Signal.realRangeY(min, max);
//
//					plot.addLinePlot("signal", Color.blue, x, y);
//					DoubleFFT_1D fft = new DoubleFFT_1D(size);
//					double[] y2 = Array.copy(y);
//					fft.realForward(y2);
//					double[] y3 = Array.copy(y2);
//					fft.realInverse(y3, true);
//					plot.addLinePlot("fourier", x, y2);
//					plot.addLinePlot("original", x, y3);
					
					
					setContentPane(plot);

				}
			});
		} catch (Exception e) {
			System.err.println("createGUI didn't complete successfully");
		}
	}
}