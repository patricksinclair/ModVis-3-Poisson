import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.JFrame;
import javax.swing.JPanel;
@SuppressWarnings("serial")
class PoissonPanel extends JPanel{

	Vector[][][] vectors;
	int N;
	double maxVal;

	public PoissonPanel(Vector[][][] vectors){
		this.vectors = vectors;
		N = vectors.length;
		setBackground(Color.BLACK);
		double tempVal = Double.MIN_VALUE; 
		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				int k = N/2;
				if(vectors[i][j][k].magnitude2D() > tempVal) tempVal = vectors[i][j][k].magnitude2D();
			}
		}
		maxVal = tempVal;
	}

	@Override
	public void paintComponent(Graphics g){
		super.paintComponent(g);

		Graphics2D g2d = (Graphics2D)g;
		//int w = getWidth()/N;
		int h = getHeight()/N;
		int scaleFactor = h;

		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				int k = N/2;
				drawArrow(g2d, i, j, k, maxVal, scaleFactor);
			}
		}

	}

	public void drawArrow(Graphics2D g, int i, int j, int k, double maxVal, int scaleFactor){
		Vector vector = vectors[i][j][k];
		int interval = scaleFactor;
		double magnitude = vector.magnitude();
		double ratio = magnitude/maxVal;

		//double magnitude2D = ((vector.magnitude()/(0.5*maxVal))*(2*scaleFactor));
		//double magnitude2D = 50.*Math.log(1/(vector.magnitude()/maxVal));
		//double ratio = (1 - Math.exp(-vector.magnitude()/maxVal));
		double magnitude2D = scaleFactor;
		double theta = vector.xyAngle();
		int x1 = i+i*interval;
		int y1 = j+j*interval;
		int x2 = x1 + (int)(magnitude2D*Math.cos(theta));
		int y2 = y1 + (int)(magnitude2D*Math.sin(theta));

		double tipL = 0.25*magnitude2D;
		double qPI = 0.25*Math.PI;
		double am = theta - qPI;
		double ap = theta + qPI;

		if(magnitude > 0.0001){
			g.setColor(new Color(255, (int)(255*(1-ratio)), 255));
			g.drawLine(x1, y1, x2, y2);
			g.drawLine(x2, y2, x2 - (int)Math.round(tipL*Math.cos(am)), y2 - (int)Math.round(tipL*Math.sin(am)));
			g.drawLine(x2, y2, x2 - (int)Math.round(tipL*Math.cos(ap)), y2 - (int)Math.round(tipL*Math.sin(ap)));
		}

	}

}


public class PoissonFrame extends JFrame{

	PoissonPanel pPan;
	Vector[][][] vectors;

	public PoissonFrame(Vector[][][] vectors){
		this.vectors = vectors;
		pPan = new PoissonPanel(vectors);
		pPan.setPreferredSize(new Dimension(700, 700));

		getContentPane().add(pPan, BorderLayout.CENTER);
		pack();

		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				e.getWindow().dispose();
				System.exit(0);
			}
		});

		setTitle("Poisson");
		setLocation(0, 0);
		setVisible(true);
		setBackground(Color.LIGHT_GRAY);
		setVisible(true);

	}









}
