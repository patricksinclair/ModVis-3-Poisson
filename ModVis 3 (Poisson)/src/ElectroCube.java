import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;

//This is the ElectroCube class, it contains a 3D lattice of doubles which hold values of various fields
public class ElectroCube {

	int N;
	double precision;
	ElectroPoint[][][] electroLattice;
	static int Np = 50;
	static double precisionVal = 1E-5;
	static Random rand = new Random();

	public ElectroCube(int N, double precision){
		this.N = N;
		this.precision = precision;
	}


	public ElectroCube(ElectroPoint[][][] electroLattice, double precision){
		this.electroLattice = electroLattice;
		this.precision = precision;
	}

	public int getN(){
		return N;
	}
	public double getPrecision(){
		return precision;
	}
	public void setPrecision(double precision){
		this.precision = precision;
	}

	public ElectroPoint[][][] getElectroLattice(){
		return electroLattice;
	}
	public void setElectroLattice(ElectroPoint[][][] electroLattice){
		this.electroLattice = electroLattice;
	}
	public ElectroPoint getElectroPoint(int i, int j, int k){
		return getElectroLattice()[i][j][k];
	}
	public double getPhiAtPoint(int i, int j, int k){
		return getElectroLattice()[i][j][k].getPhi();
	}
	public void setPhiAtPoint(int i, int j, int k, double newPhi){
		getElectroLattice()[i][j][k].setPhi(newPhi);
	}
	public double getRhoAtPoint(int i, int j, int k){
		return getElectroLattice()[i][j][k].getRho();
	}
	public double getEAtPoint(int i, int j, int k){
		return getElectroLattice()[i][j][k].getE();
	}
	public void setEAtPoint(int i, int j, int k, double newE){
		getElectroLattice()[i][j][k].setE(newE);
	}

	public void printCubeAsList(int N){
		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				for(int k = 0; k < N; k++){
					System.out.println(i + " " + j + " " + k + " " + getPhiAtPoint(i, j, k) + " " + getRhoAtPoint(i, j, k));
				}
			}
		}
	}

	//The jacobi method for calculating the potential resulting from a given charge distribution
	public void diffusePotentialJacobi(){

		int N = getElectroLattice().length;
		ElectroPoint[][][] tempELat = ElectroPoint.allzeroes(N);
		int counter = 0;

		whileloop:
			while(true){
				boolean converged = true;
				for(int i = 1; i < N-1; i++){
					for(int j = 1; j < N-1; j++){
						for(int k = 1; k < N-1; k++){

							double phiTemp = (1./6.)*(getPhiAtPoint(i+1, j, k) + getPhiAtPoint(i-1, j, k)
									+ getPhiAtPoint(i, j+1, k) + getPhiAtPoint(i, j-1, k)
									+ getPhiAtPoint(i, j, k+1) + getPhiAtPoint(i, j, k-1)
									+ getRhoAtPoint(i, j, k));

							if(Math.abs(phiTemp - getPhiAtPoint(i, j, k)) > getPrecision()) converged = false;

							tempELat[i][j][k] = new ElectroPoint(phiTemp, getRhoAtPoint(i, j, k)); 
						}
					}
				}
				if(converged) break whileloop;
				setElectroLattice(tempELat);
				System.out.println(counter);
				counter++;
			}
	}


	public void diffusePotentialGS2(){

		int N = getElectroLattice().length;
		int counter = 0;

		whileloop:
			while(true){
				boolean converged = true;
				for(int i = 1; i < N-1; i++){
					for(int j = 1; j < N-1; j++){
						for(int k = 1; k < N-1; k++){
							double newPhi = (1./6.)*(getPhiAtPoint(i+1, j, k) + getPhiAtPoint(i-1, j, k)
									+ getPhiAtPoint(i, j+1, k) + getPhiAtPoint(i, j-1, k)
									+ getPhiAtPoint(i, j, k+1) + getPhiAtPoint(i, j, k-1)
									+getRhoAtPoint(i, j, k));
							if(Math.abs(newPhi - getPhiAtPoint(i, j, k)) > getPrecision()) converged = false;
							setPhiAtPoint(i, j, k, newPhi);
						}
					}
				}
				if(converged) break whileloop;
				if(counter %50 == 0)System.out.println(counter);
				counter++;
			}
	}

	//The gauss-siedel method for calculating the potential for a given charge distribution
	public void diffusePotentialGS(){
		int N = getElectroLattice().length;
		int counter = 0;

		whileloop:
			while(true){
				boolean converged = true;
				for(int i = 1; i < N-1; i++){
					for(int j = 1; j < N-1; j++){
						for(int k = 1; k < N-1; k++){
							if((i+j+k)%2==0){
								double newPhi = (1./6.)*(getPhiAtPoint(i+1, j, k) + getPhiAtPoint(i-1, j, k)
										+ getPhiAtPoint(i, j+1, k) + getPhiAtPoint(i, j-1, k)
										+ getPhiAtPoint(i, j, k+1) + getPhiAtPoint(i, j, k-1)
										+getRhoAtPoint(i, j, k));
								if(Math.abs(newPhi - getPhiAtPoint(i, j, k)) > getPrecision()) converged = false;
								setPhiAtPoint(i, j, k, newPhi);
							}
						}
					}
				}

				for(int i = 1; i < N-1; i++){
					for(int j = 1; j < N-1; j++){
						for(int k = 1; k < N-1; k++){
							if((i+j+k)%2==1){
								double newPhi = (1./6.)*(getPhiAtPoint(i+1, j, k) + getPhiAtPoint(i-1, j, k)
										+ getPhiAtPoint(i, j+1, k) + getPhiAtPoint(i, j-1, k)
										+ getPhiAtPoint(i, j, k+1) + getPhiAtPoint(i, j, k-1)
										+getRhoAtPoint(i, j, k));
								if(Math.abs(newPhi - getPhiAtPoint(i, j, k)) > getPrecision()) converged = false;
								setPhiAtPoint(i, j, k, newPhi);
							}
						}
					}
				}
				if(converged) break whileloop;
				if(counter %50 == 0)System.out.println(counter);
				counter++;
			}
	}

	//method for calcultating the no. of iterations used for the over-relaxed gs method.
	public double relaxedIterations(double omega){

		int N = getElectroLattice().length;
		int counter = 0;

		whileloop:
			while(true){
				boolean converged = true;
				for(int i = 1; i < N-1; i++){
					for(int j = 1; j < N-1; j++){
						for(int k = 1; k < N-1; k++){

							if((i+j+k)%2==0){
								double newPhi = (1-omega)*getPhiAtPoint(i, j, k) + omega*(1./6.)*(getPhiAtPoint(i+1, j, k) + getPhiAtPoint(i-1, j, k)
										+ getPhiAtPoint(i, j+1, k) + getPhiAtPoint(i, j-1, k)
										+ getPhiAtPoint(i, j, k+1) + getPhiAtPoint(i, j, k-1)
										+getRhoAtPoint(i, j, k));

								if(Math.abs(newPhi - getPhiAtPoint(i, j, k)) > getPrecision()) converged = false;

								setPhiAtPoint(i, j, k, newPhi);
							}
						}
					}
				}
				
				for(int i = 1; i < N-1; i++){
					for(int j = 1; j < N-1; j++){
						for(int k = 1; k < N-1; k++){

							if((i+j+k)%2==1){
								double newPhi = (1-omega)*getPhiAtPoint(i, j, k) + omega*(1./6.)*(getPhiAtPoint(i+1, j, k) + getPhiAtPoint(i-1, j, k)
										+ getPhiAtPoint(i, j+1, k) + getPhiAtPoint(i, j-1, k)
										+ getPhiAtPoint(i, j, k+1) + getPhiAtPoint(i, j, k-1)
										+getRhoAtPoint(i, j, k));

								if(Math.abs(newPhi - getPhiAtPoint(i, j, k)) > getPrecision()) converged = false;

								setPhiAtPoint(i, j, k, newPhi);
							}
						}
					}
				}
				
				if(converged || counter > 5000) break whileloop;
				if(counter%20 == 0) System.out.println(counter+"  "+omega);
				counter++;
			}
		return counter;
	}

	//method for determining the electric field from the calculated potential
	public void computeElectricField(){
		int N = getElectroLattice().length;
		for(int i = 1; i < N-1; i++){
			for(int j = 1; j < N-1; j++){
				for(int k = 1; k < N-1; k++){

					double Ex = -0.5*(getPhiAtPoint(i+1, j, k) - getPhiAtPoint(i-1, j, k));
					double Ey = -0.5*(getPhiAtPoint(i, j+1, k) - getPhiAtPoint(i, j-1, k));
					double Ez = -0.5*(getPhiAtPoint(i, j, k+1) - getPhiAtPoint(i, j, k-1));
					double Etot = Math.sqrt(Ex*Ex+Ey*Ey+Ez*Ez);

					setEAtPoint(i, j, k, Etot);
				}
			}
		}
	}

	//calculates the 
	public static void pointPotentialLineGraph(Scanner keyboard){
		System.out.println("Please enter the size of the system");
		int N = keyboard.nextInt();
		//int N = 100;
		System.out.println("Please enter the desired precision");
		double precision = keyboard.nextDouble();
		//double precision = precisionVal;
		String filename = "pointChargePotentialLineGraph";
		ElectroCube eCube = ElectroCube.pointAtCentre(N, precision);

		eCube.diffusePotentialGS();

		ArrayList<Double> xData = new ArrayList<Double>();
		ArrayList<Double> yData = new ArrayList<Double>();

		for(int i = 0; i < N; i++){
			xData.add((double)i);
			yData.add(eCube.getPhiAtPoint(i, N/2, N/2));	
		}
		ElectroToolbox.printToFile(filename, xData, yData);
	}


	public static void gaussianLineGraph(Scanner keyboard){
		System.out.println("Please enter the size of the system");
		int N = keyboard.nextInt();
		//int N = 100;
		System.out.println("Please enter the desired precision");
		double precision = keyboard.nextDouble();
		//double precision = precisionVal;
		int x = N/2, y = N/2, z = N/2;
		
		String filename = "gaussianChargePotentialLineGraph";
		ElectroCube eCube = ElectroCube.gaussianPointAtPoint(x, y, z, N, precision);

		eCube.diffusePotentialGS();

		ArrayList<Double> xData = new ArrayList<Double>();
		ArrayList<Double> yData = new ArrayList<Double>();

		for(int i = 0; i < N; i++){
			xData.add((double)i);
			yData.add(eCube.getPhiAtPoint(i, N/2, N/2));
		}
		ElectroToolbox.printToFile(filename, xData, yData);
	}


	public static void ELineGraph(Scanner keyboard){
		System.out.println("Please enter the size of the system");
		int N = keyboard.nextInt();
		//int N = 100;
		System.out.println("Please enter the desired precision");
		double precision = keyboard.nextDouble();
		//double precision = precisionVal;
		int x = N/2, y = N/2, z = N/2;
		//double precision = precisionVal;
		String filename = "gaussianChargeEFieldLineGraph";
		ElectroCube eCube = ElectroCube.gaussianPointAtPoint(x, y, z, N, precision);
		//ElectroCube eCube = ElectroCube.pointAtCentre(N, precision);

		eCube.diffusePotentialGS();
		eCube.computeElectricField();

		ArrayList<Double> xData = new ArrayList<Double>();
		ArrayList<Double> yData = new ArrayList<Double>();

		for(int i = 0; i < N; i++){
			xData.add((double)i);
			yData.add(eCube.getEAtPoint(i, N/2, N/2));
		}
		ElectroToolbox.printToFile(filename, xData, yData);
	}


	public static void pointPotentialContourPlot(Scanner keyboard){
		System.out.println("Please enter the size of the system");
		int N = keyboard.nextInt();
		//int N = 100;
		System.out.println("Please enter the desired precision");
		double precision = keyboard.nextDouble();
		//double precision = precisionVal;
		
		String filename = "pointChargePotentialContourPlot";
		ElectroCube eCube = ElectroCube.pointAtCentre(N, precision);

		//eCube.diffusePotentialJacobi();
		eCube.diffusePotentialGS();

		ArrayList<Double> xData = new ArrayList<Double>();
		ArrayList<Double> yData = new ArrayList<Double>();
		ArrayList<Double> zData = new ArrayList<Double>();

		int k = N/2;
		for(int i = 0; i < N; i++){
			xData.add((double)i);
			for(int j = 0; j < N; j++){
				yData.add((double)j);
				zData.add(eCube.getPhiAtPoint(i, j , k));
			}
		}
		ElectroToolbox.contourFileWrite(filename, xData, yData, zData);
	}


	public static void quadropolePotentialContourPlot(Scanner keyboard){
		System.out.println("Please enter the size of the system");
		int N = keyboard.nextInt();
		//int N = 100;
		System.out.println("Please enter the desired precision");
		double precision = keyboard.nextDouble();
		//double precision = precisionVal;
		String filename = "gaussianQuadropolePotentialContourPlot";
		//ElectroCube eCube = ElectroCube.quadropole(N, precision);
		ElectroCube eCube = ElectroCube.gaussianQuadropole(N, precision);

		eCube.diffusePotentialGS();

		ArrayList<Double> xData = new ArrayList<Double>();
		ArrayList<Double> yData = new ArrayList<Double>();
		ArrayList<Double> zData = new ArrayList<Double>();

		int k = N/2;
		for(int i = 0; i < N; i++){
			xData.add((double)i);
			for(int j = 0; j < N; j++){
				yData.add((double)j);
				zData.add(eCube.getPhiAtPoint(i, j , k));
			}
		}
		ElectroToolbox.contourFileWrite(filename, xData, yData, zData);
	}


	public static void gaussianPotentialContourPlot(Scanner keyboard){
		System.out.println("Please enter the size of the system");
		int N = keyboard.nextInt();
		//int N = 100;
		System.out.println("Please enter the desired precision");
		double precision = keyboard.nextDouble();
		//double precision = precisionVal;
		int x = N/2;
		int y = N/2;
		int z = N/2;
	
		String filename = "gaussianChargePotentialContourPlot";
		ElectroCube eCube = ElectroCube.gaussianPointAtPoint(x, y, z, N, precision);

		//eCube.diffusePotentialJacobi();
		eCube.diffusePotentialGS();

		ArrayList<Double> xData = new ArrayList<Double>();
		ArrayList<Double> yData = new ArrayList<Double>();
		ArrayList<Double> zData = new ArrayList<Double>();

		int k = N/2;
		for(int i = 0; i < N; i++){
			xData.add((double)i);
			for(int j = 0; j < N; j++){
				yData.add((double)j);
				zData.add(eCube.getPhiAtPoint(i, j , k));
			}
		}
		ElectroToolbox.contourFileWrite(filename, xData, yData, zData);
	}


	public static void EFieldContourPlot(Scanner keyboard){
		System.out.println("Please enter the size of the system");
		int N = keyboard.nextInt();
		//int N = 100;
		System.out.println("Please enter the desired precision");
		double precision = keyboard.nextDouble();
		//double precision = precisionVal;
		int x = N/2;
		int y = N/2;
		int z = N/2;
		
		String filename = "gaussianChargeEFieldContourPlot";
		//ElectroCube eCube = ElectroCube.pointAtCentre(N, precision);
		ElectroCube eCube = ElectroCube.gaussianPointAtPoint(x, y, z, N, precision);

		//eCube.diffusePotentialJacobi();
		eCube.diffusePotentialGS();
		eCube.computeElectricField();


		ArrayList<Double> xData = new ArrayList<Double>();
		ArrayList<Double> yData = new ArrayList<Double>();
		ArrayList<Double> zData = new ArrayList<Double>();

		int k = N/2;
		for(int i = 0; i < N; i++){
			xData.add((double)i);
			for(int j = 0; j < N; j++){
				yData.add((double)j);
				zData.add(eCube.getEAtPoint(i, j, k));
			}
		}
		ElectroToolbox.contourFileWrite(filename, xData, yData, zData);
	}


	public static void quadropoleEFieldContourPlot(Scanner keyboard){
		System.out.println("Please enter the size of the system");
		int N = keyboard.nextInt();
		//int N = 100;
		System.out.println("Please enter the desired precision");
		double precision = keyboard.nextDouble();
		//double precision = precisionVal;
		String filename = "gaussianQuadropoleEFieldContourPlot";
		//ElectroCube eCube = ElectroCube.quadropole(N, precision);
		ElectroCube eCube = ElectroCube.gaussianQuadropole(N, precision);

		//eCube.diffusePotentialJacobi();
		eCube.diffusePotentialGS();
		eCube.computeElectricField();


		ArrayList<Double> xData = new ArrayList<Double>();
		ArrayList<Double> yData = new ArrayList<Double>();
		ArrayList<Double> zData = new ArrayList<Double>();

		int k = N/2;
		for(int i = 0; i < N; i++){
			xData.add((double)i);
			for(int j = 0; j < N; j++){
				yData.add((double)j);
				zData.add(eCube.getEAtPoint(i, j, k));
			}
		}
		ElectroToolbox.contourFileWrite(filename, xData, yData, zData);
	}


	public static void visualiseEField(){

		int N = 100;
		int x = N/2;
		int y = N/2;
		int z = N/2;
		double precision = precisionVal;

		//ElectroCube eCube = ElectroCube.gaussianPointAtPoint(x, y, z, N, precision);
		ElectroCube eCube = ElectroCube.gaussianQuadropole(N, precision);
		//ElectroCube eCube = ElectroCube.pointAtCentre(N, precision);
		eCube.diffusePotentialGS();
		Vector[][][] vCube = Vector.electroVectorCube(eCube);

		PoissonFrame pFrame =  new PoissonFrame(vCube);
		pFrame.setVisible(true);
	}


	public static void visualiseBField(){

		int N = 100;
		int x = N/2;
		int y = N/2;
		int z = N/2;
		double j1 = 1.0;
		double j2 = 1.0;
		double precision = precisionVal;

		//ElectroCube bCube = ElectroCube.wireThroughCentre(N, precision);
		ElectroCube bCube = ElectroCube.twoParallelWires(N, precision, j1, j2);
		bCube.diffusePotentialGS();

		Vector[][][] vCube = Vector.magnetoVectorCube(bCube);

		PoissonFrame pFrame =  new PoissonFrame(vCube);
		pFrame.setVisible(true);

	}


	public static void relaxedGS(){

		int N = 100;
		double precision = precisionVal;
		double increment = 0.01;
		String filename = "optimalOmegaQuadrupole";


		ArrayList<Double> omegaX = new ArrayList<Double>();
		ArrayList<Double> iterationsY = new ArrayList<Double>();

		for(double omega = 1.+increment; omega < 2.; omega+=increment){
			omegaX.add(omega);
			//ElectroCube eCube = ElectroCube.gaussianPointAtPoint(N/2, N/2, N/2, N, precision);
			//ElectroCube eCube = ElectroCube.pointAtCentre(N, precision);
			ElectroCube eCube = ElectroCube.quadropole(N, precision);
			double iterations = eCube.relaxedIterations(omega);
			iterationsY.add(iterations);
			System.out.println(omega);
		}

		ElectroToolbox.printToFile(filename, omegaX, iterationsY);
	}




	public static ElectroCube pointAtCentre(int N, double precision){

		ElectroPoint[][][] electroLat = new ElectroPoint[N][N][N];
		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				for(int k = 0; k < N; k++){
					electroLat[i][j][k] = new ElectroPoint(0.);
				}
			}
		}
		electroLat[N/2][N/2][N/2] = new ElectroPoint(1.);
		return new ElectroCube(electroLat, precision);
	}


	public static ElectroCube pointAtPoint(int x, int y, int z, int N, double precision){
		ElectroPoint[][][] electroLat = new ElectroPoint[N][N][N];
		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				for(int k = 0; k < N; k++){
					electroLat[i][j][k] = new ElectroPoint(0.);
				}
			}
		}
		electroLat[x][y][z] = new ElectroPoint(1.);
		return new ElectroCube(electroLat, precision);
	}


	public static ElectroCube quadropole(int N, double precision){

		ElectroPoint[][][] electroLat = new ElectroPoint[N][N][N];
		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				for(int k = 0; k < N; k++){
					electroLat[i][j][k] = new ElectroPoint(0.);
				}
			}
		}
		electroLat[N/4][N/2][N/2] = new ElectroPoint(1.);
		electroLat[N/2][3*N/4][N/2] = new ElectroPoint(-1.);
		electroLat[3*N/4][N/2][N/2] = new ElectroPoint(1.);
		electroLat[N/2][N/4][N/2] = new ElectroPoint(-1.);


		return new ElectroCube(electroLat, precision);
	}


	public static ElectroCube gaussianPointAtPoint(int x, int y, int z, int N, double precision){

		ElectroPoint[][][] electroLat = new ElectroPoint[N][N][N];
		double sigmaSq = 20;
		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				for(int k = 0; k < N; k++){
					double r = Math.sqrt((i-x)*(i-x) + (j-y)*(j-y) + (k-z)*(k-z));
					double rho = (1./Math.sqrt(2.*Math.PI))*Math.exp(-0.5*(r*r)/sigmaSq);
					electroLat[i][j][k] = new ElectroPoint(0., rho);
				}
			}
		}
		return new ElectroCube(electroLat, precision);
	}


	public static ElectroCube gaussianQuadropole(int N, double precision){

		ElectroPoint[][][] electroLat = new ElectroPoint[N][N][N];
		double sigmaSq = 20;

		int x1 = N/2-10, y1 = N/2-10;
		int x2 = N/2-10, y2 = N/2+10;
		int x3 = N/2+10, y3 = N/2+10;
		int x4 = N/2+10, y4 = N/2-10;
		int z = N/2;

		double qPos = 1.;
		double qNeg = -1.;

		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				for(int k = 0; k < N; k++){
					electroLat[i][j][k] = new ElectroPoint(0.);
				}
			}
		}

		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				for(int k = 0; k < N; k++){

					double r1 = Math.sqrt((i-x1)*(i-x1) + (j-y1)*(j-y1) + (k-z)*(k-z));
					double r2 = Math.sqrt((i-x2)*(i-x2) + (j-y2)*(j-y2) + (k-z)*(k-z));
					double r3 = Math.sqrt((i-x3)*(i-x3) + (j-y3)*(j-y3) + (k-z)*(k-z));
					double r4 = Math.sqrt((i-x4)*(i-x4) + (j-y4)*(j-y4) + (k-z)*(k-z));

					double rho1 = qPos*(1./Math.sqrt(2.*Math.PI))*Math.exp(-0.5*(r1*r1)/sigmaSq);
					double rho2 = qNeg*(1./Math.sqrt(2.*Math.PI))*Math.exp(-0.5*(r2*r2)/sigmaSq);
					double rho3 = qPos*(1./Math.sqrt(2.*Math.PI))*Math.exp(-0.5*(r3*r3)/sigmaSq);
					double rho4 = qNeg*(1./Math.sqrt(2.*Math.PI))*Math.exp(-0.5*(r4*r4)/sigmaSq);
					double rho = rho1 + rho2 + rho3 + rho4;
					electroLat[i][j][k] = new ElectroPoint(0., rho);
				}
			}
		}

		return new ElectroCube(electroLat, precision);	
	}


	public static ElectroCube wireThroughCentre(int N, double precision){

		ElectroPoint[][][] electroLat = new ElectroPoint[N][N][N];

		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				for(int k = 0; k < N; k++){
					if(i == N/2 && j == N/2) electroLat[i][j][k] = new ElectroPoint(1.);
					else electroLat[i][j][k] = new ElectroPoint(0.);
				}
			}
		}

		return new ElectroCube(electroLat, precision);
	}


	public static ElectroCube twoParallelWires(int N, double precision, double j1, double j2){

		ElectroPoint[][][] electroLat = new ElectroPoint[N][N][N];

		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				for(int k = 0; k < N; k++){
					if(i == N/4 && j == N/2) electroLat[i][j][k] = new ElectroPoint(j1);
					else if(i == 3*N/4 && j == N/2) electroLat[i][j][k] = new ElectroPoint(j2);
					else electroLat[i][j][k] = new ElectroPoint(0.);
				}
			}
		}
		return new ElectroCube(electroLat, precision);
	}









}
