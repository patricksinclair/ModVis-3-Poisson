
public class ElectroPoint {

	public double phi, rho, E;
	
	public ElectroPoint(double phi, double rho){
		this.phi = phi;
		this.rho = rho;
		this.E = 0.;
	}
	
	public ElectroPoint(double rho){
		this.rho = rho;
		this.phi = 0.;
		this.E = 0.;
	}
	
	public double getPhi(){
		return phi;
	}
	public void setPhi(double phi){
		this.phi = phi;
	}
	public double getRho(){
		return rho;
	}
	public void setRho(double rho){
		this.rho = rho;
	}
	public double getE(){
		return E;
	}
	public void setE(double E){
		this.E = E;
	}
	
	
	public static ElectroPoint[][][] allzeroes(int N){
		ElectroPoint[][][] cube = new ElectroPoint[N][N][N];
		
		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				for(int k = 0; k < N; k++){
					cube[i][j][k] = new ElectroPoint(0., 0.);
				}
			}
		}
		return cube;
	}
}
