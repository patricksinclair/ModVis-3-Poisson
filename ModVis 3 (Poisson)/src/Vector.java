
public class Vector {

	private double x, y, z;
	
	public Vector(double x, double y, double z){
		this.x = x;
		this.y = y;
		this.z = z;
	}
	
	public double getX(){
		return x;
	}
	public void setX(double x){
		this.x = x;
	}
	public double getY(){
		return y;
	}
	public void setY(double y){
		this.y = y;
	}
	public double getZ(){
		return z;
	}
	public void setZ(double z){
		this.z = z;
	}
	
	public double magnitude(){
		return Math.sqrt(getX()*getX() + getY()*getY() + getZ()*getZ());
	}
	
	public double magnitude2D(){
		return Math.sqrt(getX()*getX() + getY()*getY());
	}
	
	public double xyAngle(){
		return Math.atan2(getY(), getX());
	}
	
	
	public static Vector[][][] electroVectorCube(ElectroCube eCube){
		
		int N = eCube.getElectroLattice().length;
		Vector[][][] eVectors = Vector.zeroVectors(N);
		
		for(int i = 1; i < N-1; i++){
			for(int j = 1; j < N-1; j++){
				for(int k = 1; k < N-1; k++){
					double Ex = -0.5*(eCube.getPhiAtPoint(i+1, j, k) - eCube.getPhiAtPoint(i-1, j, k));
					double Ey = -0.5*(eCube.getPhiAtPoint(i, j+1, k) - eCube.getPhiAtPoint(i, j-1, k));
					double Ez = -0.5*(eCube.getPhiAtPoint(i, j, k+1) - eCube.getPhiAtPoint(i, j, k-1));
					eVectors[i][j][k] = new Vector(Ex, Ey, Ez);
				}
			}
		}		
		return eVectors;
	}
	
	
	public static Vector[][][] magnetoVectorCube(ElectroCube bCube){
		
		int N = bCube.getElectroLattice().length;
		Vector[][][] bVectors = Vector.zeroVectors(N);
		
		for(int i = 1; i < N-1; i++){
			for(int j = 1; j < N-1; j++){
				for(int k = 1; k < N-1; k++){
					double Bx = 0.5*(bCube.getPhiAtPoint(i, j+1, k) - bCube.getPhiAtPoint(i, j-1, k));
					double By = -0.5*(bCube.getPhiAtPoint(i+1, j, k) - bCube.getPhiAtPoint(i-1, j, k));
					double Bz = 0.;
					bVectors[i][j][k] = new Vector(Bx, By, Bz);
				}
			}
		}		
		return bVectors;
	}
	
	
	public static Vector[][][] zeroVectors(int N){
		Vector[][][] vectors = new Vector[N][N][N];
		
		for(int i = 0; i < N; i++){
			for(int j = 0; j < N; j++){
				for(int k = 0; k < N; k++){
					vectors[i][j][k] = new Vector(0., 0., 0.);
				}
			}	
		}
		return vectors;
	}
	
}
