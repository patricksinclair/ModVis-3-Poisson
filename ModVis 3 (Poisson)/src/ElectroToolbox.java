import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;


public class ElectroToolbox {

	public static void printToFile(String filename, ArrayList<Double> xData, ArrayList<Double> yData){
		try{

			File file = new File(filename+".txt");

			if(!file.exists()) file.createNewFile();

			FileWriter fw = new FileWriter(file.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);

			for(int i = 0; i < xData.size(); i ++){
				String output = String.valueOf(xData.get(i)) + ", " + String.valueOf(yData.get(i));
				bw.write(output);
				bw.newLine();
			}

			bw.close();
		}catch (IOException e){}
	}

	public static void contourFileWrite(String filename, ArrayList<Double> xData, ArrayList<Double> yData, ArrayList<Double> zData){

		try{

			File file = new File(filename+".txt");

			if(!file.exists()) file.createNewFile();

			FileWriter fw = new FileWriter(file.getAbsoluteFile());
			BufferedWriter bw = new BufferedWriter(fw);

			int n = xData.size();


			for(int i = 0; i < n; i ++){
				for(int j = i*n; j < (i+1)*n; j++){
					String output = String.valueOf(xData.get(i)) + " " + String.valueOf(yData.get(j)) +" " + String.valueOf(zData.get(j));
					bw.write(output);
					bw.newLine();
				}
				bw.newLine();
			}
			bw.close();
		}catch (IOException e){}
	}
}
