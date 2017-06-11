import java.util.InputMismatchException;
import java.util.Scanner;


public class PoissonMain {
	public static void main(String args[]){

		
	  Scanner keyboard = new Scanner(System.in);
		/*int choice;
		boolean quit = false;

		while(!quit){
			try{
				printOptions();
				choice = keyboard.nextInt();

				mainmenu:
					switch(choice){
					
					case 1:
						ElectroCube.pointPotentialLineGraph();
						break;
						
					case 2:
						ElectroCube.gaussianLineGraph();
						break;
						
					case 3:
						ElectroCube.ELineGraph();
						break;
						
					case 4:
						ElectroCube.pointPotentialContourPlot();
						break;
						
					case 5:
						ElectroCube.quadropolePotentialContourPlot();
						break;
						
					case 6:
						ElectroCube.gaussianPotentialContourPlot();
						break;
						
					case 7:
						ElectroCube.EFieldContourPlot();
						break;
						
					case 8:
						ElectroCube.quadropoleEFieldContourPlot();
						//ElectroCube.visualiseEField();
						//ElectroCube.visualiseBField();
						ElectroCube.relaxedGS();
						
					

					
					break;
					}

			}catch(InputMismatchException e){
				System.out.println("Please enter a valid number.");
				keyboard.nextLine();
			}
		}
	*/
		
		
		//ElectroCube.pointPotentialLineGraph(keyboard);
		//ElectroCube.gaussianLineGraph(keyboard);
		//ElectroCube.ELineGraph(keyboard);
		//ElectroCube.pointPotentialContourPlot(keyboard);
		//ElectroCube.quadropolePotentialContourPlot(keyboard);
		//ElectroCube.gaussianPotentialContourPlot(keyboard);
		//ElectroCube.EFieldContourPlot(keyboard);
		ElectroCube.quadropoleEFieldContourPlot(keyboard);
		//ElectroCube.visualiseEField();
		//ElectroCube.visualiseBField();
		//ElectroCube.relaxedGS();
	}

	public static void printOptions(){
		System.out.println("Please select one of the following: ");
		System.out.println("1. Line graph of the radial potential for a point charge distribution.");
		System.out.println("2. Line graph of the radial potential for a gaussian charge.");
		System.out.println("3. Line graph of the radial E field strength for a point charge.");
		System.out.println("4. Line graph of the radial E field strength for a gausian charge.");
		System.out.println("5. Contour plot of the potential for a point charge.");
		System.out.println("6. Contour plot of the potential for a gaussian charge.");
		System.out.println("7. Contour plot of the potential for a quadrapole.");
		System.out.println("8. Contour plot of the radial E field strength for a point charge.");
		System.out.println("9. Contour plot of the radial E field strength for a gaussian charge.");
		System.out.println("10. Contour plot of the radial E field strength for a quadrapole.");
		System.out.println("11. Contour plot of the radial E field strength for a gaussian quadropole.");
		System.out.println("12. Visualisation of E field.");
		System.out.println("13. Visualisaton of B field.");
		System.out.println("0. quit");
	}

	public static void printLineOptions(){
		System.out.println("Please select one of the following: ");
		System.out.println("1. Line graph of the potential of a point charge.");
		System.out.println("2. Line graph of the potential of a gaussian charge distribution.");
		System.out.println("3. Line graph of radial E field strength for a point charge.");
		System.out.println("4. Line graph of radial E field strength for a point charge.");
	}
}
