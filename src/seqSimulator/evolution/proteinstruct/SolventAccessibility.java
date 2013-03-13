package seqSimulator.evolution.proteinstruct;

import java.util.ArrayList;
import java.util.List;

public class SolventAccessibility {

	List<double[]> solventCategories;
	
	public SolventAccessibility(){
		solventCategories = new ArrayList<double[]>();
	}
	
	public double getLogProb(int firstOrderTermCategory, int codonType) {
		return solventCategories.get(firstOrderTermCategory)[codonType];
	}

	public void addCategory(double[] matrixEntries) {
		solventCategories.add(matrixEntries);
	}

	public void computeLogCategories(){
		for(int i=0; i < solventCategories.size(); i++){
			computeLogCategory(solventCategories.get(i));
		}
	}
	
	public void computeLogCategory(double[] category){
		for(int i = 0; i < category.length; i++){
			 category[i] = Math.log(category[i]);
		}
	}
	
	// for debugging only
	public void printSum(int categoryNr){
		double[] currentCategory = solventCategories.get(categoryNr);
		double sum = 0;
		for(int i = 0; i < currentCategory.length; i++){
			sum += currentCategory[i];
		}
		System.out.println(sum);
	}
	
	public void printItem(int categoryNr, int item){
		System.out.println(solventCategories.get(categoryNr)[item]);
	}
	
}
