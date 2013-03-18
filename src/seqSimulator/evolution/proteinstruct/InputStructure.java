package seqSimulator.evolution.proteinstruct;

import java.util.Arrays;

public class InputStructure {

	StructureEnv structureEnv;
	SolventAccessibility solventAccessibility;
	
	int [] firstOrderTerms; 

	public InputStructure() {}
	
	// need a sparse matrix here to store pairwise info, however, since the # of entries is much smaller than ...
	int [][] interactionTerm2EnvMap;
	
	public double getFirstOrderLogProb(int codonPosition, int codonType) {
		if(solventAccessibility==null) {
			System.out.println("solvent is null");
		}
		if(structureEnv==null) {
			System.out.println("struct is null");
		}
		int firstOrderTermCategory = firstOrderTerms[codonPosition];
	    	//System.out.println(solventAccessibility.getLogProb(firstOrderTermCategory, codonType));
		return solventAccessibility.getLogProb(firstOrderTermCategory, codonType); 
	}
	
	// needs efficiency boost
	public double getInteractionLogProb(int firstCodonPosition, int secondCodonPosition, int firstCodonType, int secondCodonType){
		int structEnvNumber = interactionTerm2EnvMap[firstCodonPosition][secondCodonPosition];
		return structureEnv.getLogProb(structEnvNumber, firstCodonType, secondCodonType);
	}

	public void setStructEnv(StructureEnv inputStructEnv) {
		structureEnv = inputStructEnv;
	}
	
	public void setSolventAccessibility(SolventAccessibility inputSolventAccessibility) {
		solventAccessibility = inputSolventAccessibility;
	}
	
	public void setFirstOrderTerms(int[] solventEntries) {
		firstOrderTerms = Arrays.copyOfRange(solventEntries, 0, solventEntries.length);
	}

	public void setInteractionTerms(int[][] matrix) {
		// TODO Auto-generated method stub
		interactionTerm2EnvMap = matrix;
	}
	
}
