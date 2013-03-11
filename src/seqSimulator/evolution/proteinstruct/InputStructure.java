package seqSimulator.evolution.proteinstruct;

public class InputStructure {

	StructureEnv structureEnv = new StructureEnv();
	SolventAccessibility solventAccessibility = new SolventAccessibility();
	
	int [] firstOrderTerms; 
	
	// need a sparse matrix here to store pairwise info, however, since the # of entries is much smaller than ...
	int [][] interactionTerm2EnvMap;
	
	public double getFirstOrderLogProb(int codonPosition, int codonType) {
		int firstOrderTermCategory = firstOrderTerms[codonPosition];
		return solventAccessibility.getLogProb(firstOrderTermCategory, codonType); 
	}
	
	// needs efficiency boost
	public double getInteractionLogProb(int firstCodonPosition, int secondCodonPosition, int firstCodonType, int secondCodonType){
		int structEnvNumber = interactionTerm2EnvMap[firstCodonPosition][secondCodonPosition];
		return structureEnv.getLogProb(structEnvNumber, firstCodonType, secondCodonType);
	}
	
}
