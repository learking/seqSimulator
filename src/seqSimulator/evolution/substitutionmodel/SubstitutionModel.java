package seqSimulator.evolution.substitutionmodel;

import java.util.List;

import seqSimulator.evolution.datatype.MutableSequence;

public interface SubstitutionModel {

	public double getSubstitutionRate(MutableSequence seqI, MutableSequence seqJ, int[] codonArrayI, int differPosition, int differCodon) throws Exception;
	
	public void parseParameters(String line);
	void parseAdditionalInfo(int sectionNr, List<String> lines);
	
	// for debugging
	public void printParameters();
	
	public void isSolventNull();
	
	public void isStructNull(); 
}
