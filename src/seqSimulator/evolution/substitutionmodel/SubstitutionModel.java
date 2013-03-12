package seqSimulator.evolution.substitutionmodel;

public interface SubstitutionModel {

	public void parseParameters(String line);
	public void parseAdditionalInfo(int sectionNr);
	
	// for debugging
	public void printParameters();
}
