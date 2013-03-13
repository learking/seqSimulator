package seqSimulator.evolution.substitutionmodel;

import java.util.List;

public interface SubstitutionModel {

	public void parseParameters(String line);
	void parseAdditionalInfo(int sectionNr, List<String> lines);
	
	// for debugging
	public void printParameters();
}
