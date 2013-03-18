package test.seqSimulator.evolution.substitutionmodel;


import org.junit.Test;

import seqSimulator.evolution.datatype.MutableSequence;
import seqSimulator.evolution.substitutionmodel.ProteinEvolutionModel;

import junit.framework.TestCase;

public class ProteinEvolutionModelTest extends TestCase {

	MutableSequence seqI, seqJ;
	
	protected void setup() {
		ProteinEvolutionModel ourModel = new ProteinEvolutionModel();

		ourModel.setKappa(1.5);
		ourModel.setInteractionRange(10);
		//ourModel.setSolventAccessibility(inputSolventAccessibility)
		
		seqI = new MutableSequence(9);
		int [] intSeqI = new int [] {1,1,1,0,2,0,3,0,1};
		seqI.setSequence(intSeqI);
		
		seqJ = new MutableSequence(9);
		int [] intSeqJ = new int [] {1,1,1,0,2,0,3,1,1};
		seqJ.setSequence(intSeqJ);
	}
	
	@Test
	public void test() {
		fail("Not yet implemented");
	}

}
