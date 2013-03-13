package seqSimulator.evolution.proteinstruct;

import java.util.List;

public class StructureEnv {

	List<double [][]> structEnv;
	
	public void addEnv(double[][] newEnv) {
		structEnv.add(newEnv);
	}
	
	public double getLogProb(int structEnvNumber, int firstCodonType,
			int secondCodonType) {
		// TODO Auto-generated method stub
		return 0;
	}

}
