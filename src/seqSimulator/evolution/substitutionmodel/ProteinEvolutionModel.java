package seqSimulator.evolution.substitutionmodel;

import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import seqSimulator.evolution.datatype.MutableSequence;
import seqSimulator.evolution.proteinstruct.InputStructure;
import seqSimulator.evolution.proteinstruct.SolventAccessibility;
import seqSimulator.evolution.proteinstruct.StructureEnv;

public class ProteinEvolutionModel implements SubstitutionModel {
	
	// define variables here
    double[] frequencies;
    double kappa;
    int interactionRange;
    
    InputStructure inputStructure;
    StructureEnv structureEnv;
    SolventAccessibility solventAccessibility;
    double scalingFactor;
    
	static Pattern parameters_pattern = Pattern.compile("(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+");
    
    public ProteinEvolutionModel(){
    	System.out.println("protein evolution model ");
    	inputStructure = new InputStructure();
    }
    
	public double getSubstitutionRate(MutableSequence seqI, MutableSequence seqJ, int[] codonArrayI, int differPosition, int differCodon) throws Exception{

		scalingFactor = getScalingFactor();
				
		double substitutionRate = 0;

		// find where these two sequences differ (both location and value)
		//int differPosition = getDifferPosition(seqI, seqJ);
		double logTAU = getLogTAU(seqI, seqJ, codonArrayI, differPosition, differCodon);

		substitutionRate = logTAU / (1 - 1/Math.exp(logTAU));
		
		if(isTransition(seqI.getNucleotide(differPosition), seqJ.getNucleotide(differPosition))){
			double k = kappa;
			substitutionRate = scalingFactor * substitutionRate * frequencies[seqJ.getNucleotide(differPosition)] * k;
		}else{
			substitutionRate = scalingFactor * substitutionRate * frequencies[seqJ.getNucleotide(differPosition)];
		}

		return substitutionRate;
	}
	
	double getScalingFactor(){
		double u;
		double freqSquareSum = 0;
		double[] freqs = frequencies;
		for (int i = 0 ; i < freqs.length ; i++) {
			freqSquareSum += freqs[i] * freqs[i];
		}
		u = 1.0 / (1.0 - freqSquareSum);
		return u; 
	}
	
    boolean isTransition(int firstNucleotide, int secondNucleotide){
    	int sum = firstNucleotide + secondNucleotide;
    	if(sum == 2 || sum == 4){
    		return true;
    	}else{
    		return false;
    	}
    }
    
    double getLogTAU(MutableSequence seqI, MutableSequence seqJ, int[] codonArrayI, int differPosition, int differCodon) throws Exception{
    	double logTAU;
    	double neutralSeqProbRatio = Math.log(frequencies[seqJ.getSequence()[differPosition]] / frequencies[seqI.getSequence()[differPosition]]);

    	int codonDifferPosition = differPosition / 3;
    	//int differCodon = getDifferCodon(seqJ, codonDifferPosition);

    	double structBasedSeqProbRatio = getStructBasedSeqProbRatio(codonArrayI, differCodon, codonDifferPosition);

    	logTAU = structBasedSeqProbRatio - neutralSeqProbRatio;
    	return logTAU;
    }
    
    // range 10 version of getStructBasedSeqProbRatio
    double getStructBasedSeqProbRatio(int[] codonArrayI, int differCodon, int codonDifferPosition) throws Exception{
    	
    	double structBasedSeqProbRatio = 0;
    	
    	double firstOrderRatio = inputStructure.getFirstOrderLogProb(codonDifferPosition, differCodon) - inputStructure.getFirstOrderLogProb(codonDifferPosition, codonArrayI[codonDifferPosition]);
    	//System.out.println("first order:" + firstOrderRatio);

    	double interactionRatio = 0;
    	
    	// here, m refers to a codon position, it should in
    	
		for (int m = getInteractionRangeLeftBound(codonDifferPosition); m < codonDifferPosition; m++) {
			interactionRatio += inputStructure.getInteractionLogProb(m, codonDifferPosition, codonArrayI[m], differCodon) - inputStructure.getInteractionLogProb(m, codonDifferPosition, codonArrayI[m], codonArrayI[codonDifferPosition]);
		}
		
		for (int n = codonDifferPosition + 1 ; n < getInteractionRangeRightBound(codonDifferPosition, codonArrayI.length); n++) {
			interactionRatio += inputStructure.getInteractionLogProb(codonDifferPosition, n, differCodon, codonArrayI[n]) - inputStructure.getInteractionLogProb(codonDifferPosition, n, codonArrayI[codonDifferPosition], codonArrayI[n]);
		}	
    	
    	structBasedSeqProbRatio = firstOrderRatio + interactionRatio;
    	//System.out.println(interactionRatio);
    	return structBasedSeqProbRatio;
    }
    
    int getInteractionRangeLeftBound(int codonDifferPosition){
    	int leftBound = -1;
    	if (codonDifferPosition <= interactionRange) {
    		leftBound = 0;
    	}else{
    		leftBound = codonDifferPosition - interactionRange;
    	}
    	return leftBound;
    }

    int getInteractionRangeRightBound(int codonDifferPosition, int codonArrayLength){
    	int rightBound = -1;
    	if (codonDifferPosition + interactionRange + 1 >= codonArrayLength) {
    		rightBound = codonArrayLength;
    	}else{
    		rightBound = codonDifferPosition + interactionRange + 1;
    	}
    	return rightBound;
    }

	@Override
	public void parseParameters(String line) {
		Matcher matcher = parameters_pattern.matcher(line);
		matcher.find();
		this.setFrequencies(Double.parseDouble(matcher.group(1)), Double.parseDouble(matcher.group(2)), Double.parseDouble(matcher.group(3)), Double.parseDouble(matcher.group(4)));
		this.setKappa(Double.parseDouble(matcher.group(5)));
		this.setInteractionRange(Integer.parseInt(matcher.group(6)));
	}
    
	void setFrequencies(double freqA, double freqC, double freqG, double freqT){
		frequencies = new double [4];
		frequencies[0] = freqA;
		frequencies[1] = freqC;
		frequencies[2] = freqG;
		frequencies[3] = freqT;
	}
	
	public void setKappa(double inputKappa){
		kappa = inputKappa;
	}
	
	public void setInteractionRange(int inputRange){
		interactionRange = inputRange;
	}

	// for Junit test only
	/*
	public void setStructEnv(StructureEnv inputStructure) {
		structureEnv = inputStructure;
	}
	
	public void setSolventAccessibility(SolventAccessibility inputSolventAccessibility) {
		solventAccessibility = inputSolventAccessibility;
	}
	*/
	
	// for debugging purpose
	@Override
	public void printParameters(){
		System.out.println("kappa:" + kappa + " range:" + interactionRange );
		System.out.println("frequencies:" + frequencies[0] + " " + frequencies[1] + " " + frequencies[2] + " " + frequencies[3]);
	}

	@Override
	public void parseAdditionalInfo(int sectionNr, List<String> lines) {
		if(sectionNr == 3){
			this.parseStructEnv(lines);
		}
		
		if(sectionNr == 4){
			this.parseSolventAccessibility(lines);
		}
		
		if(sectionNr == 5){
			this.parseFirstOrderTerms(lines);
		}

		if(sectionNr == 6){
			this.parseInteractionTerms(lines);
		}
	}

	void parseFirstOrderTerms(List<String> lines) {
		System.out.println("start parsing firstOrderTerms:");
		if(lines.size() == 1) {
			String[] solventEntriesStr = lines.get(0).split("\\s+");
			int[] solventEntries = new int[solventEntriesStr.length];
			for (int item = 0; item < solventEntries.length; item++) {
				solventEntries[item] = Integer.parseInt(solventEntriesStr[item]);
			}
			
			inputStructure.setFirstOrderTerms(solventEntries);
			
		}else {
			System.err.println("should be one line!");
		}
		System.out.println("Successfully parsing firstOrderTerms!");
	}

	void parseInteractionTerms(List<String> lines) {
		System.out.println("start parsing interactionTerms:");
		int matrixDim = lines.get(0).split("\\s+").length;
		int[][] matrix = new int[matrixDim][matrixDim];
		for (int i = 0; i < lines.size(); i++) {
			String[] matrixEntriesStr = lines.get(i).split("\\s+");
			int[] matrixEntries = new int[matrixEntriesStr.length];
			for (int item = 0; item < matrixEntries.length; item++) {
				matrixEntries[item] = Integer.parseInt(matrixEntriesStr[item]);
				//System.out.print(matrixEntries[item] + "|");
			}
			matrix[i] = matrixEntries;
		}
		inputStructure.setInteractionTerms(matrix);
		System.out.println("Successfully parsing interactionTerms!");
	}
	
	void parseStructEnv(List<String> lines) {
		System.out.println("start parsing structEnv:");
		structureEnv = new StructureEnv();
		
		for (int i = 0; i < lines.size(); i++) {
			String[] matrixEntriesStr = lines.get(i).split("\\s+");
			double[] matrixEntries = new double[matrixEntriesStr.length];
			for (int item = 0; item < matrixEntries.length; item++) {
				matrixEntries[item] = Double.parseDouble(matrixEntriesStr[item]);
				//System.out.print(matrixEntries[item] + "|");
			}
			structureEnv.addEnv(getNewEnv(matrixEntries));
		}
		
		structureEnv.calculateLogStructEnv();
		
		// debug
		//System.out.println(structureEnv.getLogProb(9, 43, 34));
		inputStructure.setStructEnv(structureEnv);
		
		System.out.println("Parsing structEnv successful!");
	}
	
	double[][] getNewEnv(double[] matrixEntries){
		double[][] newEnv = new double[61][61];
		for(int i = 0; i< 61; i++) {
			for(int j = 0; j < 61; j++) {
				newEnv[i][j] = matrixEntries[i*61 + j];
			}
		}
		return newEnv;
	}

	void parseSolventAccessibility(List<String> lines) {
		System.out.println("start parsing SolventAccessibility:");
		solventAccessibility = new SolventAccessibility();
		
		for (int i = 0; i < lines.size(); i++) {
			String[] matrixEntriesStr = lines.get(i).split("\\s+");
			double[] matrixEntries = new double[matrixEntriesStr.length];
			for (int item = 0; item < matrixEntries.length; item++) {
				matrixEntries[item] = Double.parseDouble(matrixEntriesStr[item]);
				//System.out.print(matrixEntries[item] + "|");
			}
			solventAccessibility.addCategory(matrixEntries);
		}
		// pre-compute log for each category
		solventAccessibility.computeLogCategories();
		
		//solventAccessibility.printItem(9, 60);
		inputStructure.setSolventAccessibility(solventAccessibility);
		System.out.println("Parsing SolventAccessibility successful!");
	}
	
	/*
	public void mockUpSolventAccessibility(){
		double [] codonFreq = new double[]{0.00983798,0.01745548,0.00222048,0.01443315,
				0.00844604,0.01498576,0.00190632,0.01239105,
				0.01064012,0.01887870,
				0.00469486,0.00833007,0.00688776,
				0.01592816,0.02826125,0.00359507,0.02336796,
				0.01367453,0.02426265,0.00308642,0.02006170,
				0.01722686,0.03056552,0.00388819,0.02527326,
				0.00760121,0.01348678,0.00171563,0.01115161,
				0.01574077,0.02792876,0.00355278,0.02309304,
				0.01351366,0.02397721,0.00305010,0.01982568,
				0.01702419,0.03020593,0.00384245,0.02497593,
				0.00751178,0.01332811,0.00169545,0.01102042,
				0.02525082,0.04480239,0.00569924,0.03704508,
				0.02167816,0.03846344,0.00489288,0.03180369,
				0.02730964,0.04845534,0.00616393,0.04006555,
				0.01205015,0.02138052,0.00271978,0.01767859};
		
		for(int i = 0; i < 61; i++) {
				System.out.print(codonFreq[i]+" ");
		}		
	}
	
	public double[][] mockUpMatrix(){

		double [] codonFreq = new double[]{0.00983798,0.01745548,0.00222048,0.01443315,
				0.00844604,0.01498576,0.00190632,0.01239105,
				0.01064012,0.01887870,
				0.00469486,0.00833007,0.00688776,
				0.01592816,0.02826125,0.00359507,0.02336796,
				0.01367453,0.02426265,0.00308642,0.02006170,
				0.01722686,0.03056552,0.00388819,0.02527326,
				0.00760121,0.01348678,0.00171563,0.01115161,
				0.01574077,0.02792876,0.00355278,0.02309304,
				0.01351366,0.02397721,0.00305010,0.01982568,
				0.01702419,0.03020593,0.00384245,0.02497593,
				0.00751178,0.01332811,0.00169545,0.01102042,
				0.02525082,0.04480239,0.00569924,0.03704508,
				0.02167816,0.03846344,0.00489288,0.03180369,
				0.02730964,0.04845534,0.00616393,0.04006555,
				0.01205015,0.02138052,0.00271978,0.01767859};
		
		double [][] mockupmatrix = new double [61][61];

		for (int i = 0; i < 61; i++) {
			for (int j = 0; j < 61; j++) {
				mockupmatrix[i][j] = codonFreq[i] * codonFreq[j];
			}
		}
		
		return mockupmatrix;
	}
	*/
	
	public void isSolventNull() {
		if(solventAccessibility == null) {
			System.out.println("solvent null");
		}
	}

	public void isStructNull() {
		if(structureEnv == null) {
			System.out.println("struct null");
		}
	}
}
