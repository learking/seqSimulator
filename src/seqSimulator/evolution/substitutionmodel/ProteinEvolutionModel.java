package seqSimulator.evolution.substitutionmodel;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import seqSimulator.evolution.datatype.MutableSequence;
import seqSimulator.evolution.proteinstruct.InputStructure;

public class ProteinEvolutionModel implements SubstitutionModel {
	
	// define variables here
    double[] frequencies;
    double kappa;
    int interactionRange;
    
    InputStructure inputStructure;
    double scalingFactor;
    
	static Pattern parameters_pattern = Pattern.compile("(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+");
    
    public ProteinEvolutionModel(){
    	System.out.println("protein evolution model ");
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
	
	void setKappa(double inputKappa){
		kappa = inputKappa;
	}
	
	void setInteractionRange(int inputRange){
		interactionRange = inputRange;
	}

	// for debugging purpose
	@Override
	public void printParameters(){
		System.out.println("kappa:" + kappa + " range:" + interactionRange );
		System.out.println("frequencies:" + frequencies[0] + " " + frequencies[1] + " " + frequencies[2] + " " + frequencies[3]);
	}

	@Override
	public void parseAdditionalInfo(int sectionNr) {
		// TODO Auto-generated method stub
		if(sectionNr == 3){
		}
		
		if(sectionNr == 4){
			
		}
		
		if(sectionNr == 5){
			
		}
	}
}
