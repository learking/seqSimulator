package seqSimulator.core;

import java.util.ArrayList;
import java.util.List;

import seqSimulator.evolution.datatype.CodonUtil;
import seqSimulator.evolution.datatype.MutableSequence;
import seqSimulator.evolution.substitutionmodel.SubstitutionModel;
import seqSimulator.evolution.tree.Node;
import seqSimulator.util.Randomizer;

public class CoreSimulator {

	static CodonUtil codonUtil = new CodonUtil();
	
	int seqNr, codonNr;
	Node m_tree;
	SubstitutionModel m_model;
	
	MutableSequence m_rootSeq;
	// by default, rootSeqTime = 10
	double rootSeqTime = 10;
	
	public CoreSimulator(){
	}
	
	public void init() {
		
	}
	
	public void run() throws Exception{
		/*
		MutableSequence seqI = new MutableSequence(9);
		int [] intSeqI = new int [] {1,1,1,0,2,0,3,0,1};
		seqI.setSequence(intSeqI);
		
		MutableSequence seqJ = new MutableSequence(9);
		int [] intSeqJ = new int [] {1,1,1,0,2,0,3,1,1};
		seqJ.setSequence(intSeqJ);

		CodonUtil codonUtil = new CodonUtil();
		int differCodon = codonUtil.translate(seqJ, 6);
		
		double testRate = m_model.getSubstitutionRate(seqI, seqJ, seqI.toCodonArray(), 7, differCodon);
		System.out.println("the rate is:" + testRate);
	*/
		m_rootSeq = createRootSeq();
		simulateDataSet();
	}

	MutableSequence createRootSeq() throws Exception {
		MutableSequence rootSeq = new MutableSequence(codonNr*3);
		simulateSequence(rootSeq, rootSeqTime);
		return rootSeq;
	}
	
	// core of simulation
	MutableSequence simulateSequence(MutableSequence startSeq, double timeLimit) throws Exception {
		int seqLength = startSeq.getSequence().length;
		MutableSequence endSeq = new MutableSequence(seqLength);
		double totalTime = 0;
		int totalNumOfChanges = 0;
		
		List<MutableSequence> seqList = new ArrayList<MutableSequence>();
		List<Double> timeList = new ArrayList<Double>();
		
		// tmp list, empty every time
		List<MutableSequence> allPossibleSeqs = new ArrayList<MutableSequence>();
		List<Double> seqProbs = new ArrayList<Double>();
		// tmp seq
		MutableSequence tmpSeq = new MutableSequence(seqLength);	
		MutableSequence inputSeq = new MutableSequence(seqLength);
		int[] originalNucleoSeq;
		int[] originalCodonArray;
 		
		while(totalTime < timeLimit) {
			// empty all possibleSeqs and seqProbs
			allPossibleSeqs.clear();
			seqProbs.clear();
			// initialize inputSeq
			if(seqList.size()==0){
				inputSeq.setSequence(startSeq.getSequence());
			}else{
				inputSeq.setSequence(seqList.get(seqList.size()-1).getSequence());
			}
			// initialize originalNucleoSeq and codonArray
			originalNucleoSeq = inputSeq.getSequence();
			originalCodonArray = inputSeq.toCodonArray();
			// initialize awayRate
	    	double awayRate = 0;
	    	
			// calculate R(i,j) for each possible j
			for(int site = 0; site < seqLength; site++){
				tmpSeq.setSequence(originalNucleoSeq);
				int originalNucleotide = originalNucleoSeq[site];
				int[] changableNucleotides = getChangableNucleotides(originalNucleotide);
				for (int i = 0; i < changableNucleotides.length; i++) {
					tmpSeq.mutate(site, changableNucleotides[i]);
					int mutatedCodonStartSite = site - (site%3);
					int differCodon = codonUtil.translate(tmpSeq, mutatedCodonStartSite);
	    			if(differCodon != -1) {
	    				double currentRate = m_model.getSubstitutionRate(inputSeq, tmpSeq, originalCodonArray, site, differCodon);
	    				// 3N probs
	    				seqProbs.add(currentRate);
	    				MutableSequence currentSeq = new MutableSequence(seqLength);
	    				currentSeq.setSequence(tmpSeq.getSequence());
	    				// 3N seqs
	    				allPossibleSeqs.add(currentSeq);
	    				
	    				awayRate += currentRate;
	    				
	    			}
				}
			}
			
			// pick a J based on PDF and add it to seqList
			MutableSequence chosenSeq = new MutableSequence(seqLength);
			chosenSeq.setSequence(sampleSeqFromCDF(seqProbs, allPossibleSeqs, awayRate).getSequence());
			seqList.add(chosenSeq);
			
			// draw a time from expo dist and add it to totalTime and timeList
			// sample a time from expo dist with para awayrate
			double sampledTimeInterval = Randomizer.nextExponential(awayRate);
			totalTime += sampledTimeInterval;
			
			totalNumOfChanges +=1 ;
		}
		
		if(seqList.size() > 0){
			endSeq.setSequence(seqList.get(seqList.size()-1).getSequence());
		}else{
			System.err.println("There is no substitution along this branch!");
		}
		
		System.out.println("Total number:" + totalNumOfChanges);
		System.out.println(endSeq.toString());
		System.out.println(endSeq.toACGT());
		return endSeq;
	}
	
	MutableSequence sampleSeqFromCDF(List<Double> seqProbs, List<MutableSequence> allPossibleSeqs, double awayRate){
		double[] seqProbPDF = new double[seqProbs.size()];
				
		for(int i = 0; i < seqProbs.size(); i++){
			seqProbPDF[i] = seqProbs.get(i);
		}
				
		double [] seqProbPDFNormalized = Randomizer.getNormalized(seqProbPDF);
		
		double cumulativeProb = 0;
		double [] seqProbCDF = new double[seqProbs.size()];
		for (int i = 0 ; i < seqProbPDFNormalized.length ; i++) {
			cumulativeProb += seqProbPDFNormalized[i];
			seqProbCDF[i] = cumulativeProb;
		}
		
		int chosenSeqNr = Randomizer.randomChoice(seqProbCDF);
		
		return allPossibleSeqs.get(chosenSeqNr);
	}
	
    int[] getChangableNucleotides(int originalNucleotide){
    	int [] changableNucleotides;
    	switch(originalNucleotide){
    	case 0:
    		changableNucleotides = new int[]{1,2,3};
    		break;
    	case 1:
    		changableNucleotides = new int[]{0,2,3};
    		break;
    	case 2:
    		changableNucleotides = new int[]{0,1,3};
    		break;
    	case 3:
    		changableNucleotides = new int[]{0,1,2};
    		break;
    	default:
    		changableNucleotides = new int[]{-1,-1,-1};
    	} 
    	return changableNucleotides;
    }
	
    void simulateDataSet() throws Exception{
    	// set seq for root
    	m_tree.setSequence(m_rootSeq);
    	// starting from root traverse the entire tree and simulate seq for each node
    	traverse(m_tree);
    	m_tree.printStructure();
    }
    
    void traverse(Node node) throws Exception{
    	System.out.print("Node:" + node.getNr());
    	for(int iChild = 0; iChild < 2; iChild++){
    		Node child = (iChild == 0 ? node.getLeft() : node.getRight());
    		System.out.print("Child:" + child.getNr() + " branch length:" + child.getLength());
    		if(child.getLength()==0){
    			child.setSequence(node.getSequence());
    		}else{
    		// currentseq = simulateSequence(parentseq, child.getLength);
    			MutableSequence currentSeq = simulateSequence(node.getSequence(), child.getLength());
    			child.setSequence(currentSeq);
    		}
    		
    		if(!child.isLeaf()){
    			traverse(child);
    		}
    	}
    }
    
    // traverse the entire tree
    /*
    void traverse(Node node, int[] parentSequence, int[] category, Alignment alignment) throws Exception {
        for (int iChild = 0; iChild < 2; iChild++) {
            Node child = (iChild == 0 ? node.getLeft() : node.getRight());
            for (int i = 0; i < m_categoryCount; i++) {
                getTransitionProbabilities(m_tree, child, i, m_probabilities[i]);
            }

            int[] seq = new int[m_sequenceLength];
            double[] cProb = new double[m_stateCount];
            for (int i = 0; i < m_sequenceLength; i++) {
                System.arraycopy(m_probabilities[category[i]], parentSequence[i] * m_stateCount, cProb, 0, m_stateCount);
                seq[i] = Randomizer.randomChoicePDF(cProb);
            }

            if (child.isLeaf()) {
                alignment.m_pSequences.setValue(intArray2Sequence(seq, child), alignment);
            } else {
                traverse(child, seq, category, alignment);
            }
        }
    } // traverse
    */
    
	// setter
	public void setSeqNr(int inputSeqNr) {
		seqNr = inputSeqNr;
	}
	
	public void setCodonNr(int inputCodonNr) {
		codonNr = inputCodonNr;
	}
	
	public void setTree(Node inputTree) {
		m_tree = inputTree;
	}
	
	public void setModel(SubstitutionModel inputModel) {
		m_model = inputModel;
	}
	
	public void setRootSeqTime(double inputRootSeqTime) {
		rootSeqTime = inputRootSeqTime;
	}
}
