package seqSimulator.core;

import java.util.ArrayList;
import java.util.List;

import seqSimulator.evolution.datatype.CodonUtil;
import seqSimulator.evolution.datatype.MutableSequence;
import seqSimulator.evolution.substitutionmodel.SubstitutionModel;
import seqSimulator.evolution.tree.Node;

public class CoreSimulator {

	int seqNr, codonNr;
	Node m_tree;
	SubstitutionModel m_model;
	
	MutableSequence m_rootSeq;
	// set the time it will take to get our root seq
	double rootSeqTime = 100;
	
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
	}

	MutableSequence createRootSeq() {
		MutableSequence rootSeq = new MutableSequence(codonNr*3);
		
		return rootSeq;
	}
	
	// core of simulation
	MutableSequence simulateSequence(MutableSequence startSeq, double timeLimit) {
		int seqLength = startSeq.getSequence().length;
		MutableSequence endSeq = new MutableSequence(seqLength);
		double totalTime = 0;
		
		List<MutableSequence> seqList = new ArrayList<MutableSequence>();
		List<Double> timeList = new ArrayList<Double>();
		
		while(totalTime < timeLimit) {
			// calculate R(i,j) for each possible j
			
			// pick a J based on PDF
			
			// draw a time from expo dist and add it to totalTime
			
		}
		
		return null;
	}
	
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
}
