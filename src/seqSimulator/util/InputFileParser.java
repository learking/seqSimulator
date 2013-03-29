package seqSimulator.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import seqSimulator.core.CoreSimulator;
import seqSimulator.evolution.datatype.CodonUtil;
import seqSimulator.evolution.datatype.MutableSequence;
import seqSimulator.evolution.substitutionmodel.ProteinEvolutionModel;
import seqSimulator.evolution.substitutionmodel.SubstitutionModel;
import seqSimulator.evolution.tree.Node;

public class InputFileParser {
    /**
     * default beast.tree branch length, used when that info is not in the Newick beast.tree
     */
    final static double DEFAULT_LENGTH = 0.001f;

    /**
     * labels of leafs *
     */
    List<String> m_sLabels = null;
    
	CoreSimulator m_coreSimulator;
	
	SubstitutionModel m_substitutionModel;
	
	Node m_tree;

	int seqNr;
	int codonNr;
	
	static Pattern endFile_pattern = Pattern.compile("//");
	static Pattern seqCodonNr_pattern = Pattern.compile("(\\d+)\\s+(\\d+)\\s+(\\d+)");
	static Pattern label_pattern = Pattern.compile("[a-zA-Z]+");
	static Pattern modelType = Pattern.compile("([0-9]+)\\s\\*");
	
    /**
     * used to make sure all taxa only occur once in the tree *
     */
    List<Boolean> m_bTaxonIndexInUse = new ArrayList<Boolean>();
    boolean createUnrecognizedTaxa = false;
    int m_nOffset = 1;
    double m_nThreshold = 0.0;
    Boolean m_bIsLabelledNewick = false;
    Boolean m_bAllowSingleChild = false;
    double scale = 1.0;
	
    public InputFileParser() {
    	m_sLabels = new ArrayList<String>();
    }
    
	public CoreSimulator parseFile(File inputFile) throws Exception{
		System.out.println("start of the process");

		parse(inputFile);
		
		// test whether getSubstitutionRate generates the same result as the original Model in evoprotein
		/*
		MutableSequence seqI = new MutableSequence(9);
		int [] intSeqI = new int [] {1,1,1,0,2,0,3,0,1};
		seqI.setSequence(intSeqI);
		
		MutableSequence seqJ = new MutableSequence(9);
		int [] intSeqJ = new int [] {1,1,1,0,2,0,3,1,1};
		seqJ.setSequence(intSeqJ);

		CodonUtil codonUtil = new CodonUtil();
		int differCodon = codonUtil.translate(seqJ, 6);
		
		double testRate = m_substitutionModel.getSubstitutionRate(seqI, seqJ, seqI.toCodonArray(), 7, differCodon);
		System.out.println("the rate is:" + testRate);
		*/
		// result is OK
		
		m_coreSimulator.setSeqNr(seqNr);
		m_coreSimulator.setCodonNr(codonNr);
		m_coreSimulator.setTree(m_tree);
		m_coreSimulator.setModel(m_substitutionModel);
		
        if (m_coreSimulator != null){
            return m_coreSimulator;
        }
        else {
            throw new Exception("CoreSimulator cannot be created from input file.");
        }
	}
	
	void initLabels(String line) {
	    Matcher matcher = label_pattern.matcher(line);
	    // Check all occurance
	    while (matcher.find()) {
	    	m_sLabels.add(matcher.group());
	    }
	}
	
	void parse(File inputFile) throws Exception{
		
		m_coreSimulator = new CoreSimulator();
		int sectionCount = 0;
		
		FileReader reader = new FileReader(inputFile);
		BufferedReader bufferedReader = new BufferedReader(reader);

		Matcher matcher;
		
	    String line = "";
	    while ((line = bufferedReader.readLine()) != null) {
	    	
	    	matcher = endFile_pattern.matcher(line);
	    	if(!matcher.find()){
	    		
	    		if(sectionCount == 0){
	    			matcher = seqCodonNr_pattern.matcher(line);
		    		if(matcher.find()){
		    			seqNr = Integer.parseInt(matcher.group(1));
		    			codonNr = Integer.parseInt(matcher.group(2));
		    			m_coreSimulator.setRootSeqTime((double) Integer.parseInt(matcher.group(3)));
		    			System.out.println("seq Nr:" + seqNr + " codonNr:" + codonNr);
		    		}
	    		}
	    		
	    		if(sectionCount == 1 && !line.trim().isEmpty()){
	    			System.out.println("Newick tree parse started!");
	    			line = line.trim();
	    			initLabels(line);
	    			m_tree = parseNewick(line);
	    			m_tree.printStructure();
	    			System.out.println("Newick tree parse successful!");
	    		}
	    		
	    		if(sectionCount == 2 && !line.trim().isEmpty()){
	    			// parse model type
	    			matcher = modelType.matcher(line.trim());
	    			matcher.find();
	    			m_substitutionModel = createModel(Integer.parseInt(matcher.group(1)));
	    			// parse model parameters based on model type
	    			line = bufferedReader.readLine();
	    			m_substitutionModel.parseParameters(line);
	    			m_substitutionModel.printParameters();
	    			System.out.println("Model parse successful!");
	    		}	    		

	    		if(sectionCount >= 3){
	    			// parse structure environment?
	    			List<String> lines = new ArrayList<String>();
	    			lines.add(line);
	    			
	    			while(!(line = bufferedReader.readLine()).trim().isEmpty()) {
	    				lines.add(line);
	    			}
	    			
	    			m_substitutionModel.parseAdditionalInfo(sectionCount, lines);
	    		}	
	    		
	    		if(line.trim().isEmpty()){
	    			sectionCount++;
	    		}
	    		
	    	}else{
	    		break;
	    	}
	    	
	    }
	    
	}
	
	SubstitutionModel createModel(int modelType){
		SubstitutionModel newModel;
		 switch (modelType) {
         	case 0:  newModel = new ProteinEvolutionModel();
         	break;
         	default: newModel = null;
         	break;
		 }
		 return newModel;
	}
	
    void processMetadata(Node node) throws Exception {
        if (node.isLeaf()) {
            if (m_sLabels != null) {
                node.setID(m_sLabels.get(node.getNr()));
            }
        } else {
            processMetadata(node.getLeft());
            if (node.getRight() != null) {
                processMetadata(node.getRight());
            }
        }
    }
	
    /**
     * Try to map sStr into an index. First, assume it is a number.
     * If that does not work, look in list of labels to see whether it is there.
     */
    private int getLabelIndex(String sStr) throws Exception {
        if (!m_bIsLabelledNewick && m_sLabels == null) {
            try {
                int nIndex = Integer.parseInt(sStr) - m_nOffset;
                checkTaxaIsAvailable(sStr, nIndex);
                return nIndex;
            } catch (Exception e) {
                System.out.println(e.getClass().getName() + " " + e.getMessage() + ". Perhaps taxa or taxonset is not specified?");
            }
        }
        // look it up in list of taxa
        for (int nIndex = 0; nIndex < m_sLabels.size(); nIndex++) {
            if (sStr.equals(m_sLabels.get(nIndex))) {
                checkTaxaIsAvailable(sStr, nIndex);
                return nIndex;
            }
        }

        // if createUnrecognizedTaxa==true, then do it now, otherwise labels will not be populated and
        // out of bounds error will occur in m_sLabels later.
        if (createUnrecognizedTaxa) {
            m_sLabels.add(sStr);
            int nIndex = m_sLabels.size() - 1;
            checkTaxaIsAvailable(sStr, nIndex);
            return nIndex;
        }

        // finally, check if its an integer number indicating the taxon id
        try {
            int nIndex = Integer.parseInt(sStr) - m_nOffset;
            checkTaxaIsAvailable(sStr, nIndex);
            return nIndex;
        } catch (NumberFormatException e) {
        	// apparently not a number
        }
        throw new Exception("Label '" + sStr + "' in Newick beast.tree could not be identified. Perhaps taxa or taxonset is not specified?");
    }
	
    void checkTaxaIsAvailable(String sStr, int nIndex) throws Exception {
        while (nIndex + 1 > m_bTaxonIndexInUse.size()) {
            m_bTaxonIndexInUse.add(false);
        }
        if (m_bTaxonIndexInUse.get(nIndex)) {
            throw new Exception("Duplicate taxon found: " + sStr);
        }
        m_bTaxonIndexInUse.set(nIndex, true);
    }
    
    char[] m_chars;
    int m_iTokenStart;
    int m_iTokenEnd;
    final static int COMMA = 1;
    final static int BRACE_OPEN = 3;
    final static int BRACE_CLOSE = 4;
    final static int COLON = 5;
    final static int SEMI_COLON = 8;
    final static int META_DATA = 6;
    final static int TEXT = 7;
    final static int UNKNOWN = 0;
    
    int nextToken() {
        m_iTokenStart = m_iTokenEnd;
        while (m_iTokenEnd < m_chars.length) {
            // skip spaces
            while (m_iTokenEnd < m_chars.length && (m_chars[m_iTokenEnd] == ' ' || m_chars[m_iTokenEnd] == '\t')) {
                m_iTokenStart++;
                m_iTokenEnd++;
            }
            if (m_chars[m_iTokenEnd] == '(') {
                m_iTokenEnd++;
                return BRACE_OPEN;
            }
            if (m_chars[m_iTokenEnd] == ':') {
                m_iTokenEnd++;
                return COLON;
            }
            if (m_chars[m_iTokenEnd] == ';') {
                m_iTokenEnd++;
                return SEMI_COLON;
            }
            if (m_chars[m_iTokenEnd] == ')') {
                m_iTokenEnd++;
                return BRACE_CLOSE;
            }
            if (m_chars[m_iTokenEnd] == ',') {
                m_iTokenEnd++;
                return COMMA;
            }
            if (m_chars[m_iTokenEnd] == '[') {
                m_iTokenEnd++;
                while (m_iTokenEnd < m_chars.length && m_chars[m_iTokenEnd - 1] != ']') {
                    m_iTokenEnd++;
                }
                return META_DATA;
            }
            while (m_iTokenEnd < m_chars.length && (m_chars[m_iTokenEnd] != ' ' && m_chars[m_iTokenEnd] != '\t'
                    && m_chars[m_iTokenEnd] != '(' && m_chars[m_iTokenEnd] != ')' && m_chars[m_iTokenEnd] != '['
                    && m_chars[m_iTokenEnd] != ':' && m_chars[m_iTokenEnd] != ',' && m_chars[m_iTokenEnd] != ';')) {
                m_iTokenEnd++;
            }
            return TEXT;
        }
        return UNKNOWN;
    }
    
    public Node parseNewick(String sStr) throws Exception {
        // get rid of initial and terminal spaces
        sStr = sStr.replaceAll("^\\s+", "");
        sStr = sStr.replaceAll("\\s+$", "");
        
        try {
            m_chars = sStr.toCharArray();
            if (sStr == null || sStr.length() == 0) {
                return null;
            }
            m_iTokenStart = 0;
            m_iTokenEnd = 0;
            Vector<Node> stack = new Vector<Node>();
            Vector<Boolean> isFirstChild = new Vector<Boolean>();
            stack.add(newNode());
            isFirstChild.add(true);
            stack.lastElement().setHeight(DEFAULT_LENGTH);
            boolean bIsLabel = true;
            while (m_iTokenEnd < m_chars.length) {
            	
                switch (nextToken()) {
                case BRACE_OPEN: {
                    Node node2 = newNode();
                    node2.setHeight(DEFAULT_LENGTH);
                    stack.add(node2);
                    isFirstChild.add(true);
                    bIsLabel = true;
                }
                break;
                case BRACE_CLOSE: {
                    if (isFirstChild.lastElement()) {
                        if (m_bAllowSingleChild) {
                            // process single child nodes
                            Node left = stack.lastElement();
                            stack.remove(stack.size() - 1);
                            isFirstChild.remove(isFirstChild.size() - 1);
                            Node parent = stack.lastElement();
                            parent.setLeft(left);
                            //parent.setRight(null);
                            left.setParent(parent);
                            break;
                        } else {
                            // don't know how to process single child nodes
                            throw new Exception("Node with single child found.");
                        }
                    }
                    // process multi(i.e. more than 2)-child nodes by pairwise merging.
                    while (isFirstChild.get(isFirstChild.size() - 2) == false) {
                        Node right = stack.lastElement();
                        stack.remove(stack.size() - 1);
                        isFirstChild.remove(isFirstChild.size() - 1);
                        Node left = stack.lastElement();
                        stack.remove(stack.size() - 1);
                        isFirstChild.remove(isFirstChild.size() - 1);
                        Node dummyparent = newNode();
                        dummyparent.setHeight(DEFAULT_LENGTH);
                        dummyparent.setLeft(left);
                        left.setParent(dummyparent);
                        dummyparent.setRight(right);
                        right.setParent(dummyparent);
                        stack.add(dummyparent);
                        isFirstChild.add(false);
                    }
                    // last two nodes on stack merged into single parent node
                    Node right = stack.lastElement();
                    stack.remove(stack.size() - 1);
                    isFirstChild.remove(isFirstChild.size() - 1);
                    Node left = stack.lastElement();
                    stack.remove(stack.size() - 1);
                    isFirstChild.remove(isFirstChild.size() - 1);
                    Node parent = stack.lastElement();
                    parent.setLeft(left);
                    left.setParent(parent);
                    parent.setRight(right);
                    right.setParent(parent);
                }
                break;
                case COMMA: {
                    Node node2 = newNode();
                    node2.setHeight(DEFAULT_LENGTH);
                    stack.add(node2);
                    isFirstChild.add(false);
                    bIsLabel = true;
                }
                break;
                case COLON:
                    bIsLabel = false;
                    break;
                case TEXT:
                    if (bIsLabel) {
                        String sLabel = sStr.substring(m_iTokenStart, m_iTokenEnd);
                        stack.lastElement().setNr(getLabelIndex(sLabel));
                    } else {
                        String sLength = sStr.substring(m_iTokenStart, m_iTokenEnd);
                        stack.lastElement().setHeight(Double.parseDouble(sLength));
                    }
                    break;
                    
                /*
                case META_DATA:
                    if (stack.lastElement().m_sMetaData == null) {
                        stack.lastElement().m_sMetaData = sStr.substring(m_iTokenStart + 1, m_iTokenEnd - 1);
                    } else {
                        stack.lastElement().m_sMetaData += " " + sStr.substring(m_iTokenStart + 1, m_iTokenEnd - 1);
                    }
                    break;
                */
                    
                case SEMI_COLON:
                    //System.err.println(stack.lastElement().toString());
                    Node tree = stack.lastElement();
                    tree.sort();
                    // at this stage, all heights are actually lengths
                    convertLengthToHeight(tree);
                    int n = tree.getLeafNodeCount();
                    tree.labelInternalNodes(n);
                    processMetadata(tree);
                    return stack.lastElement();
                default:
                    throw new Exception("parseNewick: unknown token");
                }         
            }
            Node tree = stack.lastElement();
            tree.sort();
            // at this stage, all heights are actually lengths
            convertLengthToHeight(tree);
            int n = tree.getLeafNodeCount();
            if (tree.getNr() == 0) {
                tree.labelInternalNodes(n);
            }
            processMetadata(tree);
            return tree;   
        } catch (Exception e) {
            System.err.println(e.getClass().toString() + "/" + e.getMessage() + ": " + sStr.substring(Math.max(0, m_iTokenStart - 100), m_iTokenStart) + " >>>" + sStr.substring(m_iTokenStart, m_iTokenEnd) + " <<< ...");
            throw new Exception(e.getMessage() + ": " + sStr.substring(Math.max(0, m_iTokenStart - 100), m_iTokenStart) + " >>>" + sStr.substring(m_iTokenStart, m_iTokenEnd) + " <<< ...");
        }   
    }
    
    Node newNode() {
        return new Node();
    }
    
    void convertLengthToHeight(Node node) {
        double fTotalHeight = convertLengthToHeight(node, 0);
        offset(node, -fTotalHeight);
    }

    double convertLengthToHeight(Node node, double fHeight) {
        double fLength = node.getHeight();
        node.setHeight((fHeight - fLength) * scale);
        if (node.isLeaf()) {
            return node.getHeight();
        } else {
            double fLeft = convertLengthToHeight(node.getLeft(), fHeight - fLength);
            if (node.getRight() == null) {
                return fLeft;
            }
            double fRight = convertLengthToHeight(node.getRight(), fHeight - fLength);
            return Math.min(fLeft, fRight);
        }
    }

    void offset(Node node, double fDelta) {
        node.setHeight(node.getHeight() + fDelta);
        if (node.isLeaf()) {
            if (node.getHeight() < m_nThreshold) {
                node.setHeight(0);
            }
        }
        if (!node.isLeaf()) {
            offset(node.getLeft(), fDelta);
            if (node.getRight() != null) {
                offset(node.getRight(), fDelta);
            }
        }
    }

}
