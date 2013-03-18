package seqSimulator.evolution.tree;

import java.util.ArrayList;
import java.util.List;

import seqSimulator.evolution.datatype.MutableSequence;

public class Node {
	
	MutableSequence m_sequence;

    /**
     * label nr of node, used mostly when this is a leaf.
     * */
    protected int m_iLabel;
	
    /**
     * height of this node.
     */
    protected double m_fHeight = Double.MAX_VALUE;
	
    /**
     * list of children of this node *
     * Don't use m_left and m_right directly
     * Use getChildCount() and getChild(x) or getChildren() instead
     */
    List<Node> children = new ArrayList<Node>();
    
//    @Deprecated
//	private Node m_left;
//    @Deprecated
//	private Node m_right;

    /**
     * parent node in the beast.tree, null if root *
     */
    Node m_Parent = null;
    
    public Node() {}
    
    // identifiable
    protected String m_sID;

    public String getID() {
        return m_sID;
    }
    
    public void setID(String sID) {
        m_sID = sID;
    }
    
    /**
     * @return true if current node is a leaf node *
     */
    public boolean isLeaf() {
    	return children.size() == 0;
        //return getLeft() == null && getRight() == null;
    }
    
    public void setHeight(double fHeight) {
        m_fHeight = fHeight;
    }
    
    public double getHeight() {
        return m_fHeight;
    }
    
    public int getNr() {
        return m_iLabel;
    }

    public int getLeafNodeCount() {
        if (isLeaf()) {
            return 1;
        }
        int nodes = 0;
        for (Node child : children) {
       		nodes += child.getLeafNodeCount();
        }
        return nodes;
    }
    
    public void setNr(int iLabel) {
        m_iLabel = iLabel;
    }
    
    public void setParent(Node parent) {
        if (m_Parent != parent) {
            m_Parent = parent;
        }
    }

	public Node getLeft() {
		if (children.size() == 0) {
			return null;
		}
		return children.get(0);
	}

	public Node getRight() {
		if (children.size() <= 1) {
			return null;
		}
		return children.get(1);
	}
    
	public void setLeft(Node m_left) {
		if (children.size() == 0) {
	    	children.add(m_left);
		} else {
			children.set(0, m_left);
    	}
	}

	public void setRight(Node m_right) {
		switch (children.size()) {
		case 0:
	    	children.add(null);
		case 1:
	    	children.add(m_right);
	    	break;
		default:
			children.set(1, m_right);
	    	break;
    	}
	}
	
	// set sequence
	public void setSequence(MutableSequence inputSequence){
		m_sequence = inputSequence;
	}
	
	public MutableSequence getSequence(){
		return m_sequence.copy();
	}
	
    /**
     * @return length of branch in the beast.tree *
     */
    public double getLength() {
        if (isRoot()) {
            return 0;
        } else {
            return getParent().m_fHeight - m_fHeight;
        }
    }
	
    /**
     * @return parent node, or null if this is root *
     */
    public Node getParent() {
        return m_Parent;
    }
    
    /**
     * @return true if current node is root node *
     */
    public boolean isRoot() {
        return m_Parent == null;
    }
    
    /**
     * sorts nodes in children according to lowest numbered label in subtree
     *
     * @return
     */
    public int sort() {
        if (getLeft() != null) {
            int iChild1 = getLeft().sort();
            if (getRight() != null) {
                int iChild2 = getRight().sort();
                if (iChild1 > iChild2) {
                    Node tmp = getLeft();
                    setLeft(getRight());
                    setRight(tmp);
                    return iChild2;
                }
            }
            return iChild1;
        }
        // this is a leaf node, just return the label nr
        return m_iLabel;
    } // sort
	
    /**
     * during parsing, leaf nodes are numbered 0...m_nNrOfLabels-1
     * but internal nodes are left to zero. After labeling internal
     * nodes, m_iLabel uniquely identifies a node in a beast.tree.
     *
     * @param iLabel
     * @return
     */
    public int labelInternalNodes(int iLabel) {
        if (isLeaf()) {
            return iLabel;
        } else {
            iLabel = getLeft().labelInternalNodes(iLabel);
            if (getRight() != null) {
                iLabel = getRight().labelInternalNodes(iLabel);
            }
            m_iLabel = iLabel++;
        }
        return iLabel;
    } // labelInternalNodes
    
    // temp
    public void printStructure() {
    	if(!this.isLeaf()) {
    	Node left = this.getLeft();
    	Node right = this.getRight();
    	System.out.print("Node:" + this.getNr() + " Left:" + left.getNr() + " length:" + left.getLength() + " right:" + right.getNr() + " length:" + right.getLength());
    	if(left.isLeaf()) {
    		System.out.print(" left ID:" + left.getID());
    	}
    	if(right.isLeaf()) {
    		System.out.print(" right ID:" + right.getID());
    	}
    	System.out.println();
    	left.printStructure();
    	right.printStructure();
    	}
    }
    
}
