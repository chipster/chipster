package fi.csc.microarray.cluster;

import fi.csc.microarray.util.LookaheadStringReader;

/**
 * Parser for parenthesis formatted trees. Parses String and generates 
 * a corresponding tree of ClusterBranchNodes and ClusterLeafNodes.
 * 
 * @author Aleksi Kallio
 *
 */
public class ClusterParser {

	private static final String SEPARATE_PARTS = ":";
	private static final String CLOSE_TREE = ";";
	private static final String SEPARATE_BRANCHES = ",";
	private static final String CLOSE_BRANCH = ")";
	private static final String OPEN_BRANCH = "(";
	
	private LookaheadStringReader reader;

	public ClusterParser(String tree) {
		// remove whitespace first
		this.reader = new LookaheadStringReader(tree.replace("\n", "").replace(" ", ""));
	}

	public ClusterBranchNode getTree() throws TreeParseException {
		if (reader == null) {
			throw new IllegalStateException("Tree already parsed, cannot reparse");
		}
		
		// parse root branch 
		//   (it is different (no length), so it is parsed here 
		//   instead of using parseBranch)
		readAndCheck(OPEN_BRANCH);
		ClusterNode leftBranch = parseBranchOrLeaf();
		readAndCheck(SEPARATE_BRANCHES);
		ClusterNode rightBranch = parseBranchOrLeaf();
		readAndCheck(CLOSE_BRANCH);
		readAndCheck(CLOSE_TREE);
		
		// disable parsing
		reader = null;
		
		// parse AST
		ClusterBranchNode rootBranch = new ClusterBranchNode();
		rootBranch.setLeftBranch(leftBranch);
		rootBranch.setRightBranch(rightBranch);
		return rootBranch;
	}

	private ClusterBranchNode parseBranch() throws TreeParseException {

		// parse
		readAndCheck(OPEN_BRANCH);
		ClusterNode leftBranch = parseBranchOrLeaf();
		readAndCheck(SEPARATE_BRANCHES);
		ClusterNode rightBranch = parseBranchOrLeaf();
		readAndCheck(CLOSE_BRANCH);
		String length = parseLength();
		
		// construct AST
		ClusterBranchNode branch = new ClusterBranchNode();
		branch.setLeftBranch(leftBranch);
		branch.setRightBranch(rightBranch);
		branch.setLength(length);
		return branch;
	}

	private ClusterNode parseBranchOrLeaf() throws TreeParseException {
		ClusterNode node;
		if (reader.lookahead().equals(OPEN_BRANCH)) {
			node = parseBranch();
		} else {
			node = parseLeaf();
		}
		return node;
	}
	
	private ClusterLeafNode parseLeaf() throws TreeParseException {
		
		// parse
		String gene = parseGene();
		String length = parseLength();
		
		// construct AST
		ClusterLeafNode leaf = new ClusterLeafNode();
		leaf.setGene(gene);
		leaf.setLength(length);
		return leaf;
	}

	private String parseGene() {
		String gene = reader.readTo(SEPARATE_PARTS);
		return gene;
	}

	private String parseLength() throws TreeParseException {
		readAndCheck(SEPARATE_PARTS);
		String length = reader.readTo(SEPARATE_BRANCHES, CLOSE_BRANCH);
		return length;
		
	}

	private void readAndCheck(String s) throws TreeParseException {
		try {
			String r = reader.read();
			if (!r.equals(s)) {
				throw new TreeParseException("expected \"" + s + "\" but found \"" + r + "\"", reader);
			}
		} catch (StringIndexOutOfBoundsException se) {
			throw new TreeParseException("unexpected end of string", reader);			
		}
	}

}
