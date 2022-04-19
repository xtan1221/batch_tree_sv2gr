/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package phylo.tree.reader;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import basic.Pair;
import basic.Triple;

/**
 * key word: 
 * Branch label: including edge length or bootstrap value or none or both
 * Node label: normally the name of the node, can be empty
 * 
 * for simple newick tree, only leaf node can have label
 * @author tanxu
 */
public class SimpleNewickParserUtils {
    /**
     * delimit first level nodes in the same branch
     */
    private final static String COMMA = ",";
    /**
     * tag for start of an internal node/branch
     */
    private final static String LEFT_PARENTHESIS = "(";
    /**
     * tag for end of an internal node/branch
     */
    private final static String RIGHT_PARENTHESIS = ")";
    /**
     * end of newick tree string
     */
    private final static String SEMICOLON=";";
    /**
     * branch length tag
     */
    private final static String COLON=":";
    
    /**
     * 
     */
    private final static String NON_RESERVED_CHARACTERS = "[^\\,\\(\\)\\;\\:]";
    /**
     * tag for start of bootstrap value
     */
    private final static String LEFT_SQUARE_BRACKET="[";
    /**
     * tag for end of bootstrap value
     */
    private final static String RIGHT_SQUARE_BRACKET="]";
    
    private final static String WHITE_SPACES="\\s+";
    
    private final static String EMPTY="";
    
    private final static String UNDERSCORE="_";
    
    
    
    /**
     * non-empty node label string and branch label string;
     * 
     * leaf_node_label:leaf_branch_label
     * 
     * group 1 = leaf_node_label; group 2 = leaf_branch_label
     */
    private final static Pattern FULL_LEAF = Pattern.compile(
                "([^\\,\\(\\)\\;\\:]+)" //group 1
                +"\\:"
                +"([^\\,\\(\\)\\;\\:]+)" //group 2
    );
    
    /**
     * non-empty node label string and empty branch label string
     * leaf_node_label
     * 
     * group 0 = leaf_node_label
     */
    private final static Pattern LEAF_NODE_LABEL = Pattern.compile("[^\\,\\(\\)\\;\\:]+");
    
    /**
     * empty node label string and non-empty branch label string
     * :leaf_branch_label
     * group 1 = leaf_branch_label
     */
    private final static Pattern LEAF_BRANCH_LABEL = Pattern.compile("\\:([^\\,\\(\\)\\;\\:]+)");
    
    //////////////////////////
    /**
     * (raccoon:19.19959,bear:6.80041)mammal:0.84600[50]
     * (cn1,cn2,...)node_label:branch_label
     * note that branch_label could contains branch length and/or bootstrap
     */
    private final static Pattern SIMPLE_NEWICK_1_FULL_INTERNAL_NODE = Pattern.compile(
                "\\((.+)\\)"   //children newick strings with outer parenthsis removed = group 1
                + "([^\\,\\(\\)\\;\\:]+)"  //node_label = group 2
                + "\\:"  //colon
                + "([^\\,\\(\\)\\;\\:]+)" //branch_label = group 3
    );
    
    
    /**
     * (cn1,cn2,...)node_label
     */
    private final static Pattern SIMPLE_NEWICK_1_INTERNAL_NODE_LABEL_ONLY = Pattern.compile(
                "\\((.+)\\)"   //children newick strings with outer parenthsis removed = group 1
                + "([^\\,\\(\\)\\;\\:]+)"  //node_label = group 2
    );
    
    /**
     * (cn1,cn2,...):branch_label
     */
    private final static Pattern SIMPLE_NEWICK_1_INTERNAL_BRANCH_LABEL_ONLY = Pattern.compile(
                "\\((.+)\\)"   //children newick strings with outer parenthesis removed = group 1
                + "\\:"  //colon
                + "([^\\,\\(\\)\\;\\:]+)" //branch_label = group 2
    );
    
    /////////////////////////////////////////////////////
    ///for simple newick format 2;
    ///note that simple newick format 2 cannot have internal node label
    /**
     * (raccoon:19.19959,bear:6.80041)90:0.84600
     * (cn1,cn2,...)branch_label
     */
    private final static Pattern SIMPLE_NEWICK_2_INTERNAL_BRANCH_LABEL_ONLY = Pattern.compile(
                "\\((.+)\\)"   //children newick strings with outer parenthsis removed = group 1
                + "([^\\,\\(\\)\\;]+)" //branch_label = group 2
    );
    
    
    /**
     * not containing any node label nor branch label;
     * (cn1,cn2,...)
     */
    private final static Pattern EMPTY_INTERNAL_NODE = Pattern.compile(
                "\\((.+)\\)"   //children newick strings with outer parenthsis removed = group 1
    );
    
    
    /**
     * check if the given newick tree string is in format '(...)'
     * @param newickTreeString
     * @return
     */
    public static boolean isEnclosedByAPairOfParenthesis(String newickTreeString) {
    	String[] splits = newickTreeString.split("");
    	
    	int layer=0;
    	
    	for(int i=0;i<splits.length;i++) {
    		if(splits[i].equals("(")) {
    			layer++;
    		}else if(splits[i].equals(")")) {
    			layer--;
    		}
    		
    		if(layer==0) {//layer can only reach 0 again when the last character is reached;
    			if(i==splits.length-1) {
    				return true;
    			}else {
    				return false;
    			}
    		}
    	}
    	
    	if(layer!=0)
    		throw new IllegalArgumentException("given newick tree string is invalid!");
    	
    	return true;
    }
    
    /**
     * pre-process the newick string read from a newick format file;
     * 1. remove ending semicolon;
     * 2. if implicit root, add a pair of parenthesis to cover the full newick string
     * 
     * @param rawString
     * @return 
     */
    public static String preprocessNewickStringFromFile(String rawString){
        
        if(rawString.endsWith(SEMICOLON)){
            rawString = rawString.substring(0,rawString.length()-1);
        }
        
        //
        if(!isEnclosedByAPairOfParenthesis(rawString)){
            rawString = LEFT_PARENTHESIS.concat(rawString).concat(RIGHT_PARENTHESIS);
        }
        
        //add any processing if necessary
        
        
        return rawString;
        
    }
    
    /**
     * 
     * @param newickString
     * @param formatType
     * @return
     */
    public static Triple<String, String, String> extractChildrenNodeStringNodeLabelStringAndBranchLabelString(String newickString, NewickFileFormatType formatType){
    	if(formatType.equals(NewickFileFormatType.SIMPLE_NEWICK_1)) {
			return extract_SIMPLE_NEWICK_1_ChildrenNodeStringNodeLabelStringAndBranchLabelString(newickString);
		}else if(formatType.equals(NewickFileFormatType.SIMPLE_NEWICK_2)) {
			return extract_SIMPLE_NEWICK_2_ChildrenNodeStringNodeLabelStringAndBranchLabelString(newickString);
		}else {
			throw new UnsupportedOperationException("given VfTreeDataFileFormatType is not supported yet!");
		}
    }
    
    /**
     * {@link NewickFileFormatType#SIMPLE_NEWICK_1} format specific
     * 
     * 
     * 
     * extract a newick string of a node N into children node string, node label string of N, branch label string between N and its parent node;
     * the input string must be a VALID newick string except that the node string is an empty leaf without node and branch label strings (trivial case); 
     * 
     * if no children nodes exist(leaf), return a null value;
     * if no node or branch label string exist, return an empty string;
     * 
     * @param newickString
     * @return 
     */
    protected static Triple<String, String, String> extract_SIMPLE_NEWICK_1_ChildrenNodeStringNodeLabelStringAndBranchLabelString(String newickString){
//        System.out.println(newickString);
        String childrenNodesString = null;
        String nodeLabelString = "";
        String branchLabelString = "";
        
        Matcher m;
        if(newickString.startsWith(LEFT_PARENTHESIS)){//internal node
            if((m = SIMPLE_NEWICK_1_FULL_INTERNAL_NODE.matcher(newickString)).matches()){
                childrenNodesString = m.group(1);
                nodeLabelString = m.group(2);
                branchLabelString = m.group(3);
            }else if((m = SIMPLE_NEWICK_1_INTERNAL_NODE_LABEL_ONLY.matcher(newickString)).matches()){
                childrenNodesString = m.group(1);
                nodeLabelString = m.group(2);
            }else if((m = SIMPLE_NEWICK_1_INTERNAL_BRANCH_LABEL_ONLY.matcher(newickString)).matches()){
                childrenNodesString = m.group(1);
                branchLabelString = m.group(2);
            }else if((m = EMPTY_INTERNAL_NODE.matcher(newickString)).matches()){// the internal node has no node label string nor branch label string
                childrenNodesString = m.group(1);
            }else{
                throw new IllegalArgumentException("unrecognized internal node string:"+newickString);
            }
            
        }else{//leaf
            if((m = FULL_LEAF.matcher(newickString)).matches()){
                nodeLabelString = m.group(1);
                branchLabelString = m.group(2);
            }else if((m = LEAF_NODE_LABEL.matcher(newickString)).matches()){
                nodeLabelString = m.group(0);
            }else if((m = LEAF_BRANCH_LABEL.matcher(newickString)).matches()){
                branchLabelString = m.group(1);
            }else if(newickString.isEmpty()){//empty leaf node == input string is an empty string
                //do nothing
            }else{
                throw new IllegalArgumentException("unrecognized internal leaf string:"+newickString);
            }
        }
        
        
        return new Triple<>(childrenNodesString, nodeLabelString, branchLabelString);
    }
    
    /**
     * {@link NewickFileFormatType#SIMPLE_NEWICK_1} format specific
     * 
     * 
     * 
     * extract a newick string of a node N into children node string, node label string of N, branch label string between N and its parent node;
     * the input string must be a VALID newick string except that the node string is an empty leaf without node and branch label strings (trivial case); 
     * 
     * if no children nodes exist(leaf), return a null value;
     * if no node or branch label string exist, return an empty string;
     * 
     * @param newickString
     * @return 
     */
    protected static Triple<String, String, String> extract_SIMPLE_NEWICK_2_ChildrenNodeStringNodeLabelStringAndBranchLabelString(String newickString){
        String childrenNodesString = null;
        String nodeLabelString = "";
        String branchLabelString = "";
        
        Matcher m;
        if(newickString.startsWith(LEFT_PARENTHESIS)){//internal node
            if((m = SIMPLE_NEWICK_2_INTERNAL_BRANCH_LABEL_ONLY.matcher(newickString)).matches()){
                childrenNodesString = m.group(1);
                branchLabelString = m.group(2);
            }else if((m = EMPTY_INTERNAL_NODE.matcher(newickString)).matches()){// the internal node has no node label string nor branch label string
                childrenNodesString = m.group(1);
            }else{
                throw new IllegalArgumentException("unrecognized internal node string:"+newickString);
            }
            
        }else{//leaf
            if((m = FULL_LEAF.matcher(newickString)).matches()){
                nodeLabelString = m.group(1);
                branchLabelString = m.group(2);
            }else if((m = LEAF_NODE_LABEL.matcher(newickString)).matches()){
                nodeLabelString = m.group(0);
            }else if((m = LEAF_BRANCH_LABEL.matcher(newickString)).matches()){
                branchLabelString = m.group(1);
            }else if(newickString.isEmpty()){//empty leaf node == input string is an empty string
                //do nothing
            }else{
                throw new IllegalArgumentException("unrecognized internal leaf string:"+newickString);
            }
            if(!branchLabelString.isEmpty()) {
            	branchLabelString=":".concat(branchLabelString);
            }
        }
        
        
        return new Triple<>(childrenNodesString, nodeLabelString, branchLabelString);
    }
    
    
    /**
     * parse the given internal node string containing a list of children nodes into a list of children node string;
     * the input node's labels and the outer parenthesis must be removed before feed to this method!
     * (c1,c2,c3,....)node_label:branch_label   =====>  c1,c2,c3,... (this is the correct input string);
     * for example: (raccoon:19.19959,bear:6.80041):0.84600[50] is not acceptable as input string, 
     * need to be raccoon:19.19959,bear:6.80041 with the node/branch labels and outer parenthesis of the node parsed and removed;
     * @param internalNodeString
     * @return 
     */
    public static List<String> splitNakedInternalNodeStringIntoChildrenNodeStrings(String internalNodeString){
        List<String> ret = new ArrayList<>();
        String[] splits = internalNodeString.split(EMPTY); //split newick string into single characters
        
        int layer=0;
        int previousUnParsedPos = 0; //index of the first unparsed char
        
        for(int i=0;i<splits.length;i++){
            
            if(splits[i].equals(LEFT_PARENTHESIS)){
                layer++;
            }
            
            if(splits[i].equals(RIGHT_PARENTHESIS)){
                layer--;
            }
//            System.out.print(layer);
            
            
            if(splits[i].equals(COMMA)){
                if(layer == 0){ //first layer comma indicate end of a sub-branch
                    
                    String branchString = internalNodeString.substring(previousUnParsedPos, i);
                    ret.add(branchString);
                    
                    previousUnParsedPos = i+1;
                }
            }   
            
            if(i == splits.length-1){//end of string
                if(layer != 0){
                    throw new IllegalArgumentException("Given newick string is invalid with unequal number of left and right parenthesis:"+internalNodeString);
                }


                String branchString = internalNodeString.substring(previousUnParsedPos);

                ret.add(branchString);


            }   
            
        }
        
        
        return ret;
    }
    

    /**
     * 
     * @param edgeLabelString
     * @param formatType
     * @return
     */
    public static Pair<Double, Double> parseEdgeLabelStringForLengthAndBootstrap(String edgeLabelString,
			NewickFileFormatType formatType) {
		if(formatType.equals(NewickFileFormatType.SIMPLE_NEWICK_1)) {
			return parse_NEWICK_1_EdgeLabelStringForLengthAndBootstrap(edgeLabelString);
		}else if(formatType.equals(NewickFileFormatType.SIMPLE_NEWICK_2)) {
			return parse_NEWICK_2_EdgeLabelStringForLengthAndBootstrap(edgeLabelString);
		}else {
			throw new UnsupportedOperationException("given VfTreeDataFileFormatType is not supported yet!");
		}
	}
	
    /**
     * parse the edge label string for the possible length and bootstrap with respect to {@link NewickFileFormatType#SIMPLE_NEWICK_1} format:
     * length[bootstrap]  ===== 0.05818[100]
     * length ==== 0.15420
     * [bootstrap] ===== [92]
     * empty string ==== no length nor bootstrap
     * 
     * 
     * 
     * @param branchLabelString
     * @return 
     */
    protected static Pair<Double,Double> parse_NEWICK_1_EdgeLabelStringForLengthAndBootstrap(String branchLabelString){
        
        Double length = null;
        Double bootstrap = null;
        
        if(branchLabelString.trim().isEmpty()){
        	//both length and bootstrap are empty
        }else{
            if(branchLabelString.contains("[")&&branchLabelString.contains("]")){//has bootstrap;
                Pattern p = Pattern.compile("(.*)\\[(.+)\\]");
                Matcher m = p.matcher(branchLabelString);
                if(m.matches()){
                    if(m.group(1).trim().isEmpty()){
                        
                    }else{
                        length = Double.parseDouble(m.group(1));
                    }
                    
                    bootstrap = Double.parseDouble(m.group(2).trim());
                }else{
                    throw new IllegalArgumentException("input branchLabelString is invalid");
                }
                
            }else{//no bootstrap;
                length = Double.parseDouble(branchLabelString.trim());
            }
            
            
        }
        
        //if length is not null and negative, set it to 0;
        length = length==null?null:length<0?0:length;
        
        Pair<Double,Double> ret = new Pair<>(length, bootstrap);
        
        return ret;
    }
    
	
    /**
     * parse the edge label string for the possible length and bootstrap with respect to {@link NewickFileFormatType#SIMPLE_NEWICK_2} format:
     * bootstrap:length  ===== 50:0.84600
     * :length ==== :25.46154
     * bootstrap ===== 92
     * empty string ==== no length nor bootstrap
     * @param branchLabelString
     * @return 
     */
    protected static Pair<Double,Double> parse_NEWICK_2_EdgeLabelStringForLengthAndBootstrap(String branchLabelString){
        
        Double length = null;
        Double bootstrap = null;
        
        if(branchLabelString.trim().isEmpty()){
        	//both length and bootstrap are empty
        }else{
            if(branchLabelString.contains(":")){//has length;
            	String[] splits = branchLabelString.split(":");
//                Pattern p = Pattern.compile("(.*)\\:(.+)");
//                Matcher m = p.matcher(branchLabelString);
//                if(m.matches()){
                    if(splits[0].isEmpty()){//bootstrap is empty
                        
                    }else{
                    	bootstrap = Double.parseDouble(splits[0]);
                    }
                    
                    length = Double.parseDouble(splits[1].trim());
//                }else{
//                    throw new IllegalArgumentException("input branchLabelString is invalid");
//                }
                
            }else{//no length;
            	bootstrap = Double.parseDouble(branchLabelString.trim());
            }
            
        }
        
        //if length is not null and negative, set it to 0;
        length = length==null?null:length<0?0:length;
        
        Pair<Double,Double> ret = new Pair<>(length, bootstrap);
        
        
        return ret;
    }
    
    
    public static void main(String[] args){
//        String s = "a:1sldkjf";
//        String ss = "(raccoon:19.19959,bear:6.80041)mammal:0.84600";
//        Matcher m = FULL_INTERNAL_NODE.matcher(ss);
//        if(m.matches()){
//            for(int i = 0;i<m.groupCount()+1;i++){
//                System.out.println(i+ ":" + m.group(i));
//            }
//        }
//        
        
        
//        String nakedInternalNodeString1 = "(raccoon, bear),((sea_lion,seal),((monkey,cat), weasel)),dog";
//        String nakedInternalNodeString2 = "(raccoon:19.19959,bear:6.80041):0.84600[50],((sea_lion:11.99700, seal:12.00300):7.52973[100],((monkey:100.85930,cat:47.14069):20.59201[80], weasel:18.87953):2.09460[75]):3.87382[50],dog:25.46154";
//        String nakedInternalNodeString3 = "(sea_lion:11.99700, seal:12.00300):7.52973[100],((monkey:100.85930,cat:47.14069):20.59201[80], weasel:18.87953):2.09460[75]";
//        List<String> clist = splitNakedInternalNodeStringIntoChildrenNodeStrings(nakedInternalNodeString3);
//        for(String c:clist){
//           System.out.println(c);
//        }
        
        
        
//        String tree1 = "((raccoon, bear),((sea_lion,seal),((monkey,cat), weasel)),dog);";
//        
//        Triple<String,String,String> splits = extractChildrenNodeStringNodeLabelStringAndBranchLabelString(tree1);
//        System.out.println("left:"+splits.getLeft());
//        System.out.println("middle:"+splits.getMiddle());
//        System.out.println("right:"+splits.getRight());
        


        String branchLabelString = " [232]";
        Pair<Double,Double> ret = parse_NEWICK_1_EdgeLabelStringForLengthAndBootstrap(branchLabelString);
        
        System.out.println("length=="+ret.getFirst());
        System.out.println("bootstrap=="+ret.getSecond());
        
    }
    
    
    
    
}
