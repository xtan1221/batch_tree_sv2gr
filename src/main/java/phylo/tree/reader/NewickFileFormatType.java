package phylo.tree.reader;



/**
 * visframe supported tree file format that can be imported as VfTreeDataMetadata
 * @author tanxu
 * 
 */
public enum NewickFileFormatType{
	//internal node label is allowed
	SIMPLE_NEWICK_1("SIMPLE_NEWICK_1", "bootstrap value (if exist) is in squared bracket after branch length like in \n\t ((raccoon:19.19959,bear:6.80041):0.84600[50],((sea_lion:11.99700, seal:12.00300):7.52973[100],((monkey:100.85930,cat:47.14069):20.59201[80], weasel:18.87953):2.09460[75]):3.87382[50],dog:25.46154);"),  
	//this type of newick seems not allow internal node labels since the edge length takes the position of the internal node label of SIMPLE_NEWICK_1?
	//thus parsing of this type regarding the 
	SIMPLE_NEWICK_2("SIMPLE_NEWICK_2", "bootstrap value (if exist) is before branch length with a colon like in \n\t ((raccoon:19.19959,bear:6.80041)50:0.84600,((sea_lion:11.99700, seal:12.00300)100:7.52973,((monkey:100.85930,cat:47.14069)80:20.59201, weasel:18.87953)75:2.09460)50:3.87382,dog:25.46154);");
//	EXTENDED_NEWICK("Extended_NEWICK", ""),
//	RICH_NEWICK("Rich newick", ""),
//	NHX("New Hampshire X", ""),
//	NEXUS("NEXUS", ""), //contains a newick tree data string with other supportive features;
//	PhyloXML("PhyloXML", "");//do not contain any newick tree data string
	
	
	///////////////////////////////////////////////
	private final String fullName;
	private final String notes;
	
	
	/**
	 * constructor
	 * @param fullName
	 * @param notes
	 */
	NewickFileFormatType(String fullName, String notes){
		this.fullName = fullName;
		this.notes = notes;
	}
	
	
	public String getFullName() {
		return fullName;
	}


	/**
	 * @return the notes
	 */
	public String getNotes() {
		return notes;
	}

}
