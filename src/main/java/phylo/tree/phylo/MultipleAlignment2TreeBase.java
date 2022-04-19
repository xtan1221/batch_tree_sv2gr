package phylo.tree.phylo;

import java.nio.file.Path;

import phylo.alignment.MultipleAlignment;
import phylo.exception.NewickTreeFileNotCreatedException;
import phylo.tree.reader.NewickFileFormatType;
import phylo.tree.reader.Tree;

/**
 * base class for building a Tree object from a MultipleAlignment
 * 
 * specifically,
 * 
 * 
 * 
 * @author tanxu
 *
 */
public abstract class MultipleAlignment2TreeBase {
//	/**
//	 * the name of the data; facilitate building output file names;
//	 */
//	private final String dataName;
//	/**
//	 * 
//	 */
//	private final MultipleAlignment msa;
//
//	
//	/**
//	 * directory for all produced temporary files for this tree only 
//	 */
//	private final Path tmpDir;
	////////////////////////////
	protected Process process; //the Process that runs the bash code
	protected Tree tree;
	
	/**
	 * 
	 * @param dataName
	 * @param msa
	 * @param tmpDir
	 */
	MultipleAlignment2TreeBase(){
	}
	
	
	/**
	 * 
	 * @param dataName
	 * @param msa
	 * @param tmpDir
	 * @throws NewickTreeFileNotCreatedException if the newick tree file is not successfully built (due to tree building programming errors or input msa data errors)
	 */
	public abstract void build(String dataName, MultipleAlignment msa, Path tmpDir) throws NewickTreeFileNotCreatedException;
	
	
	/**
	 * return the type of newick format of the generated tree file
	 * @return
	 */
	abstract NewickFileFormatType getNewickType();
	
	
	/////////////////////////////////////
//	/**
//	 * @return the msa
//	 */
//	public MultipleAlignment getMsa() {
//		return msa;
//	}
//
//	/**
//	 * @return the tmpDir
//	 */
//	public Path getTmpDir() {
//		return tmpDir;
//	}

	/**
	 * @return the tree
	 */
	public Tree getTree() {
		return tree;
	}
	
	public Process getProcess() {
		return this.process;
	}

	/**
	 * @param process the process to set
	 */
	public void setProcess(Process process) {
		this.process = process;
	}

//	/**
//	 * @return the dataName
//	 */
//	public String getDataName() {
//		return dataName;
//	}

}
