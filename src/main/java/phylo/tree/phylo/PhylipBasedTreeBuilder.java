package phylo.tree.phylo;

import java.nio.file.Path;

import phylo.alignment.MultipleAlignment;
import phylo.alignment.SequentialPhylipAlignmentFileWriter;
import phylo.exception.NewickTreeFileNotCreatedException;
import phylo.tree.phylo.PhylipBashRunner.DistanceModel;
import phylo.tree.reader.NewickFileFormatType;
import phylo.tree.reader.Tree;

public class PhylipBasedTreeBuilder extends MultipleAlignment2TreeBase{
	private final int outgroupIndex;
	private final DistanceModel distanceModel;
	
	/**
	 * constructor
	 * @param uniqueID
	 * @param msa
	 * @param tmpDir
	 * @param outgroupIndex
	 */
	PhylipBasedTreeBuilder(int outgroupIndex, DistanceModel distanceModel) {
		
		this.outgroupIndex=outgroupIndex;
		this.distanceModel=distanceModel;
	}
	
	/**
	 * do not specify the outgroup index, use the default one (1)
	 * @param dataName
	 * @param msa
	 * @param tmpDir
	 * @param distanceModel
	 */
	PhylipBasedTreeBuilder(DistanceModel distanceModel) {
		this(1, distanceModel);
	}
	
	
	/**
	 * major steps
	 * 
	 * 1. write the MultipleAlignment to a phylip file
	 * 2. run the phylip bash file to generate the newick tree file
	 * 3. read the newick tree file and create a Tree object
	 * 
	 * 4. post-processing
	 * 		deleting temporary files and folders
	 */
	@Override
	public void build(String dataName, MultipleAlignment msa, Path tmpDir) throws NewickTreeFileNotCreatedException{
		//1
		SequentialPhylipAlignmentFileWriter writer = 
				SequentialPhylipAlignmentFileWriter.makeWithOriginalSeqName(msa, tmpDir, dataName);
		
		writer.write();
		
		//2
//		String phylipFilePathString, String workingDirString, String outputDirString, DistanceModel distanceModel, int outgroupIndex
		Path workingDir = tmpDir;
		Path outputDir = tmpDir;
		
		PhylipBashRunner.runNJTreeBash(writer.getOutputFile().toString(), workingDir.toString(), outputDir.toString(), this.distanceModel, this.outgroupIndex);
		
		//${file_basename}.${distance_model}.dist.nj.tree.nwk
		//must be consistent with the bash file
		String newickTreeFileName = writer.getOutFileBaseName().concat(".").concat(this.distanceModel.getValue()).concat(".dist.nj.tree.nwk");
		Path outputNewickTreeFile = Path.of(outputDir.toString(), newickTreeFileName);
		
		if(!outputNewickTreeFile.toFile().exists()) {
			throw new NewickTreeFileNotCreatedException("tree newick file is not successfully created for:"+dataName);
		}
		
		//3
		this.tree = Tree.fromNewickFile(outputNewickTreeFile, this.getNewickType());
	}

	/**
	 * TODO
	 * need to verify for bootstrapped tree
	 */
	@Override
	NewickFileFormatType getNewickType() {
		return NewickFileFormatType.SIMPLE_NEWICK_1;
	}
	
}
