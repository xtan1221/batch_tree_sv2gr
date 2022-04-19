package phylo.tree.phylo;

import java.nio.file.Path;

import htsjdk.samtools.util.Log;
import phylo.alignment.FastaAlignmentFileWriter;
import phylo.alignment.MultipleAlignment;
import phylo.exception.NewickTreeFileNotCreatedException;
import phylo.tree.reader.NewickFileFormatType;
import phylo.tree.reader.Tree;

public class MegaXBasedTreeBuilder extends MultipleAlignment2TreeBase{
	public static final Log log=Log.getInstance(MegaXBasedTreeBuilder.class);
	
	/**
	 * the file containing all settings/configurations of running MegaX from command line;
	 */
	private final Path maoFile;
	
	/////////////////////////////////
	/**
	 * constructor
	 * @param uniqueID
	 * @param msa
	 * @param tmpDir
	 * @param outgroupIndex
	 */
	public MegaXBasedTreeBuilder(Path maoFile) {
		if(maoFile==null||!maoFile.toFile().exists())
			throw new IllegalArgumentException("given mao file cannot be null or does not exist!");
		
		this.maoFile=maoFile;
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
		log.info("start build tree with mega10 for tree id=="+dataName);
		//1
		FastaAlignmentFileWriter writer = 
				new FastaAlignmentFileWriter(msa, tmpDir, dataName);
		
		writer.write();
		
		//2
		try {
			MegaXBashRunner.runBash(this, this.maoFile.toString(), writer.getOutputFile().toString(), tmpDir.toString(), dataName);
//			this.process.waitFor();
		} catch (Exception e) {
//			e.printStackTrace();
			log.error("Exception from MegaXBashRunner "+e.getMessage()+" for tree id=="+dataName);
		}
		
		log.info("mega10 bash run is done for tree id=="+dataName);
		
		//${output_file_basename}.nwk
		//must be consistent with the bash file
		String newickTreeFileName = writer.getOutFileBaseName().concat(".nwk");
		Path outputNewickTreeFile = Path.of(tmpDir.toString(), newickTreeFileName);
		
		if(!outputNewickTreeFile.toFile().exists()) {
			throw new NewickTreeFileNotCreatedException("tree newick file is not successfully created for tree id=="+dataName);
		}
		
		log.info("tree newick file is successfully created for tree id=="+dataName);
		
		//3
		this.tree = Tree.fromNewickFile(outputNewickTreeFile, this.getNewickType());
		
		log.info("tree object is successfully created for tree id=="+dataName);
	}
	
	@Override
	NewickFileFormatType getNewickType() {
		return NewickFileFormatType.SIMPLE_NEWICK_2;
	}
	
}
