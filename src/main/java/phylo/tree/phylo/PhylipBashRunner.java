package phylo.tree.phylo;

import java.nio.file.Path;

import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProcessExecutor;

/**
 * utility class
 * @author tanxu
 *
 */
public class PhylipBashRunner {
	public static final Log log = Log.getInstance(PhylipBashRunner.class);
	
	/**
	 * the path to the phylip_dist_nj_tree_command_line_parameter.sh file
	 */
	public static final Path NJ_TREE_BASH_FILE=
			Path.of("/home/tanxu/phylogeny/phylip/phylip_dist_nj_tree_command_line_parameter.sh");
	
	/**
	 * run the neighbor-joining tree bash that first calculate distance matrix with a phylip file using dnadist module then use it to build a nj tree with neighbor module in phylip
	 * for more detail, see the phylip_dist_nj_tree_command_line_parameter.sh
	 * 
	 * @param phylipFilePathString
	 * @param workingDirString
	 * @param outputDirString
	 * @param distanceModel
	 * @param outgroupIndex
	 * @return the exit value of the process
	 */
	public static int runNJTreeBash(String phylipFilePathString, String workingDirString, String outputDirString, DistanceModel distanceModel, int outgroupIndex) {
		String commandLineString=NJ_TREE_BASH_FILE.toString().concat(" ").concat(phylipFilePathString).concat(" ")
				.concat(workingDirString).concat(" ")
				.concat(outputDirString).concat(" ")
				.concat(distanceModel.getValue()).concat(" ")
				.concat(Integer.toString(outgroupIndex));
		
		log.info("run command to build tree: '".concat(commandLineString).concat("'"));
		
		//
		return ProcessExecutor.execute(commandLineString);
		
//		try {
//			Process process=Runtime.getRuntime().exec(commandLineString);
//			
//			BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
//		    String line = "";
//		    while ((line = reader.readLine()) != null) {
//		        System.out.println(line);
//		    }
//		    
//		    return process.waitFor(); //if directly invoke exitValue() method, 'java.lang.IllegalThreadStateException: process hasn't exited' will be thrown
////		    return process.exitValue();
//		    
//		} catch (IOException | InterruptedException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		
//		return 1;//not terminated normally
	}
	
	/**
	 * the path to the phylip_dist_nj_tree_command_line_parameter.sh file
	 */
	public static final Path DIST_MATRIX_TO_NJ_TREE_BASH=
			Path.of("/home/tanxu/phylogeny/phylip/phylip_dist_matrix_2_nj_tree_command_line_parameter.sh");
	
	/**
	 * 
	 * @param distMatrixFilePathString
	 * @param outputDirString
	 * @param outputFileBaseName
	 * @param outgroupIndex
	 * @return the file Path of the built newick tree file
	 */
	public static Path runDistMatrixToNJTreeBash(String distMatrixFilePathString, String outputDirString, String outputFileBaseName, int outgroupIndex) {
		String commandLineString=
				DIST_MATRIX_TO_NJ_TREE_BASH.toString().concat(" ")
				.concat(distMatrixFilePathString).concat(" ")
				.concat(outputDirString).concat(" ")
//				.concat(outputDirString).concat(" ")
				.concat(outputFileBaseName).concat(" ")
				.concat(Integer.toString(outgroupIndex));
		
		log.info("run command to build tree: '".concat(commandLineString).concat("'"));
		
		//
		ProcessExecutor.execute(commandLineString);
		
		return Path.of(outputDirString, outputFileBaseName.concat(".nj.nwk"));
	}
	
	
	/////////////////////////////////////////
	/**
	 * distance models used by distance matrix calculation by dnadist in phylip
	 * @author tanxu
	 *
	 */
	public static enum DistanceModel{
		F84("F84"), 
		Kimura("Kimura"), 
		Jukes_Cantor("Jukes-Cantor"), 
		LogDet("LogDet");
		/**
		 * 
		 */
		private final String value;
		DistanceModel(String value){
			this.value=value;
		}
		/**
		 * @return the value
		 */
		public String getValue() {
			return value;
		}
	}
	
	public static void main(String[] args) {
		String phylipFilePathString = "/scratch/tanxu/phylogenetic/phylip/propinquum.bicolor.first100000.phy";
		String workingDirString = "/scratch/tanxu/reseq/test/batch_local_tree/working";
		String outputDirString = "/scratch/tanxu/reseq/test/batch_local_tree/output";
		DistanceModel distanceModel = DistanceModel.Kimura;
		int outgroupIndex = 7;
		
		PhylipBashRunner.runNJTreeBash(phylipFilePathString, workingDirString, outputDirString, distanceModel, outgroupIndex);
	}
}
