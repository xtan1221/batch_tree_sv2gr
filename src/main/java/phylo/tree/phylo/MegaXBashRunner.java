package phylo.tree.phylo;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Path;

import htsjdk.samtools.util.Log;

public class MegaXBashRunner {
	public static final Log log = Log.getInstance(MegaXBashRunner.class);
	
	/**
	 * the path to the phylip_dist_nj_tree_command_line_parameter.sh file
	 */
	public static final Path MEGAX_BASH_FILE= //TODO
			Path.of("/home/tanxu/phylogeny/megaX/megax_tree_with_bootstrapping_command_line_parameter.sh");
//			"sv2gr/tree/phylip_dist_nj_tree_command_line_parameter.sh";
	
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
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void runBash(MegaXBasedTreeBuilder builder, String maoFilePathString, String multipleAlignmentFastaFile, String outputDirString, String dataName) throws IOException, InterruptedException {
		String commandLineString=
				MEGAX_BASH_FILE.toString().concat(" ")
				.concat(maoFilePathString).concat(" ")
				.concat(multipleAlignmentFastaFile).concat(" ")
				.concat(outputDirString).concat(" ");
		
		log.info("run command to build tree with mega10: '".concat(commandLineString).concat("'"));
		
		//
//		return ProcessExecutor.execute(commandLineString);
		
//		try {
		Process process=Runtime.getRuntime().exec(commandLineString);
		builder.setProcess(process);
		BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
	    String line;
	    while ((line = reader.readLine()) != null) {
	    	log.info("tree id=="+dataName+"  Mega10 stardand out >>>>>> "+line);
	    }
	    process.waitFor(); //wait until the process is terminated; if somehow it does not finish, will stuck here!!!!!
	    log.info("Process: ".concat(process.toString()));
//		    return process; //if directly invoke exitValue() method, 'java.lang.IllegalThreadStateException: process hasn't exited' will be thrown
//		    return process.exitValue();
		    
//		} catch (IOException | InterruptedException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
		
//		return 1;//not terminated normally
	}
}
