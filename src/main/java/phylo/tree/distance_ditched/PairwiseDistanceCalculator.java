package phylo.tree.distance_ditched;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.concurrent.Callable;

import htsjdk.samtools.util.Log;

/**
 * calculate 
 * @author tanxu
 *
 */
public class PairwiseDistanceCalculator implements Callable<Integer>{
	public static final Log log=Log.getInstance(PairwiseDistanceCalculator.class);
	
	
//	private final BatchPairwiseTreeDistCalculatorManager batchPairwiseTreeDistCalculatorManager;
	/**
	 * FileWriter to the file where all pairwise distance are to be gathered and stored calculated in the same batch run;
	 * provided by {@link BatchPairwiseTreeDistCalculatorManager}
	 */
	private final FileWriter batchRunFileWriter;
	
	/**
	 * 
	 */
	private final Path tmpOutFile;
	
	private final String tree1ID;
	private final String tree1NewickString;
	private final String tree2ID;
	private final String tree2NewickString;
	
	PairwiseDistanceCalculator(FileWriter fileWriter, Path tmpOutFile, String tree1ID, String tree1NewickString, String tree2ID, String tree2NewickString){
		
//		this.batchPairwiseTreeDistCalculatorManager=batchPairwiseTreeDistCalculatorManager;
		this.batchRunFileWriter = fileWriter;
		
		this.tmpOutFile=tmpOutFile;
		this.tree1ID=tree1ID;
		this.tree1NewickString=tree1NewickString;
		this.tree2ID=tree2ID;
		this.tree2NewickString=tree2NewickString;
	}
	
	/**
	 * invoke the r script to calculate the distance and output to the tmpOutFile
	 * read the calculated distance from the tmpOutFile and delete tmpOutFile
	 * write to the gathering file using fileWriter;
	 * 
	 * return 0 if succeed; 1 otherwise
	 */
	@Override
	public Integer call() throws Exception {
		int ret=0;
		//run r script and output the calculated distance to the tmpOutFile
		TreeDistanceRscriptRunner.run(tree1ID, tree1NewickString, tree2ID, tree2NewickString, tmpOutFile.toString());
		
		//
		if(!tmpOutFile.toFile().exists()) {
			log.error("tmpOutFile of calculated distance is not found!");
			ret=1;
		}else {
			//read the distance from tmpOutFile and write to the final output file with fileWriter
			try {
				BufferedReader lineReader = new BufferedReader(new FileReader(this.tmpOutFile.toFile()));
				
				double dist = Double.parseDouble(lineReader.readLine());
				
				lineReader.close();
				
				this.tmpOutFile.toFile().delete();
				
				
				/////////////////////
				String line = this.tree1ID.concat("\t").concat(this.tree2ID).concat("\t").concat(Double.toString(dist));
				
				synchronized(batchRunFileWriter) {
					BufferedWriter writer = new BufferedWriter(this.batchRunFileWriter);
		            // synchronizing the outputTreeFileWriter object
		            try {
						writer.append(line);
						writer.newLine();
						writer.flush();
					} catch (IOException e) {
						e.printStackTrace();
						ret=1;
					}
				}
				
			} catch (IOException ex) {
				ex.printStackTrace();
			    ret=1;
			} catch (NumberFormatException ne) {
				ne.printStackTrace();
				log.error("no distance is calculated between "+this.tree1ID+" and "+this.tree2ID);
				ret=1;
			}
		}
		
		//update 
//		synchronized(this.batchPairwiseTreeDistCalculatorManager) {
//			this.batchPairwiseTreeDistCalculatorManager.udpate(ret);
//		}
		
		return ret;
	}
}
