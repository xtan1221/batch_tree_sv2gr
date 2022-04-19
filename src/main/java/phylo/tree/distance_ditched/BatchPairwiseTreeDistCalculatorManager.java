package phylo.tree.distance_ditched;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.apache.commons.io.FileUtils;

import basic.Pair;
import htsjdk.samtools.util.Log;
import phylo.concurrency.TimeUtils;

public class BatchPairwiseTreeDistCalculatorManager {
	public static Log log=Log.getInstance(BatchPairwiseTreeDistCalculatorManager.class);
	
	//////////////////////////////////////
	/**
	 * empty space delimited file 
	 * each line contains the full set of information for a single tree;
	 * mandatory columns: 
	 * 		treeID column, which is the unique among all trees in the same file;
	 * 		full newick tree string column;
	 */
	private final Path treeDataFile;
	
	/**
	 * column index (first is 0) of the column containing the tree id string;
	 */
	private final int treeIDColumnIndex;
	/**
	 * column index of the column containing the full newick tree string;
	 */
	private final int newickTreeStringColumnIndex;
	
	/**
	 * the final output file containing all calculated pairwise tree distances for all trees in the {@link #treeDataFile}
	 */
	private final Path outputDistanceFile;
	
	/**
	 * directory to put all intermediate temporary files
	 */
	private final Path tmpDir;
	
	/**
	 * 
	 */
	private final int threadNum;
	
	
	/////////////////////////
	private List<Pair<String,String>> treeIDNewickStringPairList;
	/**
	 * FileWriter to the final output file for all calcualted pairwise tree distances;
	 */
	private FileWriter outputDistanceFileWriter;
	
	///////////////for monitoring all running jobs
//	private int totalJobNum=0;
//	private int finishedJobNum=0;
//	private Map<Integer,Integer> finishedJobReturnedStatusNumMap;
	
	
	
	BatchPairwiseTreeDistCalculatorManager(Path treeDataFile, int treeIDColumnIndex, int newickTreeStringColumnIndex, Path outputDistanceFile, Path tmpDir, int threadNum){
		this.treeDataFile=treeDataFile;
		this.treeIDColumnIndex = treeIDColumnIndex;
		this.newickTreeStringColumnIndex = newickTreeStringColumnIndex;
		this.outputDistanceFile=outputDistanceFile;
		this.tmpDir=tmpDir;
		this.threadNum=threadNum;
	}
	
//	/**
//	 * invoked by a job that is about to finish;
//	 * 
//	 * @param returnedStatus
//	 */
//	synchronized void udpate(int returnedStatus) {
//		this.finishedJobNum++;
//		if(!this.finishedJobReturnedStatusNumMap.containsKey(returnedStatus)) {
//			this.finishedJobReturnedStatusNumMap.put(returnedStatus, 0);
//		}
//		this.finishedJobReturnedStatusNumMap.put(returnedStatus, this.finishedJobReturnedStatusNumMap.get(returnedStatus)+1);
//		
//		///all jobs are finished
//		if(this.finishedJobNum==this.totalJobNum) {
//			log.info("all done!");
//			
//			//summarize the termination status of each submitted job
//			log.info("summarizing termination status of submitted jobs...");
//			List<Integer> terminationStatusList = new ArrayList<>();
//			terminationStatusList.addAll(this.finishedJobReturnedStatusNumMap.keySet());
//			//print to standard output
//			Collections.sort(terminationStatusList);
//			log.info("In total "+totalJobNum+" submitted jobs:");
//			for(int terminationStatus:terminationStatusList) {
//				log.info(finishedJobReturnedStatusNumMap.get(terminationStatus)+" have termination status = "+terminationStatus);
//			}
//			System.exit(0);
//		}
//	}
	
	/**
	 * step 1. read in all trees into a list
	 * step 2. for each pair of trees, calculate the distance
	 */
	void run() {
		/////
		this.treeIDNewickStringPairList = new ArrayList<>();
		
		log.info("start reading tree file "+this.treeDataFile.toString());
		try {
			BufferedReader reader = new BufferedReader(new FileReader(this.treeDataFile.toFile()));
			String line;
			while ((line = reader.readLine()) != null) {
				String[] splits = line.split("\\s+");
				String treeID=splits[this.treeIDColumnIndex];
				String treeNewickString=splits[this.newickTreeStringColumnIndex];
				
				if(treeID.isEmpty() || treeNewickString.isEmpty()) {
					log.error("empty tree ID or newick string is found! "+line);
					System.exit(1);
				}
				
				this.treeIDNewickStringPairList.add(new Pair<>(treeID, treeNewickString));
			}
			reader.close();
			
		}catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
		log.info("in total "+this.treeIDNewickStringPairList.size()+" trees are read in from tree file: "+this.treeDataFile.toFile());
		
		//initialize the outputDistanceFileWriter for all threads
		try {
			this.outputDistanceFileWriter = new FileWriter(this.outputDistanceFile.toFile(), true);
		} catch (IOException e1) {
			e1.printStackTrace();
			log.error("error!");
			System.exit(1);
		}
		
		//////////////////
		ExecutorService executorService = Executors.newFixedThreadPool(this.threadNum); 
		
		List<Future<Integer>> futureList = new ArrayList<>();
		
		///
//		this.finishedJobNum=0;
		int totalJobNum=0;
//		this.finishedJobReturnedStatusNumMap=new HashMap<>();
		////
		for(int i=0;i<this.treeIDNewickStringPairList.size()-1;i++) {
			String tree1ID = this.treeIDNewickStringPairList.get(i).getFirst();
			String tree1NewickTreeString = this.treeIDNewickStringPairList.get(i).getSecond();
			for(int j=i+1;j<this.treeIDNewickStringPairList.size();j++) {
				String tree2ID = this.treeIDNewickStringPairList.get(j).getFirst();
				String tree2NewickTreeString = this.treeIDNewickStringPairList.get(j).getSecond();
				
				String tmpOutFileName = tree1ID.concat("_").concat(tree2ID).concat(".dist");
				
				Path tmpOutFile = Path.of(this.tmpDir.toString(), tmpOutFileName);
				
				//
				PairwiseDistanceCalculator calculator = new PairwiseDistanceCalculator(this.outputDistanceFileWriter, tmpOutFile, tree1ID, tree1NewickTreeString, tree2ID, tree2NewickTreeString);
				//
				futureList.add(executorService.submit(calculator));
				
				totalJobNum++;
			}
		}
		
		log.info("total job num is: "+totalJobNum);
		
		//shut down ExecutorService so that program can terminate
		executorService.shutdown();
		////////////////////////////////////////
		//monitor all submitted jobs until all done (due to normal termination, an exception, or cancellation)
		boolean allDone=false;
		
		try {
			while(!allDone) {
				Thread.sleep(100000);
				log.info("still running...");
				allDone=true;
				for(Future<Integer> future:futureList) {
					allDone=allDone&&future.isDone();
					if(!allDone)
						break;
				}
			}
			log.info("all done!");
			
			//summarize the termination status of each submitted job
			log.info("summarizing termination status of submitted jobs...");
			List<Integer> terminationStatusList = new ArrayList<>();
			Map<Integer, Integer> terminationStatusJobNumberMap = new LinkedHashMap<>();
			for(Future<Integer> future:futureList) {
				int status=future.get();
				if(!terminationStatusJobNumberMap.containsKey(status)) {
					terminationStatusJobNumberMap.put(status, 0);
					terminationStatusList.add(status);
				}
				terminationStatusJobNumberMap.put(status, terminationStatusJobNumberMap.get(status)+1);
			}
			
			//print to standard output
			Collections.sort(terminationStatusList);
			log.info("In total "+futureList.size()+" submitted jobs:");
			for(int terminationStatus:terminationStatusList) {
				log.info(terminationStatusJobNumberMap.get(terminationStatus)+" have termination status = "+terminationStatus);
			}
		} catch (InterruptedException | ExecutionException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
	}
	
	
	public static void main(String[] args) {
		if(args.length!=6) {
			System.out.println("java treeDataFile treeIDColumnIndex newickTreeStringColumnIndex outputDistanceFile tmpDir threadNum");
			System.exit(1);
		}
		
		
		Path treeDataFile = Path.of(args[0]);
		int treeIDColumnIndex = Integer.parseInt(args[1]); 
		int newickTreeStringColumnIndex = Integer.parseInt(args[2]); 
		Path outputDistanceFile = Path.of(args[3]);
		Path tmpDir = Path.of(args[4]);
		int threadNum = Integer.parseInt(args[5]);
		
		//check treeDataFile file
		if(!treeDataFile.toFile().exists()) {
			log.error("treeDataFile is not found:"+treeDataFile.toString());
			System.exit(1);
		}
		
		//
		if(treeIDColumnIndex<0) {
			log.error("given treeIDColumnIndex: "+treeIDColumnIndex+" must be non-negative integer!");
			System.exit(1);
		}
		//
		if(newickTreeStringColumnIndex<0) {
			log.error("given newickTreeStringColumnIndex: "+newickTreeStringColumnIndex+" must be non-negative integer!");
			System.exit(1);
		}
		
		//check outputTreeFile
		if(outputDistanceFile.toFile().exists()) {
			log.error("existing file with the same name with outputDistanceFile is found!"+outputDistanceFile.toString());
			System.exit(1);
		}
		
		//check root tmp dir
		if(tmpDir.toFile().exists()) {
			if(tmpDir.toFile().isFile()) {
				log.error("existing file with the same path of the tmp directory is found!");
				System.exit(1);
			}else {
				try {
					FileUtils.cleanDirectory(tmpDir.toFile());
				} catch (IOException e) {
					e.printStackTrace();
					log.error("root tmp dir cannot be cleaned! "+tmpDir.toString());
					System.exit(1);
				}
			}
		}else {
			if(!tmpDir.toFile().mkdir()) {
				log.error("root tmp directory cannot be created!");
				System.exit(1);
			}
		}
		
		//check 
		if(threadNum<=0) {
			log.error("given threadNum: "+threadNum+" is not positive integer!");
			System.exit(1);
		}
		
		//////////////////////////////////
		
		BatchPairwiseTreeDistCalculatorManager manager = 
				new BatchPairwiseTreeDistCalculatorManager(
						treeDataFile,
						treeIDColumnIndex,
						newickTreeStringColumnIndex,
						outputDistanceFile,
						tmpDir,
						threadNum
						);
		
		///////////////////////
		long startTime=System.nanoTime();
		manager.run();
		log.info("elpased time:"+TimeUtils.getReadableTimeFromNanoTime(System.nanoTime() - startTime));
		//
		System.exit(0);
	}
}
