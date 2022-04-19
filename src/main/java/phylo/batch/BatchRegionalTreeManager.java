package phylo.batch;

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

import htsjdk.samtools.util.Log;
import htsjdk.variant.vcf.VCFFileReader;
import phylo.concurrency.TimeUtils;
import phylo.tree.phylo.MegaXBasedTreeBuilderFactory;
import phylo.tree.phylo.MultipleAlignment2TreeFactory;
import phylo.vcf.Vcf2AlignmentSeqBuilderFactory;
import population.vcf.utils.GenotypeRecoderFactory;
import population.vcf.utils.VariantContextFilterFactory;

/**
 * when running in linux environment 
 * 1. java version 1.8 or higher
 * 2. dependencies:
 * 		bcftools, phylip/mega10
 * 
 * ==================================
 * major input files
 * 1. bcf file containing SNPs (and maybe indels)
 * 2. region file each line of which contains a list of regions to build a tree based on the regional SNP extracted from the bcf file;
 * 		see -r option of bcftools
 * 
 * major output file
 * 		a file each line containing the data of a tree
 * 
 * ==================================
 * this class does the following
 * for each line/record/region in the bed file
 * 	1. extract the SNPs based on provided VariantContextFilter (filter for SNP site)
 * 	2. build a sequence alignment for all individuals based on given GenotypeRecoder (recode genotype to bases) and output to phylip format
 * 	3. build a tree based with the phylip file based on given tree inference strategy;
 * 
 * simple concurrency is supported to run multiple regions in parallel in a predefined thread pool;
 * 
 * @author tanxu
 * 
 */
public class BatchRegionalTreeManager {
	public static Log log=Log.getInstance(BatchRegionalTreeManager.class);
	
	/////////////////////
	/**
	 * 
	 */
	private final Bcf2SingleTreeFactory bcf2SingleTree2Factory;
	
	/**
	 * each line contains two columns delimited by empty space
	 * 
	 * 		column 1 is the unique id of the list of regions
	 * 		column 2 is a region list consistent with the -r option of bcftools
	 *
	 * SNP within the regions are to be used to build a tree
	 */
	private final Path regionListFile;
	
	/**
	 * the file to store all generated trees;
	 * specifically, each line contains the information for one tree in three columns delimited by empty spaces;
	 * 		column 1 is the unique id of the list of regions
	 * 		column 2 is a region list consistent with the -r option of bcftools
	 * 		column 3 is the full newick tree string;
	 */
	private final Path outputTreeFile;
	
	/**
	 * the root tmp directory where to put all temp folders and files of all trees;
	 * 
	 * note that for each tree, a separate folder will be created under this directory;
	 */
	private final Path rootTmpDir;
	/**
	 * 
	 */
	private final int threadNum;
	
	
	/////////////////////
	/**
	 * FileWriter to the final output file for all constructed trees;
	 */
	private FileWriter outputTreeFileWriter;
	/**
	 * constructor
	 * @param bcfFile
	 * @param outgroupIndex
	 * @param bedFile
	 * @param workdingDir
	 * @param outputDir
	 * @param threadNum
	 */
	BatchRegionalTreeManager(Bcf2SingleTreeFactory bcf2SingleTree2Factory, Path regionListFile, Path outputTreeFile, Path rootTmpDir, int threadNum){
		//TODO validate
		
		
		/////////////////////////////////
		this.bcf2SingleTree2Factory = bcf2SingleTree2Factory;
		this.regionListFile = regionListFile;
		this.outputTreeFile = outputTreeFile;
		this.rootTmpDir=rootTmpDir;
		this.threadNum = threadNum;
	}
	
	
	/**
	 * 
	 */
	void run() {
		try {
			this.outputTreeFileWriter = new FileWriter(this.outputTreeFile.toFile(), true);
		} catch (IOException e1) {
			e1.printStackTrace();
			log.error("error!");
			return;
		}
		
		//////////////////
		ExecutorService executorService = Executors.newFixedThreadPool(this.threadNum); 
		
//		List<Future<Integer>> futureList = new ArrayList<>();
		Map<Future<Integer>, Bcf2SingleTree> allFutureJobMap = new LinkedHashMap<>();
		try {
			BufferedReader reader = new BufferedReader(new FileReader(this.regionListFile.toFile()));
			String line;
			while ((line = reader.readLine()) != null) {
				String[] splits = line.split("\\s+");
				String regionID=splits[0];
				String regionString=splits[1];
				
				//specify a temporary directory under the rootTmpDir for the tree
				Path treeTempDir = Path.of(this.rootTmpDir.toString(), regionID);
				
				assert !treeTempDir.toFile().exists();
				
				log.info("create temp dir for region:"+regionID);
				if(!treeTempDir.toFile().mkdir()) {
					log.error("temp dir for region:"+regionID+" cannot be created!");
					continue;
				}
				
				//
				Bcf2SingleTree regionalTree =
						this.bcf2SingleTree2Factory.makeNew(regionString, regionID, this.outputTreeFileWriter, treeTempDir);
				
				allFutureJobMap.put(executorService.submit(regionalTree), regionalTree);
			}
			
			reader.close();
			
			////////////////////////////////////////
			//monitor all submitted jobs until all done (due to normal termination, an exception, or cancellation)
			boolean allDone=false;
			
			long maxFullRunTimeForSuccessfullyFinishedJob=0;
			double flexibility=3;
			double ratio = 0.5; //when the finished job num/total job num > this value, start to restart unfinished jobs with running time > maxFullRunTimeForSuccessfullyFinishedJob * flexibility
			double minSuccessfullyDoneJobNum=100; //
			//////////////////////////////
			while(!allDone) {
				Thread.sleep(60000);
				log.info("update job status ...");
				allDone=true;
				
				Map<Future<Integer>, Bcf2SingleTree> runningJobFutureMap = new LinkedHashMap<>();
				
				int totalSuccessfullyDoneJobNum=0;
				//check status of all jobs
				for(Future<Integer> future:allFutureJobMap.keySet()) {
					Bcf2SingleTree job = allFutureJobMap.get(future);
//					allDone=allDone&&future.isDone();
					if(!future.isDone()) {
						allDone=false;
						
						runningJobFutureMap.put(future, job);
					}else {//the job is done, check if successfully done
						if(job.isSuccessfullyFinished()) {//update the max success run full time length
							totalSuccessfullyDoneJobNum++;
							if(job.getFullSuccessfulRunTime()>maxFullRunTimeForSuccessfullyFinishedJob) {
								maxFullRunTimeForSuccessfullyFinishedJob = job.getFullSuccessfulRunTime();
							}
						}
					}
				}
				
				//
				if(allDone)
					continue;
				
				//
				log.info("max full run time for successfully finished job:"+TimeUtils.getReadableTimeFromNanoTime(maxFullRunTimeForSuccessfullyFinishedJob));
				
				//print out infor for all running jobs
				StringBuilder sb=new StringBuilder();
				int runningJobNum=0;
				for(Future<Integer> future:runningJobFutureMap.keySet()) {
					Bcf2SingleTree job = runningJobFutureMap.get(future);
					if(job.started()) {
						runningJobNum++;
						if(!sb.toString().isEmpty()) {
							sb.append(";");
						}
						sb.append("id(").append(job.getUniqueRegionIdentifier()).append(")");
						
						sb.append(":run time(").append(TimeUtils.getReadableTimeFromNanoTime(job.hasRunFor())).append(")");
					}
				}
				
				log.info("total jobs num:"+allFutureJobMap.size()+
						"; total undone jobs num:"+runningJobFutureMap.size()+
						"; running jobs:"+runningJobNum+" (started by not finished) ===> ("+sb.toString()+")");
				
				//restart jobs that run too long if the finished job (both success or unsuccess) num/total job num > ratio or 
				if((double)(allFutureJobMap.size()-runningJobFutureMap.size())/allFutureJobMap.size()>ratio||totalSuccessfullyDoneJobNum>minSuccessfullyDoneJobNum) {
					log.info("finished job num/total job num >"+ratio+" or totalSuccessfullyDoneJobNum >"+minSuccessfullyDoneJobNum+"; start checking if running jobs need to be restarted...");
					//check each running jobs, if its running time exceed
					StringBuilder cancledAndRestartedJobInfo = new StringBuilder();
					for(Future<Integer> future:runningJobFutureMap.keySet()) {
						Bcf2SingleTree job = runningJobFutureMap.get(future);
						if(job.started() && job.hasRunFor()>maxFullRunTimeForSuccessfullyFinishedJob*flexibility) {
							log.warn("job of tree id=="+job.getUniqueRegionIdentifier()+" has run for too long, cancel it...");
							//cancel job
							boolean s = future.cancel(true);
//							future.isCancelled();
							
							if(!future.isCancelled())
								log.error("job cancellation of tree id=="+job.getUniqueRegionIdentifier()+" is failed!"+s);
							
							//rest before re-submit
							job.reset();
							
							//remove from map
							allFutureJobMap.remove(future);
							//restart the job
							allFutureJobMap.put(executorService.submit(job), job);
							
							if(!cancledAndRestartedJobInfo.toString().isEmpty()) {
								cancledAndRestartedJobInfo.append(",");
							}
							cancledAndRestartedJobInfo.append(job.getUniqueRegionIdentifier());
						}
					}
					
					if(!cancledAndRestartedJobInfo.toString().isEmpty())
						log.warn("canceled and restarted jobs:"+ cancledAndRestartedJobInfo.toString());
				}
				
				////
			}
			
			//shut down ExecutorService so that program can terminate
			executorService.shutdown();
			/////////////////
			log.info("all done!");
			
			//summarize the termination status of each submitted job
			log.info("summarizing termination status of submitted jobs...");
			List<Integer> terminationStatusList = new ArrayList<>();
			Map<Integer, Integer> terminationStatusJobNumberMap = new LinkedHashMap<>();
			for(Future<Integer> future:allFutureJobMap.keySet()) {
				int status=future.get();
				if(!terminationStatusJobNumberMap.containsKey(status)) {
					terminationStatusJobNumberMap.put(status, 0);
					terminationStatusList.add(status);
				}
				terminationStatusJobNumberMap.put(status, terminationStatusJobNumberMap.get(status)+1);
			}
			
			//print to standard output
			Collections.sort(terminationStatusList);
			log.info("In total "+allFutureJobMap.size()+" submitted jobs:");
			for(int terminationStatus:terminationStatusList) {
				log.info(terminationStatusJobNumberMap.get(terminationStatus)+" have termination status = "+terminationStatus);
			}
			
		} catch (IOException | InterruptedException | ExecutionException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
	}
	
	/**
	 * to run this in linux, need to load java 8 or above, bcftools and phylip
	 * @param args
	 */
	public static void main(String[] args) {
		if(args.length!=10) {
			System.out.println("java bcfFile regionFile outputTreeFile outgroupName minAlignmentLen tmpDir threadNum maxMissingGenotype minQUALForVariantSites megaxMaoFile");
			System.exit(1);
		}
		
		Path bcfFile = Path.of(args[0]);
		Path regionFile=Path.of(args[1]);
		Path outputTreeFile=Path.of(args[2]);
		String outgroupName = args[3];
		int minAlignmentLen = Integer.parseInt(args[4]);
		Path rootTmpDir=Path.of(args[5]);
		int threadNum=Integer.parseInt(args[6]);
		
		//additional parameters
		int maxMissingGenotype = Integer.parseInt(args[7]); //20
		int minQUALForVariantSites = Integer.parseInt(args[8]); //30
		Path megaxMaoFile=Path.of(args[9]);
		/////////////////////////////////////////
		//check bcf file
		if(!bcfFile.toString().endsWith(".bcf")) {
			log.error("bcf file name is not ended with .bcf:"+bcfFile.toString());
			System.exit(1);
		}
		if(!bcfFile.toFile().exists()) {
			log.error("bcf file is not found:"+bcfFile.toString());
			System.exit(1);
		}
		
		//index the input bcf file if not indexed yet
		Path bcfFileIndex=Path.of(bcfFile.toString().concat(".csi"));
		if(!bcfFileIndex.toFile().exists()) {
			log.info("index bcf file...");
			BcfUtils.indexBcfFile(bcfFile);
		}else {
			log.info("csi index file of bcf file is found!");
		}
		
		//check region file
		if(!regionFile.toFile().exists()) {
			log.error("region file is not found:"+regionFile.toString());
			System.exit(1);
		}
		
		//check outputTreeFile
		if(outputTreeFile.toFile().exists()) {
			log.error("existing file with the same name with output tree file is found!"+outputTreeFile.toString());
			System.exit(1);
		}
		
		//check root tmp dir
		if(rootTmpDir.toFile().exists()) {
			if(rootTmpDir.toFile().isFile()) {
				log.error("existing file with the same path of the root tmp directory is found!");
				System.exit(1);
			}else {
				try {
					FileUtils.cleanDirectory(rootTmpDir.toFile());
				} catch (IOException e) {
					e.printStackTrace();
					log.error("root tmp dir cannot be cleaned! "+rootTmpDir.toString());
					System.exit(1);
				}
			}
		}else {
			if(!rootTmpDir.toFile().mkdir()) {
				log.error("root tmp directory cannot be created!");
				System.exit(1);
			}
		}
		
		
		//check if the outgroupName exist in the bcf file
		//first extract header section of bcf file into a vcf file (VCFFileReader cannot read bcf file encoded by current specification!)
//		File headerSectionVcfFile = new File(rootTmpDir.toFile().getAbsolutePath() + File.separator + "header.vcf");
		Path headerSectionVcfFile = Path.of(rootTmpDir.toString(), "header.vcf");
		BcfUtils.queryHeaderRegionFromBcfFile(bcfFile.toString(), headerSectionVcfFile.toString());
		VCFFileReader reader = new VCFFileReader(headerSectionVcfFile, false);
		if(!reader.getFileHeader().getSampleNamesInOrder().contains(outgroupName)) {
			log.error("given outgroupName: "+outgroupName+" is not found in the header line of the given bcf file!");
			headerSectionVcfFile.toFile().delete();
			System.exit(1);
		}
		headerSectionVcfFile.toFile().delete();
		reader.close();
		
		//
		if(minAlignmentLen<=0) {
			log.error("given minAlignmentLen: "+threadNum+" is not positive integer!");
			System.exit(1);
		}
		
		
		//check 
		if(threadNum<=0) {
			log.error("given threadNum: "+threadNum+" is not positive integer!");
			System.exit(1);
		}
		
		//check
		if(!megaxMaoFile.toFile().exists()) {
			log.error("mega10 mao file is not found:"+megaxMaoFile.toString());
			System.exit(1);
		}
		
		
		//////////////////////////////additional parameters for bcf record filtering and recorder from genotype to sequence
//		int maxMissingGenotype = 20;
//		int minQUALForVariantSites = 30;
		Vcf2AlignmentSeqBuilderFactory msaBuilderFactory = new Vcf2AlignmentSeqBuilderFactory(
//				VariantContextFilterFactory.onlySNPSites().and(VariantContextFilterFactory.maxMissingCount(maxMissingGenotype)), 
				VariantContextFilterFactory.nonIndelSite().and(VariantContextFilterFactory.nonMixedTypeSite()).and(VariantContextFilterFactory.filterVariantSitesByQual(minQUALForVariantSites)).and(VariantContextFilterFactory.allIndividualsAreGenericHomozygous()).and(VariantContextFilterFactory.maxMissingCount(maxMissingGenotype)),
//				GenotypeRecoderFactory.biBaseSeqRecoder()
				GenotypeRecoderFactory.singleBaseSeqRecoder()
				);
		
//		Path maoFile=Path.of("/home/tanxu/phylogeny/megaX/infer_NJ_nucleotide_pairwise_deletion.mao");
		MultipleAlignment2TreeFactory<?> msa2treeFactory = new MegaXBasedTreeBuilderFactory(megaxMaoFile);
		
		Bcf2SingleTreeFactory bcf2SingleTree2Factory = 
				new Bcf2SingleTreeFactory(
						bcfFile, msaBuilderFactory, msa2treeFactory, outgroupName, minAlignmentLen);
		
		////////////////////////////
		BatchRegionalTreeManager manager = 
				new BatchRegionalTreeManager(bcf2SingleTree2Factory, regionFile, outputTreeFile, rootTmpDir, threadNum);
		
		///////////////////////
		long startTime=System.nanoTime();
		manager.run();
		log.info("elpased time:"+TimeUtils.getReadableTimeFromNanoTime(System.nanoTime() - startTime));
		//
		System.exit(0);
	}
	
}
