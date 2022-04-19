package phylo.tree.dist;

import java.io.BufferedWriter;
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
import phylo.batch.BcfUtils;
import phylo.concurrency.TimeUtils;
import phylo.ref.Region;
import phylo.ref.RegionUtils;
import phylo.tree.dist.diff.PairwiseDiffMatrixFromVcf;
import population.vcf.utils.GenotypeRecoderFactory;
import population.vcf.utils.VariantContextFilterFactory;

/**
 * calculate the nucleotide no. difference matrix for all samples in a bcf files for each region partitioned with a specific length;
 * 
 * @author tanxu
 * 
 */
public class BatchRegionalPairwiseDiffMatrixFromVcfManager{
	public static Log log=Log.getInstance(BatchRegionalPairwiseDiffMatrixFromVcfManager.class);
	//////////////////
	
	///////////////////////////////////
	/**
	 * bcf file containing all the variant and/or non-variant sites of all individuals to be processed
	 */
	private final Path bcfFile;
	
	/**
	 * the bed file in which all regions are to be calculated together
	 */
	private final Path regionBedFile;
	
	/**
	 * length of region to parition the whole genome in {@link #regionBedFile}
	 */
	private final int regionLen;
	
	
	/**
	 * minimal QUAL value for variant site to be included as informative site to build the tree
	 */
	private final int minQUALForVariantSites;
	
	/**
	 * 
	 */
	private final int maxMissingGenotype;
	
	/**
	 * root output directory
	 */
	private final Path rootOutDir;
	
	/**
	 * thread number
	 */
	private final int threadNum;
	
	
	/////////////////////////////////////////
	/**
	 * the directory to put all pairwise difference matrix for each region
	 */
	private Path regionDiffMatrixOutDir;
	
	/**
	 * 
	 */
	private Path rootTmpDir;
	
	/**
	 * 
	 */
	private List<String> orderedSampleNameList;
	
	
	/**
	 * map from the index of a region to the region string in the format of -r parameter of bcftools
	 */
	private Map<Integer, String> regionIndexStringMap;
	
	/**
	 * path to the file for each region in the {@link #regionBedFile} with unique integer index
	 * 	four columns: chrom	start	end	index
	 */
	private Path outputRegionIndexFile;
	
	
	/**
	 * 
	 * @param bcfFile
	 * @param regionBedFile
	 * @param outputPath
	 * @param outputFileBaseName
	 * @param threadNum
	 * @throws IOException 
	 */
	public BatchRegionalPairwiseDiffMatrixFromVcfManager(
			Path bcfFile, Path regionBedFile, int regionLen, 
			int minQUALForVariantSites, int maxMissingGenotype,
			Path outputPath, int threadNum) throws IOException{
		
		/////////////////////////////////////////validations
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
		if(!regionBedFile.toFile().exists()) {
			log.error("region bed file is not found:"+regionBedFile.toString());
			System.exit(1);
		}
		
		
		//check 
		if(regionLen<=0) {
			log.error("given region length: "+regionLen+" is not positive integer!");
			System.exit(1);
		}
		
		
		
		//check output path
		if(outputPath.toFile().exists()) {
			if(outputPath.toFile().isFile()) {
				log.error("existing file with the same path of the outputPath is found!");
				System.exit(1);
			}else {
				try {
					FileUtils.cleanDirectory(outputPath.toFile());
				} catch (IOException e) {
					e.printStackTrace();
					log.error("root tmp dir cannot be cleaned! "+outputPath.toString());
					System.exit(1);
				}
			}
		}else {
			if(!outputPath.toFile().mkdir()) {
				log.error("outputPath cannot be created!");
				System.exit(1);
			}
		}
		
		
		//check 
		if(threadNum<=0) {
			log.error("given threadNum: "+threadNum+" is not positive integer!");
			System.exit(1);
		}
		
		
		/////////////////////////////////
		this.bcfFile=bcfFile;
		this.regionBedFile=regionBedFile;
		this.regionLen=regionLen;
		this.minQUALForVariantSites=minQUALForVariantSites;
		this.maxMissingGenotype=maxMissingGenotype;
		this.rootOutDir=outputPath;
		this.threadNum=threadNum;
		
//		try {
			this.preprocess();
//		} catch (IOException e) {
//			log.error("Exception thrown when preprocessing..."+e.getMessage());
//		}
	}
	
	//////////////////////////////////////
	/**
	 * @throws IOException 
	 * 
	 */
	private void preprocess() throws IOException {
		//
		this.regionDiffMatrixOutDir = Path.of(this.rootOutDir.toString(), "region_diff_matrix_out");
		if(this.regionDiffMatrixOutDir.toFile().exists()) {
			if(this.regionDiffMatrixOutDir.toFile().isFile()) {
				log.error("existing file with the same path of the region diff matrix output dir is found!");
				System.exit(1);
			}else {
				try {
					FileUtils.cleanDirectory(this.regionDiffMatrixOutDir.toFile());
				} catch (IOException e) {
					e.printStackTrace();
					log.error("region diff matrix output dir cannot be cleaned! "+this.regionDiffMatrixOutDir.toString());
					System.exit(1);
				}
			}
		}else {
			if(!this.regionDiffMatrixOutDir.toFile().mkdir()) {
				log.error("region diff matrix output dir cannot be created!");
				System.exit(1);
			}
		}
		
		
		//create and check rootTmpDir
		this.rootTmpDir=Path.of(this.rootOutDir.toString(), "regional_vcf_files");
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
		
		
		////create VcfFileReader and create the ordered list of sample names
		this.orderedSampleNameList=new ArrayList<>();
		Path headerSectionVcfFile = Path.of(rootTmpDir.toString(), "header.vcf");
		BcfUtils.queryHeaderRegionFromBcfFile(bcfFile.toString(), headerSectionVcfFile.toString());
		VCFFileReader reader = new VCFFileReader(headerSectionVcfFile, false);
		reader.getFileHeader().getSampleNamesInOrder().forEach(n->{
			orderedSampleNameList.add(n);
		});
		headerSectionVcfFile.toFile().delete();
		reader.close();
		
		//check and read the regionBedFile and divide into regions
		this.regionIndexStringMap = RegionUtils.readBedFileIntoRegionString(this.regionBedFile, regionLen);
		
		
		//output to the outputRegionIndexFile
		this.outputRegionIndexFile = Path.of(this.rootOutDir.toString(),"region_index.txt");
	    BufferedWriter writer = new BufferedWriter(new FileWriter(this.outputRegionIndexFile.toString(), true));
	    for(int index:this.regionIndexStringMap.keySet()) {
	    	Region r = Region.fromString(this.regionIndexStringMap.get(index));
	    	writer.append(
	    			r.getReferenceName().concat("\t")
	    			.concat(Integer.toString(r.getStart()).concat("\t")
	    					.concat(Integer.toString(r.getEnd()).concat("\t").concat(Integer.toString(index)))));
	    
	    	writer.newLine();
	    }
	    writer.flush();
	    writer.close();
		
		
	}
	
	
	/**
	 * 1. divide the full region into smaller regions
	 * @throws ExecutionException 
	 * @throws InterruptedException 
	 */
	public void run() throws InterruptedException, ExecutionException {
		//////////////////
		ExecutorService executorService = Executors.newFixedThreadPool(this.threadNum); 
		
		Map<Future<Integer>, PairwiseDiffMatrixFromVcf> allFutureJobMap = new LinkedHashMap<>();
		
		for(int index:this.regionIndexStringMap.keySet()) {
			String regionString = this.regionIndexStringMap.get(index);
			
			PairwiseDiffMatrixFromVcf collector = 
					new PairwiseDiffMatrixFromVcf(
							this.orderedSampleNameList,//<String> orderedSampleNameList, 
							this.bcfFile,// bcfFile, 
							regionString,//String regionListString, 
							Integer.toString(index),// String uniqueRegionIdentifier, 
							VariantContextFilterFactory.nonIndelSite().and(VariantContextFilterFactory.nonMixedTypeSite()).and(VariantContextFilterFactory.filterVariantSitesByQual(minQUALForVariantSites)).and(VariantContextFilterFactory.allIndividualsAreGenericHomozygous()).and(VariantContextFilterFactory.maxMissingCount(maxMissingGenotype)),//VariantContextFilter variantFilter, 
							GenotypeRecoderFactory.singleBaseSeqRecoder(),//GenotypeRecoder genotypeRecoder,
							this.regionDiffMatrixOutDir,
							this.rootTmpDir//Path tmpDir, 
							);
			
			allFutureJobMap.put(executorService.submit(collector), collector);
		}
		
		
		//monitor all submitted jobs until all done (due to normal termination, an exception, or cancellation)
		boolean allDone=false;
		while(!allDone) {
			Thread.sleep(60000);
			
			
			log.info("update job status ...");
			allDone=true;
			
			/////////
			//check status of all jobs
			for(Future<Integer> future:allFutureJobMap.keySet()) {
				if(!future.isDone()) {
					allDone=false;
//						runningJobFutureMap.put(future, job);
				}
			}
			
			//
			if(allDone)
				continue;
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
		
	}

	/**
	 * @return the regionDiffMatrixOutDir
	 */
	public Path getRegionDiffMatrixOutDir() {
		return regionDiffMatrixOutDir;
	}

	/**
	 * @return the orderedSampleNameList
	 */
	public List<String> getOrderedSampleNameList() {
		return orderedSampleNameList;
	}

	/**
	 * @return the outputRegionIndexFile
	 */
	public Path getOutputRegionIndexFile() {
		return outputRegionIndexFile;
	}

	////////////////////////////////
	/**
	 * to run this in linux, need to load java 8 or above, bcftools and phylip
	 * @param args
	 * @param outputTreeFile 
	 */
	public static void main(String[] args) {
		if(args.length!=7) {
			System.out.println("java this bcfFile regionBedFile regionLen minQUALForVariantSites maxMissingGenotype outputPath threadNum");
			System.exit(1);
		}
		Path bcfFile = Path.of(args[0]);
		Path regionBedFile=Path.of(args[1]);
		int regionLen = Integer.parseInt(args[2]);
		int minQUALForVariantSites=Integer.parseInt(args[3]);
		int maxMissingGenotype=Integer.parseInt(args[4]);
		Path outputPath=Path.of(args[5]);
		int threadNum=Integer.parseInt(args[6]);
		
		////////////////////////////
		
		try {
			BatchRegionalPairwiseDiffMatrixFromVcfManager manager = new BatchRegionalPairwiseDiffMatrixFromVcfManager(
					bcfFile, 
					regionBedFile,
					regionLen,
					minQUALForVariantSites,
					maxMissingGenotype,
					outputPath,
					threadNum
					);
			long startTime=System.nanoTime();
			manager.run();
			log.info("elpased time:"+TimeUtils.getReadableTimeFromNanoTime(System.nanoTime() - startTime));
			//
			System.exit(0);
		} catch (IOException | InterruptedException | ExecutionException e) {
			e.printStackTrace();
			log.error("Exception thrown :"+e.getMessage());
			System.exit(1);
		}
	}
	
	
	
}
