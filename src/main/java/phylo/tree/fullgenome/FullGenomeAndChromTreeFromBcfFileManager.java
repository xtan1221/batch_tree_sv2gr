package phylo.tree.fullgenome;

import java.io.IOException;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.util.function.Predicate;

import org.apache.commons.io.FileUtils;

import htsjdk.samtools.util.Log;
import phylo.batch.BcfUtils;
import phylo.concurrency.TimeUtils;
import phylo.ref.Region;
import phylo.tree.dist.BatchRegionalPairwiseDiffMatrixFromVcfManager;
import phylo.tree.dist.diff.PairwiseDiffMatrixFromVcf;
import phylo.tree.dist.diff.merge.MergerUtils;
import phylo.tree.dist.nj.DiffMatrixToNJTree;
/**
 * build tree for the whole genome and each chromosome from data in a bcf file
 * 
 * must be running in linux with following packages loaded
 * 1. java
 * 2. bcftools (bcf region query)
 * 3. phylip (tree building)
 * 
 * the core class is {@link PairwiseDiffMatrixFromVcf}
 * 
 * @author tanxu
 */
public class FullGenomeAndChromTreeFromBcfFileManager {
	public static Log log=Log.getInstance(FullGenomeAndChromTreeFromBcfFileManager.class);
	////////////////////////////
	/**
	 * bcf file containing all the variant and/or non-variant sites of all individuals to be processed
	 */
	private final Path bcfFile;
	
	/**
	 * the bed file in which all regions are to be calculated together
	 */
	private final Path regionBedFile;
	
	/**
	 * root output directory
	 */
	private final Path rootOutDir;
	
	/**
	 * 
	 */
	private final String outgroupName;
	
	/**
	 * length of region to partition the whole genome in {@link #regionBedFile}
	 * 
	 * for each region, a pairwise nucleotide diff matrix will be calcualted and all regions' matrix will be added up to build the matrix for whole chrom and genome
	 * then these matrices will be used to build p-distance matrix or other distance matrices
	 * 
	 * see {@link PairwiseDiffMatrixFromVcf} for more details
	 */
	private final int regionLen;
	
	/**
	 * minimal QUAL value for variant site to be included as informative site to build the tree
	 */
	private final int minQUALForVariantSites;
	
	/**
	 * max number of samples with missing genotype so that the site will be included in the tree building
	 */
	private final int maxMissingGenotype;
	/**
	 * thread number
	 */
	private final int threadNum;
	
	////////////////////////
	private BatchRegionalPairwiseDiffMatrixFromVcfManager batchRegionalPairwiseDiffMatrixFromVcfManager;
	private DiffMatrixToNJTree diffMatrixToNJTree;
	
	
	
	
	FullGenomeAndChromTreeFromBcfFileManager(
			Path bcfFile, Path regionBedFile, Path rootOutDir, String outgroupName, 
			int regionLen, 
			int minQUALForVariantSites, int maxMissingGenotype,
			int threadNum){
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
		if(rootOutDir.toFile().exists()) {
			if(rootOutDir.toFile().isFile()) {
				log.error("existing file with the same path of the rootOutDir is found!");
				System.exit(1);
			}else {
				try {
					FileUtils.cleanDirectory(rootOutDir.toFile());
				} catch (IOException e) {
					e.printStackTrace();
					log.error("root tmp dir cannot be cleaned! "+rootOutDir.toString());
					System.exit(1);
				}
			}
		}else {
			if(!rootOutDir.toFile().mkdir()) {
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
		this.rootOutDir=rootOutDir;
		this.outgroupName=outgroupName;
		this.regionLen=regionLen;
		this.minQUALForVariantSites=minQUALForVariantSites;
		this.maxMissingGenotype=maxMissingGenotype;
		this.threadNum=threadNum;
	}
	
	
	
	void run() {
		log.info("start partition genome into regions and calculate regional diff matrix...");
		try {
			this.partitionGenomeAndBuildRegionalDiffMatrix();
		} catch (IOException | InterruptedException | ExecutionException e) {
			log.error("Exception thrown when partition genome into regions and calculate regional diff matrix:"+e.getMessage());
			System.exit(1);
		}
		log.info("partition genome into regions and calculating regional diff matrix is done");
		
		
		log.info("start building tree for whole genome and each chromosome...");
		try {
			this.buildChromAndGenomeTree();
		} catch (IOException e) {
			log.error("Exception thrown when building tree for whole genome and each chromosome:"+e.getMessage());
			System.exit(1);
		}
		log.info("building tree for whole genome and each chromosome is done");
	}
	
	/**
	 * build regional diff matrix
	 * major steps
	 * 1. partition full genome into regions
	 * 2. calculate regional diff matrix for each region in concurrency
	 * @throws IOException 
	 * @throws ExecutionException 
	 * @throws InterruptedException 
	 * 
	 */
	private void partitionGenomeAndBuildRegionalDiffMatrix() throws IOException, InterruptedException, ExecutionException {
		batchRegionalPairwiseDiffMatrixFromVcfManager = 
				new BatchRegionalPairwiseDiffMatrixFromVcfManager(
						bcfFile, regionBedFile, regionLen, 
						minQUALForVariantSites, maxMissingGenotype,
						rootOutDir, threadNum
						);
		
		batchRegionalPairwiseDiffMatrixFromVcfManager.run();
	}
	
	/**
	 * 
	 * merge single region diff matrix into full chromosome and whole genome diff matrix
	 * major steps 
	 * see {@link DiffMatrixToNJTree#run()}
	 * @throws IOException 
	 */
	private void buildChromAndGenomeTree() throws IOException {
		Path diffMatrixOutDir = this.batchRegionalPairwiseDiffMatrixFromVcfManager.getRegionDiffMatrixOutDir();
		Path regionIndexFile = this.batchRegionalPairwiseDiffMatrixFromVcfManager.getOutputRegionIndexFile();
		Map<String,Predicate<Region>> outputDataNameRegionFilterMap = new LinkedHashMap<>();
		outputDataNameRegionFilterMap.put("all_chrom", e->{return true;});
		outputDataNameRegionFilterMap.putAll(MergerUtils.makeChromNamePredicateMap(regionIndexFile, 0));
		
		int seqNum = this.batchRegionalPairwiseDiffMatrixFromVcfManager.getOrderedSampleNameList().size();
		Path outputDir = Path.of(this.rootOutDir.toString(),"genome_chrom_tree");
		
		if(!outputDir.toFile().exists()) {
			outputDir.toFile().mkdir();
		}
		
		diffMatrixToNJTree = 
				new DiffMatrixToNJTree(
						diffMatrixOutDir, regionIndexFile, 
						outputDataNameRegionFilterMap,
						seqNum, outgroupName, outputDir);
		
		this.diffMatrixToNJTree.run();
	}
	
	
	
	
	/**
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		if(args.length!=8) {
			System.out.println("java this bcfFile regionBedFile rootOutDir outgroupName regionLen minQUALForVariantSites maxMissingGenotype threadNum");
			System.exit(1);
		}
		
		Path bcfFile = Path.of(args[0]);
		Path regionBedFile=Path.of(args[1]);
		Path rootOutDir=Path.of(args[2]);
		String outgroupName = args[3];
		int regionLen=Integer.parseInt(args[4]);
		int minQUALForVariantSites=Integer.parseInt(args[5]);
		int maxMissingGenotype=Integer.parseInt(args[6]);
		
		int threadNum=Integer.parseInt(args[7]);
		
		
		////////////////////////
		FullGenomeAndChromTreeFromBcfFileManager manager = 
				new FullGenomeAndChromTreeFromBcfFileManager(
				bcfFile, regionBedFile, rootOutDir, outgroupName, regionLen, minQUALForVariantSites, maxMissingGenotype, threadNum
				);
		long startTime=System.nanoTime();
		manager.run();
		log.info("elpased time:"+TimeUtils.getReadableTimeFromNanoTime(System.nanoTime() - startTime));
		//
		System.exit(0);
	}
}
