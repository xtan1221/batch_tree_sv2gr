package phylo.tree.dist.nj;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.function.Predicate;

import htsjdk.samtools.util.Log;
import phylo.concurrency.TimeUtils;
import phylo.ref.Region;
import phylo.tree.dist.DistMatrixUtils;
import phylo.tree.dist.diff.merge.MergerUtils;
import phylo.tree.dist.diff.merge.RegionDiffMatrixAndNonMissingSiteNumMatrixFilesMerger;
import phylo.tree.phylo.PhylipBashRunner;
import phylo.tree.reader.NewickFileFormatType;
import phylo.tree.reader.Tree;

/**
 * 
 * @author tanxu
 *
 */
public class DiffMatrixToNJTree {
	public static Log log=Log.getInstance(DiffMatrixToNJTree.class);
	
	////////////////////////
	/**
	 * the directory where all the regional matrix file are stored;
	 * 
	 * including both non-missing sites num matrix and nucleotide diff no. matrix for each region in the same batch run
	 * 
	 */
	private final Path regionalMatrixFilesOutDir;
	/**
	 * four column file each line containing information of a region
	 * col 1 is the chrom name; col 2 is the start of the region, col 3 is the end of the region, col 4 is the index of the region
	 * 
	 * note that the diff matrix file of a region with index i has name {i}{diffMatrixFileSuffix} in the diffMatrixOutDir
	 * 		for example, if diffMatrixFileSuffix is "_diff_matrix.txt" and region index is 5, the file name is "5_diff_matrix.txt"
	 */
	private final Path regionIndexFile;
	
	/**
	 * the map from the data name to the filter of regions to be merged to the data
	 */
	private final Map<String,Predicate<Region>> outputDataNameRegionFilterMap;
	
	///////////////////
	/**
	 * number of sequence in each matrix
	 */
	private final int seqNum;
	
	////////////////
	private final String outgroupName;
	private final Path outputDir;
	
	////////////////
	private Map<String, Tree> dataNameRerootedTreeMap;
	/**
	 * file to put all built tree with three columns
	 * col1 data name
	 * col2 regionString
	 * col3 full newick tree string
	 */
	private Path mergedAllTreeFile;
	
	
	public DiffMatrixToNJTree(
		Path diffMatrixOutDir, Path regionIndexFile, 
		Map<String,Predicate<Region>> outputDataNameRegionFilterMap,
		int seqNum, String outgroupName, Path outputDir) {
		//TODO validations
		this.regionalMatrixFilesOutDir = diffMatrixOutDir;
		this.regionIndexFile = regionIndexFile;
		this.outputDataNameRegionFilterMap = outputDataNameRegionFilterMap;
		this.seqNum=seqNum;
		/////////////
		this.outgroupName=outgroupName;
		this.outputDir=outputDir;
		
	}
	
	
	public void run() throws NumberFormatException, IOException {
		log.info("start merge diff matrix for each of the targeted dataset...");
		//merge the non-missing site num matrix and nuc diff no matrix of all regions
		RegionDiffMatrixAndNonMissingSiteNumMatrixFilesMerger merger = 
				new RegionDiffMatrixAndNonMissingSiteNumMatrixFilesMerger(
						this.regionalMatrixFilesOutDir, 
						this.regionIndexFile, 
						this.outputDataNameRegionFilterMap, 
						this.seqNum);
		merger.run();
		
		log.info("merging diff matrix for each of the targeted dataset is done...");
		
		
		//for each data
		dataNameRerootedTreeMap = new LinkedHashMap<>();
		
		log.info("start processing merged diff matrix of each targeted data set and build tree...");
		for(String dn:merger.getDataNameSummedPairwiseNonMissingSitesNumMatrixMap().keySet()) {
			log.info("calculate Jukes-Cantor distance matrix for data =="+dn);
			//calculate the jc distance matrix
			int[][] diffMatrix = merger.getDataNameSummedPairwiseDiffMatrixMap().get(dn);
			long[][] nonMissingSiteNumMatrix=merger.getDataNameSummedPairwiseNonMissingSitesNumMatrixMap().get(dn);
			
			//calculate the p-distance matrix
			double[][] pDistMatrix = DistMatrixUtils.calculatePDistaFromDiffMatrix(nonMissingSiteNumMatrix, diffMatrix);
			//Calculate the JC model based distance matrix
			double[][] jcDistMatrix = DistMatrixUtils.calcualteJCDistMatrixFromPDistMatrix(pDistMatrix);
			
			log.info("write Jukes-Cantor distance matrix to file for data =="+dn);
			//output directory for current data
			Path dataOutDir = Path.of(this.outputDir.toString(),dn);
			if(!dataOutDir.toFile().exists()) {
				dataOutDir.toFile().mkdir();
			}
			
			//write diff matrix to file
			Path diffMatrixFile = Path.of(dataOutDir.toString(),"diff.no.matrix.txt");
			if(diffMatrixFile.toFile().exists()) {
				diffMatrixFile.toFile().delete();
			}
			DistMatrixUtils.writeMatrixToFile(diffMatrixFile, merger.getSeqNameList(), diffMatrix);
			
			//write non missing sites number matrix to file
			Path nonMissingSiteNumMatrixFile = Path.of(dataOutDir.toString(),"all.non.missing.sites.no.matrix.txt");
			if(nonMissingSiteNumMatrixFile.toFile().exists()) {
				nonMissingSiteNumMatrixFile.toFile().delete();
			}
			DistMatrixUtils.writeMatrixToFile(nonMissingSiteNumMatrixFile, merger.getSeqNameList(), nonMissingSiteNumMatrix);
			
			//output file for p-dist matrix
			Path pDistanceMatrixPhylipFile = Path.of(dataOutDir.toString(),"p.dist.phylip.matrix.txt");
			if(pDistanceMatrixPhylipFile.toFile().exists()) {
				pDistanceMatrixPhylipFile.toFile().delete();
			}
			DistMatrixUtils.writeToPhylipDistMatrixFile(pDistanceMatrixPhylipFile, merger.getSeqNameList(), DistMatrixUtils.toStringMatrix(pDistMatrix));
			
			//output file for jc distance matrix in phylip format
			Path njDistanceMatrixPhylipFile = Path.of(dataOutDir.toString(),"jc.dist.phylip.matrix.txt");
			if(njDistanceMatrixPhylipFile.toFile().exists()) {
				njDistanceMatrixPhylipFile.toFile().delete();
			}
			DistMatrixUtils.writeToPhylipDistMatrixFile(njDistanceMatrixPhylipFile, merger.getSeqNameList(), DistMatrixUtils.toStringMatrix(jcDistMatrix));
			
			
			///////////////////////////////below tested in linux
			log.info("build neighbor joining tree with Jukes-Cantor distance matrix for data =="+dn);
			Path njNewickTreeFile = PhylipBashRunner.runDistMatrixToNJTreeBash(njDistanceMatrixPhylipFile.toString(), dataOutDir.toString(), dn, 1);
			
			log.info("reroot tree with outgroup=="+this.outgroupName+" for data =="+dn);
			//reroot tree
			Tree tree = Tree.fromNewickFile(njNewickTreeFile, NewickFileFormatType.SIMPLE_NEWICK_1);
			Tree rerootedTree = tree.reroot(this.outgroupName);
			
			dataNameRerootedTreeMap.put(dn, rerootedTree);
			//
			log.info("write rerooted tree to newick file for data =="+dn);
			Path rerootedTreeFile = Path.of(dataOutDir.toString(),dn.concat(".reroot.nj.nwk"));
			BufferedWriter writer = new BufferedWriter(new FileWriter(rerootedTreeFile.toFile(), false));
		    writer.write(rerootedTree.toFullNewickString(NewickFileFormatType.SIMPLE_NEWICK_2));
		    
		    writer.close();
		    
		    log.info("diff matrix to NJ tree for data =="+dn+" is done!");
		}
		
		///write all trees to mergedAllTreeFile
		log.info("write all trees to a single file ...");
		this.mergedAllTreeFile = Path.of(this.outputDir.toString(),"all.trees.txt");
		if(this.mergedAllTreeFile.toFile().exists()) {
			this.mergedAllTreeFile.toFile().delete();
		}
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(this.mergedAllTreeFile.toFile(), true));
		for(String dn:this.dataNameRerootedTreeMap.keySet()) {
			writer.append(dn.concat("\t").concat(dn).concat("\t").concat(this.dataNameRerootedTreeMap.get(dn).toFullNewickString(NewickFileFormatType.SIMPLE_NEWICK_2)));
			writer.newLine();
		}
	    writer.flush();
	    writer.close();
	    log.info("write all trees to a single file is done!");
	}
	
	
	public static void main(String[] args) throws IOException {
		if(args.length!=5) {
			System.out.println("java diffMatrixOutDir regionIndexFile seqNum outgroupName outputDir");
			System.exit(1);
		}
		
		Path diffMatrixOutDir = Path.of(args[0]);
		Path regionIndexFile = Path.of(args[1]);
		
		int seqNum = Integer.parseInt(args[2]);
		String outgroupName = args[3];
		Path outputDir = Path.of(args[4]);
		
		//////////////////////
		// validations
		if(!diffMatrixOutDir.toFile().exists()) {
			log.error("diffMatrixOutDir is not found:"+diffMatrixOutDir.toString());
			System.exit(1);
		}
		if(!regionIndexFile.toFile().exists()) {
			log.error("regionIndexFile is not found:"+regionIndexFile.toString());
			System.exit(1);
		}
		if(seqNum<=1) {
			log.error("seqNum must be no less than 2!");
		}
		if(!outputDir.toFile().exists()) {
			log.error("outputDir is not found:"+outputDir.toString());
			System.exit(1);
		}
		
		///////////////////////
		Map<String,Predicate<Region>> outputDataNameRegionFilterMap = new LinkedHashMap<>();
		outputDataNameRegionFilterMap.put("all_chrom", e->{return true;});
		outputDataNameRegionFilterMap.putAll(MergerUtils.makeChromNamePredicateMap(regionIndexFile, 0));
		
		
		
		//////////////////////
		DiffMatrixToNJTree diffMatrixToNJTree = new DiffMatrixToNJTree(
				diffMatrixOutDir, regionIndexFile, 
				outputDataNameRegionFilterMap,
				seqNum,
				outgroupName, outputDir);
		
		
		//////////////////
		long startTime=System.nanoTime();
		diffMatrixToNJTree.run();
		log.info("elpased time:"+TimeUtils.getReadableTimeFromNanoTime(System.nanoTime() - startTime));
		//
		System.exit(0);
	}
}
