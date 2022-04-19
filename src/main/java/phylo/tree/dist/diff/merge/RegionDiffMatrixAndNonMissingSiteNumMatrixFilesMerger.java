package phylo.tree.dist.diff.merge;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;

import basic.Pair;
import htsjdk.samtools.util.Log;
import phylo.ref.Region;
import phylo.tree.dist.DistMatrixUtils;
import phylo.tree.dist.diff.PairwiseDiffMatrixFromVcf;
import phylo.tree.dist.nj.DiffMatrixToNJTree;

/**
 * for all regions in the same batch run, 
 * 1. this class read the non-missing site num matrix and nuc diff no. matrix for each region
 * 2. build a summed non-missing site num matrix and a summed nuc diff no. matrix for all regions that can be used to calculate the p-distance matrix
 * 
 * 		see {@link #dataNameSummedPairwiseNonMissingSitesNumMatrixMap}
 * 		see {@link #dataNameSummedPairwiseDiffMatrixMap}
 * ======================
 * 
 * this class is utilized by {@link DiffMatrixToNJTree}
 * 
 * @author tanxu
 * 
 */
public class RegionDiffMatrixAndNonMissingSiteNumMatrixFilesMerger{
	public static Log log=Log.getInstance(RegionDiffMatrixAndNonMissingSiteNumMatrixFilesMerger.class);

	/////////////////////////
	/**
	 * the directory where all the regional diff matrix file and non-missing total sites matrix files are stored;
	 */
	private final Path regionMatrixFilesOutDir;
//	/**
//	 * the suffix of regional diff matrix file names
//	 */
//	private final String diffMatrixFileSuffix;
	/**
	 * four column file each line containing information of a region
	 * col 1 is the chrom name; col 2 is the start of the region, col 3 is the end of the region, col 4 is the index of the region
	 * 
	 * for the matrix file name for non-missing sites number 
	 * see {@link PairwiseDiffMatrixFromVcf#buildNonMissingSiteNumMatrixFileName(String)} 
	 * for the matrix file name for the diff nuc num
	 * see {@link PairwiseDiffMatrixFromVcf#buildRegionDiffNumMatrixFileName(String)}
	 * 
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
	
	/////////////////////////////////
	/**
	 * list of sequence name in the same order as in each diff matrix file (note that the order of seq names in all regional diff matrix file should be the same!)
	 */
	private List<String> seqNameList;
//	/**
//	 * 
//	 */
//	private Map<String,Integer> dataNameSeqLenMap;
	
	
	/**
	 * the summed pairwise nucleotide site number with both samples with non-missing data
	 */
	private Map<String,long[][]> dataNameSummedPairwiseNonMissingSitesNumMatrixMap;
	/**
	 * the summed pairwise nucleotide no. difference
	 */
	private Map<String,int[][]> dataNameSummedPairwiseDiffMatrixMap;
	
	
	public RegionDiffMatrixAndNonMissingSiteNumMatrixFilesMerger(
			Path diffMatrixOutDir, 
//			String diffMatrixFileSuffix, 
			Path regionIndexFile, 
			Map<String,Predicate<Region>> outputDataNameRegionFilterMap,
			int seqNum) {
		//TODO validations
		this.regionMatrixFilesOutDir = diffMatrixOutDir;
//		this.diffMatrixFileSuffix = diffMatrixFileSuffix;
		this.regionIndexFile = regionIndexFile;
		this.outputDataNameRegionFilterMap = outputDataNameRegionFilterMap;
		this.seqNum=seqNum;
		/////////////
		this.initialize();
	}
	
	
	private void initialize() {
		
		this.seqNameList=new ArrayList<>();
		
		//////////////////
//		this.dataNameSeqLenMap=new LinkedHashMap<>();
		this.dataNameSummedPairwiseNonMissingSitesNumMatrixMap=new LinkedHashMap<>();
		this.dataNameSummedPairwiseDiffMatrixMap = new LinkedHashMap<>();
		
		this.outputDataNameRegionFilterMap.keySet().forEach(dn->{
//			this.dataNameSeqLenMap.put(dn, 0);
			long[][] summedPairwiseNonMissingSitesNumMatrix=new long[seqNum][seqNum];
			int[][] summedPairwiseDiffMatrix=new int[seqNum][seqNum];
			
			for(int i=0;i<this.seqNum;i++) {
				for(int j=0;j<this.seqNum;j++) {
					summedPairwiseNonMissingSitesNumMatrix[i][j]=0;
					summedPairwiseDiffMatrix[i][j]=0;
				}
			}
			
			this.dataNameSummedPairwiseNonMissingSitesNumMatrixMap.put(dn, summedPairwiseNonMissingSitesNumMatrix);
			this.dataNameSummedPairwiseDiffMatrixMap.put(dn, summedPairwiseDiffMatrix);
		});
		
	}
	
	
	/**
	 * @throws IOException 
	 * @throws NumberFormatException 
	 * 
	 */
	public void run() throws NumberFormatException, IOException {
	    BufferedReader lineReader = new BufferedReader(new FileReader(this.regionIndexFile.toFile()));
	    String line = null;
	    
	    while ((line = lineReader.readLine()) != null) {
	    	String[] splits = line.split("\\s+");
	    	Region region = new Region(splits[0], Integer.parseInt(splits[1]), Integer.parseInt(splits[2]));
	    	String regionIndex = splits[3];
	    	log.info("process region with index="+regionIndex+" and region="+region.toString());
	    	
	    	//////////////////////////add the region non missing site num matrix to the summed matrix
	    	Path nonMissingSitesNumMatrixFile = Path.of(
	    			this.regionMatrixFilesOutDir.toString(),
	    			PairwiseDiffMatrixFromVcf.buildNonMissingSiteNumMatrixFileName(regionIndex));
	    	
	    	log.info("read non-missing sites num matrix from file "+nonMissingSitesNumMatrixFile.toString()+" for region with index="+regionIndex+" and region="+region.toString());
	    	
	    	Pair<List<String>, int[][]> seqNameListNonMissingSiteNumMatrixPair = DistMatrixUtils.readFromMatrixFile(nonMissingSitesNumMatrixFile);
	    	
	    	//seq num in the current diff matrix file is not the same with the expected seqNum
	    	if(seqNameListNonMissingSiteNumMatrixPair.getFirst().size()!=this.seqNum) {
	    		lineReader.close();
	    		throw new IllegalArgumentException("inconsistent sequence num in regional diff matrix file with the given one!!");
	    	}
	    	//
//	    	DistMatrixUtils.printMatrix(DistMatrixUtils.toStringMatrix(seqNameListLenDiffMatrixTriple.getRight()));
	    	//
	    	if(this.seqNameList.isEmpty()) {
	    		this.seqNameList.addAll(seqNameListNonMissingSiteNumMatrixPair.getFirst());
	    	}else {
	    		if(!this.seqNameList.equals(seqNameListNonMissingSiteNumMatrixPair.getFirst())) {
	    			lineReader.close();
	    			throw new IllegalArgumentException("inconsistent sequence names in regional diff matrix file!");
	    		}
	    	}
	    	
	    	this.outputDataNameRegionFilterMap.forEach((dn, filter)->{
	    		if(filter.test(region)) {
	    			//update the seq len
//	    			this.dataNameSeqLenMap.put(dn, this.dataNameSeqLenMap.get(dn)+seqNameListDiffMatrixPair.getMiddle());
	    			
	    			//update the diff matrix
	    			long[][] currentMatrix = this.dataNameSummedPairwiseNonMissingSitesNumMatrixMap.get(dn);
//	    			DistMatrixUtils.printMatrix(DistMatrixUtils.toStringMatrix(diffMatrix));
	    			
	    			///////////////
	    			for(int i=0;i<seqNameListNonMissingSiteNumMatrixPair.getSecond().length;i++) {
	    				for(int j=0;j<seqNameListNonMissingSiteNumMatrixPair.getSecond()[i].length;j++) {
	    					currentMatrix[i][j]=seqNameListNonMissingSiteNumMatrixPair.getSecond()[i][j]+currentMatrix[i][j];
	    				}
	    			}
	    			
	    			this.dataNameSummedPairwiseNonMissingSitesNumMatrixMap.put(dn, currentMatrix);
//	    			DistMatrixUtils.printMatrix(DistMatrixUtils.toStringMatrix(diffMatrix));
	    		}
	    	});
	    		
	    	//////////////////////////add the region nuc different no. matrix to the summed matrix
	    	Path diffMatrixFile = Path.of(
	    			this.regionMatrixFilesOutDir.toString(),
	    			PairwiseDiffMatrixFromVcf.buildRegionDiffNumMatrixFileName(regionIndex));
	    	
	    	log.info("read diff matrix from file "+diffMatrixFile.toString()+" for region with index="+regionIndex+" and region="+region.toString());
	    	
	    	//////////////////////////
	    	Pair<List<String>, int[][]> seqNameListDiffMatrixPair = DistMatrixUtils.readFromMatrixFile(diffMatrixFile);
	    	
	    	//seq num in the current diff matrix file is not the same with the expected seqNum
	    	if(seqNameListDiffMatrixPair.getFirst().size()!=this.seqNum) {
	    		lineReader.close();
	    		throw new IllegalArgumentException("inconsistent sequence num in regional diff matrix file with the given one!!");
	    	}
	    	//
//	    	DistMatrixUtils.printMatrix(DistMatrixUtils.toStringMatrix(seqNameListLenDiffMatrixTriple.getRight()));
	    	//
	    	if(this.seqNameList.isEmpty()) {
	    		this.seqNameList.addAll(seqNameListDiffMatrixPair.getFirst());
	    	}else {
	    		if(!this.seqNameList.equals(seqNameListDiffMatrixPair.getFirst())) {
	    			lineReader.close();
	    			throw new IllegalArgumentException("inconsistent sequence names in regional diff matrix file!");
	    		}
	    	}
	    	
	    	this.outputDataNameRegionFilterMap.forEach((dn, filter)->{
	    		if(filter.test(region)) {
	    			//update the seq len
//	    			this.dataNameSeqLenMap.put(dn, this.dataNameSeqLenMap.get(dn)+seqNameListDiffMatrixPair.getMiddle());
	    			
	    			//update the diff matrix
	    			int[][] currentMatrix = this.dataNameSummedPairwiseDiffMatrixMap.get(dn);
//	    			DistMatrixUtils.printMatrix(DistMatrixUtils.toStringMatrix(diffMatrix));
	    			
	    			///////////////
	    			for(int i=0;i<seqNameListDiffMatrixPair.getSecond().length;i++) {
	    				for(int j=0;j<seqNameListDiffMatrixPair.getSecond()[i].length;j++) {
	    					currentMatrix[i][j]=seqNameListDiffMatrixPair.getSecond()[i][j]+currentMatrix[i][j];
	    				}
	    			}
	    			
	    			this.dataNameSummedPairwiseDiffMatrixMap.put(dn, currentMatrix);
//	    			DistMatrixUtils.printMatrix(DistMatrixUtils.toStringMatrix(diffMatrix));
	    		}
	    	});
	    	
	    }
	    
	    lineReader.close();
	}
	

	/**
	 * @return the seqNameList
	 */
	public List<String> getSeqNameList() {
		return seqNameList;
	}


//	/**
//	 * @return the dataNameSeqLenMap
//	 */
//	public Map<String, Integer> getDataNameSeqLenMap() {
//		return dataNameSeqLenMap;
//	}


	/**
	 * @return the dataNameSummedPairwiseNonMissingSitesNumMatrixMap
	 */
	public Map<String, long[][]> getDataNameSummedPairwiseNonMissingSitesNumMatrixMap() {
		return dataNameSummedPairwiseNonMissingSitesNumMatrixMap;
	}

	
	/**
	 * @return the dataNameSummedPairwiseDiffMatrixMap
	 */
	public Map<String, int[][]> getDataNameSummedPairwiseDiffMatrixMap() {
		return dataNameSummedPairwiseDiffMatrixMap;
	}


	
}
