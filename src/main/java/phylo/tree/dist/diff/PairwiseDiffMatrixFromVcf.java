package phylo.tree.dist.diff;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Stream;

import basic.Pair;
import htsjdk.samtools.util.Log;
import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import phylo.batch.BcfUtils;
import phylo.tree.dist.DistMatrixUtils;
import population.vcf.utils.GenotypeRecoder;
import population.vcf.utils.GenotypeRecoderFactory;
import population.vcf.utils.VariantContextFilter;

/**
 * collector for pairwise nucleotide no. difference from all sites in a vcf file as well as the total seq length (based on the VariantContextFilter)
 * 
 * for a pair of samples, a different nucleotide of a specific site is counted if both samples have non-missing alleles that are different from each other
 * 		in cases one or both samples have missing data (./. genotype), the site is not counted in the nucleotide no. difference 
 * 
 * ==================================
 * Distinguish between the difference matrix, p-distance matrix and j-c/kimura-2/etc model based distance matrix
 * 		difference matrix contains the pairwise nucleotide no. difference
 * 		p-distance matrix contains the pairwise nucleotide proportion difference
 * 	see https://www.megasoftware.net/mega1_manual/Distance.html
 * 
 * @author tanxu
 * 
 */
public class PairwiseDiffMatrixFromVcf implements Callable<Integer>{
	public static Log log=Log.getInstance(PairwiseDiffMatrixFromVcf.class);
	
	////////////////////////
	/**
	 * matrix file contains number of different nucleotides of each pair of samples
	 */
	public static String diffNucNumMatrixFileSuffix="_diff_nuc_num_Matrix.txt";
	/**
	 * matrix file contains the number of sites with both samples with non-missing data of each pair of samples
	 */
	public static String nonMissingSiteNumMatrixFileSuffix="_non_missing_sites_num_Matrix.txt";
	
	///////////////////////////////////////////////////
	/**
	 * 
	 */
	private final List<String> orderedSampleNameList;
	/**
	 * 
	 */
	private final Path bcfFile;
	
	
	/**
	 * a list of disjoint regions SNP within which are to be used to build regional tree;
	 * see -r of bcftools
	 */
	private final String regionListString;
	
	/**
	 * a unique string identifier of the regionList that distinguish it from other BcfRegionalSNP2SingleTree running in the same batch;
	 * help to create temporary folder and name output files
	 */
	private final String uniqueRegionIdentifier;
	
	private final VariantContextFilter variantFilter;
	
	private final GenotypeRecoder genotypeRecoder;
	
	/**
	 * directory where to write all the matrix files of all regions in the same batch
	 */
	private final Path regionMatrixOutDir;
	
	/**
	 * the directory where to put the temporary files and folders
	 */
	private final Path tmpOutDir;
	
	
	//////////////////
	private Path queriedVcfFile;
	
	
	/////////////////////////////
	/**
	 * the matrix for pairwise nucleotide sites num with non-missing data
	 */
	protected int [][] nonMissingSitesNumMatrix;
//	protected int totalSites;
	
	/**
	 * the matrix for pairwise sequence difference;
	 */
	private int[][] diffMatrix;
	
	
	/**
	 * path of the file to write the matrix of pairwise total sites number with both samples with non-missing data
	 */
	private Path outNonMissingSiteNumMatrixFile;
	
	
	/**
	 * path of the file to write the pairwise nucleotide diff no. Matrix
	 */
	private Path outDiffMatrixFile;
	

	
	/**
	 * 
	 * @param reader
	 */
	public PairwiseDiffMatrixFromVcf(
			List<String> orderedSampleNameList, Path bcfFile, String regionListString, String uniqueRegionIdentifier, 
			VariantContextFilter variantFilter, GenotypeRecoder genotypeRecoder,
			Path regionMatrixOutDir, Path tmpOutDir
			) {
		// TODO Auto-generated constructor stub
		
		
		
		
		
		this.orderedSampleNameList=orderedSampleNameList;
		this.bcfFile=bcfFile;
		this.regionListString=regionListString;
		this.uniqueRegionIdentifier=uniqueRegionIdentifier;
		this.variantFilter=variantFilter;
		this.genotypeRecoder=genotypeRecoder;
		this.tmpOutDir=tmpOutDir;
		this.regionMatrixOutDir=regionMatrixOutDir;
	}
	
	
	/**
	 */
	@Override
	public Integer call() throws Exception {
		
		if(this.queryVCF()!=0) {
			return 1;
		}
		
		if(this.parseVcfFile()!=0) {
			return 2;
		}
		
//		this.update();
		
		if(this.writeToFile()!=0) {
			return 3;
		}
		
		return 0;
	}
	
	/**
	 * return 0 if succeed; return 1 otherwise
	 */
	private int queryVCF() {
		log.info("start querying regional VCF from bcf file"+" for region id=="+this.uniqueRegionIdentifier);
		this.queriedVcfFile = Path.of(this.tmpOutDir.toString(), this.uniqueRegionIdentifier.concat(".vcf"));
		
		BcfUtils.queryRegionIntoVcfFileFromIndexedBcfFile(this.bcfFile.toString(), this.regionListString, this.queriedVcfFile.toString());
		
		if(this.queriedVcfFile.toFile().exists()) {
			log.info("querying regional VCF from bcf file is successfully finished "+" for region id=="+this.uniqueRegionIdentifier);
			return 0;
		}else {
			log.error("queried regional vcf file is not found"+" for region id=="+this.uniqueRegionIdentifier);
			return 1;
		}
	}
	
	/**
	 * return 0 if succeed, 1 otherwise
	 * @return
	 */
	private int parseVcfFile() {
		log.info("start parsing queried vcf file"+" for region id=="+this.uniqueRegionIdentifier);
		///
		this.diffMatrix=new int[this.orderedSampleNameList.size()][this.orderedSampleNameList.size()];
		this.nonMissingSitesNumMatrix=new int[this.orderedSampleNameList.size()][this.orderedSampleNameList.size()];
		
		for(int i=0;i<this.orderedSampleNameList.size();i++) {
			for(int j=0;j<this.orderedSampleNameList.size();j++) {
				this.diffMatrix[i][j]=0;
				this.nonMissingSitesNumMatrix[i][j]=0;
			}
		}
		
		try {
		////
			VCFFileReader reader = new VCFFileReader(this.queriedVcfFile, false);
			
			Stream<VariantContext> stream=reader.iterator().stream();
			
			stream=stream.filter(this.variantFilter); //pass the given filter
			
			stream.forEach(v->{
//				if(v.getStart()%100000==0)
//					log.info("vcf file progress:"+v.getContig()+" "+v.getStart());
				
				
				//
				List<String> sampleSeqList = new ArrayList<>();
				
				for(String sample:this.orderedSampleNameList) {
					Genotype g = v.getGenotype(sample);
					
					sampleSeqList.add(this.genotypeRecoder.apply(new Pair<>(v,g)));
				}
				
				for(int i=0;i<this.orderedSampleNameList.size();i++) {
					for(int j=0;j<this.orderedSampleNameList.size();j++) {
						//
						if(sampleSeqList.get(i).equals(GenotypeRecoderFactory.MISSING_BASE) //missing data
								||sampleSeqList.get(j).equals(GenotypeRecoderFactory.MISSING_BASE)  //missing data
								) { 
							//skip sites with one or two samples of the pair with missing data
						}else {
							//both samples of the sample pair has non-missing data
//							this.nonMissingSitesNumMatrix[i][j]=this.nonMissingSitesNumMatrix[i][j]+v.getLengthOnReference();
							this.nonMissingSitesNumMatrix[i][j]=this.nonMissingSitesNumMatrix[i][j]+1;
							
							if(sampleSeqList.get(i).equals(sampleSeqList.get(j))) {
								//
							}else {
								this.diffMatrix[i][j] = this.diffMatrix[i][j]+1; //TODO this assume there is only one nucleotide 
							}
						}
					}
				}
			});
			
			
			reader.close();
		}catch(TribbleException e) {
			log.error("Exception is thrown when parsing queried vcf file for region id=="+this.uniqueRegionIdentifier+" :" +e.getMessage());
			return 1;
		}
		
		
		log.info("parsing queried vcf file"+" for region id=="+this.uniqueRegionIdentifier+" is successfully done");
		return 0;
	}
	
	/**
	 * build and return the file name (not including the path) of matrix file containing the pairwise total non-missing sites for the region with the given identifier
	 * @param uniqueRegionIdentifier
	 * @return
	 */
	public static String buildNonMissingSiteNumMatrixFileName(String uniqueRegionIdentifier) {
		return uniqueRegionIdentifier.concat(nonMissingSiteNumMatrixFileSuffix);
	}
	
	/**
	 * build and return the file name (not including the path) of matrix file containing the pairwise diff nucleotide num for the region with the given identifier
	 * @param uniqueRegionIdentifier
	 * @return
	 */
	public static String buildRegionDiffNumMatrixFileName(String uniqueRegionIdentifier) {
		return uniqueRegionIdentifier.concat(diffNucNumMatrixFileSuffix);
	}
	
	private int writeToFile() {
		log.info("start write to output matrix file for region id=="+this.uniqueRegionIdentifier);
		
		this.outNonMissingSiteNumMatrixFile=Path.of(this.regionMatrixOutDir.toString(),buildNonMissingSiteNumMatrixFileName(uniqueRegionIdentifier));
		
		this.outDiffMatrixFile = Path.of(this.regionMatrixOutDir.toString(),buildRegionDiffNumMatrixFileName(this.uniqueRegionIdentifier));
		
		if(this.outDiffMatrixFile.toFile().exists()) {
			this.outDiffMatrixFile.toFile().delete();
		}
		
		try {
			//write to the 
			DistMatrixUtils.writeMatrixToFile(
					this.outNonMissingSiteNumMatrixFile, 
					this.orderedSampleNameList, 
//					this.totalSites,
					this.nonMissingSitesNumMatrix);
			
			DistMatrixUtils.writeMatrixToFile(
					this.outDiffMatrixFile, 
					this.orderedSampleNameList, 
//					this.totalSites, 
					this.diffMatrix);
			
		} catch (IOException e) {
			log.error("Exception is thrown when writing to output matrix file for region id=="+this.uniqueRegionIdentifier+"  :"+e.getMessage());
			return 1;
		}
		
		log.info("writing to output matrix file for region id=="+this.uniqueRegionIdentifier+" is successfully done!");
		return 0;
	}
}
