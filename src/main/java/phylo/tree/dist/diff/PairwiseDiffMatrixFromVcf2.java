package phylo.tree.dist.diff;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Stream;

import basic.Pair;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import phylo.batch.BcfUtils;
import phylo.tree.dist.DistMatrixUtils;
import population.vcf.utils.GenotypeRecoder;
import population.vcf.utils.VariantContextFilter;

/**
 * collector for pairwise difference that can be used to calculate the distance based on Jukes-Cantor model
 * 
 * specifically, collect the number of positions at which two sequences differ in additional to the total number of sites(nucleotides) or sequence length
 * 
 * 
 * @author tanxu
 * 
 */
public class PairwiseDiffMatrixFromVcf2 implements Callable<Integer>{
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
	
	private final Path regionDiffMatrixOutDir;
	/**
	 * the directory where to put the temporary files and folders
	 */
	private final Path tmpOutDir;
	
	/**
	 * 
	 */
//	private final JCModelBasedDiffMatrixMerger JCModelBasedDiffMatrixMerger;
	//////////////////
	private Log log=Log.getInstance(PairwiseDiffMatrixFromVcf2.class);

	//////////////////
	private Path queriedVcfFile;
	
	
	/////////////////////////////
	
	protected int totalSites;
	
	/**
	 * the matrix for pairwise sequence difference;
	 * 
	 */
	private int[][] diffMatrix;
	
	/**
	 * path of the file to output the totalSites and diffMatrix
	 */
	private Path outDiffMatrixFile;
	
	/**
	 * 
	 * @param reader
	 */
	public PairwiseDiffMatrixFromVcf2(
			List<String> orderedSampleNameList, Path bcfFile, String regionListString, String uniqueRegionIdentifier, 
			VariantContextFilter variantFilter, GenotypeRecoder genotypeRecoder,
			Path regionDiffMatrixOutDir, Path tmpOutDir
			) {
		
		// TODO Auto-generated constructor stub
		
		this.orderedSampleNameList=orderedSampleNameList;
		this.bcfFile=bcfFile;
		this.regionListString=regionListString;
		this.uniqueRegionIdentifier=uniqueRegionIdentifier;
		this.variantFilter=variantFilter;
		this.genotypeRecoder=genotypeRecoder;
		this.tmpOutDir=tmpOutDir;
		this.regionDiffMatrixOutDir=regionDiffMatrixOutDir;
//		this.JCModelBasedDiffMatrixMerger=JCModelBasedDiffMatrixMerger;
	}
	
	
	/**
	 */
	@Override
	public Integer call() throws Exception {
		
//		if(this.queryVCF()!=0) {
//			return 1;
//		}
		
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
		log.info("start querying regional VCF from bcf file"+" for tree id=="+this.uniqueRegionIdentifier);
		this.queriedVcfFile = Path.of(this.tmpOutDir.toString(), this.uniqueRegionIdentifier.concat(".vcf"));
		
		BcfUtils.queryRegionIntoVcfFileFromIndexedBcfFile(this.bcfFile.toString(), this.regionListString, this.queriedVcfFile.toString());
		
		if(this.queriedVcfFile.toFile().exists()) {
			log.info("querying regional VCF from bcf file is successfully finished "+" for tree id=="+this.uniqueRegionIdentifier);
			return 0;
		}else {
			log.error("queried regional vcf file is not found"+" for tree id=="+this.uniqueRegionIdentifier);
			return 1;
		}
	}
	
	/**
	 * return 0 if succeed, 1 otherwise
	 * @return
	 */
	private int parseVcfFile() {
		log.info("start parsing queried vcf file"+" for tree id=="+this.uniqueRegionIdentifier);
		
		///
		this.diffMatrix=new int[this.orderedSampleNameList.size()][this.orderedSampleNameList.size()];
		for(int i=0;i<this.orderedSampleNameList.size();i++) {
			for(int j=0;j<this.orderedSampleNameList.size();j++) {
				diffMatrix[i][j]=0;
			}
		}
		////
		VCFFileReader reader = new VCFFileReader(this.bcfFile, false);
		
		Stream<VariantContext> stream=reader.iterator().stream();
		
		stream=stream.filter(this.variantFilter); //pass the given filter
		
		stream.forEach(v->{
			if(v.getStart()%100000==0)
				log.info("vcf file progress:"+v.getContig()+" "+v.getStart());
			
			//TODO
			this.totalSites = this.totalSites+v.getLengthOnReference();
			
			//
			List<String> sampleSeqList = new ArrayList<>();
			
			for(String sample:this.orderedSampleNameList) {
				Genotype g = v.getGenotype(sample);
				
				sampleSeqList.add(this.genotypeRecoder.apply(new Pair<>(v,g)));
			}
			
			for(int i=0;i<this.orderedSampleNameList.size();i++) {
				for(int j=0;j<this.orderedSampleNameList.size();j++) {
					if(sampleSeqList.get(i).equals(sampleSeqList.get(j))) {
						//
					}else {
						this.diffMatrix[i][j] = this.diffMatrix[i][j]+1; //TODO this assume there is only one nucleotide 
					}
				}
			}
		});
		
		reader.close();
		
		return 0;
	}
	
	
//	private int update() {
//		this.JCModelBasedDiffMatrixMerger.merge(totalSites, diffMatrix);
//		
//		return 0;
//	}
//	
	
	
	private int writeToFile() {
		this.outDiffMatrixFile = Path.of(this.regionDiffMatrixOutDir.toString(),this.uniqueRegionIdentifier.concat("_totalSites_diffMatrix.txt"));
		
		if(this.outDiffMatrixFile.toFile().exists()) {
			this.outDiffMatrixFile.toFile().delete();
		}
		try {
			DistMatrixUtils.writeToPhylipDistMatrixFile(
					this.outDiffMatrixFile, 
					this.orderedSampleNameList, 
//					this.totalSites, 
					DistMatrixUtils.toStringMatrix(
							DistMatrixUtils.calcualteJCDistMatrixFromPDistMatrix(
									DistMatrixUtils.calculatePDistaFromDiffMatrix(totalSites, diffMatrix))));
//			PhylipDistMatrixUtils.writeDiffMatrixAndSeqLenToFile(
//					this.outDiffMatrixFile, 
//					this.orderedSampleNameList, 
//					this.totalSites, 
//					this.diffMatrix);
		} catch (IOException e) {
			e.printStackTrace();
			return 1;
		}
		return 0;
	}
}
