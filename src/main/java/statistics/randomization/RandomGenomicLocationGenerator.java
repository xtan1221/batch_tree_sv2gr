package statistics.randomization;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import basic.Bin;
import genomics.chrom.ChromLenReader;
import genomics.utils.SimpleGenomicRegion;


/**
 * generate a set of random genomic location 
 * 
 * @author tanxu
 * 
 */
public class RandomGenomicLocationGenerator {
	private final ChromLenReader chromLenReader;
	
	/**
	 * a list of genomic region sizes from which size of a random genomic region should be uniformally drawn
	 */
	private final List<Integer> genomicRegionSizes;
	
	/**
	 * number of random genomic regions to generate
	 */
	private final int replicates;
	
	
	////////////////////////
	private int totalLen;
	private List<String> orderedChroms; //store a specific ordering of chroms for lookup; the specific ordering strategy is trival
	private Map<String, Bin> chromCumulativeProportionBinMap; //map from chrom name to the cumulative proportions of total length 
	
	/////////////////////
	private Map<String, List<SimpleGenomicRegion>> chromGeneratedGenomicRegionsMap;
	private List<SimpleGenomicRegion> generatedGenomicRegions;
	
	/**
	 * 
	 * @param chromLenReader
	 * @param genomicRegionSizes
	 * @param replicates
	 */
	public RandomGenomicLocationGenerator(ChromLenReader chromLenReader, List<Integer> genomicRegionSizes, int replicates) {
		super();
		this.chromLenReader = chromLenReader;
		this.genomicRegionSizes=genomicRegionSizes;
		this.replicates = replicates;
		
		this.prepare();
		this.run();
	}
	

	void prepare() {
		this.totalLen=0;
		this.orderedChroms=new ArrayList<>();
		for(String chrom:this.chromLenReader.getChromNameLengthMap().keySet()) {
			this.totalLen+=this.chromLenReader.getChromNameLengthMap().get(chrom);
			this.orderedChroms.add(chrom);
		}
		
		///////////////////////////
		this.chromCumulativeProportionBinMap=new LinkedHashMap<>();
		double previousChromUpper=0d;
		double upperProportion=0d;
		for(String chrom:this.orderedChroms) {
			upperProportion=(double)this.chromLenReader.getChromNameLengthMap().get(chrom)/this.totalLen+previousChromUpper;
			
			this.chromCumulativeProportionBinMap.put(chrom, new Bin(previousChromUpper, true, upperProportion, true));
			previousChromUpper=upperProportion;
		}
		
		//////////////////
	}
	
	
	void run() {
		Random rand=new Random();
		this.chromGeneratedGenomicRegionsMap=new HashMap<>();
		int i=0; //number of successfully generated genomic regions
		while(i<this.replicates) {
			double value=rand.nextDouble();
			
			String chrom=null;
			Integer start=null;
			for(String c:this.orderedChroms) {
				Bin bin=this.chromCumulativeProportionBinMap.get(c);
				if(bin.inRange(value)) {
					chrom=c;
					start=(int)((value-bin.getMin())/bin.getBinSize()*this.chromLenReader.getChromNameLengthMap().get(c));
					break;
				}
			}
			
			//generate a random length
			int len=this.genomicRegionSizes.get(rand.nextInt(this.genomicRegionSizes.size()));
			//check if the region is out of the chromosome or not
			if(start+len<=this.chromLenReader.getChromNameLengthMap().get(chrom)) {//the region is fully within the chromosome
				if(!this.chromGeneratedGenomicRegionsMap.containsKey(chrom)) {
					this.chromGeneratedGenomicRegionsMap.put(chrom, new ArrayList<>());
				}
				this.chromGeneratedGenomicRegionsMap.get(chrom).add(new SimpleGenomicRegion(chrom, start, start+len));
				i++;
			}
		}
	}
	

	/**
	 * @return the chromGenomicRegionsMap
	 */
	public Map<String, List<SimpleGenomicRegion>> getChromGenomicRegionsMap() {
		return chromGeneratedGenomicRegionsMap;
	}


	/**
	 * @return the generatedGenomicRegions
	 */
	public List<SimpleGenomicRegion> getGeneratedRandomGenomicRegions() {
		if(this.generatedGenomicRegions==null) {
			this.generatedGenomicRegions=new ArrayList<>();
			
			for(String chrom:this.chromGeneratedGenomicRegionsMap.keySet()) {
				this.generatedGenomicRegions.addAll(this.chromGeneratedGenomicRegionsMap.get(chrom));
			}
		}
		return generatedGenomicRegions;
	}
	
	
}
