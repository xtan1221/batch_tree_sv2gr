package population.sv.utils;

import java.util.Map;

import genomics.utils.Strand;
import htsjdk.variant.variantcontext.VariantContext;
import population.utils.Genotype;
import population.utils.SingleLocusBase;

/**
 * a locus in a population with one or more samples containing the SV allele
 * 
 * @author tanxu
 * 
 */
public class SimpleSVLocus extends SingleLocusBase{
	/**
	 * 
	 */
	private final SimpleSVType type;
	/**
	 * 
	 */
	private final int start;
	/**
	 * 
	 */
	private final int end;
	
	/////////////////////////////////////
	/**
	 * note that the first sample has index 1 rather than 0!!!!!!!!!!!!!
	 * ??
	 */
	private Map<Integer, Genotype> sampleIndexGenotypeMap;
	/**
	 * 
	 */
	private VariantContext vc;
	
	/**
	 * 
	 * @param type
	 * @param chrom
	 * @param start
	 * @param end
	 * @param sampleIndexGenotypeMap sample index starting from 1
	 */
	public SimpleSVLocus(SimpleSVType type,String chrom, int start, int end, Map<Integer, Genotype> sampleIndexGenotypeMap) {
		super(chrom);
		this.type = type;
		this.start = start;
		this.end = end;
		this.sampleIndexGenotypeMap = sampleIndexGenotypeMap;
	}
	
	public void setVariantContext(VariantContext vc) {
		this.vc=vc;
	}
	
	public VariantContext getVariantContext() {
		return this.vc;
	}
	
	/**
	 * @return the type
	 */
	public SimpleSVType getType() {
		return type;
	}

	/**
	 * @return the start
	 * 
	 */
	@Override
	public int getStart() {
		return start;
	}
	
	/**
	 * @return the end
	 */
	@Override
	public int getEnd() {
		return end;
	}
	
	/**
	 * return the size of the SV, equal to the {@link #getLen()}
	 * @return
	 */
	public int getSize() {
		return this.getLen();
	}
	
	/**
	 * recode the genotype of the given sample index with the given {@link Genotype}
	 * @param index
	 * @param gt
	 */
	public void recodeSampleGenotype(int index, Genotype gt) {
		this.sampleIndexGenotypeMap.put(index, gt);
	}
	
	////////////////////////////////////////
	/**
	 * return the map from sample index to the genotype;
	 * note that that sample index starts from 1 rather than 0
	 */
	@Override
	public Map<Integer, Genotype> getSampleIndexGenotypeMap() {
		return this.sampleIndexGenotypeMap;
	}
	
	
	
	/////////////////////////////////unused methods
	/**
	 * strand is irrelevant for SV locus
	 */
	@Override
	public Strand getStrand() {
		throw new UnsupportedOperationException("strand of sv is not supported!");
	}
	
	
	//////////////////////////////////////////////
	@Override
	public String toString() {
		return "SimpleSVLocus [type=" + type + ", start=" + start + ", end=" + end + ", getChrom()="
				+ getChrom() + "]";
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = super.hashCode();
		result = prime * result + end;
		result = prime * result + start;
		result = prime * result + ((type == null) ? 0 : type.hashCode());
		return result;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (!super.equals(obj))
			return false;
		if (!(obj instanceof SimpleSVLocus))
			return false;
		SimpleSVLocus other = (SimpleSVLocus) obj;
		if (end != other.end)
			return false;
		if (start != other.start)
			return false;
		if (type != other.type)
			return false;
		return true;
	}
	
	
}
