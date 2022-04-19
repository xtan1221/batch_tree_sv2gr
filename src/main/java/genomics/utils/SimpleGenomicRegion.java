package genomics.utils;

/**
 * base class for a simple genomic region not corresponding to a called locus of variant allele for a sample or population
 * 
 * see {@link SingleLocusBase} for the genomic region to a called locus of variant allele for a sample or population
 * 
 * @author tanxu
 *
 */
public class SimpleGenomicRegion extends AbstractGenomicRegion{
	/**
	 * 1-based inclusive
	 */
	private int start;
	/**
	 * 1-based inclusive
	 */
	private int end;
	
	/**
	 * 
	 * @param chrom
	 * @param start
	 * @param end
	 * @param strand
	 */
	public SimpleGenomicRegion(String chrom, int start, int end, Strand strand) {
		super(chrom, strand);
		this.start = start;
		this.end = end;
	}
	
	/**
	 * 
	 * @param chrom
	 * @param start
	 * @param end
	 */
	public SimpleGenomicRegion(String chrom, int start, int end) {
		super(chrom, null);
		this.start = start;
		this.end = end;
	}
	
	
	/**
	 * @param start the start to set
	 */
	public void setStart(int start) {
		this.start = start;
	}

	/**
	 * @param end the end to set
	 */
	public void setEnd(int end) {
		this.end = end;
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


	
	//////////////////////////////////////////
	@Override
	public String toString() {
		return "SimpleGenomicRegion [start=" + start + ", end=" + end + ", strand=" + strand + ", getChrom()="
				+ getChrom() + "]";
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = super.hashCode();
		result = prime * result + end;
		result = prime * result + start;
		return result;
	}
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (!super.equals(obj))
			return false;
		if (!(obj instanceof SimpleGenomicRegion))
			return false;
		SimpleGenomicRegion other = (SimpleGenomicRegion) obj;
		if (end != other.end)
			return false;
		if (start != other.start)
			return false;
		return true;
	}
	
}
