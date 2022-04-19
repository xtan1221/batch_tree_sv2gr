package genomics.utils;


public class DomainGenomicRegion extends SimpleGenomicRegion{
	private final String domain;

	public DomainGenomicRegion(String chrom, int start, int end, Strand strand, String domain) {
		super(chrom, start, end, strand);
		this.domain = domain;
	}

	/**
	 * @return the domain
	 */
	public String getDomain() {
		return domain;
	}
	
	
	/**
	 * return a new DomainGenomicRegion that contains the region of the given start and end pos;
	 * if the given region is not overlapping with this DomainGenomicRegion, return null;
	 * @param start
	 * @param end
	 * @return
	 */
	public DomainGenomicRegion getOverlappedRegion(int start, int end) {
		if(start>end) {
			throw new IllegalArgumentException("given start cannot be larger than end!");
		}
		if(end<this.getStart() || this.getEnd()<start) {//no overlapping
			return null;
		}
		
		return new DomainGenomicRegion(
				this.getChrom(),
				this.getStart()<start?start:this.getStart(),
				this.getEnd()>end?end:this.getEnd(),
				this.getStrand(),
				this.getDomain()
				);
		
	}
	
}
