package genomics.utils;

public class CDS extends SimpleGenomicRegion{
	private final Integer phase;
	
	public CDS(String chrom, int start, int end, Strand strand, Integer phase) {
		super(chrom, start, end, strand);
		// TODO Auto-generated constructor stub
		
		this.phase=phase;
	}

	/**
	 * @return the phase
	 */
	public Integer getPhase() {
		return phase;
	}
	
}
