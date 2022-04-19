package genomics.utils;

/**
 * abstract class for a genomic region
 * 
 * @author tanxu
 *
 */
public abstract class AbstractGenomicRegion {
	/**
	 * chrom/contig name
	 */
	private String chrom;
	
	/**
	 * strand of the chrom of this {@link AbstractGenomicRegion}
	 */
	protected Strand strand;
	
	/**
	 * 
	 */
	public AbstractGenomicRegion() {
		//
	}
	
	/**
	 * 
	 * @param chrom
	 * @param strand
	 */
	public AbstractGenomicRegion(String chrom, Strand strand) {
		super();
		this.chrom = chrom;
		this.strand = strand;
	}
	
	/**
	 * return whether this {@link AbstractGenomicRegion} covers the given position
	 * @param pos
	 * @return
	 */
	public boolean covers(int pos) {
		return this.getStart()<=pos && this.getEnd()>=pos;
	}
	
	/**
	 * build and return a data line of bed file format of this {@link AbstractGenomicRegion};
	 * 
	 * for subclass with additional information, this method should be overridden
	 * 
	 * @return
	 */
	public String toBedLine() {
		StringBuilder sb =new StringBuilder();
		sb.append(this.chrom).append("\t")
		.append(this.getStart()-1).append("\t") //start of bed file!!!!
		.append(this.getEnd()).append("\t")
		.append(".").append("\t") //name
		.append(".").append("\t") //score
		.append(this.strand);
		
		return sb.toString();
	}
	
	/**
	 * 
	 * @return
	 */
	public Position getStartPosition() {
		return new Position(this.getChrom(), this.getStart());
	}
	/**
	 * 
	 * @return
	 */
	public Position getEndPosition() {
		return new Position(this.getChrom(), this.getEnd());
	}
	
	//////////////////////
	/**
	 * return the start position (1-based inclusive)
	 * @return the start
	 */
	public abstract int getStart();

	/**
	 * return the end position (1-based inclusive)
	 * 
	 * @return the end
	 */
	public abstract int getEnd();
	
	///////////////////////////////
	/**
	 * @param chrom the chrom to set
	 */
	public void setChrom(String chrom) {
		this.chrom = chrom;
	}
	/**
	 * @param strand the strand to set
	 */
	public void setStrand(Strand strand) {
		this.strand = strand;
	}
	/**
	 * @return the chrom
	 */
	public String getChrom() {
		return chrom;
	}
	
	/**
	 * @return the strand
	 */
	public Strand getStrand() {
		return strand;
	}
	
	/**
	 * return the length
	 * @return
	 */
	public int getLen() {
		return this.getEnd()-this.getStart()+1;
	}
	
	/**
	 * return the center position
	 * @return
	 */
	public int getCenter() {
		return (this.getStart()+this.getEnd())/2;
	}

	

	
	//////////////////////////////////////

	@Override
	public String toString() {
		return "AbstractGenomicRegion [chrom=" + chrom + ", strand=" + strand + "]";
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((chrom == null) ? 0 : chrom.hashCode());
		result = prime * result + ((strand == null) ? 0 : strand.hashCode());
		return result;
	}
	
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (!(obj instanceof AbstractGenomicRegion))
			return false;
		AbstractGenomicRegion other = (AbstractGenomicRegion) obj;
		if (chrom == null) {
			if (other.chrom != null)
				return false;
		} else if (!chrom.equals(other.chrom))
			return false;
		if (strand != other.strand)
			return false;
		return true;
	}
	
}
