package population.sv.analysis.adjacentRegion.polymophismSiteProportion;

public class Window {
	private final String chrom;
	private final int start;
	private final int end;
	
	/**
	 * 
	 */
	private int totalSiteNum;
	private int totalPolymophSiteNum;
	
	public Window(String chrom, int start, int end) {
		super();
		this.chrom = chrom;
		this.start = start;
		this.end = end;
		
		this.totalPolymophSiteNum=0;
		this.totalSiteNum=0;
	}

	public void addOneToTotalSite() {
		this.totalSiteNum++;
	}
	public void addOneToPolymorphSite() {
		this.totalPolymophSiteNum++;
	}
	
	/**
	 * @return the chrom
	 */
	public String getChrom() {
		return chrom;
	}


	/**
	 * @return the start
	 */
	public int getStart() {
		return start;
	}


	/**
	 * @return the end
	 */
	public int getEnd() {
		return end;
	}
	
	/**
	 * @return the totalSiteNum
	 */
	public int getTotalSiteNum() {
		return totalSiteNum;
	}

	/**
	 * @return the totalPolymophSiteNum
	 */
	public int getTotalPolymophSiteNum() {
		return totalPolymophSiteNum;
	}

	
	@Override
	public String toString() {
		return "Window [chrom=" + chrom + ", start=" + start + ", end=" + end + ", totalSiteNum=" + totalSiteNum
				+ ", totalPolymophSiteNum=" + totalPolymophSiteNum + "]";
	}

	
	///////////////////////////////////
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((chrom == null) ? 0 : chrom.hashCode());
		result = prime * result + end;
		result = prime * result + start;
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (!(obj instanceof Window))
			return false;
		Window other = (Window) obj;
		if (chrom == null) {
			if (other.chrom != null)
				return false;
		} else if (!chrom.equals(other.chrom))
			return false;
		if (end != other.end)
			return false;
		if (start != other.start)
			return false;
		return true;
	}
}
