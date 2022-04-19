package genomics.utils;

public class Position extends SimpleGenomicRegion{
	private final String chrom;
	private final int coordinate;
	
	/**
	 * 
	 * @param chrom
	 * @param coordinate
	 */
	public Position(String chrom, int coordinate) {
		super(chrom, coordinate, coordinate, null);
		this.chrom = chrom;
		this.coordinate = coordinate;
	}
	

	/**
	 * @return the coordinate
	 */
	public int getCoordinate() {
		return coordinate;
	}


	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((chrom == null) ? 0 : chrom.hashCode());
		result = prime * result + coordinate;
		return result;
	}


	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (!(obj instanceof Position))
			return false;
		Position other = (Position) obj;
		if (chrom == null) {
			if (other.chrom != null)
				return false;
		} else if (!chrom.equals(other.chrom))
			return false;
		if (coordinate != other.coordinate)
			return false;
		return true;
	}
	
	
	
}
