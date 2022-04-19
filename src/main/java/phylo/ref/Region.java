package phylo.ref;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * a consecutive region on a chrom
 * 1-based coordinate
 * 
 * 
 * consistent with -r parameter of bcftools
 * 		chr:beg-end
 * 
 * 
 * @author tanxu
 * 
 */
public class Region implements Comparable<Region>{
	private final String chrom;
	private final int start;
	private final int end;
	
	public Region(String chrom, int start, int end){
		if(chrom==null||chrom.isEmpty()) {
			throw new IllegalArgumentException("given chrom cannot be null or empty!");
		}
		
		if(start<=0)
			throw new IllegalArgumentException("start must be positive integer!");
		
		if(end<=0)
			throw new IllegalArgumentException("end must be positive integer!");
		
		if(start>end) {
			throw new IllegalArgumentException("given start cannot be larger than end!");
		}
		
		this.chrom=chrom;
		this.start=start;
		this.end=end;
	}
	
	@Override
	public int compareTo(Region o) {
		if(this.chrom.equals(o.getReferenceName())) {
			return start-o.getStart();
		}else {
			return this.chrom.compareTo(o.chrom);
		}
	}
	
	/**
	 * @return the chrom
	 */
	public String getReferenceName() {
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
	 * return the string representation of this Region consistent with the -r parameter of bcftools
	 * 		chr:beg-end
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(this.chrom).append(":").append(this.start).append("-").append(this.end);
		return sb.toString();
	}
	
	/**
	 * 
	 * @param regionString
	 * @return
	 */
	public static Region fromString(String regionString) {
		String[] splits = regionString.split(":");
		String[] posSplits = splits[1].split("-");
		
		return new Region(splits[0], Integer.parseInt(posSplits[0]), Integer.parseInt(posSplits[1]));
	}
	
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
		if (!(obj instanceof Region))
			return false;
		Region other = (Region) obj;
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
	
	public static void main(String[] args) {
		Region r1 = new Region("Chr03", 4444, 222222);
		Region r2 = new Region("Chr01", 345, 222222);
		
		Region r3 = new Region("Chr01", 1, 1000);
		
		List<Region> list= new ArrayList<>();
		list.add(r1);
		list.add(r2);
		list.add(r3);
		
		
		Collections.sort(list);
		
		System.out.println(list);
	}

	
	
}
