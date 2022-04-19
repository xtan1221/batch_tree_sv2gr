package basic;

public class Bin {
	private final double min;
	private final boolean minInclusive;
	private final double max;
	private final boolean maxInclusive;
	
	public Bin(double min, boolean minInclusive, double max, boolean maxInclusive) {
		if(min>max) {
			throw new IllegalArgumentException("min cannot be larger than max!");
		}
		this.min = min;
		this.minInclusive = minInclusive;
		this.max = max;
		this.maxInclusive = maxInclusive;
	}



	/**
	 * 
	 * @param value
	 * @return
	 */
	public boolean inRange(double value) {
		if(minInclusive) {
			if(value<min) 
				return false;
		}else {
			if(value<=min)
				return false;
		}
		
		if(this.maxInclusive) {
			if(value>max)
				return false;
		}else {
			if(value>=max)
				return false;
		}
		
		return true;
	}

	

	/**
	 * @return the min
	 */
	public double getMin() {
		return min;
	}



	/**
	 * @return the max
	 */
	public double getMax() {
		return max;
	}

	public double getBinSize() {
		return this.max-this.min;
	}

	@Override
	public String toString() {
		StringBuilder sb=new StringBuilder();
		sb.append(this.minInclusive?"[":"(");
		sb.append(this.min);
		sb.append(",");
		sb.append(this.max);
		sb.append(this.maxInclusive?"]":")");
		
		return sb.toString();
//		return "Bin [min=" + min + ", minInclusive=" + minInclusive + ", max=" + max + ", maxInclusive=" + maxInclusive
//				+ "]";
	}
	
	public static void main(String[] args) {
		
		Bin bin=new Bin(0, false, 3, true);
		
		System.out.println(bin.inRange(4));
	}
}
