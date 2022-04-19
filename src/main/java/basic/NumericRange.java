package basic;

import java.util.ArrayList;
import java.util.List;

/**
 * a numeric range that can be used to check if a given numeric value is included by this one
 * 
 * @author tanxu
 *
 */
public class NumericRange{
	
	/**
	 * min value; if null, skip checking
	 * inclusive
	 */
	private final Double min;//inclusive
	/**
	 * max value; if null, skip checking
	 * inclusive
	 */
	private final Double max; //inclusive
	/**
	 * whether or not null value should be considered as covered by this range
	 */
	private final boolean toIncludeNullValue;
	
	
	public NumericRange(Double min, Double max, boolean toIncludeNullValue) {
		super();
		if(min==null&&max==null)
			throw new IllegalArgumentException("given min and max cannot both be null!");
		
		this.min = min;
		this.max = max;
		this.toIncludeNullValue=toIncludeNullValue;
	}
	
	public NumericRange(Integer min, Integer max, boolean toIncludeNullValue) {
		super();
		if(min==null&&max==null)
			throw new IllegalArgumentException("given min and max cannot both be null!");
		
		this.min = min==null?null:(double)min;
		this.max = max==null?null:(double)max;
		
		this.toIncludeNullValue=toIncludeNullValue;
	}
	
	/**
	 * 
	 * @param value
	 * @return
	 */
	public boolean covers(Double value) {
		if(value==null) {
			return this.toIncludeNullValue;
		}
		
		if(this.min==null)
			return value <=max;
		if(this.max==null)
			return value>=this.min;
		
		return value>=min && value <=max;
	}
	
	/**
	 * 
	 * @param value
	 * @return
	 */
	public boolean covers(Integer value) {
		if(value==null) {
			return this.toIncludeNullValue;
		}
		
		if(this.min==null)
			return value <=max;
		if(this.max==null)
			return value>=this.min;
		
//		boolean ret=value>=min && value <=max;
		return value>=min && value <=max;
	}

	/**
	 * @return the min
	 */
	public final Double getMin() {
		return min;
	}
	
	/**
	 * @return the max
	 */
	public final Double getMax() {
		return max;
	}

	/**
	 * @return the toIncludeNullValue
	 */
	public final boolean isToIncludeNullValue() {
		return toIncludeNullValue;
	}
	
	
	/**
	 * build and return a list of {@link NumericRange}s
	 * @param type
	 * @param gap
	 * @param num
	 * @param toIncludeNullValue
	 * @return
	 */
	public static List<NumericRange> buildBatchRanges(int firstStart, int gap, int num, boolean toIncludeNullValue) {
		List<NumericRange> ret = new ArrayList<>();
		for(int i=0;i<num;i++) {
			int start=i*gap+firstStart;
			int end=start+gap-1;
			ret.add(new NumericRange(start,end, toIncludeNullValue));
		}
		return ret;
	}

	@Override
	public String toString() {
		return "NumericRange [min=" + min + ", max=" + max + ", toIncludeNullValue=" + toIncludeNullValue + "]";
	}
	
	
}
