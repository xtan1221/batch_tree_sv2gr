package population.sv.analysis.assignSiteToClosestSV.postprocess;

import java.util.ArrayList;
import java.util.List;

import basic.NumericRange;
import population.sv.utils.SimpleSVType;

public class SVAdjacentNucSiteNumericRangeUtils {
	public static final int MAX_SV_SIZE=20000;
	
	/**
	 * build and return the list of all {@link SimpleSVType}s
	 * 
	 * note that only DEL, INV, DUP types were included in the analysis
	 * @return
	 */
	public static List<SimpleSVType> getAllSVTypes(){
		List<SimpleSVType> ret = new ArrayList<>();
		ret.add(SimpleSVType.DEL);
		ret.add(SimpleSVType.INV);
		ret.add(SimpleSVType.DUP);
		ret.add(SimpleSVType.INS);
		return ret;
	}
	
	/////////////////////////////////
	/**
	 * return the {@link NumericRange} for all included SVs
	 * @return
	 */
	public static NumericRange getSVFullSizeRanges(){
		return new NumericRange(1, MAX_SV_SIZE, false);
	}
	
	
	/**
	 * 
	 * @return
	 */
	public static List<NumericRange> getSVSizeRanges1To1000WithSize100(){
		List<NumericRange> ret = new ArrayList<>();
		ret.addAll(NumericRange.buildBatchRanges(1, 100, 10, false));
		return ret;
	}
	
	/**
	 * 
	 * @return
	 */
	public static List<NumericRange> getSVSizeRanges1To10000WithSize1000(){
		List<NumericRange> ret = new ArrayList<>();
		ret.addAll(NumericRange.buildBatchRanges(1, 1000, 10, false));
		return ret;
	}
	
	///////////////////////////////////
	
	/**
	 * return the {@link NumericRange} for all SVs regardless of the value of derived allele prop (either null or not)
	 * @return
	 */
	public static NumericRange getSVFullDerivedAllelePropRange(){
		return new NumericRange(0d, 1d, true); 
	}
	
	
	/**
	 * return the list of numeric ranges for derived allele prop; only SVs with non-null DAP will be included
	 * @return
	 */
	public static List<NumericRange> getSVDerivedAllelePropRanges(){
		List<NumericRange> ret = new ArrayList<>();
		ret.add(new NumericRange(0d, 0.1, false)); //
		ret.add(new NumericRange(0.1, 0.5, false)); //
		ret.add(new NumericRange(0.5, 0.9, false)); //
		ret.add(new NumericRange(0.9, 1d, false)); //
		return ret;
	}
	
	
	///////////////////////////////////////
	/**
	 * 
	 * @param windSize
	 * @param windNum
	 * @return
	 */
	public static List<NumericRange> getOrderedDistToClosestSVRanges(int windSize, int windNum){
		List<NumericRange> ret = new ArrayList<>();
		ret.addAll(NumericRange.buildBatchRanges(1, windSize, windNum, false));
		return ret;
	}
}
