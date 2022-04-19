package population.sv.analysis.adjacentRegion.postprocess;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import basic.NumericRange;
import population.sv.analysis.adjacentRegion.postprocess.CollectedWindowedDataFileReader.SVRecord;
import population.sv.utils.SimpleSVType;

/**
 * 
 * @author tanxu
 * 
 */
public class GroupStatSummariesBuilder {
	/**
	 * all SVs of each type
	 */
	private final Map<SimpleSVType, List<SVRecord>> svTypeSVRecordsMap;
	
	/**
	 * a list of size range for each of which to group the SVs
	 * 
	 * each {@link SVRecord} can only be assigned to at most one of these size ranges
	 * 
	 * cannot be null or empty;
	 */
	private final List<NumericRange> sizeRanges;
	
	/**
	 * a list of ranges for derived allele proportion for each of which to group the SVs
	 * 
	 * each {@link SVRecord} can only be assigned to at most one of these DAP ranges
	 * 
	 * cannot be null or empty;
	 */
	private final List<NumericRange> derivedAlleleProportionRanges;
	
	/**
	 * 
	 */
	private final int windowNum;
	
	///////////////////////
	private Map<SimpleSVType, Map<NumericRange, Map<NumericRange,List<SVRecord>>>> svTypeSizeRangeDAPRangeSVRecordsMapMapMap;
	/**
	 * map from the sv type to the 
	 * 		map from sv size range to the 
	 * 			map from derived allele proportion range to the {@link GroupStatSummary}
	 */
	private Map<SimpleSVType, Map<NumericRange, Map<NumericRange,GroupStatSummary>>> svTypeSizeRangeDAPRangeGroupStatSummaryMapMapMap;
	
	/**
	 * 
	 * @param svTypeSVRecordsMap
	 * @param sizeRanges
	 * @param derivedAlleleProportionRanges
	 * @param windowNum
	 */
	public GroupStatSummariesBuilder(
			Map<SimpleSVType, List<SVRecord>> svTypeSVRecordsMap,
			List<NumericRange> sizeRanges, 
			List<NumericRange> derivedAlleleProportionRanges,
			int windowNum) {
		super();
		this.svTypeSVRecordsMap = svTypeSVRecordsMap;
		this.sizeRanges = sizeRanges;
		this.derivedAlleleProportionRanges=derivedAlleleProportionRanges;
		this.windowNum = windowNum;
		
		
		this.build();
	}


	
	void build() {
		//initialize
		this.svTypeSizeRangeDAPRangeSVRecordsMapMapMap=new HashMap<>();
		this.svTypeSizeRangeDAPRangeGroupStatSummaryMapMapMap=new HashMap<>();
		
		for(SimpleSVType svType:this.svTypeSVRecordsMap.keySet()) {
			this.svTypeSizeRangeDAPRangeSVRecordsMapMapMap.put(svType, new HashMap<>());
			this.svTypeSizeRangeDAPRangeGroupStatSummaryMapMapMap.put(svType, new HashMap<>());
			for(NumericRange sr:sizeRanges) {
				this.svTypeSizeRangeDAPRangeSVRecordsMapMapMap.get(svType).put(sr, new HashMap<>());
				this.svTypeSizeRangeDAPRangeGroupStatSummaryMapMapMap.get(svType).put(sr, new HashMap<>());
				
				for(NumericRange dr:this.derivedAlleleProportionRanges) {
					this.svTypeSizeRangeDAPRangeSVRecordsMapMapMap.get(svType).get(sr).put(dr, new ArrayList<>());
				}
			}
		}
		
		//assign SVRecords to each group
		for(SimpleSVType svType:this.svTypeSVRecordsMap.keySet()) {
			for(SVRecord sv:this.svTypeSVRecordsMap.get(svType)) {
				NumericRange sizeRange=null;
				for(NumericRange sr:this.sizeRanges) {
					if(sr.covers(sv.getSize())) {
						sizeRange=sr;
						break;
					}
				}
				
				NumericRange dapRange=null;
				for(NumericRange dr:this.derivedAlleleProportionRanges) {
					if(dr.covers(sv.getDerivedAlleleProportion())) {
						dapRange=dr;
						break;
					}
				}
				if(sizeRange!=null && dapRange!=null) {
					this.svTypeSizeRangeDAPRangeSVRecordsMapMapMap.get(svType).get(sizeRange).get(dapRange).add(sv);
				}
			}
		}
		
		
		///build GroupStatSummaries
		for(SimpleSVType svType:this.svTypeSizeRangeDAPRangeSVRecordsMapMapMap.keySet()) {
			for(NumericRange sizeRange:this.svTypeSizeRangeDAPRangeSVRecordsMapMapMap.get(svType).keySet()) {
				for(NumericRange dapRange:this.svTypeSizeRangeDAPRangeSVRecordsMapMapMap.get(svType).get(sizeRange).keySet()) {
					List<SVRecord> svRecords=this.svTypeSizeRangeDAPRangeSVRecordsMapMapMap.get(svType).get(sizeRange).get(dapRange);
					if(!svRecords.isEmpty()) {
						this.svTypeSizeRangeDAPRangeGroupStatSummaryMapMapMap.get(svType).get(sizeRange)
						.put(
								dapRange, 
								new GroupStatSummary(svType, sizeRange, dapRange, svRecords, this.windowNum));
					}
				}
			}
		}
	}
	
	
	/**
	 * @return the svTypeSizeRangeDAPRangeGroupStatSummaryMapMapMap
	 */
	public final Map<SimpleSVType, Map<NumericRange, Map<NumericRange, GroupStatSummary>>> getSvTypeSizeRangeDAPRangeGroupStatSummaryMapMapMap() {
		return svTypeSizeRangeDAPRangeGroupStatSummaryMapMapMap;
	}


	
	

}
