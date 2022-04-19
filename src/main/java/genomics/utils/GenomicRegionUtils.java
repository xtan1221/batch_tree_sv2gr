package genomics.utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


public class GenomicRegionUtils {
	
	
	/**
	 * return the sorter by chrom and start position (asc order)
	 * @return
	 */
	public static Comparator<? super AbstractGenomicRegion> sorterByChromAndStartPos(){
		return (a,b)->{
			if(!a.getChrom().equals(b.getChrom())) {
				return a.getChrom().compareTo(b.getChrom());
			}else {
				return a.getStart()-b.getStart();
			}
		};
	}
	
	/**
	 * merge the overlapping SimpleGenomicRegions in the given list into a list of non-overlapping ones
	 * 
	 * the strand is ignored;
	 * 
	 * @param regions
	 * @return
	 */
	public static List<SimpleGenomicRegion> mergeOverlappingRegions(List<SimpleGenomicRegion> regions){
		List<SimpleGenomicRegion> ret = new ArrayList<>();
		
		Collections.sort(regions, sorterByChromAndStartPos());
		
		String currentChrom=null;
		Integer currentRegionStart=null;
		Integer currentRegionEnd=null;
		for(SimpleGenomicRegion r:regions) {
			if(currentChrom==null||!currentChrom.equals(r.getChrom())) {
				if(currentChrom!=null&&!currentChrom.equals(r.getChrom())) {
					ret.add(new SimpleGenomicRegion(currentChrom, currentRegionStart, currentRegionEnd));
				}
				//update
				currentChrom=r.getChrom();
				currentRegionStart=r.getStart();
				currentRegionEnd=r.getEnd();
			}else {//same chrom
				if(r.getStart()<=currentRegionEnd) {//overlapping with current region
					if(r.getEnd()>currentRegionEnd) {//update end if needed
						currentRegionEnd=r.getEnd();
					}
				}else {//not overlapping with current region
					ret.add(new SimpleGenomicRegion(currentChrom, currentRegionStart, currentRegionEnd));
					//update
					currentChrom=r.getChrom();
					currentRegionStart=r.getStart();
					currentRegionEnd=r.getEnd();
				}
			}
		}
		
		//add the last region
		ret.add(new SimpleGenomicRegion(currentChrom, currentRegionStart, currentRegionEnd));
		
		////
		return ret;
	}
	
	/**
	 * merge the given GenomicRegions into groups by the given distance;
	 * 
	 * strand is not considered;
	 * 
	 * if strand should be considered, pre-process the GenomicRegions before feeding to this method
	 * @param regions
	 * @param minDist
	 * @param maxDist
	 * @return
	 */
	public static <T extends SimpleGenomicRegion> Set<List<T>> mergeByDistance(List<T> regions, int maxDist){
		
		Set<List<T>> ret = new HashSet<>();
		
		Collections.sort(regions, sorterByChromAndStartPos());
		
		Integer currentGroupEnd=null;
		List<T> currentGroup=new ArrayList<>();
		for(T gr:regions) {
			if(currentGroupEnd==null ||currentGroup.isEmpty()) {//
				currentGroupEnd=gr.getEnd();
			}else {//
				if(currentGroup.get(0).getChrom().equals(gr.getChrom())) {//same chrom
					if(gr.getStart()-currentGroupEnd>maxDist) {//new group
						ret.add(currentGroup);
						currentGroup=new ArrayList<>();
						currentGroupEnd=gr.getEnd();
					}else {//within the maxDist
						if(gr.getEnd()>currentGroupEnd) //update the currentGroupEnd
							currentGroupEnd=gr.getEnd();
					}
				}else {//different chrom
					ret.add(currentGroup);
					currentGroup=new ArrayList<>();
					currentGroupEnd=gr.getEnd();
				}
			}
			//
			currentGroup.add(gr);
		}
		
		//add the last group
		if(!currentGroup.isEmpty())
			ret.add(currentGroup);
		
		return ret;
	}
	
	
	/**
	 * 
	 * @param gr1
	 * @param gr2
	 * @return
	 */
	public static boolean overlapping(AbstractGenomicRegion gr1, AbstractGenomicRegion gr2) {
		if(gr1.getChrom().equals(gr2.getChrom())) {
			return gr1.covers(gr2.getStart()) || gr1.covers(gr2.getEnd()) || gr2.covers(gr1.getStart()) ||gr2.covers(gr1.getEnd());
		}else {
			return false;
		}
	}
	
	/**
	 * 
	 * @param gr1
	 * @param gr2
	 * @return
	 */
	public static boolean region1FullyCoversRegion2(AbstractGenomicRegion gr1, AbstractGenomicRegion gr2) {
		if(gr1.getChrom().equals(gr2.getChrom())) {
			return gr1.getStart()<=gr2.getStart() && gr1.getEnd() >=gr2.getEnd();
		}else {
			return false;
		}
		
	}
}
