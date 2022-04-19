package population.sv.analysis.adjacentRegion.polymophismSiteProportion;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import genomics.chrom.ChromLenReader;
import population.sv.utils.SimpleSVLocus;


/**
 * 
 * @author tanxu
 *
 */
public class SVAdjacentRegionWindowBuilder {
	/**
	 * 
	 */
	private final List<SimpleSVLocus> svLocusList;
	/**
	 * 
	 */
	private final int windowSize;
	/**
	 * number of windows to calculate at each side of each SV
	 */
	private final int windowNum;
	
	/**
	 * only SVs on chroms included in this reader are to be included to build windows!
	 */
	private final ChromLenReader chromLenReader;
	
	////////////////////////
	/**
	 * left windows, the one with index = 0 is the one closest to the start of the sv
	 */
	private Map<SimpleSVLocus, List<Window>> svLeftSideWindowsMap;
	/**
	 * right windows, the one with the index=0 is the one closest to the end of the sv
	 */
	private Map<SimpleSVLocus, List<Window>> svRightSideWindowsMap;
	
	////////////
	private Map<String,List<Window>> chromBuiltWindowsMap;
	
	
	public SVAdjacentRegionWindowBuilder(
			List<SimpleSVLocus> svLocusList, Integer windowSize, Integer windowNum, ChromLenReader chromLenReader) {
		super();
		if(windowSize==null||windowNum==null)
			throw new IllegalArgumentException("given windowSize and windowNum cannot be null!");
		
		this.svLocusList = svLocusList;
		this.windowSize = windowSize;
		this.windowNum = windowNum;
		this.chromLenReader=chromLenReader;
		
		
		////////////////////////
		this.run();
	}
	
	
	/**
	 * 
	 */
	void run() {
		this.svLeftSideWindowsMap=new HashMap<>();
		this.svRightSideWindowsMap=new HashMap<>();
		
		this.chromBuiltWindowsMap = new HashMap<>();
		
		for(SimpleSVLocus locus:this.svLocusList) {
			if(!this.chromLenReader.getChromNameLengthMap().containsKey(locus.getChrom())) {//only include SVs on the chrom of the given chrom length file
				continue;
			}
			////////////////////
			this.svLeftSideWindowsMap.put(locus, new ArrayList<>());
			this.svRightSideWindowsMap.put(locus, new ArrayList<>());
			
			if(!this.chromBuiltWindowsMap.containsKey(locus.getChrom())) {
				this.chromBuiltWindowsMap.put(locus.getChrom(), new ArrayList<>());
			}
			
			for(int i=0;i<this.windowNum;i++) {
				///left window
				int start=locus.getStart()-(i+1)*this.windowSize;
				int end=start+this.windowSize-1;
				if(start>0) {
					Window leftWindow=new Window(locus.getChrom(), start, end);
					this.svLeftSideWindowsMap.get(locus).add(leftWindow);
					this.chromBuiltWindowsMap.get(locus.getChrom()).add(leftWindow);
				}else {
					this.svLeftSideWindowsMap.get(locus).add(null);
				}
				
				///right window
				start=locus.getEnd()+i*this.windowSize;
				end=start+this.windowSize-1;
				if(this.chromLenReader.getChromNameLengthMap().get(locus.getChrom())>=end) {
					Window rightWindow=new Window(locus.getChrom(), locus.getEnd()+i*this.windowSize, locus.getEnd()+(i+1)*this.windowSize);
					this.svRightSideWindowsMap.get(locus).add(rightWindow);
					this.chromBuiltWindowsMap.get(locus.getChrom()).add(rightWindow);
				}else {
					this.svRightSideWindowsMap.get(locus).add(null);
				}
			}
		}
	}
	
	
	/**
	 * @return the windowNum
	 */
	public int getWindowNum() {
		return windowNum;
	}


	/**
	 * @return the chromBuiltWindowsMap
	 */
	public Map<String, List<Window>> getChromBuiltWindowsMap() {
		return chromBuiltWindowsMap;
	}


	/**
	 * @return the svLeftSideWindowsMap
	 */
	public Map<SimpleSVLocus, List<Window>> getSvLeftSideWindowsMap() {
		return svLeftSideWindowsMap;
	}


	/**
	 * @return the svRightSideWindowsMap
	 */
	public Map<SimpleSVLocus, List<Window>> getSvRightSideWindowsMap() {
		return svRightSideWindowsMap;
	}


	/**
	 * @return the windowSize
	 */
	public int getWindowSize() {
		return windowSize;
	}

	
}
