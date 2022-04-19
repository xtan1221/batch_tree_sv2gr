package statistics;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import genomics.utils.AbstractGenomicRegion;
import genomics.utils.GenomicRegionUtils;

/**
 * generic class to find out the set of target {@link GenomicRegion}s that are overlapping with one or more query {@link GenomicRegion}s
 * 
 * @author tanxu
 * 
 */
public class VariantOverlappedGenomicRegionQuerier<T extends AbstractGenomicRegion,Q extends AbstractGenomicRegion> {
	private final List<T> targetGenomicRegions;
	
	private final List<Q> queryGenomicRegions;
	
	/**
	 * if true, only target regions fully covered by query genomic regions will be extracted;
	 * if false, all target regions overlapping with query genomic region will be extracted;
	 */
	private final boolean fullyCovered;
	
	///////////////////////////
	private Map<String, List<T>> chromTargetGenomicRegionsMap;
	private Map<String, List<Q>> chromQueryGenomicRegionsMap;
	
	/**
	 * 
	 */
	private List<T> queriedTargetGenomicRegions;
	
	
	public VariantOverlappedGenomicRegionQuerier(
			List<T> targetGenomicRegions,
			List<Q> queryGenomicRegions, 
			boolean fullyCovered) {
		super();
		this.targetGenomicRegions = targetGenomicRegions;
		this.queryGenomicRegions = queryGenomicRegions;
		this.fullyCovered = fullyCovered;
		
		this.prepare();
		this.run();
	}

	
	void prepare() {
		this.chromTargetGenomicRegionsMap=new HashMap<>();
		
		for(T r:this.targetGenomicRegions) {
			if(!this.chromTargetGenomicRegionsMap.containsKey(r.getChrom())) {
				this.chromTargetGenomicRegionsMap.put(r.getChrom(), new ArrayList<>());
			}
			this.chromTargetGenomicRegionsMap.get(r.getChrom()).add(r);
		}
		
		for(String chrom:this.chromTargetGenomicRegionsMap.keySet()) {
			Collections.sort(this.chromTargetGenomicRegionsMap.get(chrom), GenomicRegionUtils.sorterByChromAndStartPos());
		}
		
		////////////////////////
		this.chromQueryGenomicRegionsMap=new HashMap<>();
		
		for(Q r:this.queryGenomicRegions) {
			if(!this.chromQueryGenomicRegionsMap.containsKey(r.getChrom())) {
				this.chromQueryGenomicRegionsMap.put(r.getChrom(), new ArrayList<>());
			}
			this.chromQueryGenomicRegionsMap.get(r.getChrom()).add(r);
		}
		
		for(String chrom:this.chromQueryGenomicRegionsMap.keySet()) {
			Collections.sort(this.chromQueryGenomicRegionsMap.get(chrom), GenomicRegionUtils.sorterByChromAndStartPos());
		}
		
	}
	
	
	void run() {
		this.queriedTargetGenomicRegions=new ArrayList<>();
		
		for(String chrom:this.chromTargetGenomicRegionsMap.keySet()) {
			if(this.chromQueryGenomicRegionsMap.containsKey(chrom)) {
				for(T target:this.chromTargetGenomicRegionsMap.get(chrom)) {
					
					for(Q query:this.chromQueryGenomicRegionsMap.get(chrom)) {
						if(this.fullyCovered) {
							if(GenomicRegionUtils.region1FullyCoversRegion2(query, target)) {
								this.queriedTargetGenomicRegions.add(target);
								break;
							}
						}else {
							if(GenomicRegionUtils.overlapping(query, target)) {
								this.queriedTargetGenomicRegions.add(target);
								break;
							}
						}
					}
				}
			}
		}
	}
	
	
	/**
	 * @return the queriedTargetGenomicRegions
	 */
	public List<T> getQueriedTargetGenomicRegions() {
		return queriedTargetGenomicRegions;
	}
	
	
}
