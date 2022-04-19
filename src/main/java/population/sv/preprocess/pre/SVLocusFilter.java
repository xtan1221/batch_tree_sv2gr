package population.sv.preprocess.pre;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import population.sv.reader.DellyJointGenotypeVCFReader;
import population.sv.reader.LumpyJointGenotypeVCFReader;
import population.sv.utils.SimpleSVLocus;
import population.sv.utils.SimpleSVType;


/**
 * 
 * find out the {@link SimpleSVLocus} in a population that are detected by both delly and lumpy 
 * 
 * Strategy:
 * 
 * for every {@link SimpleSVLocus} from {@link DellyJointGenotypeVCFReader}, check if there is a {@link SimpleSVLocus} from LumpyJointGenotypeVCFReader that 
 * 1. of the same sv type
 * 2. regions are close to each other
 * 		1. if overlapping, true
 * 		2. if non-overlapping but distance is within a threshold, true; see {@link #maxDist}
 * 
 * if true for both, the {@link SimpleSVLocus} from {@link DellyJointGenotypeVCFReader} will be kept, otherwise, it will be filtered out
 * 
 * 
 * @author tanxu
 *
 */
public class SVLocusFilter {
	/**
	 * 
	 */
	private final DellyJointGenotypeVCFReader dellyJointGenotypeVCFReader;
	/**
	 * 
	 */
	private final LumpyJointGenotypeVCFReader lumpyJointGenotypeVCFReader;
	
	/**
	 * one SV locus 'd' detected by Delly is considered also detected by lumpy if 
	 * 1. there is a SV locus 'l' detected by lump with the same SV type and 
	 * 2. the region of 'd' is not away from 'l' more than this distance
	 */
	private final int maxDist;

	
	//////////////////
	private Map<SimpleSVType, List<SimpleSVLocus>> svTypeFilteredLocusListMap; 
	
	
	public SVLocusFilter(DellyJointGenotypeVCFReader dellyJointGenotypeVCFReader,
			LumpyJointGenotypeVCFReader lumpyJointGenotypeVCFReader, int maxDist) {
		super();
		this.dellyJointGenotypeVCFReader = dellyJointGenotypeVCFReader;
		this.lumpyJointGenotypeVCFReader = lumpyJointGenotypeVCFReader;
		this.maxDist = maxDist;
		
		
		this.run();
	}
	

	void run() {
		this.svTypeFilteredLocusListMap=new HashMap<>();
		
		System.out.println("start");
		for(SimpleSVType svType:this.dellyJointGenotypeVCFReader.getSvTypeChromLocusListMapMap().keySet()) {
			
			for(String chrom:this.dellyJointGenotypeVCFReader.getSvTypeChromLocusListMapMap().get(svType).keySet()) {
				System.out.println(svType+"\t"+chrom);
				for(SimpleSVLocus dellySV:this.dellyJointGenotypeVCFReader.getSvTypeChromLocusListMapMap().get(svType).get(chrom)) {
					
					boolean sameSVDetectedByLumpy=false;
					//
					if(this.lumpyJointGenotypeVCFReader.getSvTypeChromLocusListMapMap().containsKey(svType) && 
							this.lumpyJointGenotypeVCFReader.getSvTypeChromLocusListMapMap().get(svType).containsKey(chrom)) {
						//
						for(SimpleSVLocus lumpySV:this.lumpyJointGenotypeVCFReader.getSvTypeChromLocusListMapMap().get(svType).get(chrom)) {
							//lumpy sv is at left side of delly sv within max dist, stop
							if(lumpySV.getEnd()<dellySV.getStart() && dellySV.getStart()-lumpySV.getEnd()<=this.maxDist) {
								sameSVDetectedByLumpy=true;
							}
//							//lumpy sv is at left side of delly sv within max dist, stop
							if(lumpySV.getStart()>dellySV.getEnd() && lumpySV.getStart()-dellySV.getEnd()<=this.maxDist) {
								sameSVDetectedByLumpy=true;
							}
							
							//check if overlapping
							if(regionCovering(dellySV.getStart(), dellySV.getEnd(), lumpySV.getStart())
									||regionCovering(dellySV.getStart(), dellySV.getEnd(), lumpySV.getEnd())
									||regionCovering(lumpySV.getStart(), lumpySV.getEnd(), dellySV.getStart())
									||regionCovering(lumpySV.getStart(), lumpySV.getEnd(), dellySV.getEnd())
									) {
								sameSVDetectedByLumpy=true;
							}
							
							
							//stop when lumpy SV is at right side of delly sv with distance > maxDist
							if(lumpySV.getStart()-dellySV.getEnd()>this.maxDist || sameSVDetectedByLumpy) {
								break;
							}
						}
					}
					
					if(sameSVDetectedByLumpy) {
						if(!this.svTypeFilteredLocusListMap.containsKey(svType))
							this.svTypeFilteredLocusListMap.put(svType, new ArrayList<>());
						
						this.svTypeFilteredLocusListMap.get(svType).add(dellySV);
					}
				}
			}
		}
	}

	private boolean regionCovering(int regionStart, int regionEnd, int pos) {
		return regionStart<=pos&&regionEnd>=pos;
	}
	
	
	/**
	 * @return the svTypeFilteredLocusListMap
	 */
	public Map<SimpleSVType, List<SimpleSVLocus>> getSvTypeFilteredLocusListMap() {
		return svTypeFilteredLocusListMap;
	}
	
	public List<SimpleSVLocus> getAllQualifiedSVLoci(){
		List<SimpleSVLocus> ret = new ArrayList<>();
		
		for(SimpleSVType type:this.svTypeFilteredLocusListMap.keySet()) {
			ret.addAll(this.svTypeFilteredLocusListMap.get(type));
		}
		
		return ret;
	}
	
	/**
	 * return the set of SimpleSVTypes with at least one SimpleSVLocus that passed filter
	 * @return
	 */
	public Set<SimpleSVType> getSimpleSVTypes(){
		Set<SimpleSVType> ret = new HashSet<>();
		
		for(SimpleSVType type:this.svTypeFilteredLocusListMap.keySet()) {
			if(!this.svTypeFilteredLocusListMap.get(type).isEmpty()) {
				ret.add(type);
			}
		}
		return ret;
	}
	
}
