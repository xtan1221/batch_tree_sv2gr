package population.sv.analysis.assignSiteToClosestSV.postprocess;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import basic.NumericRange;
import basic.Pair;
import population.sv.analysis.assignSiteToClosestSV.postprocess.utils.WindowedNucSitesCollectorAndWriter;
import population.sv.utils.SimpleSVType;

/**
 * build a batch of {@link WindowedNucSitesCollectorAndWriter}s based on the given conditions
 * 
 * 
 * @author tanxu
 *
 */
public class WindowedNucSitesCollectorAndWriterBatchBuilder {
	private final boolean toWriteNucRecordToFile;
	/**
	 * output dir for the files to write the record data for each nuc site of the corresponding {@link WindowedNucSitesCollectorAndWriter}
	 */
	private final Path outDir;
	
	/**
	 * target SV types
	 */
	private final List<SimpleSVType> svTypes;
	
	/**
	 * ranges for SV size 
	 */
	private final List<NumericRange> svSizeRanges;
	
	/**
	 * 
	 */
	private final List<NumericRange> closestSVDerivedAllelePropRanges;
	
	/**
	 * ordered ranges for distance to closest sv
	 * the window indices are built based on the order of these ranges
	 */
	private final List<NumericRange> orderedDistToClosestSVRanges;
	
	
	////////////////////////////////////////////
	/**
	 * map from sv type to 
	 * 		map from sv size range to 
	 * 			map from closest sv derived allele proportion range to
	 * 				map from dist to closest sv range to the pair of 
	 * 					{@link WindowedNucSitesCollectorAndWriter} of nuc sites at the left side of closest sv
	 * 					{@link WindowedNucSitesCollectorAndWriter} of nuc sites at the right side of closest sv
	 */
	private Map<SimpleSVType, Map<NumericRange,Map<NumericRange, Map<NumericRange, Pair<WindowedNucSitesCollectorAndWriter,WindowedNucSitesCollectorAndWriter>>>>> groupMap;
	/**
	 * all {@link WindowedNucSitesCollectorAndWriter}s in the {@link #groupMap}
	 */
	private List<WindowedNucSitesCollectorAndWriter> allGroups;
	
	/**
	 * 
	 * @param outDir
	 * @param svTypes
	 * @param svSizeRanges
	 * @param closestSVDerivedAllelePropRanges
	 * @param orderedDistToClosestSVRanges
	 * @throws IOException
	 */
	public WindowedNucSitesCollectorAndWriterBatchBuilder(
			boolean toWriteNucRecordToFile, Path outDir, 
			List<SimpleSVType> svTypes, List<NumericRange> svSizeRanges,
			List<NumericRange> closestSVDerivedAllelePropRanges, List<NumericRange> orderedDistToClosestSVRanges
			) throws IOException {
		super();
		this.toWriteNucRecordToFile=toWriteNucRecordToFile;
		this.outDir = outDir;
		this.svTypes = svTypes;
		this.svSizeRanges = svSizeRanges;
		this.closestSVDerivedAllelePropRanges = closestSVDerivedAllelePropRanges;
		this.orderedDistToClosestSVRanges = orderedDistToClosestSVRanges;
		
		////////////////////////
		this.run();
	}



	void run() throws IOException {
		this.groupMap=new HashMap<>();
		this.allGroups=new ArrayList<>();
		
		for(SimpleSVType svType:this.svTypes) {
			if(!this.groupMap.containsKey(svType)) {
				this.groupMap.put(svType, new HashMap<>());
			}
			for(NumericRange svSizeRange:this.svSizeRanges) {
				if(!this.groupMap.get(svType).containsKey(svSizeRange)) {
					this.groupMap.get(svType).put(svSizeRange, new HashMap<>());
				}
				for(NumericRange svDerivedAlleleProRange:this.closestSVDerivedAllelePropRanges) {
					if(!this.groupMap.get(svType).get(svSizeRange).containsKey(svDerivedAlleleProRange)) {
						this.groupMap.get(svType).get(svSizeRange).put(svDerivedAlleleProRange, new HashMap<>());
					}
					
					for(int i=0;i<this.orderedDistToClosestSVRanges.size();i++) {
						NumericRange distToSVRange=this.orderedDistToClosestSVRanges.get(i);
						//
//						System.out.println(svType+"\t"+svSizeRange.toString()+"\t"+svDerivedAlleleProRange.toString()+"\t"+distToSVRange);
						
						double leftWindowIndex = -(0.5+i);
						double rightWindowIndex = 0.5+i;
						
						WindowedNucSitesCollectorAndWriter left =
								new WindowedNucSitesCollectorAndWriter(
										this.toWriteNucRecordToFile,
										this.outDir,//Path outDir, 
										svType,//SimpleSVType targetSVType, 
										svSizeRange,//NumericRange svSizeRange, 
										svDerivedAlleleProRange, //NumericRange svDAPRange, 
										distToSVRange, //NumericRange distToClosestSVRange, 
										leftWindowIndex,//double windowIndex, 
										true//boolean atLeftSideOfClosestSV
										);
						WindowedNucSitesCollectorAndWriter right =
								new WindowedNucSitesCollectorAndWriter(
										this.toWriteNucRecordToFile,
										this.outDir,//Path outDir, 
										svType,//SimpleSVType targetSVType, 
										svSizeRange,//NumericRange svSizeRange, 
										svDerivedAlleleProRange, //NumericRange svDAPRange, 
										distToSVRange, //NumericRange distToClosestSVRange, 
										rightWindowIndex,//double windowIndex, 
										false//boolean atLeftSideOfClosestSV
										);
						
						this.groupMap.get(svType).get(svSizeRange).get(svDerivedAlleleProRange).put(distToSVRange, new Pair<>(left,right));
						this.allGroups.add(left);
						this.allGroups.add(right);
					}
				}
			}
		}
	}



	//////////////////////////////////////////
	/**
	 * return the 
	 * 	map from sv type to 
	 * 		map from sv size range to 
	 * 			map from closest sv derived allele proportion range to
	 * 				map from dist to closest sv range to the pair of 
	 * 					{@link WindowedNucSitesCollectorAndWriter} of nuc sites at the left side of closest sv
	 * 					{@link WindowedNucSitesCollectorAndWriter} of nuc sites at the right side of closest sv
	 *
	 * @return the groupMap
	 */
	public final Map<SimpleSVType, Map<NumericRange, Map<NumericRange, Map<NumericRange, Pair<WindowedNucSitesCollectorAndWriter, WindowedNucSitesCollectorAndWriter>>>>> getGroupMap() {
		return groupMap;
	}
	
	/**
	 * 
	 * @return
	 */
	public final List<WindowedNucSitesCollectorAndWriter> getAllGroups() {
		return allGroups;
	}
	
	
	
}
