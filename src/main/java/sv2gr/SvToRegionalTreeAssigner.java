package sv2gr;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import population.sv.utils.SimpleSVLocus;
import population.sv.utils.SimpleSVType;
import population.utils.Genotype;
import sv2gr.tree.RegionalTree;
import genomics.utils.GenomicRegionUtils;

/**
 * utility class to assign each {@link SimpleSVLocus} to all {@link RegionalTree}s that covering the SV
 * 
 * note that a {@link SimpleSVLocus} may be assigned to multiple {@link RegionalTree}s
 * 
 * also note that DUP and CNV type SVs are not included in GR event inference
 * 
 * @author tanxu
 *
 */
public class SvToRegionalTreeAssigner {
	/**
	 * all regional trees
	 */
	private final List<RegionalTree> regionalTrees;
	/**
	 * all SVs
	 */
	private final List<SimpleSVLocus> allSVLocusList;
	
//	/**
//	 * information for ingroup samples
//	 */
//	private final PopulationStructureFileReader populationStructureFileReader;
	
	/**
	 * indices of ingroup samples
	 */
	private final List<Integer> ingroupSampleIndices;
	
	
	/**
	 * the minimal covered proportion of the SV by the genomic region of a regional tree so that the SV can be assigned to the tree
	 */
	private final double minCoveredPropertion;
	
	/**
	 * if true, all SV with one or more ingroup samples with HETER genotype will be excluded from GR event inference
	 * if false, they will be included
	 */
	private final boolean toIgnoreSVWithOneOrMoreIngroupSampleWithHeterGenotype;

	
	/////////////////////////////////
	/**
	 * the list of {@link SimpleSVLocus} with none of the ingroup samples contains '0/1' genotype
	 */
	private Map<String, List<SimpleSVLocus>> chromQualifiedSvLocusListMap;
	
	private Map<String, List<RegionalTree>> chromRegionalTreesMap;
	
	public SvToRegionalTreeAssigner(
			List<RegionalTree> regionalTrees, 
			List<SimpleSVLocus> allSVLocusList,
//			PopulationStructureFileReader populationStructureFileReader, 
			List<Integer> ingroupSampleIndices,
			double minCoveredPropertion,
			boolean toIgnoreSVWithOneOrMoreIngroupSampleWithHeterGenotype
			) {
		super();
		this.regionalTrees = regionalTrees;
		this.allSVLocusList = allSVLocusList;
//		this.populationStructureFileReader = populationStructureFileReader;
		this.ingroupSampleIndices=ingroupSampleIndices;
		this.minCoveredPropertion = minCoveredPropertion;
		this.toIgnoreSVWithOneOrMoreIngroupSampleWithHeterGenotype=toIgnoreSVWithOneOrMoreIngroupSampleWithHeterGenotype;
//		this.toTreatHeterGenotypeAsPresence=toTreatHeterGenotypeAsPresence;
		
		//////////////////////////////
		this.processSimpleSVLocusByIngroupSampleGenotype();
		
		this.sortRegionalTreesAndQualifedSVs();
		
//		this.assignSVLocusToRegionalTree();
		this.assignSVLocusToRegionalTree2();
	}

	/**
	 * 
	 */
	void processSimpleSVLocusByIngroupSampleGenotype() {
		this.chromQualifiedSvLocusListMap = new HashMap<>();
		
		for(SimpleSVLocus svLocus:this.allSVLocusList) {
			if(svLocus.getType().equals(SimpleSVType.CNV) || svLocus.getType().equals(SimpleSVType.DUP)) {//skip CNV and DUP type SVs since their reversed GR events are difficult to interpret
				continue;
			}
			
			
			if(this.toIgnoreSVWithOneOrMoreIngroupSampleWithHeterGenotype) {
				boolean ingroupSampleWithHeterGenotypeFound=false;
				
				for(int ingroupSampleIndex:this.ingroupSampleIndices) {
					if(svLocus.getSampleIndexGenotypeMap().get(ingroupSampleIndex).equals(Genotype.HETER)) {
						ingroupSampleWithHeterGenotypeFound=true;
					}
				}
				
				
				if(!ingroupSampleWithHeterGenotypeFound) {
					if(!this.chromQualifiedSvLocusListMap.containsKey(svLocus.getChrom())) {
						this.chromQualifiedSvLocusListMap.put(svLocus.getChrom(), new ArrayList<>());
					}
					
					this.chromQualifiedSvLocusListMap.get(svLocus.getChrom()).add(svLocus);
				}
			}else {
				if(!this.chromQualifiedSvLocusListMap.containsKey(svLocus.getChrom())) {
					this.chromQualifiedSvLocusListMap.put(svLocus.getChrom(), new ArrayList<>());
				}
				this.chromQualifiedSvLocusListMap.get(svLocus.getChrom()).add(svLocus);
			}
		}
	}
	
	
	void sortRegionalTreesAndQualifedSVs() {
		////////////
		this.chromRegionalTreesMap=new HashMap<>();
		
		for(RegionalTree tree:this.regionalTrees) {
			if(!this.chromRegionalTreesMap.containsKey(tree.getChrom())) {
				this.chromRegionalTreesMap.put(tree.getChrom(), new ArrayList<>());
			}
			this.chromRegionalTreesMap.get(tree.getChrom()).add(tree);
		}
		//////////////
		for(String chrom:this.chromRegionalTreesMap.keySet()) {
			Collections.sort(
					this.chromRegionalTreesMap.get(chrom),
					(a,b)->{
						if(a.getChrom().equals(b.getChrom())) {
							return a.getStart()-b.getStart();
						}else {
							return a.getChrom().compareTo(b.getChrom());
						}
					}
					);
			
		}
		
		///////////////
		for(String chrom:this.chromQualifiedSvLocusListMap.keySet()) {
			Collections.sort(
					this.chromQualifiedSvLocusListMap.get(chrom),
					GenomicRegionUtils.sorterByChromAndStartPos()
//					(a,b)->{
//						if(a.getChrom().equals(b.getChrom())) {
//							return a.getStart()-b.getStart();
//						}else {
//							return a.getChrom().compareTo(b.getChrom());
//						}
//					}
					);
		}
		
	}
	
	
	
	void assignSVLocusToRegionalTree() {
		for(String chrom:this.chromQualifiedSvLocusListMap.keySet()) {
			if(this.chromRegionalTreesMap.containsKey(chrom)) {
				for(SimpleSVLocus sv:this.chromQualifiedSvLocusListMap.get(chrom)) {
					RegionalTree coveringTree=null;
					Double coveredSVProportion =null;
					for(RegionalTree tree:this.chromRegionalTreesMap.get(chrom)) {
						if(tree.getStart()>sv.getEnd()) {
							break;
						}else if(tree.getEnd()<sv.getStart()) {
							//
						}else {
							int coveredLen=this.getCoveredRegionLength(tree.getStart(), tree.getEnd(), sv.getStart(), sv.getEnd());
							double coveredProportion=(double)coveredLen/sv.getSize();
							if(coveredProportion>=this.minCoveredPropertion) {
								if(coveringTree==null || coveredProportion>coveredSVProportion) {//update if a regional tree covers more proportion of the current sv
									coveringTree=tree;
									coveredSVProportion=coveredProportion;
								}
							}
						}
					}
					if(coveringTree!=null) {
						coveringTree.addCoveredSVLocus(sv);
					}
				}
			}
		}
	}
	
	/**
	 * assign each {@link SimpleSVLocus} to all {@link RegionalTree}s that covering proportion of the {@link SimpleSVLocus} more than {@link #minCoveredPropertion}
	 * TODO testing
	 */
	void assignSVLocusToRegionalTree2() {
		for(String chrom:this.chromQualifiedSvLocusListMap.keySet()) {
			if(this.chromRegionalTreesMap.containsKey(chrom)) {
				for(SimpleSVLocus sv:this.chromQualifiedSvLocusListMap.get(chrom)) {
					for(RegionalTree tree:this.chromRegionalTreesMap.get(chrom)) {
						if(tree.getStart()>sv.getEnd()) {
							break;
						}else if(tree.getEnd()<sv.getStart()) {
							//
						}else {//overlapping
							int coveredLen=this.getCoveredRegionLength(tree.getStart(), tree.getEnd(), sv.getStart(), sv.getEnd());
							double coveredProportion=(double)coveredLen/sv.getSize();
							if(coveredProportion>=this.minCoveredPropertion) {
								tree.addCoveredSVLocus(sv);
							}
						}
					}
				}
			}
		}
	}
	
	
	/**
	 * return the length of the region 2 covered by region 1
	 * @param start1
	 * @param end1
	 * @param start2
	 * @param end2
	 * @return
	 */
	private int getCoveredRegionLength(int start1, int end1, int start2, int end2) {
		if(end2<start1||start2>end1) {//not overlapping
			return 0;
		}else {//overlapping
			if(start1<=start2 && end1>=end2) {//region 1 fully covers region2
				return end2-start2+1;
			}else {//region 1 covers a part of region2
				if(end1>=end2) {
					return end2-start1+1;
				}else {
					return end1-start2+1;
				}
			}
		}
	}
	
}
