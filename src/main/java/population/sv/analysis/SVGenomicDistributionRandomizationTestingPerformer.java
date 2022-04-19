package population.sv.analysis;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import genomics.chrom.ChromLenReader;
import genomics.utils.AbstractGenomicRegion;
import genomics.utils.SimpleGenomicRegion;
import population.sv.utils.SimpleSVLocus;
import population.sv.utils.SimpleSVType;
import statistics.VariantOverlappedGenomicRegionQuerier;
import statistics.randomization.RandomGenomicLocationGenerator;

/**
 * perform randomization testing for each type of SVs on whether there is preference for SV to overlap with a set of target genomic regions
 * 
 * for example, the target genomic regions could be the genes or exons/CDS/introns/UTRs, etc
 * 
 * @author tanxu
 * @param <R> the type of genomic region
 */
public class SVGenomicDistributionRandomizationTestingPerformer<R extends AbstractGenomicRegion> {
	private final List<SimpleSVLocus> targetSVs;
	private final ChromLenReader chromLenReader;
	private final List<R> targetGenomicRegions;
	private final int iterations;
	private final boolean fullyCovered;
	private final Class<R> clazz;
	private final Path outDir;
	
	///////////////////////////////
	/**
	 * map from sv type to the list of sv locus in the {@link #targetSVs}
	 */
	private Map<SimpleSVType, List<SimpleSVLocus>> svTypeListMap;
	///////////////////
	/**
	 * number of genomic regions in {@link #targetGenomicRegions} that are overlapping with one or more of the {@link #targetSVs}
	 */
	private int targetGenomicRegionNumOverlappingWithTargetSV;
	/**
	 * map from the sv type to the number of genomic regions in {@link #targetGenomicRegions} that are overlapping with one or more of the {@link #targetSVs} of that specific type
	 */
	private Map<SimpleSVType, Integer> svTypeTargetGenomicRegionNumOverlappingWithTargetSVMap;
	//////////////////////
	/**
	 * map from sv type to the list of sv sizes of that type
	 */
	private Map<SimpleSVType, List<Integer>> svTypeSizeListMap;
	private List<Integer> allTypeSVSizeList;
	////////////////
	private List<Integer> allSVTypeRandomOverlappingGeneNumList;
	private Map<SimpleSVType, List<Integer>> svTypeRandomOverlappingGeneNumListMap;
	
	//////////////
	/**
	 * output file to store the {@link #targetGenomicRegionNumOverlappingWithTargetSV} and {@link #svTypeTargetGenomicRegionNumOverlappingWithTargetSVMap}
	 */
	private Path targetGenomicRegionNumOverlappingWithTargetSVOutputFile;
	/**
	 * output file to store the number of genomic regions in {@link #targetGenomicRegions} that are overlapping with one or more of the simulated SVs of each iteration of randomization testing one a single line
	 * the simulated SVs include all sv types in the given {@link #targetSVs}
	 */
	private Path allSVTypeTestResultFile;
	/**
	 * output file to store simulation result for each sv type
	 */
	private Map<SimpleSVType, Path> svTypeTestResultOutputFileMap;
	
	
	public SVGenomicDistributionRandomizationTestingPerformer(
			List<SimpleSVLocus> targetSVs, ChromLenReader chromLenReader,
			List<R> targetGenomicRegions, Class<R> clazz,
			int iterations, boolean fullyCovered,
			Path outDir) {
		super();
		this.targetSVs = targetSVs;
		this.chromLenReader = chromLenReader;
		this.targetGenomicRegions = targetGenomicRegions;
		this.iterations = iterations;
		this.fullyCovered = fullyCovered;
		this.clazz=clazz;
		this.outDir = outDir;
		
		////////////////
		this.prepare();
		this.findoutTargetGenomicRegionNumOverlappingWithTargetSV();
		
		/////
		this.generateAllSVTypeRandomGeneNumList();
		this.generateOverlappingGeneNumListForEachSVType();
		
		//////
		this.outputTargetGenomicRegionNumOverlappingWithTargetSVOutputFile();
		this.outputAllSVTypeGeneNumListToFile();
		this.outputGeneNumListOfEachSVTypeToFiles();
	}
	
	void prepare() {
		this.svTypeListMap=new HashMap<>();
		
		for(SimpleSVLocus sv:this.targetSVs) {
			if(!this.svTypeListMap.containsKey(sv.getType())) {
				this.svTypeListMap.put(sv.getType(), new ArrayList<>());
			}
			this.svTypeListMap.get(sv.getType()).add(sv);
		}
		
		this.svTypeListMap.forEach((k,v)->{
			System.out.println(k+"\t"+v.size());
		});
		
		/////////////////////
		this.allTypeSVSizeList=new ArrayList<>();
		svTypeSizeListMap = new HashMap<>();
		for(SimpleSVLocus sv:this.targetSVs) {
			this.allTypeSVSizeList.add(sv.getSize());
			if(!svTypeSizeListMap.containsKey(sv.getType())) {
				svTypeSizeListMap.put(sv.getType(), new ArrayList<>());
			}
			svTypeSizeListMap.get(sv.getType()).add(sv.getSize());
		}
		
		/////////////
		this.targetGenomicRegionNumOverlappingWithTargetSVOutputFile=Path.of(this.outDir.toString(), this.clazz.getSimpleName().concat(".target.genomic.region.num.overlapping.with.target.sv.txt"));
		if(targetGenomicRegionNumOverlappingWithTargetSVOutputFile.toFile().exists()) {
			System.out.println("targetGenomicRegionNumOverlappingWithTargetSVOutputFile exists, delete it...");
			targetGenomicRegionNumOverlappingWithTargetSVOutputFile.toFile().delete();
		}
		///
		this.allSVTypeTestResultFile=Path.of(this.outDir.toString(), this.clazz.getSimpleName().concat(".all.sv.type.").concat("random.test.result.tsv"));
		if(allSVTypeTestResultFile.toFile().exists()) {
			System.out.println("allSVTypeTestResultFile exists, delete it...");
			allSVTypeTestResultFile.toFile().delete();
		}
		////
		this.svTypeTestResultOutputFileMap=new HashMap<>();
		for(SimpleSVType type:svTypeListMap.keySet()) {
			Path outFile=Path.of(this.outDir.toString(), this.clazz.getSimpleName().concat(".").concat(type.toString()).concat(".type.".concat("random.test.result.tsv")));
			if(outFile.toFile().exists()) {
				System.out.println("sv type outFile exists, delete it...");
				outFile.toFile().delete();
			}
			
			this.svTypeTestResultOutputFileMap.put(type, outFile);
		}
		
	}
	
	
	void findoutTargetGenomicRegionNumOverlappingWithTargetSV() {
		VariantOverlappedGenomicRegionQuerier<R, SimpleSVLocus> querier = 
				new VariantOverlappedGenomicRegionQuerier<>(
						targetGenomicRegions,//List<T> targetGenomicRegions,
						this.targetSVs,//List<Q> queryGenomicRegions, 
						fullyCovered//boolean fullyCovered
						);
		
		this.targetGenomicRegionNumOverlappingWithTargetSV=querier.getQueriedTargetGenomicRegions().size();
		
		
		////////////
		this.svTypeTargetGenomicRegionNumOverlappingWithTargetSVMap=new HashMap<>();
		
		for(SimpleSVType type:this.svTypeListMap.keySet()) {
			VariantOverlappedGenomicRegionQuerier<R, SimpleSVLocus> querier2 = 
					new VariantOverlappedGenomicRegionQuerier<>(
							targetGenomicRegions,//List<T> targetGenomicRegions,
							this.svTypeListMap.get(type),//List<Q> queryGenomicRegions, 
							fullyCovered//boolean fullyCovered
							);
			
			this.svTypeTargetGenomicRegionNumOverlappingWithTargetSVMap.put(type, querier2.getQueriedTargetGenomicRegions().size());
		}
	}
	
	
	void generateAllSVTypeRandomGeneNumList() {
		this.allSVTypeRandomOverlappingGeneNumList=new ArrayList<>();
		
		for(int i=0;i<iterations;i++) {
			RandomGenomicLocationGenerator generator=
					new RandomGenomicLocationGenerator(
							chromLenReader,//ChromLenReader chromLenReader, 
							this.allTypeSVSizeList,//List<Integer> genomicRegionSizes, 
							this.targetSVs.size()//int replicates
							);
			
			VariantOverlappedGenomicRegionQuerier<R, SimpleGenomicRegion> querier = 
					new VariantOverlappedGenomicRegionQuerier<>(
							targetGenomicRegions,//List<T> targetGenomicRegions,
							generator.getGeneratedRandomGenomicRegions(),//List<Q> queryGenomicRegions, 
							fullyCovered//boolean fullyCovered
							);
			
			this.allSVTypeRandomOverlappingGeneNumList.add(querier.getQueriedTargetGenomicRegions().size());
			System.out.println(i+":"+querier.getQueriedTargetGenomicRegions().size());
		}
		
		/////////////////
		Collections.sort(this.allSVTypeRandomOverlappingGeneNumList);
	}
	
	void generateOverlappingGeneNumListForEachSVType() {
		//////////////randomization
		svTypeRandomOverlappingGeneNumListMap=new HashMap<>();
		for(SimpleSVType svType:this.svTypeListMap.keySet()) {
			List<Integer> geneNumList=new ArrayList<>();
			for(int i=0;i<iterations;i++) {
				RandomGenomicLocationGenerator generator=
						new RandomGenomicLocationGenerator(
								chromLenReader,//ChromLenReader chromLenReader, 
								this.svTypeSizeListMap.get(svType),//List<Integer> genomicRegionSizes, 
								this.svTypeListMap.get(svType).size()//int replicates
								);
				
				VariantOverlappedGenomicRegionQuerier<R, SimpleGenomicRegion> querier = 
						new VariantOverlappedGenomicRegionQuerier<>(
								targetGenomicRegions,//List<T> targetGenomicRegions,
								generator.getGeneratedRandomGenomicRegions(),//List<Q> queryGenomicRegions, 
								fullyCovered//boolean fullyCovered
								);
				
				geneNumList.add(querier.getQueriedTargetGenomicRegions().size());
				System.out.println(i+":"+querier.getQueriedTargetGenomicRegions().size());
			}
			
			/////////////////
			Collections.sort(geneNumList);
			this.svTypeRandomOverlappingGeneNumListMap.put(svType, geneNumList);
		}
	}
	
	void outputTargetGenomicRegionNumOverlappingWithTargetSVOutputFile() {
		try {
			////////////////////////
			BufferedWriter writer = new BufferedWriter(new FileWriter(this.targetGenomicRegionNumOverlappingWithTargetSVOutputFile.toString()));
			writer.append("all SV types:\t").append(Integer.toString(this.targetGenomicRegionNumOverlappingWithTargetSV));
			writer.newLine();
			/////////////
			for(SimpleSVType type:this.svTypeTargetGenomicRegionNumOverlappingWithTargetSVMap.keySet()) {
				writer.append(type.toString()).append(":\t").append(Integer.toString(this.svTypeTargetGenomicRegionNumOverlappingWithTargetSVMap.get(type)));
				writer.newLine();
			}
			
			
			//////////////
			writer.flush();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	void outputAllSVTypeGeneNumListToFile() {
		try {
			////////////////////////
			BufferedWriter writer = new BufferedWriter(new FileWriter(this.allSVTypeTestResultFile.toString()));
			for(int num:this.allSVTypeRandomOverlappingGeneNumList) {
				writer.append(Integer.toString(num));
				writer.newLine();
			}
			writer.flush();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	void outputGeneNumListOfEachSVTypeToFiles() {
		try {
			////////////////////////
			for(SimpleSVType type:this.svTypeRandomOverlappingGeneNumListMap.keySet()) {
				BufferedWriter writer = new BufferedWriter(new FileWriter(this.svTypeTestResultOutputFileMap.get(type).toString()));
				for(int num:this.svTypeRandomOverlappingGeneNumListMap.get(type)) {
					writer.append(Integer.toString(num));
					writer.newLine();
				}
				writer.flush();
				writer.close();
			}
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
}
