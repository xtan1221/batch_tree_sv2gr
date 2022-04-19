package population.sv.analysis;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import genomics.geneAnnotation.GeneExonAnnotationGff3Reader;
import genomics.utils.Gene;
import genomics.utils.GenomicRegionUtils;
import population.sv.utils.SimpleSVLocus;
import population.sv.utils.SimpleSVType;

/**
 * identify the genes that are overlapping with at least one of the given {@link SimpleSVLocus}
 * 
 * also write the detected genes into an output file
 * 
 * @author tanxu
 * 
 */
public class GeneOverlappingWithSVIdentifier {
	/**
	 * 
	 */
	private final List<SimpleSVLocus> targetSimpleSVLocusList;
	
	/**
	 * 
	 */
	private final GeneExonAnnotationGff3Reader geneExonAnnotationGff3Reader;
	
	/**
	 * whether to include only genes that are fully covered by one or more SV or not
	 */
	private final boolean onlyIncludeFullyCoveredGene;
	
	/**
	 * information of target SVs to facilitate make name of output file that distinguishes from other files of genes overlapping with SVs;
	 * for example, 
	 * 		for target SVs that is found in ingroup samples, use 'ingroup.SV'
	 * 		for target SVs with derived allele frequency in ingroup population between 0.8-1, use 'derived.allele.freq.0.8-1'
	 */
	private final String SVInfor;
	
	
	///////////////////////////////////////
	/**
	 * map from gene to the list of overlapping SVs;
	 */
	private Map<Gene, Map<SimpleSVType,Set<SimpleSVLocus>>> geneSVTypeOverlappingSVsMap;
	
	/**
	 * map from the sv type to the set of genes overlapping with one or more sv of that type;
	 * 
	 * note that some genes may be overlapping with multiple types of SVs
	 */
	private Map<SimpleSVType, Set<Gene>> svTypeOverlappingGenesMap;
	
	/**
	 * map from the sv type to the set of SVs of that type that is overlapping with one or more genes;
	 */
	private Map<SimpleSVType, List<SimpleSVLocus>> svTypeSVsOverlappingWithGeneMap;
	
	/**
	 * the full set of genes that are overlapping with one or more SVs
	 */
	private List<Gene> genesOverlappingWithSV;
	
	/**
	 * list of SV types with one or more overlapping genes
	 */
	private List<SimpleSVType> svTypeList;
	
	/**
	 * output file
	 */
	private Path allSVTypeGeneListFile;
	

	public GeneOverlappingWithSVIdentifier(
			List<SimpleSVLocus> targetSimpleSVLocusList,
			GeneExonAnnotationGff3Reader geneExonAnnotationGff3Reader, boolean onlyIncludeFullyCoveredGene,
			String SVInfor
			) {
		super();
		this.targetSimpleSVLocusList = targetSimpleSVLocusList;
		this.geneExonAnnotationGff3Reader = geneExonAnnotationGff3Reader;
		this.onlyIncludeFullyCoveredGene = onlyIncludeFullyCoveredGene;
		this.SVInfor=SVInfor;
		
		this.run();
	}
	
	
	void run() {
		this.geneSVTypeOverlappingSVsMap=new HashMap<>();
		
		this.svTypeOverlappingGenesMap = new HashMap<>();
		this.svTypeSVsOverlappingWithGeneMap=new HashMap<>();
		this.genesOverlappingWithSV=new ArrayList<>();
		
		for(SimpleSVLocus sv: this.targetSimpleSVLocusList) {
			List<Gene> overlappingGenes = 
					this.geneExonAnnotationGff3Reader.queryGene(
							sv.getChrom(), sv.getStart(), sv.getEnd(), onlyIncludeFullyCoveredGene
							);
			
			if(!overlappingGenes.isEmpty()) {
				if(!this.svTypeOverlappingGenesMap.containsKey(sv.getType())) {
					this.svTypeOverlappingGenesMap.put(sv.getType(), new HashSet<>());
				}
				this.svTypeOverlappingGenesMap.get(sv.getType()).addAll(overlappingGenes);
				
				for(Gene gene:overlappingGenes) {
					if(!this.genesOverlappingWithSV.contains(gene)) {
						this.genesOverlappingWithSV.add(gene);
					}
				}
				
				////////////////
				for(Gene gene:overlappingGenes) {
					if(!this.geneSVTypeOverlappingSVsMap.containsKey(gene)) {
						this.geneSVTypeOverlappingSVsMap.put(gene, new HashMap<>());
					}
					
					if(!this.geneSVTypeOverlappingSVsMap.get(gene).containsKey(sv.getType())) {
						this.geneSVTypeOverlappingSVsMap.get(gene).put(sv.getType(), new HashSet<>());
					}
					this.geneSVTypeOverlappingSVsMap.get(gene).get(sv.getType()).add(sv);
				}
				
				///////////////
				if(!this.svTypeSVsOverlappingWithGeneMap.containsKey(sv.getType())) {
					this.svTypeSVsOverlappingWithGeneMap.put(sv.getType(), new ArrayList<>());
				}
				this.svTypeSVsOverlappingWithGeneMap.get(sv.getType()).add(sv);
			}
		}
		
		Collections.sort(this.genesOverlappingWithSV, GenomicRegionUtils.sorterByChromAndStartPos());
	}

	
	/**
	 * 
	 * @param outDir
	 */
	public void writeToFile(Path outDir) {
		this.svTypeList=new ArrayList<>();
		this.svTypeList.addAll(this.svTypeOverlappingGenesMap.keySet());
		
		/////////////////
		this.allSVTypeGeneListFile=Path.of(outDir.toString(), "gene.list.overlapping.".concat(this.SVInfor).concat(".tsv"));
		if(this.allSVTypeGeneListFile.toFile().exists()) {
			System.out.println("output gene list file exists, delete it...");
			this.allSVTypeGeneListFile.toFile().delete();
		}
		
		
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(this.allSVTypeGeneListFile.toString()));

		    writer.append(this.getOutputFileHeaderLine());
		    writer.newLine();
		    
		    for(Gene gene:this.genesOverlappingWithSV) {
		    	StringBuilder sb=new StringBuilder();
		    	sb.append(gene.getName());
		    	
		    	int totalNum=0;
		    	for(SimpleSVType svType:this.svTypeList) {
		    		if(this.geneSVTypeOverlappingSVsMap.get(gene).containsKey(svType)) {
		    			sb.append("\t").append(this.geneSVTypeOverlappingSVsMap.get(gene).get(svType).size());
		    			totalNum+=this.geneSVTypeOverlappingSVsMap.get(gene).get(svType).size();
		    		}else {
		    			sb.append("\t").append(0);
		    		}
		    	}
		    	sb.append("\t").append(totalNum);
		    	
		    	writer.append(sb.toString());
		    	writer.newLine();
		    }
		    
		    writer.flush();
		    writer.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}
	
	
	private String getOutputFileHeaderLine() {
		StringBuilder sb=new StringBuilder();
		
		sb.append("#gene");
		
		for(SimpleSVType type:this.svTypeList) {
			sb.append("\t").append(type.toString());
		}
		
		sb.append("\t").append("total");
		
		return sb.toString();
	}
	
	/**
	 * @return the svTypeOverlappingGenesMap note that some genes may be overlapping with multiple type of SVs
	 */
	public Map<SimpleSVType, Set<Gene>> getSvTypeOverlappingGenesMap() {
		return svTypeOverlappingGenesMap;
	}
	
	
	/**
	 * 
	 * @return
	 */
	public List<Gene> getAllGenesOverlappingWithSV(){
		return this.genesOverlappingWithSV;
	}

	
	/**
	 * @return the svTypeSVsOverlappingWithGeneMap
	 */
	public Map<SimpleSVType, List<SimpleSVLocus>> getSvTypeSVsOverlappingWithGeneMap() {
		return svTypeSVsOverlappingWithGeneMap;
	}
	

	
	/**
	 * @return the allSVTypeGeneListFile
	 */
	public final Path getAllSVTypeGeneListFile() {
		return allSVTypeGeneListFile;
	}

}
