package population.popHierarchy;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * reader and lookup for population information of samples of the same taxa under study;
 * 
 * @author tanxu
 * 
 */
public class PopulationStructureFileReader {
	/**
	 * four column space delimited table file containing the population relationship of the samples being studied
	 * 
	 * col1 = sample index in the related files
	 * 		the samples are ordered by the index
	 * 		the first sample index can either be 0 or 1 based on the specific data source
	 * 			for example, for {@link VariantContext} type data, the first index should be 0
	 * 			for {@link SingleLocusBase} type data, the first index should be 1
	 * 
	 * col2 = level or order, 
	 * 		0==> ingroup samples, for example, sorghum bicolor samples/accessions
	 * 		1, 2, 3, ... => outgroup, where the samples with the same level are thought to be from same outgroup taxon; sample with higher level is more distantly related than those with smaller level
	 * 				for example, sorghum propinquum has level = 1 because it is more closely related with sorghum bicolor;
	 * 						sorghum versicolor and sorghum timonrense have level =2 because they are from the same outgroup taxa and are more distantly related with sorghum bicolor than sorghum propinquum;
	 * 
	 * col3 = sample name
	 * 
	 * col4 = taxon
	 * 		a more detailed classification of each sample in addition to the population level; for example, 'indica', 'japonica' for rice samples with pop level=0
	 */
	private final Path samplePopulationInfoTableFile;
	
	/**
	 * the first sample index can either be 0 or 1 based on the specific data source
	 * 			for example, for {@link VariantContext} type data, the first index should be 0
	 * 			for {@link SingleLocusBase} type data, the first index should be 1
	 * 
	 * if true, the smallest sample index should be 0;
	 * if false, the smallest sample index should be 1;
	 */
	private final boolean firstIndexIs0;
	
	/////////////////////////////
	private Map<Integer, Sample> sampleIndexMap;
	
	////////////////////////////////
	/**
	 * 
	 */
	private Map<Integer, String> sampleIndexSampeNameMap;
	/**
	 * 
	 */
	private Map<Integer, Integer> sampleIndexPopulationLevelMap;
	
	/**
	 * 
	 */
	private Map<Integer, List<Integer>> populationLevelSampleIndicesMap;
	
	/**
	 * 
	 */
	private List<Integer> orderedSampleIndices;
	
	private Map<Integer, Integer> outgroupSampleIndexPopLevelMap;
	
	private Map<String, List<Integer>> taxonNameSampleIndicesMap;
	private Map<String, List<String>> taxonNameSampleNamesMap;
	
	public PopulationStructureFileReader(Path samplePopulationInfoTableFile, boolean firstIndexIs0) {
		super();
		this.samplePopulationInfoTableFile = samplePopulationInfoTableFile;
		this.firstIndexIs0=firstIndexIs0;
		
		///////////////////////
		this.readFile();
		//check the sample index
		if(this.firstIndexIs0 && !this.sampleIndexMap.containsKey(0)) {
			throw new IllegalArgumentException("no sample with index 0 is found!");
		}
	}

	

	void readFile() {
		this.sampleIndexMap=new HashMap<>();
		this.sampleIndexSampeNameMap = new HashMap<>();
		this.sampleIndexPopulationLevelMap = new HashMap<>();
		this.populationLevelSampleIndicesMap = new HashMap<>();
		this.orderedSampleIndices = new ArrayList<>();
		
		try {
			////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(this.samplePopulationInfoTableFile.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty() || line.startsWith("#")) { //header line starts with #
					continue;
				}
				
				//#sample_index	population_level	sample_name
				//1		1	propinquum_4
				String[] splits=line.split("\\s+");
				
				int sampleIndex = Integer.parseInt(splits[0]);
				int populationLevel=Integer.parseInt(splits[1]);
				String sampleName=splits[2];
				String taxon=splits[3];
				
				Sample sample=new Sample(sampleIndex, populationLevel, sampleName, taxon);
				this.sampleIndexMap.put(sampleIndex, sample);
				
				if(this.sampleIndexSampeNameMap.containsKey(sampleIndex)) {
					lineReader.close();
					throw new IllegalArgumentException("duplicate sample index found:"+sampleIndex);
				}
				if(this.sampleIndexSampeNameMap.containsValue(sampleName)) {
					lineReader.close();
					throw new IllegalArgumentException("duplicate sample name found:"+sampleName);
				}
				
				this.sampleIndexSampeNameMap.put(sampleIndex, sampleName);
				this.sampleIndexPopulationLevelMap.put(sampleIndex, populationLevel);
				
				if(!this.populationLevelSampleIndicesMap.containsKey(populationLevel)) {
					this.populationLevelSampleIndicesMap.put(populationLevel, new ArrayList<>());
				}
				this.populationLevelSampleIndicesMap.get(populationLevel).add(sampleIndex);
				
				this.orderedSampleIndices.add(sampleIndex);
			}
			
			lineReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		Collections.sort(this.orderedSampleIndices);
	}


	//////////////////////////////
	/**
	 * lookup and return the sample name of the given sample index
	 * @return the sampleIndexSampeNameMap
	 */
	public String lookupSampeNameMap(int sampleIndex) {
		return this.sampleIndexSampeNameMap.get(sampleIndex);
	}
	
	
	/**
	 * lookup and return the population level of the given sample index
	 * 
	 * 
	 * @return the populationLevelSampleIndicesMap
	 */
	public int lookupPopulationLevel(int sampleIndex) {
		return this.sampleIndexPopulationLevelMap.get(sampleIndex);
	}
	
	
	/**
	 * @return the sampleIndexPopulationLevelMap
	 */
	public Map<Integer, Integer> getSampleIndexPopulationLevelMap() {
		return sampleIndexPopulationLevelMap;
	}



	/**
	 * 
	 * @return
	 */
	public List<Integer> getAllIngroupSampleIndices(){
		return this.populationLevelSampleIndicesMap.get(0);
	}
	
	/**
	 * 
	 * @return
	 */
	public List<Integer> getAllOutgroupSampleIndices(){
		List<Integer> ret = new ArrayList<>();
		
		for(int level:this.populationLevelSampleIndicesMap.keySet()) {
			if(level>0){
				ret.addAll(this.populationLevelSampleIndicesMap.get(level));
			}
		}
		
		return ret;
	}
	
	/**
	 * 
	 * @param populationLevel
	 * @return
	 */
	public List<Integer> getAllSampleIndicesOfPopulationLevel(int populationLevel){
		return this.populationLevelSampleIndicesMap.get(populationLevel);
	}
	
	
	/**
	 * return the ordered list of index of all samples
	 * @return the orderedSampleIndices
	 */
	public List<Integer> getOrderedSampleIndices() {
		return orderedSampleIndices;
	}
	
	/**
	 * return the list of sample names ordered by the index
	 * @return
	 */
	public List<String> getOrderedSampleNameList(){
		List<String> ret=new ArrayList<>();
		
		for(int i:this.getOrderedSampleIndices()) {
			ret.add(this.sampleIndexMap.get(i).getName());
		}
		
		return ret;
	}
	
	public List<String> getOutgroupSampleNames(){
		List<String> ret=new ArrayList<>();
		
		for(int i:this.sampleIndexMap.keySet()) {
			if(this.sampleIndexMap.get(i).getPopulationLevel()>0) {
				ret.add(this.sampleIndexMap.get(i).getName());
			}
		}
		
		return ret;
	}

	
	public List<String> getIngroupSampleNames(){
		List<String>  ret = new ArrayList<>();
		
		for(String sample:this.getSampleNameIndexMap().keySet()) {
			if(this.sampleIndexPopulationLevelMap.get(this.getSampleNameIndexMap().get(sample))==0) {
				ret.add(sample);
			}
		}
		
		return ret;
	}
	
	/**
	 * 
	 * @return
	 */
	public Map<String, Integer> getSampleNameIndexMap(){
		Map<String, Integer> ret = new HashMap<>();
		
		for(int index:this.sampleIndexSampeNameMap.keySet()) {
			ret.put(this.sampleIndexSampeNameMap.get(index), index);
		}
		
		return ret;
	}
	
	/**
	 * 
	 * @return
	 */
	public Map<Integer, Integer> getOutgroupSampleIndexPopLevelMap(){
		if(this.outgroupSampleIndexPopLevelMap==null) {
			this.outgroupSampleIndexPopLevelMap=new HashMap<>();
			
			for(int index:this.getAllOutgroupSampleIndices()) {
				this.outgroupSampleIndexPopLevelMap.put(index, this.sampleIndexPopulationLevelMap.get(index));
			}
		}
		
		return this.outgroupSampleIndexPopLevelMap;
	}
	
	/**
	 * 
	 * @param taxon
	 * @return
	 */
	public List<Integer> getSampleIndicesByTaxon(String taxon){
		if(this.taxonNameSampleIndicesMap==null) {
			this.taxonNameSampleIndicesMap=new HashMap<>();
			
			for(int index:this.sampleIndexMap.keySet()) {
				Sample sample=this.sampleIndexMap.get(index);
				
				if(!this.taxonNameSampleIndicesMap.containsKey(sample.getTaxon())) {
					this.taxonNameSampleIndicesMap.put(sample.getTaxon(), new ArrayList<>());
				}
				
				this.taxonNameSampleIndicesMap.get(sample.getTaxon()).add(sample.getIndex());
			}
		}
		
		return this.taxonNameSampleIndicesMap.get(taxon);
	}
	
	/**
	 * 
	 * @param taxon
	 * @return
	 */
	public List<String> getSampleNamesByTaxon(String taxon){
		if(this.taxonNameSampleNamesMap==null) {
			this.taxonNameSampleNamesMap=new HashMap<>();
			
			for(int index:this.sampleIndexMap.keySet()) {
				Sample sample=this.sampleIndexMap.get(index);
				
				if(!this.taxonNameSampleNamesMap.containsKey(sample.getTaxon())) {
					this.taxonNameSampleNamesMap.put(sample.getTaxon(), new ArrayList<>());
				}
				
				this.taxonNameSampleNamesMap.get(sample.getTaxon()).add(sample.getName());
			}
		}
		
		return this.taxonNameSampleNamesMap.get(taxon);
	}
	
	/**
	 * 
	 * @param sampleIndices
	 * @return
	 */
	public List<String> getSampleNames(List<Integer> sampleIndices){
		List<String> ret = new ArrayList<>();
		
		for(int index:sampleIndices) {
			ret.add(this.sampleIndexMap.get(index).getName());
		}
		
		return ret;
	}
}
