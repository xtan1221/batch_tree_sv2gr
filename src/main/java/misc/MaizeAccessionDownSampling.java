package misc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

/**
 * basic:
 * 		each bio sample is from an individual belonging to certain genotype/taxon;
 * 		each accession is from a bio sample;
 * 		========================
 * 		multiple accessions may come from the same bio sample;
 * 		multiple bio samples may come from the same genotype/taxon
 * 			individuals of the same genotype/taxon
 * 
 * ===============================
 * strategy:
 * 
 * downsampling from the HapMap3.2.1 dataset
 * 
 * 1. for each genotype, select the accession with the largest base num;
 * 
 * 2. 
 * 
 * 
 * @author tanxu
 *
 */
public class MaizeAccessionDownSampling {
	/**
	 * tsv file
	 * first line is header line starting with #
	 */
	private final Path accessionMetadataFile;
	/**
	 * 
	 */
	private final Path outputTSVFile;
	/**
	 * column index of dataset name (dataset)
	 * starting from 0
	 */
	private final int dataSetColIndex;
	/**
	 * column index of accession name (Run)
	 * starting from 0
	 */
	private final int accessionColIndex;
	/**
	 * column index of total number of resequenced bases (Bases)
	 * starting from 0
	 */
	private final int baseNumColIndex;
	
	/**
	 * column index of bio sample name (BioSample)
	 * starting from 0
	 */
	private final int bioSampleColIndex;
	
	/**
	 * column index of genotype (Genotype)
	 * starting from 0
	 */
	private final int genotypeColIndex;
	
	/**
	 * column index of library name (Library Name)
	 * starting from 0
	 */
	private final int libraryNameColIndex;
	
	
	
	////////////////////////////
	private String headerLine;
	/**
	 * 
	 */
	private Map<String, Accession> accessionMap;
	private Map<String, Set<Accession>> sampleNameAccessionSetMap; 
	private Map<String, Set<String>> genotypeNameSampleNameSetMap;
	/**
	 * map from genotype to the Accession with the largest base num
	 */
	private Map<String, Accession> genotypeLargestBaseNumAccession;
	//////////////////////////////
	/**
	 * map from sample name with one single accession to the accession object
	 */
	private Map<String, Accession> singleAccessionSampleNameAccessionMap;
	
	

	
	public MaizeAccessionDownSampling(
			Path accessionMetadataFile, Path outputTSVFile,
			int dataSetColIndex, int accessionColIndex,
			int baseNumColIndex, int bioSampleColIndex, 
			int genotypeColIndex, int libraryNameColIndex) {
		super();
		this.accessionMetadataFile = accessionMetadataFile;
		this.outputTSVFile=outputTSVFile;
		this.dataSetColIndex = dataSetColIndex;
		this.accessionColIndex = accessionColIndex;
		this.baseNumColIndex = baseNumColIndex;
		this.bioSampleColIndex = bioSampleColIndex;
		this.genotypeColIndex = genotypeColIndex;
		this.libraryNameColIndex = libraryNameColIndex;
	}


	/**
	 * map from the genotype to the accession with largest base num among all the accessions of the single-accession samples of the genotype;
	 * only genotypes with at least one sample with single accession will be included in this map;
	 */
	private Map<String, Accession> genotypeLargestBaseNumAccessionOfSingleAccessionSample;
	
	
	void readFile() {
		this.accessionMap=new LinkedHashMap<>();
		this.sampleNameAccessionSetMap = new LinkedHashMap<>();
		this.genotypeNameSampleNameSetMap=new LinkedHashMap<>();
		this.genotypeLargestBaseNumAccession=new LinkedHashMap<>();
		try {
		    BufferedReader lineReader = new BufferedReader(new FileReader(this.accessionMetadataFile.toFile()));
		    String line = null;
		    String[] splits;
		    while ((line = lineReader.readLine()) != null) {
//		        System.out.println(line);
		        if(line.startsWith("#")) {
		        	this.headerLine=line;
		        	continue;
		        }
		        
		        splits=line.split("\t");
		        
		    	/**
		    	 *  dataset name (dataset)
		    	 */
		    	String dataSet = splits[this.dataSetColIndex];
		    	/**
		    	 *  accession name (Run)
		    	 */
		    	String accession = splits[this.accessionColIndex];
		    	/**
		    	 *  total number of resequenced bases (Bases)
		    	 */
		    	long baseNum = Double.valueOf(splits[this.baseNumColIndex]).longValue();

		    	/**
		    	 *  bio sample name (BioSample)
		    	 */
		    	String bioSample = splits[this.bioSampleColIndex];
		    	
		    	/**
		    	 *  genotype (Genotype)
		    	 */
		    	String genotype = splits[this.genotypeColIndex];
		    	
		    	/**
		    	 *  library name (Library Name)
		    	 */
		    	String libraryName = splits[this.libraryNameColIndex];
		    	
		    	
		    	Accession acc = new Accession(line);
		    	acc.setAccession(accession);
		        acc.setBaseNum(baseNum);
		        acc.setBioSample(bioSample);
		        acc.setDataSet(dataSet);
		        acc.setGenotype(genotype);
		        acc.setLibraryName(libraryName);
		        
		       
		        this.accessionMap.put(accession, acc);
		        
		        
		        if(!this.sampleNameAccessionSetMap.containsKey(bioSample)) {
		        	this.sampleNameAccessionSetMap.put(bioSample, new LinkedHashSet<>());
		        }
		        
		        this.sampleNameAccessionSetMap.get(bioSample).add(acc);
		        
		        if(!this.genotypeNameSampleNameSetMap.containsKey(genotype)) {
		        	this.genotypeNameSampleNameSetMap.put(genotype, new LinkedHashSet<>());
		        }
		        this.genotypeNameSampleNameSetMap.get(genotype).add(bioSample);
		        
		        
		        if(!this.genotypeLargestBaseNumAccession.containsKey(acc.getGenotype())) {
		        	this.genotypeLargestBaseNumAccession.put(acc.getGenotype(), acc);
		        }else {
		        	if(this.genotypeLargestBaseNumAccession.get(acc.getGenotype()).getBaseNum()<acc.getBaseNum()) {
		        		this.genotypeLargestBaseNumAccession.put(acc.getGenotype(), acc);
		        	}
		        }
		    }
		    
		    
		    lineReader.close();
		} catch (IOException ex) {
			
		    System.err.println(ex);
		}
	}
	
	
	
	void process() {
		this.singleAccessionSampleNameAccessionMap = new LinkedHashMap<>();
		this.genotypeLargestBaseNumAccessionOfSingleAccessionSample=new LinkedHashMap<>();
		
		//first identify bio sample with one single accessions
		this.sampleNameAccessionSetMap.forEach((sample, accSet)->{
			if(accSet.size()==1) {//only process the sample with one single accession
				Accession acc=accSet.iterator().next();
				this.singleAccessionSampleNameAccessionMap.put(sample, acc);
				
				if(!this.genotypeLargestBaseNumAccessionOfSingleAccessionSample.containsKey(acc.getGenotype())) {
					this.genotypeLargestBaseNumAccessionOfSingleAccessionSample.put(acc.getGenotype(), acc);
				}else {//if the genotype has already assigned an accession, check its base num with the current one
					if(this.genotypeLargestBaseNumAccessionOfSingleAccessionSample.get(acc.getGenotype()).getBaseNum()<acc.getBaseNum()) {
						this.genotypeLargestBaseNumAccessionOfSingleAccessionSample.put(acc.getGenotype(), acc);
					}
				}
			}
		});
		
		
		//
	}
	
	
	void print() {
		System.out.println("hello");
//		this.genotypeLargestBaseNumAccessionOfSingleAccessionSample.forEach((gt, acc)->{
//			System.out.println(acc.getOriginalDataLine());
//		});
		System.out.println(this.genotypeLargestBaseNumAccession.size());
		this.genotypeLargestBaseNumAccession.forEach((gt, acc)->{
			System.out.println(acc.getOriginalDataLine());
		});
	}
	
	void outputToFile() {
		String output = this.headerLine;
		
		for(String genotype:this.genotypeLargestBaseNumAccession.keySet()){
			output=output.concat("\n").concat(this.genotypeLargestBaseNumAccession.get(genotype).getOriginalDataLine());
		}
		
		
	    try {
	    	BufferedWriter writer = new BufferedWriter(new FileWriter(this.outputTSVFile.toFile()));
			writer.write(output);
			 writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	   
	}
}
