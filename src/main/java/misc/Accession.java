package misc;

public class Accession {
	private final String originalDataLine;
	
	////////////////////////
	/**
	 *  dataset name (dataset)
	 */
	private String dataSet;
	/**
	 *  accession name (Run)
	 */
	private String accession;
	/**
	 *  total number of resequenced bases (Bases)
	 */
	private long baseNum;
	
	/**
	 *  bio sample name (BioSample)
	 */
	private String bioSample;
	
	/**
	 *  genotype (Genotype)
	 */
	private String genotype;
	
	/**
	 *  library name (Library Name)
	 */
	private String libraryName;
	
	
	Accession(String originalDataLine){
		this.originalDataLine=originalDataLine;
	}

	////////////////////////////////////
	/**
	 * @return the dataSet
	 */
	public String getDataSet() {
		return dataSet;
	}


	/**
	 * @param dataSet the dataSet to set
	 */
	public void setDataSet(String dataSet) {
		this.dataSet = dataSet;
	}


	/**
	 * @return the accessionColIndex
	 */
	public String getAccession() {
		return accession;
	}


	/**
	 * @param accession the accessionColIndex to set
	 */
	public void setAccession(String accession) {
		this.accession = accession;
	}


	/**
	 * @return the baseNum
	 */
	public long getBaseNum() {
		return baseNum;
	}


	/**
	 * @param baseNum the baseNum to set
	 */
	public void setBaseNum(long baseNum) {
		this.baseNum = baseNum;
	}


	/**
	 * @return the bioSample
	 */
	public String getBioSample() {
		return bioSample;
	}


	/**
	 * @param bioSample the bioSample to set
	 */
	public void setBioSample(String bioSample) {
		this.bioSample = bioSample;
	}


	/**
	 * @return the genotype
	 */
	public String getGenotype() {
		return genotype;
	}


	/**
	 * @param genotype the genotype to set
	 */
	public void setGenotype(String genotype) {
		this.genotype = genotype;
	}


	/**
	 * @return the libraryName
	 */
	public String getLibraryName() {
		return libraryName;
	}


	/**
	 * @param libraryName the libraryName to set
	 */
	public void setLibraryName(String libraryName) {
		this.libraryName = libraryName;
	}


	/**
	 * @return the originalDataLine
	 */
	public String getOriginalDataLine() {
		return originalDataLine;
	}
	
}
