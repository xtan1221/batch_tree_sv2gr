package population.utils;

/**
 * genotype of a locus;
 * 
 * note that 1 indicate the presence of a r gene while 0 indicate absence of a r gene
 * 
 * @author tanxu
 *
 */
public enum Genotype {

	MISSING ("./."),
	/**
	 * for SNP and SV, this simply indicates homozygous refernece genotype; 
	 * for signature based variant locus such as MEI and R gene, this indicates the total absence of the target sequence regardless of the reference genotype!
	 */
	ABSENCE ("0/0"), 
	/**
	 * 
	 */
	HETER ("0/1"),
	/**
	 * for SNP and SV, this simply indicates homozygous refernece genotype, this simply indicates homozygous non-reference genotype
	 * for signature based variant locus such as MEI and R gene, this indicates the homozygous genotype of the presence of the target sequence regardless of the reference genotype!
	 */
	PRESENCE("1/1");
	
	private final String stringValue;
	
	Genotype(String stringValue){
		this.stringValue=stringValue;
	}

	/**
	 * @return the stringValue
	 */
	public String getStringValue() {
		return stringValue;
	}
		
}
