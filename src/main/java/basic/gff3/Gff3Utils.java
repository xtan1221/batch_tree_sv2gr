package basic.gff3;

public class Gff3Utils {
	/**
	 * Fields must be tab-separated for GFF3 file
	 */
	public static final String COLUMN_DELIMITER = "\t";
	
	/**
	 * 
	 */
	public static final String NULL_VALUE_STRING=".";
	
	/**
	 * delimiter for attributes in the attribute column;
	 * for example, the ';' in following line:
	 * 		ID=Sobic.001G000100.v3.2;Name=Sobic.001G000100;ancestorIdentifier=Sobic.001G000100.v2.1
	 */
	public static final String ATTRIBUTE_COLUMN_DELIMETER=";";
	
	/**
	 * linker string for attribute name and string value;
	 * for example, the "=" in following
	 * 		ID=Sobic.001G000100.v3.2
	 */
	public static final String SINGLE_ATTRIBUTE_NAME_VALUE_LINKER="=";
	
	//////////////////////////////////
	public static final String GENE_TYPE="gene";
	public static final String TRANSCRIPT_TYPE="mRNA";
	public static final String EXON_TYPE="exon";
	public static final String FIVE_PRIME_UTR_TYPE="five_prime_UTR";
	public static final String THREE_PRIME_UTR_TYPE="three_prime_UTR";
	public static final String CDS_TYPE="CDS";
}
