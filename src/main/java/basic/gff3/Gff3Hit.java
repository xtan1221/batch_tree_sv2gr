package basic.gff3;

import java.util.LinkedHashMap;
import java.util.Map;

import genomics.utils.SimpleGenomicRegion;
import genomics.utils.Strand;

/**
 * one record from the gff3 file containing the gene annotation of a reference genome;
 * 
 * the details about gff3 file format are found in 
 * 		https://m.ensembl.org/info/website/upload/gff3.html
 * 
 * @author tanxu
 *
 */
public class Gff3Hit extends SimpleGenomicRegion{
//	/**
//	 * name of the chromosome or scaffold;
//	 * column 1 (starting from 1)
//	 * chromosome names can be given with or without the 'chr' prefix. Important note: the seq ID must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
//	 */
//	private final String ref;
	
	/**
	 * type of feature. Must be a term or accession from the SOFA sequence ontology
	 * column 3 (starting from 1)
	 */
	private final String type;
	
//	/**
//	 * Start position of the feature, with sequence numbering starting at 1.
//	 * column 4 (starting from 1)
//	 */
//	private final int start;
//	
//	/**
//	 * End position of the feature, with sequence numbering starting at 1.
//	 * column 5 (starting from 1)
//	 */
//	private final int end;
	
	/**
	 * A floating point value.
	 * column 6 (starting from 1)
	 * can be null
	 */
	private final Double score;
	
//	/**
//	 * defined as + (forward) or - (reverse).
//	 * column 7 (starting from 1)
//	 */
//	private final String strand;

	/**
	 * One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
	 * column 8 (starting from 1)
	 * can be null
	 */
	private final Integer phase;
	
	/**
	 * 
	 */
	private final Map<String, String> attributeNameStringValueMap;
	
	/**
	 * full constructor
	 * @param chrom
	 * @param start
	 * @param end
	 * @param strand
	 * @param type
	 * @param score
	 * @param phase
	 * @param attributeNameStringValueMap
	 */
	public Gff3Hit(
			String chrom, int start, int end, Strand strand, 
			String type, Double score, Integer phase,
			Map<String, String> attributeNameStringValueMap) {
		super(chrom, start, end, strand);
		
		this.type = type;
		this.score = score;
		this.phase = phase;
		this.attributeNameStringValueMap = attributeNameStringValueMap;
	}


	
	/**
	 * factory method to build a Gff3Hit from the given data line from a gff3 file
	 * @param dataLine
	 * @return
	 */
	public static Gff3Hit build(String dataLine) {
		String[] splits=dataLine.trim().split(Gff3Utils.COLUMN_DELIMITER);
		String ref=splits[0];
		String type=splits[2];
		int start=Integer.parseInt(splits[3]);
		int end=Integer.parseInt(splits[4]);
		Double score=splits[5].equals(Gff3Utils.NULL_VALUE_STRING)?null:Double.parseDouble(splits[5]);
		Strand strand=splits[6].equals(Gff3Utils.NULL_VALUE_STRING)?Strand.UNKNOWN:Strand.find(splits[6]);
		Integer phase=splits[7].equals(Gff3Utils.NULL_VALUE_STRING)?null:Integer.parseInt(splits[7]);
		
		//
		Map<String, String> attributeNameStringValueMap=new LinkedHashMap<>();
		
		String attributesString=splits[8];
		String[] attributes=attributesString.split(Gff3Utils.ATTRIBUTE_COLUMN_DELIMETER);
		
		for(String attribute:attributes) {
			String[] attributeNameStringValueSplits=attribute.split(Gff3Utils.SINGLE_ATTRIBUTE_NAME_VALUE_LINKER);
			attributeNameStringValueMap.put(attributeNameStringValueSplits[0], attributeNameStringValueSplits[1]);
		}
		
		return new Gff3Hit(ref, start, end, strand,
				type, score, phase,
				attributeNameStringValueMap
				);
	}
	
	
	/**
	 * @return name of the chromosome or scaffold;
	 * @Override
	 */
	public String getChrom() {
		return super.getChrom();
	}
	
	
	/**
	 * @return the type
	 */
	public String getType() {
		return type;
	}

	
	/**
	 * @return the start
	 * @Override
	 */
	public int getStart() {
		return super.getStart();
	}


	/**
	 * @return the end
	 * @Override
	 */
	public int getEnd() {
		return super.getEnd();
	}

	
	/**
	 * @return the score
	 */
	public Double getScore() {
		return score;
	}


	/**
	 * @return the strand
	 * @override
	 */
	public Strand getStrand() {
		return super.getStrand();
	}


	/**
	 * @return the phase
	 */
	public Integer getPhase() {
		return phase;
	}


	/**
	 * @return the attributeNameStringValueMap
	 */
	public Map<String, String> getAttributeNameStringValueMap() {
		return attributeNameStringValueMap;
	}
	
	
}
