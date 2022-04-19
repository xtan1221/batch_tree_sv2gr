package population.unfoldedSFS;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * class to store locus number of for all variant types with derived allele proportion in a specific range
 * 
 * the derived allele frequency is the number of derived alleles at a locus, which is ranged from [0, 2n], where n is the number of ingroup samples
 * 
 * 1. note that when 0, there is no derived allele ==> outgroup and ingroup samples are all ancient allele ==> non-variant site/locus and not included in SFS
 * 2. since some sample may have missing data (./.) for some loci, for such loci, the maximal derived allele frequency is 2*m, where m is the the number of ingroup samples with called genotype
 * 3. 
 * 
 * for a population with n ingroup samples, the 2*n bins are 
 * 		(0, 1/2n], (1/2n,2/2n], ......, ((2n-1)/2n, 1] 
 * 
 * the 0 is not included in the first bin because it is non-variant sites/loci, which is not included in SFS construction;
 * 
 * @author tanxu
 */
public class CompositeDerivedAlleleProportionBin{
	/**
	 * the derived allele count of this bin when the full set of ingroup samples have non-missing data;
	 */
	private final int derivedAlleleCountWhenAllIngroupSamplesHaveNonMissingData;
	
	/**
	 * start of the bin,exclusive 
	 * the min value of derived allele proportion for a locus to be included in this bin
	 */
	private final double start;
	/**
	 * end of the bin, inclusive
	 * the max value of derived allele proportion for a locus to be included in this bin
	 */
	private final double end;
	
	/**
	 * map from the type to the number of locus with derived allele frequency falling into the range of this bin;
	 * note that this map only stores locus number for variant type with at least one locus, if none, this map will not contain the key for that type;
	 */
	private final Map<String, Integer> typeLocusNumMap;
	
	/**
	 * 
	 */
	private final Map<String, Integer> typeAncestralAlleleStateOfPresenceTypeLocusNumMap;
	/**
	 * 
	 */
	private final Map<String, Integer> typeAncestralAlleleStateOfAbsenceTypeLocusNumMap;
	
	/**
	 * 
	 * @param start
	 * @param end
	 */
	public CompositeDerivedAlleleProportionBin(int derivedAlleleCountWhenAllIngroupSamplesHaveNonMissingData, double start, double end) {
		super();
		this.derivedAlleleCountWhenAllIngroupSamplesHaveNonMissingData=derivedAlleleCountWhenAllIngroupSamplesHaveNonMissingData;
		this.start = start;
		this.end = end;
		
		this.typeLocusNumMap=new LinkedHashMap<>();
		this.typeAncestralAlleleStateOfPresenceTypeLocusNumMap=new LinkedHashMap<>();
		this.typeAncestralAlleleStateOfAbsenceTypeLocusNumMap=new LinkedHashMap<>();
	}
	
	/**
	 * 
	 * @param type
	 */
	public void addOneToType(String type, boolean ancestralAlleleStateOfPresenceType) {
		if(!this.typeLocusNumMap.containsKey(type)) {
			this.typeLocusNumMap.put(type, 0);
			this.typeAncestralAlleleStateOfPresenceTypeLocusNumMap.put(type, 0);
			this.typeAncestralAlleleStateOfAbsenceTypeLocusNumMap.put(type, 0);
		}
		
		this.typeLocusNumMap.put(type, this.typeLocusNumMap.get(type)+1);
		
		if(ancestralAlleleStateOfPresenceType) {
			this.typeAncestralAlleleStateOfPresenceTypeLocusNumMap.put(type, this.typeAncestralAlleleStateOfPresenceTypeLocusNumMap.get(type)+1);
		}else {
			this.typeAncestralAlleleStateOfAbsenceTypeLocusNumMap.put(type, this.typeAncestralAlleleStateOfAbsenceTypeLocusNumMap.get(type)+1);
		}
		
	}


	/**
	 * @return the typeLocusNumMap
	 */
	public Map<String, Integer> getTypeLocusNumMap() {
		return typeLocusNumMap;
	}
	
	
	
	/**
	 * 
	 */
	public int getLocusNum(String type) {
		if(this.typeLocusNumMap.containsKey(type)) {
			return this.typeLocusNumMap.get(type);
		}else {//this is the case where there is no locus of a specific type with derived allele proportion of this bin is identified, thus the type will never be added to the map key
			return 0;
		}
	}
	
	/**
	 * 
	 * @param type
	 * @return
	 */
	public int getLocusNumWithPresenceTypeAncestralAllele(String type) {
		if(this.typeAncestralAlleleStateOfPresenceTypeLocusNumMap.containsKey(type)) {
			return this.typeAncestralAlleleStateOfPresenceTypeLocusNumMap.get(type);
		}else {
			return 0;
		}
	}
	/**
	 * 
	 * @param type
	 * @return
	 */
	public int getLocusNumWithAbsenceTypeAncestralAllele(String type) {
		if(this.typeAncestralAlleleStateOfAbsenceTypeLocusNumMap.containsKey(type)) {
			return this.typeAncestralAlleleStateOfAbsenceTypeLocusNumMap.get(type);
		}else {
			return 0;
		}
	}
	
	/**
	 * check if the given derived allele proportion is covered by this bin;
	 * 
	 * note that the {@link #start} is exclusive while the {@link #end} is inclusive
	 * @param p
	 * @return
	 */
	public boolean isCovering(double p) {
		return p>this.start&&p<=this.end;
	}
	
	
	
	/**
	 * @return the derivedAlleleCountWithAllIngroupSamplesOfNonMissingData
	 */
	public final int getDerivedAlleleCountWhenAllIngroupSamplesHaveNonMissingData() {
		return derivedAlleleCountWhenAllIngroupSamplesHaveNonMissingData;
	}

	/**
	 * @return the start
	 */
	public double getStart() {
		return start;
	}


	/**
	 * @return the end
	 */
	public double getEnd() {
		return end;
	}
}

