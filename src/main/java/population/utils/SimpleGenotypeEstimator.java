package population.utils;

/**
 * 
 * strategy class to genotype an individual of a locus ({@link PairedSignature}) called by PoPoolationTE2 based on the calculated proportion of read coverage supporting the target sequence such as TE and R gene
 * 
 * @author tanxu
 * 
 */
public class SimpleGenotypeEstimator {
	/**
	 * if the proportion of coverage of a group of signatures supporting the target sequence (a TE or r gene) of a called variant locus for an individual equal or less than this value, 
	 * the signature group will be genotyped as {@link Genotype#ABSENCE}
	 */
	private final double maxAbsenceSupportingCovFrequency;
	
	/**
	 * if the proportion of coverage of a group of signatures supporting the target sequence (a TE or r gene) of a called variant locus for an individual equal or larger than this value, 
	 * the signature group will be genotyped as {@link Genotype#PRESENCE}
	 */
	private final double minPresenceSupportingCovFrequency;
	
	/**
	 * 
	 * @param maxAbsenceSupportingCovFrequency
	 * @param minPresenceSupportingCovFrequency
	 */
	public SimpleGenotypeEstimator(double maxAbsenceSupportingCovFrequency, double minPresenceSupportingCovFrequency) {
		super();
		this.maxAbsenceSupportingCovFrequency = maxAbsenceSupportingCovFrequency;
		this.minPresenceSupportingCovFrequency = minPresenceSupportingCovFrequency;
	}
	
	/**
	 * estimate and return the genotype of the given proportion of reads coverage supporting a specific target sequence
	 * 
	 * @param targetSeqSupportingCovProportion
	 * @return
	 */
	public Genotype estimateGenotype(double targetSeqSupportingCovProportion) {
		if(targetSeqSupportingCovProportion<=this.maxAbsenceSupportingCovFrequency) {
			return Genotype.ABSENCE;
		}else if(targetSeqSupportingCovProportion<this.minPresenceSupportingCovFrequency) {
			return Genotype.HETER;
		}else {
			return Genotype.PRESENCE;
		}
	}

	/**
	 * @return the maxAbsenceSupportingCovFrequency
	 */
	public double getMaxAbsenceSupportingCovFrequency() {
		return maxAbsenceSupportingCovFrequency;
	}

	/**
	 * @return the minPresenceSupportingCovFrequency
	 */
	public double getMinPresenceSupportingCovFrequency() {
		return minPresenceSupportingCovFrequency;
	}
	
	
	
}
