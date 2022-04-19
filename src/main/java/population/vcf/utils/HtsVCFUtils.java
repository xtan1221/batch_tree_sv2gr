package population.vcf.utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;

public class HtsVCFUtils {
	/**
	 * VCF 4.3 specification
	 * 
	 * The '*' allele is reserved to indicate that the allele is missing due to a an overlapping deletion.
	 * 
	 * @param gt
	 * @return
	 */
	public static boolean isMissingDueToOverlapWithDeletion(Genotype gt) {
		return gt.getGenotypeString().equals("*/*")||gt.getGenotypeString().equals("*|*");
	}
	
	public static boolean isMissing(Genotype gt) {
		return isMissingDueToOverlapWithDeletion(gt)||gt.isNoCall();
	}
	
	/**
	 * return whether the given Allele is missing data or not
	 * @param a
	 * @return
	 */
	public static boolean isMissing(Allele a) {
		return a.isNoCall() || a.getBaseString().equals("*");
	}
}
