package population.vcf.utils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;

/**
 * utilities for {@link Allele} in {@link Genotype} in htsjdk
 * 
 * @author tanxu
 *
 */
public class GenotypeAndAlleleUtils {
	
	/**
	 * return whether the given allele is a single nucleotide or not
	 * @param a
	 * @return
	 */
	public static boolean isSingleNucleotide(Allele a) {
		return a.equals(Allele.ALT_A) || a.equals(Allele.ALT_C) || a.equals(Allele.ALT_G)||a.equals(Allele.ALT_T)
				||a.equals(Allele.REF_A) || a.equals(Allele.REF_C) || a.equals(Allele.REF_G)||a.equals(Allele.REF_T);
	}
	
	
	
}
