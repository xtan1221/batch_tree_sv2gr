package basic.vcf;

import java.util.HashSet;
import java.util.Set;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public class VariantContextUtils {
	
	/**
	 * check if the given variant is biallelic or not;
	 * 
	 * this method is only implemented for SNP and indel type variant;
	 * 
	 * note that for sv type variant, they are always biallelic ?? (further check)
	 * @param vc
	 * @param toIgnoreMissingDueToOverlappingDeletion
	 * @return
	 */
	public static boolean isBiallelic(VariantContext vc, boolean toIgnoreMissingDueToOverlappingDeletion) {
		if(vc.isSNP()) {
			return isBiallelicSNP(vc, toIgnoreMissingDueToOverlappingDeletion);
		}else if(vc.isIndel()) {
			return isBiallelicIndel(vc, toIgnoreMissingDueToOverlappingDeletion);
		}else {
			return false;
//			throw new UnsupportedOperationException("given VariantContext type is not supported!");
		}
	}
	
	
	/**
	 * check and return whether the given VariantContext contains exactly two types of alleles, one is the reference, the other is the alternative 
	 * 
	 * ========
	 * note that the '*' allele is reserved to indicate that the allele is missing due to a an overlapping deletion (from VCF format 4.3).
	 * 
	 * note that the {@link VariantContext#isBiallelic()} method will take the missing alternative allele '*' as a valid allele, 
	 * thus for a site with reference being A, alternative as [T, *], the site will not be considered as biallelic
	 * 
	 * =========
	 * this method has only been tested for SNP site
	 * 
	 * @param vc
	 * @param toIgnoreMissingDueToOverlappingDeletion if true, ignore the '*' type allele, false otherwise
	 * @return
	 */
	public static boolean isBiallelicSNP(VariantContext vc, boolean toIgnoreMissingDueToOverlappingDeletion) {
		if(!vc.isSNP()) {
			throw new IllegalArgumentException("given variant is not SNP!");
		}
		
		if(toIgnoreMissingDueToOverlappingDeletion) {
			Set<Allele> alleles= new HashSet<>();
			for(Allele allele:vc.getAlternateAlleles()) {
				if(allele.getBaseString().equals("*")) {
					
				}else {
					alleles.add(allele);
				}
			}
			
			alleles.add(vc.getReference());
			
			return alleles.size()==2;
			
		}else {
			return vc.isBiallelic();
		}
		
	}
	
	/**
	 * check and return whether the given VariantContext contains exactly two types of alleles, one is the reference, the other is the alternative 
	 * 
	 * ========
	 * note that the '*' allele is reserved to indicate that the allele is missing due to a an overlapping deletion (from VCF format 4.3).
	 * 
	 * note that the {@link VariantContext#isBiallelic()} method will take the missing alternative allele '*' as a valid allele, 
	 * thus for a site with reference being A, alternative as [T, *], the site will not be considered as biallelic
	 * 
	 * =========
	 * this method has only been tested for SNP site
	 * 
	 * @param vc
	 * @param toIgnoreMissingDueToOverlappingDeletion if true, ignore the '*' type allele, false otherwise
	 * @return
	 */
	public static boolean isBiallelicIndel(VariantContext vc, boolean toIgnoreMissingDueToOverlappingDeletion) {
		if(!vc.isIndel()) {
			throw new IllegalArgumentException("given variant is not indel!");
		}
		
		if(toIgnoreMissingDueToOverlappingDeletion) {
			Set<Allele> alleles= new HashSet<>();
			for(Allele allele:vc.getAlternateAlleles()) {
				if(allele.getBaseString().equals("*")) {
					
				}else {
					alleles.add(allele);
				}
			}
			
			alleles.add(vc.getReference());
			
			return alleles.size()==2;
			
		}else {
			return vc.isBiallelic();
		}
		
	}
}
