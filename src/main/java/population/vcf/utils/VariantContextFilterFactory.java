package population.vcf.utils;

import htsjdk.variant.variantcontext.Genotype;

public class VariantContextFilterFactory {
	
	/**
	 * check if the site is variant site but not biallelic
	 * @return
	 */
	public static VariantContextFilter multiAllelicSite() {
		return e->{
			return e.isVariant() && !e.isBiallelic();
		};
	}
	
	public static VariantContextFilter nonVariantSite() {
		return e->{
			return !e.isVariant();
		};
	}
	
	/**
	 * 
	 * @return
	 */
	public static VariantContextFilter onlySNPSites() {
		return e->{
			return e.isSNP();
		};
	}
	
	/**
	 * return true if the 'FILTER' column does NOT contain 'LowQual'
	 * 
	 * @return
	 */
	public static VariantContextFilter notLowQual() {
		return e->{
			return !e.getFilters().contains("LowQual");
		};
	}
	

	public static VariantContextFilter notBiallelic() {
		return e->{
			return !e.isBiallelic();
		};
	}
	
	
	
	public static VariantContextFilter biallelicSite() {
		return e->{
			return e.isBiallelic();
		};
	}
	
	/**
	 * filter that only affect variant sites;
	 * if the site is not variant (all individual are either homo ref or missing), return true;
	 * else, only return true if the QUAL is not less than the given minQual
	 * @param minQual
	 * @return
	 */
	public static VariantContextFilter filterVariantSitesByQual(int minQual) {
		return e->{
			if(e.isVariant()) {
//				boolean a = e.getPhredScaledQual()>=minQual;
				return e.getPhredScaledQual()>=minQual;
			}else {
				return true;
			}
		};
	}
	
	/**
	 * only allow sites with genotype of every individual being either missing data (no call) or homozygous (reference or alternative)
	 * 
	 * @return
	 */
	public static VariantContextFilter allIndividualsAreGenericHomozygous() {
		return e->{
			for(Genotype gt:e.getGenotypesOrderedByName()) {
				if(!gt.isNoCall()&&!gt.isHom()) {
					return false;
				}
			}
			return true;
		};
	}
	
	
	/**
	 * check if genotype of all individuals are not indel;
	 * 
	 * note that it is allowed for some individuals to have missing genotype!
	 * 
	 * @return
	 */
	public static VariantContextFilter nonIndelSite() {
		return e->{
//			if(e.getStart()==154479) {
//				
//				Type t=e.getType();
//				
//				System.out.println(t);
//			}
			return !e.isIndel();
		};
	}
	
//	static int count=0;
	/**
	 * check if the variant site is not mixed type (containing both SNP and INDEL)
	 * 
	 * for example
	 * 		Chr1	42263	.	C	CT,*
	 * 		Chr1	48086	.	GCCCTT	ACCCTT,G
	 * 
	 * @return
	 */
	public static VariantContextFilter nonMixedTypeSite() {
		return e->{
//			if(e.isMixed()) {
//				count++;
//				System.out.println(count);
//			}
			return !e.isMixed();
		};
	}
	
	/**
	 * return a VariantContextFilter that will filter out VariantContext with number of sample with missing genotype exceeding the given value;
	 * 
	 * @param maxMissing
	 * @return
	 */
	public static VariantContextFilter maxMissingCount(int maxMissing) {
		return e->{
			int missingGenotypes=0;
			for(Genotype gt:e.getGenotypesOrderedByName()) {
				if(HtsVCFUtils.isMissing(gt)) {
					missingGenotypes++;
					if(missingGenotypes>maxMissing) {
						return false;
					}
				}
			}
			
			return true;
//			return e.isBiallelic();
		};
	}
	
	/**
	 * only allow sites that are homozygous for ALL individuals whose genotypes are not missing
	 * Equivalently, all sites with at least one individual's genotype being heterozygous will be filtered out
	 * 
	 * UNTESTED
	 * @return
	 */
	public static VariantContextFilter onlyAllowSitesWithHomozygousIndividuals() {
		return e->{
			for(Genotype gt:e.getGenotypesOrderedByName()) {
				if(HtsVCFUtils.isMissing(gt)) {//missing genotype, skip
					
				}else {
					if(!gt.isHom()) {
						return false;
					}
				}
			}
			
			return true;
		};
	}
	
	
	/**
	 * first check whether number of individuals with missing genotype exceed the given maxMissing value;
	 * then check whether all individuals with non-missing genotype are homozygous or not
	 * @param maxMissing
	 * @return
	 */
	public static VariantContextFilter maxMissingAndHomozygousFilter(int maxMissing) {
		return maxMissingCount(maxMissing).and(onlyAllowSitesWithHomozygousIndividuals());
	}
}
