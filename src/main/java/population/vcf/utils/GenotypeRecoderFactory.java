package population.vcf.utils;

import htsjdk.variant.variantcontext.Allele;

public class GenotypeRecoderFactory {
	public static final String GAP="-";
	public static final String MISSING_BASE="?"; //consistent with the MegaX
	
	/**
	 * genotype recoder that recodes a genotype to a sequence composed of each allele base consecutively;
	 * applicable for any ploidy;
	 * for example, if the genotype is A/G, it will be recoded to 'AG'
	 * @return
	 */
	public static GenotypeRecoder biBaseSeqRecoder() {
		return p->{
			String ret="";
			for(Allele a:p.getSecond().getAlleles()) {
				if(HtsVCFUtils.isMissing(a)) {
					ret=ret.concat(MISSING_BASE);
				}else {
					ret=ret.concat(a.getBaseString());
				}
			}
			return ret;
		};
	}
	
	/**
	 * recode a genotype
	 * 1. if genotype is missing, encode to the string of MISSING_BASE with the same length as the reference seq
	 * 2. if genotype is not missing
	 * 		1. if homozygous, trivial
	 * 		2. if the genotype is heterozygous, a random allele will be selected (normally the first one)
	 * 
	 * 
	 * @return
	 */
	public static GenotypeRecoder singleBaseSeqRecoder() {
		return p ->{
			return HtsVCFUtils.isMissing(p.getSecond().getAllele(0))?makeMissingBaseSeq(p.getFirst().getReference().getBaseString().length()):p.getSecond().getAllele(0).getBaseString();
		};
	}
	
	static String makeMissingBaseSeq(int len) {
		if(len<=0)
			throw new IllegalArgumentException("given len must be positive integer!");
		
		String ret = "";
		for(int i=0;i<len;i++)
			ret=ret.concat(MISSING_BASE);
		
		return ret;
	}
	
}
