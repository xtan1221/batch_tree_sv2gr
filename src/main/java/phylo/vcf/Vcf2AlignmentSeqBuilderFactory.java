package phylo.vcf;

import population.vcf.utils.GenotypeRecoder;
import population.vcf.utils.VariantContextFilter;


/**
 * sequence builder that transform records in a vcf file for each individual/sample into a sequence that can be used for phylogenetic tree inference for the group of individuals/samples
 * 
 * whether a specific site should be included or not is fully decided by the given {@link @VariantContextFilter};
 * 
 * note that the VCF file can contain non-variant sites where individuals' genotype are either homozygous reference or missing;
 * 
 * @author tanxu
 *
 */
public class Vcf2AlignmentSeqBuilderFactory {
	//////////////
	private final GenotypeRecoder genotypeRecoder;
	
	private final VariantContextFilter variantFilter;
	/**
	 * 
	 * @param genotypeRecoder
	 * @param variantFilter
	 */
	public Vcf2AlignmentSeqBuilderFactory(VariantContextFilter variantFilter, GenotypeRecoder genotypeRecoder){
		//
		
		
		this.variantFilter = variantFilter;
		this.genotypeRecoder = genotypeRecoder;
		
	}
	
	public Vcf2AlignmentSeqBuilder make() {
		return new Vcf2AlignmentSeqBuilder(this.variantFilter, this.genotypeRecoder);
	}
}
