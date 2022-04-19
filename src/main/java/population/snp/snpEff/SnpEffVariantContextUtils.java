package population.snp.snpEff;

import java.util.ArrayList;
import java.util.List;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * 
 * @author tanxu
 * 
 */
public class SnpEffVariantContextUtils {
	//a list of possible annotation for SNP by SnpEff
	//	downstream_gene_variant
	//	initiator_codon_variant
	//	initiator_codon_variant&splice_region_variant
	//	intergenic_region
	//	intron_variant
	//	missense_variant =======================================================NON-Synonymous
	//	missense_variant&splice_region_variant=======================================================NON-Synonymous
	//	non_coding_transcript_exon_variant
	//	splice_acceptor_variant&intron_variant
	//	splice_acceptor_variant&splice_region_variant&intron_variant
	//	splice_donor_variant&intron_variant
	//	splice_donor_variant&splice_region_variant&intron_variant
	//	splice_region_variant
	//	splice_region_variant&intron_variant
	//	splice_region_variant&non_coding_transcript_exon_variant
	//	splice_region_variant&stop_retained_variant
	//	splice_region_variant&synonymous_variant
	//	start_lost
	//	start_lost&splice_region_variant
	//	stop_gained
	//	stop_gained&splice_region_variant
	//	stop_lost
	//	stop_lost&splice_region_variant
	//	stop_retained_variant
	//	synonymous_variant =======================
	//	upstream_gene_variant
	//	3_prime_UTR_variant
	//	5_prime_UTR_premature_start_codon_gain_variant
	//	5_prime_UTR_variant

	
	/**
	 * retrieve and return the annotated type of the given SNP site
	 * 
	 * ============================
	 * first retrieve the 'ANN' attribute value from the INFO column
	 * the types of the SNP is at the second entry ('upstream_gene_variant' in the following example)
	 * ;ANN=T|upstream_gene_variant|MODIFIER|SORBI_3001G000100|SORBI_3001G000100|transcript|EER90453|protein_coding||c.-720C>T|||||720|,T|intergenic_region|MODIFIER|CHR_START-SORBI_3001G000100|CHR_START-SORBI_3001G000100|intergenic_region|CHR_START-SORBI_3001G000100|||n.1231C>T||||||
	 * 
	 * sometimes a SNP is annotated as multiple types as following
	 * 		splice_acceptor_variant&splice_region_variant&intron_variant
	 * 
	 * @return
	 */
	public static List<String> getSNPAnnotatedTypes(VariantContext vc) {
		
//		String s = vc.getCommonInfo()..toString();
		String annotationString = vc.getCommonInfo().getAttribute("ANN").toString();
		String[] splits=annotationString.split("\\|");
		String[] splits2 = splits[1].split("&");
		
		List<String> ret = new ArrayList<>();
		
		for(String type:splits2) {
			if(!type.trim().isEmpty())
				ret.add(type);
		}
		
		return ret;
	}
	
	/**
	 * return whether the given SNP site is a synonymous mutation site or not 
	 * 
	 * annotated with 'synonymous_variant' by SnpEff
	 * =====================
	 * 
	 * @return
	 */
	public static boolean isSynSite(VariantContext vc) {
		return getSNPAnnotatedTypes(vc).contains("synonymous_variant");
	}
	
	
	/**
	 * return whether the given SNP site is a non-synonymous mutation site or not 
	 * annotated with 'missense_variant' by SnpEff
	 * @param vc
	 * @return
	 */
	public static boolean isNonSynSite(VariantContext vc) {
		return getSNPAnnotatedTypes(vc).contains("missense_variant");
	}
	
	
	
	
}
