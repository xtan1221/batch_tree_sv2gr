package population.unfoldedSFS;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.function.Predicate;

import basic.Pair;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import population.utils.SingleLocusBase;
import population.vcf.utils.VariantContextUtils;

/**
 * abstract class to build unfolded SFS for a set of biallelic SNP sites of {@link VariantContext} type
 * 
 * all unfolded SFS calculation should be subclass of this class for the following variation types:
 * 		1. SNP and indel (including snpEff annotated SNPs)
 * 		
 * 
 * @author tanxu
 * 
 * @param <T> type of the variant locus, could be a {@link VariantContext} of htsjdk or a subtype of {@link SingleLocusBase} of this package
 */
public class VariantContextBiallelicSNPUnfoldedSFSBuilder extends AbstractUnfoldedSFSBuilder<VariantContext>{
	
	private Allele currentLocusPresenceTypeAllele;
	private Allele currentLocusAbsenceTypeAllele;
	
	
	public VariantContextBiallelicSNPUnfoldedSFSBuilder(
			Iterator<VariantContext> locusIterator,
			Predicate<VariantContext> locusFilter,
			Map<String, Predicate<VariantContext>> locusTypeNamePredicateMap,
			Function<VariantContext, List<String>> locusTypeNamesProducer,
			Map<Integer, Integer> selectedOutgroupSampleIndexPopLevelMap, int minOutgroupPopLevels,
			List<Integer> ingroupSampleIndices, int minIngroupSampleNumWithNonMissingDataToInclude) {
		super(locusIterator, locusFilter, locusTypeNamePredicateMap, locusTypeNamesProducer, selectedOutgroupSampleIndexPopLevelMap,
				minOutgroupPopLevels, ingroupSampleIndices, minIngroupSampleNumWithNonMissingDataToInclude,
				true,//boolean onlyIncludeBiallelicLoci
				false, //boolean toStoreCalculatedSingleLocusData
				false //boolean toBuildSeparateSFSForAncestralAlleleOfPresenceAndAbsenceType
				);
		// TODO Auto-generated constructor stub
	}
	
	
	/**
	 * the given locus should contains the {@link #currentLocusAbsenceTypeAllele} and/or {@link #currentLocusPresenceTypeAllele} in the selected outgroup samples
	 * 
	 */
	@Override
	protected Pair<Boolean, Boolean> checkAlleleTypeOfGenotypeCalledOutgroupSamples(VariantContext locus) {
		//check if outgroup sample with called genotype are all homozygous of the same type;
		boolean allOutgroupSampleWithAbsenceGenotype = true;
		boolean allOutgroupSampleWithPresenceGenotype = true;
		
		for(int outgroupSampleIndex:this.getSelectedOutgroupSampleIndexPopLevelMap().keySet()) {
			Genotype gt=locus.getGenotype(outgroupSampleIndex);
			if(gt.isCalled()) {
				if(gt.isHom()) {//homozygous genotype
					if(gt.getAllele(0).equals(this.currentLocusAbsenceTypeAllele)) {
						allOutgroupSampleWithPresenceGenotype=false;
					}else if(gt.getAllele(0).equals(this.currentLocusPresenceTypeAllele)){
						allOutgroupSampleWithAbsenceGenotype=false;
					}else {
						throw new IllegalArgumentException("given locus contains allele type different from the assigned 'presence' and 'absence' types!");
					}
				}else {//HETER
					allOutgroupSampleWithAbsenceGenotype=false;
					allOutgroupSampleWithPresenceGenotype=false;
				}
			}
		}
		
		return new Pair<>(allOutgroupSampleWithAbsenceGenotype, allOutgroupSampleWithPresenceGenotype);
	}
	
	@Override
	protected List<Integer> getSampleIndicesWithNonMissingData(VariantContext locus,
			Collection<Integer> targetSampleIndices) {
		List<Integer> ret=new ArrayList<>();
		
		for(int index:targetSampleIndices) {
			if(locus.getGenotype(index).isCalled()) {//TODO??? confirmation is needed
				ret.add(index);
			}
		}
		
		return ret;
	}
	
	
	/**
	 * this method will set the {@link #currentLocusAbsenceTypeAllele} to one of the called allele of the given biallelic locus and {@link #currentLocusPresenceTypeAllele} to the other called allele
	 * 
	 * @param locus
	 * @param targetSampleIndices
	 * @return
	 */
	private void setPresenceAbsenceTypeAllele(VariantContext biallelicLocus) {
		List<Allele> alleles=new ArrayList<>(VariantContextUtils.getCalledAlleles(biallelicLocus, getAllSampleIndices()));
		
		if(alleles.size()!=2) {
			throw new IllegalArgumentException("given locus is not biallelic!");
		}
		
		Collections.sort(alleles);
		
		this.currentLocusPresenceTypeAllele=alleles.get(0);
		this.currentLocusAbsenceTypeAllele=alleles.get(1);
	}
	
	/**
	 * for SNP locus, 
	 * 		
	 * 
	 */
	@Override
	protected int getPresenceTypeAlleleCount(VariantContext locus, Collection<Integer> targetSampleIndices) {
		return VariantContextUtils.getAlleleCount(locus, targetSampleIndices, this.currentLocusPresenceTypeAllele);
	}
	
	@Override
	protected int getAbsenceTypeAlleleCount(VariantContext locus, Collection<Integer> targetSampleIndices) {
		return VariantContextUtils.getAlleleCount(locus, targetSampleIndices, this.currentLocusAbsenceTypeAllele);
	}
	
	@Override
	protected boolean locusIsVariant(VariantContext locus) {
		return VariantContextUtils.isVariant(locus, getAllSampleIndices());
	}

	
	@Override
	protected boolean locusIsBiallelic(VariantContext locus) {
		boolean ret=VariantContextUtils.isBiallelic(locus, getAllSampleIndices());
		if(ret) {
			//set the presence and absence type allele
			this.setPresenceAbsenceTypeAllele(locus);
		}
		return ret;
	}
	
}
