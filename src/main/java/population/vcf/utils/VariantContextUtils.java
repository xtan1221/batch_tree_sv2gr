package population.vcf.utils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;


public class VariantContextUtils {
	
	/**
	 * return the number of allele of the given allele type in the given samples in the given VariantContext
	 * @param vc
	 * @param targetSampleIndices
	 * @param allele
	 * @return
	 */
	public static int getAlleleCount(VariantContext vc, Collection<Integer> targetSampleIndices, Allele allele) {
		int ret = 0;
		for(int index:targetSampleIndices) {
			for(Allele a:vc.getGenotype(index).getAlleles()) {
				if(a.isCalled()) {
					if(a.equals(allele)) {
						ret++;
					}
				}
			}
		}
		return ret;
	}
	
	/**
	 * extract and return the set of called alleles of the given VariantContext in the given samples
	 * 
	 * @param vc
	 * @param targetSampleIndices
	 * @return
	 */
	public static Set<Allele> getCalledAlleles(VariantContext vc, Collection<Integer> targetSampleIndices){
		Set<Allele> ret = new HashSet<>();
		for(int index:targetSampleIndices) {
			for(Allele a:vc.getGenotype(index).getAlleles()) {
				if(a.isCalled()) {
					ret.add(a);
				}
			}
		}
		return ret;
	}
	
	/**
	 * check whether there are exactly two types (biallelic) of called alleles of given {@link VariantContext} in terms of the given set of sample indices;
	 * 
	 * all non-called alleles ({@link Allele#NO_CALL}) will be skipped
	 * 
	 * @param original
	 * @param targetSamples
	 * @return
	 */
	public static boolean isBiallelic(VariantContext vc, Collection<Integer> targetSampleIndices) {
		Set<Allele> alleles=new HashSet<>();
		
		for(int index:targetSampleIndices) {
			Genotype gt=vc.getGenotype(index);
			for(Allele a:gt.getAlleles()) {
				if(a.isCalled()) { //there is an called allele not a single nucleotide
					alleles.add(a);
				}
			}
		}
		
		return alleles.size()==2;
	}
	
	/**
	 * check if the given locus is variant in terms of the given samples
	 * 
	 * for a locus to be variant, it must have two or more types of called alleles
	 * 
	 * all non-called alleles ({@link Allele#NO_CALL}) will be skipped
	 * 
	 * @param original
	 * @param targetSamples
	 * @return
	 */
	public static boolean isVariant(VariantContext vc, Collection<Integer> targetSampleIndices) {
		Set<Allele> alleles=new HashSet<>();
		
		for(int index:targetSampleIndices) {
			Genotype gt=vc.getGenotype(index);
			for(Allele a:gt.getAlleles()) {
				if(a.isCalled()) { //there is an called allele not a single nucleotide
					alleles.add(a);
				}
			}
		}
		
		return alleles.size()>=2;
	}
	
	/**
	 * check if the called alleles of the given locus are all single nucleotide or not 
	 * note that all non-called alleles ({@link Allele#NO_CALL}) will be skipped
	 * 
	 * @param original
	 * @param targetSamples
	 * @return
	 */
	public static boolean isSingleNucleotideSite(VariantContext vc, Collection<Integer> targetSampleIndices) {
		for(int index:targetSampleIndices) {
			Genotype gt=vc.getGenotype(index);
			for(Allele a:gt.getAlleles()) {
				if(a.isCalled() && !GenotypeAndAlleleUtils.isSingleNucleotide(a)) { //there is an called allele not a single nucleotide
					return false;
				}
			}
		}
		
		return true;
	}
	
	
	/**
	 * build and return a new {@link VariantContext} containing the genotypes of the given set of samples that is a subset of the samples of the given {@link VariantContext};
	 * return null if none of the given samples have called non-reference alleles
	 * 
	 * @param original
	 * @param sampleNames
	 * @return
	 */
	public static VariantContext subsetGenotypes(VariantContext original, List<String> targetSamples) {
		Set<String> targetSampleSet=new HashSet<>();
		targetSampleSet.addAll(targetSamples);
		if(targetSampleSet.size()!=targetSamples.size()) {
			throw new IllegalArgumentException("given targetSampleSet contains duplicate samples!");
		}
		
		VariantContextBuilder vcb=new VariantContextBuilder(original);
		
		//the genotypes containing only the given samples
		List<Genotype> genotypes = new ArrayList<>();
		
		boolean nonRefAlleleFound=false;
		Set<String> addedSamples=new HashSet<>();
		for(int i=0;i<original.getGenotypes().size();i++) {
			Genotype gt=original.getGenotypes().get(i);
			if(targetSamples.contains(gt.getSampleName())) {
				if(gt.isCalled() && !gt.isHomRef()) {
					nonRefAlleleFound=true; //called non-ref allele is found
				}
				genotypes.add(gt);
				addedSamples.add(gt.getSampleName());
			}
		}
		
		if(!addedSamples.equals(targetSampleSet)) {
			throw new IllegalArgumentException("at least one of the given samples is not found in the given VariantContext!");
		}
		
		if(nonRefAlleleFound) {
			//set the genotypes
			vcb=vcb.genotypes(genotypes);
			return vcb.make();
		}else {
			return null;
		}
	}
	
	
	/**
	 * return true if the site containing the given genotypes of samples all are reference allele
	 * 
	 * @param genotypes
	 * @return
	 */
	public static boolean nonVariantSite(List<Genotype> genotypes, boolean allowingNonCalledGenotype) {
		for(Genotype gt:genotypes) {
			if(gt.isNoCall() && !allowingNonCalledGenotype) {
				return false;
			}
			if(!gt.isNoCall()) {
				if(!gt.isHomRef()) {
					return false;
				}
			}
		}
		
		return true;
	}
	
	
	
	/**
	 * return whether the site containing the given samples (genotypes) is a biallelic site or not
	 * @param genotypes
	 * @param allowingNonCalledGenotype
	 * @return
	 */
	public static boolean biallelicSNPSite(List<Genotype> genotypes, boolean allowingNonCalledGenotype) {
		Set<Allele> calledSNPAlleles=new HashSet<>();
		
		for(Genotype gt:genotypes) {
			if(gt.isNoCall() && !allowingNonCalledGenotype) {
				return false;
			}
			
			if(!gt.isNoCall()) {
				for(Allele a:gt.getAlleles()) {
					if(isSNP(a)) {
						calledSNPAlleles.add(a);
					}else {
						return false;
					}
				}
			}
		}
		return calledSNPAlleles.size()==2;
	}
	
	
	static boolean isSNP(Allele a) {
		return a.getBaseString().equals("A") 
				|| a.getBaseString().equals("T") 
				|| a.getBaseString().equals("G") 
				|| a.getBaseString().equals("C") ;
	}
	
	
	/**
	 * return whether the site containing the given samples (genotypes) is a multi-allelic site or not
	 * @param genotypes
	 * @param allowingNonCalledGenotype
	 * @return
	 */
	public static boolean multiAllelicSNPSite(List<Genotype> genotypes, boolean allowingNonCalledGenotype) {
		Set<Allele> calledSNPAlleles=new HashSet<>();
		
		for(Genotype gt:genotypes) {
			if(gt.isNoCall() && !allowingNonCalledGenotype) {
				return false;
			}
			
			if(!gt.isNoCall()) {
				for(Allele a:gt.getAlleles()) {
					if(isSNP(a)) {
						calledSNPAlleles.add(a);
					}else {
						return false;
					}
				}
			}
		}
		return calledSNPAlleles.size()>2;
	}
	
	
	
	public static Set<Allele> retrieveAllAlleles(VariantContext vc, List<Integer> targetSampleIndexList, boolean sampleIndex1Based) {
		Set<Allele> alleles = new HashSet<>();
		
		for(int index:targetSampleIndexList) {
			Genotype gt = vc.getGenotype(sampleIndex1Based?index-1:index);
			
			alleles.addAll(gt.getAlleles());
			
		}
		
		return alleles;
	}
	
	
	
	
	/**
	 * 
	 * @param alleles
	 * @return
	 */
	public static boolean allRefOrNonCalled(Set<Allele> alleles) {
		
		for(Allele a:alleles) {
			if(!a.isReference() && !a.isNoCall()) {//
				
				
				
				return false;
			}
		}
		
		return true;
	}
	
	
	
	/**
	 * return the number of samples with non-missing data of the given {@link VariantContext}
	 * 
	 * @param targetSampleIndices
	 * @param sampleIndex1Based if true, the given sample index starts from 1, should be adjusted to 0-based (VariantContext genotype sample index)
	 * @return
	 */
	public static int getTotalSampleNumWithNonMissingData(VariantContext vc, Collection<Integer> targetSampleIndices, boolean sampleIndex1Based) {
		Set<Integer> indices = new HashSet<>();
		indices.addAll(targetSampleIndices);
		int ret=0;
		for(int sampleIndex:indices) {
			Genotype gt;
			if(sampleIndex1Based) {
				gt = vc.getGenotype(sampleIndex-1);
			}else {
				gt=vc.getGenotype(sampleIndex);
			}
			
			if(!gt.isNoCall()) {
				ret++;
			}
		}
		
		return ret;
	}
	
	/**
	 * 
	 * @param vc
	 * @param sampleIndexPopulationHierarchyLevelMap
	 * @param targetSampleIndices
	 * @param sampleIndex1Based
	 * @return
	 */
	public static int getTotalPopulationHierarchyLevelsOfSampleWithNonMissingData(
			VariantContext vc,
			Map<Integer, Integer> sampleIndexPopulationHierarchyLevelMap,
			Collection<Integer> targetSampleIndices, 
			boolean sampleIndex1Based) {
		Set<Integer> PopulationHierarchyLevels = new HashSet<>();
		
		Set<Integer> indices = new HashSet<>();
		indices.addAll(targetSampleIndices);
		for(int sampleIndex:indices) {
			Genotype gt;
			if(sampleIndex1Based) {
				gt = vc.getGenotype(sampleIndex-1);
			}else {
				gt=vc.getGenotype(sampleIndex);
			}
			
			if(!gt.isNoCall()) {
				PopulationHierarchyLevels.add(sampleIndexPopulationHierarchyLevelMap.get(sampleIndex));
			}
		}
		
		return PopulationHierarchyLevels.size();
	}
	
}
