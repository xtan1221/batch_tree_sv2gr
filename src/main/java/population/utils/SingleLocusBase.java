package population.utils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import basic.Pair;
import genomics.utils.AbstractGenomicRegion;

/**
 * base class for a called variant allele containing locus in a population with each sample assigned a {@link Genotype}
 * 
 * for example, a locus for a SV, MEI or R gene insert in a population of a group of samples
 * 
 * @author tanxu
 * 
 */
public abstract class SingleLocusBase extends AbstractGenomicRegion{
	
	public SingleLocusBase() {
		super();
		// TODO Auto-generated constructor stub
	}
	
	public SingleLocusBase(String chrom) {
		super(chrom, null);
		// TODO Auto-generated constructor stub
	}
	
	/**
	 * note that the first sample has index 1 rather than 0!!!!!!!!!!!!!
	 * 
	 * @return the sampleIndexGenotypeMap
	 */
	public abstract Map<Integer, Genotype> getSampleIndexGenotypeMap();
	
	/**
	 * return the number of samples whose genotype is not {@link Genotype#MISSING} for the given set of sample index
	 * @return
	 */
	public final int getTotalSampleNumWithNonMissingData(Collection<Integer> targetSampleIndices) {
		return this.getSampleIndicesWithNonMissingData(targetSampleIndices).size();
	}
	
	/**
	 * return the number of samples whose genotype is not {@link Genotype#MISSING} for the given set of sample index
	 * @return
	 */
	public final List<Integer> getSampleIndicesWithNonMissingData(Collection<Integer> targetSampleIndices) {
		List<Integer> ret=new ArrayList<>();
		
		for(int sampleIndex:this.getSampleIndexGenotypeMap().keySet()) {
			if(targetSampleIndices.contains(sampleIndex)) {
				if(!this.getSampleIndexGenotypeMap().get(sampleIndex).equals(Genotype.MISSING)) {
					ret.add(sampleIndex);
				}
			}
		}
		
		return ret;
	}
	
	/**
	 * return the total allele counts of the given samples with non-missing genotype ;
	 * 
	 * which is equal to 2 * number of given sample with non-missing data
	 * @return
	 */
	public final int getTotalAlleleCountOfSampleWithNonMissingData(Collection<Integer> targetSampleIndices) {
		return this.getTotalSampleNumWithNonMissingData(targetSampleIndices)*2;
	}
	
	
	/**
	 * return the total allele counts representing presence of the genomic feature at this locus, for the given set of samples
	 * which is equal to the 2*(sample num with {@link Genotype#PRESENCE} genotype) + (sample num with {@link Genotype#HETER} genotype);
	 * 
	 * @param targetSampleIndices
	 * @return
	 */
	public final int getTotalPresenceAlleleCount(Collection<Integer> targetSampleIndices) {
		int ret = 0;
		
		for(int sampleIndex:this.getSampleIndexGenotypeMap().keySet()) {
			if(targetSampleIndices.contains(sampleIndex)) {
				if(this.getSampleIndexGenotypeMap().get(sampleIndex).equals(Genotype.PRESENCE)) {
					ret+=2;
				}else if(this.getSampleIndexGenotypeMap().get(sampleIndex).equals(Genotype.HETER)) {
					ret+=1;
				}
			}
		}
		
		return ret;
	}
	
	/**
	 * return the total allele counts representing absence of the genomic feature at this locus, for the given set of samples
	 * which is equal to the 2*(sample num with {@link Genotype#PRESENCE} genotype) + (sample num with {@link Genotype#HETER} genotype);
	 * 
	 * @param targetSampleIndices
	 * @return
	 */
	public final int getTotalAbsenceAlleleCount(Collection<Integer> targetSampleIndices) {
		int ret = 0;
		
		for(int sampleIndex:this.getSampleIndexGenotypeMap().keySet()) {
			if(targetSampleIndices.contains(sampleIndex)) {
				if(this.getSampleIndexGenotypeMap().get(sampleIndex).equals(Genotype.ABSENCE)) {
					ret+=2;
				}else if(this.getSampleIndexGenotypeMap().get(sampleIndex).equals(Genotype.HETER)) {
					ret+=1;
				}
			}
		}
		
		return ret;
	}
	
	/**
	 * return the total allele counts of all samples with non-missing genotype that is estimated as {@link Genotype#PRESENCE},
	 * which is equal to the 2*(sample num with {@link Genotype#PRESENCE} genotype) + (sample num with {@link Genotype#HETER} genotype);
	 * 
	 * @return
	 */
	public final int getTotalPresenceAlleleCount() {
		return this.getTotalPresenceAlleleCount(this.getSampleIndexGenotypeMap().keySet());
	}
	
	/**
	 * return the number of samples whose genotype is not {@link Genotype#MISSING}
	 * @return
	 */
	public final int getTotalSampleNumWithNonMissingData() {
		return this.getTotalSampleNumWithNonMissingData(this.getSampleIndexGenotypeMap().keySet());
	}
	
	/**
	 * return the total allele counts of all samples with non-missing genotype;
	 * which is equal to 2*{@link #getTotalSampleNumWithNonMissingData()} because of diploid
	 * @return
	 */
	public final int getTotalAlleleCountOfSampleWithNonMissingData() {
		return this.getTotalSampleNumWithNonMissingData()*2;
	}
	
	
	/**
	 * check whether this locus is variant or not in terms of the given set of samples;
	 * 
	 * note that a locus may be variant for a set of samples while non-variant for a different set of samples!
	 * 
	 * @param locus
	 * @param outgroupSampleIndices
	 * @param ingroupSampleIndices
	 * @return
	 */
	public boolean isVariant(Collection<Integer> sampleIndices) {
		//no sample has called genotype
		if(this.getSampleIndicesWithNonMissingData(sampleIndices).size()==0) {
			return false;
		}
		
		//check if all samples with called genotype are homozygous of the same allele
		int presenceTypeAlleleCount=this.getTotalPresenceAlleleCount(sampleIndices);
		int absenceTypeAlleleCount=this.getTotalAbsenceAlleleCount(sampleIndices);
		
		
		if(presenceTypeAlleleCount==0||absenceTypeAlleleCount==0) {
			return false;
		}
		
		/////////////
		return true;
	}
	
	/**
	 * return whether this locus is biallelic in terms of the given set of samples
	 * 
	 * note that currently all variant locus of {@link SingleLocusBase} subtype are biallelic 
	 * 		sv, mei, rgene
	 * ==========================
	 * this may not be the case in future version?
	 */
	public boolean isBiallelic(Collection<Integer> sampleIndices) {
		// TODO Auto-generated method stub
		return this.isVariant(sampleIndices);
	}
	
	/**
	 * check and return the state of alleles of the given set of outgroup samples;
	 * 		the result is used to facilitate inference of the ancient allele state of this locus in terms of the given outgroup samples
	 * 		
	 * @param locus
	 * @return first boolean value indicates whether all given outgroup samples with called non-missing genotype are homozygous absence of variant in this locus; second boolean value indicates whether all given outgroup samples with called non-missing genotype are homozygous presence of variant in this locus
	 * 
	 */
	public Pair<Boolean, Boolean> checkAlleleTypeOfGenotypeCalledOutgroupSamples(Collection<Integer> outgroupSampleIndices) {
		//check if outgroup sample with called genotype are all homozygous of the same type;
		boolean allOutgroupSampleWithAbsenceGenotype = true;
		boolean allOutgroupSampleWithPresenceGenotype = true;
		
		for(int outgroupSampleIndex:outgroupSampleIndices) {
			Genotype gt=this.getSampleIndexGenotypeMap().get(outgroupSampleIndex);
			if(gt.equals(Genotype.ABSENCE)) {
				allOutgroupSampleWithPresenceGenotype=false;
			}else if(gt.equals(Genotype.PRESENCE)) {
				allOutgroupSampleWithAbsenceGenotype=false;
			}else {//HETER
				allOutgroupSampleWithAbsenceGenotype=false;
				allOutgroupSampleWithPresenceGenotype=false;
			}
		}
		
		return new Pair<>(allOutgroupSampleWithAbsenceGenotype, allOutgroupSampleWithPresenceGenotype);
	}
	
	
}
