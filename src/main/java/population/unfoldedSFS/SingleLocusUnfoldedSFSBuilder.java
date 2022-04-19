package population.unfoldedSFS;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.function.Predicate;

import basic.Pair;
import htsjdk.variant.variantcontext.VariantContext;
import population.utils.SingleLocusBase;

/**
 * class to build unfolded SFS for a set of variant loci of {@link SingleLocusBase} type
 * 
 * all unfolded SFS calculation should be subclass of this class for the following variation types:
 * 		1. SV
 * 		2. MEI
 * 		3. R gene presence/absence
 * 
 * @author tanxu
 * 
 * @param <T> type of the variant locus, could be a {@link VariantContext} of htsjdk or a subtype of {@link SingleLocusBase} of this package
 */
public class SingleLocusUnfoldedSFSBuilder<T extends SingleLocusBase> extends AbstractUnfoldedSFSBuilder<T>{
	
	/**
	 * 
	 * @param locusIterator
	 * @param locusFilter
	 * @param locusTypeNamePredicateMap
	 * @param locusTypeNamesProducer
	 * @param selectedOutgroupSampleIndexPopLevelMap
	 * @param minOutgroupPopLevels
	 * @param ingroupSampleIndices
	 * @param minIngroupSampleNumWithNonMissingDataToInclude
	 */
	public SingleLocusUnfoldedSFSBuilder(
			Iterator<T> locusIterator, Predicate<T> locusFilter,
			Map<String, Predicate<T>> locusTypeNamePredicateMap,
			Function<T, List<String>> locusTypeNamesProducer, 
			Map<Integer, Integer> selectedOutgroupSampleIndexPopLevelMap,
			int minOutgroupPopLevels, List<Integer> ingroupSampleIndices,
			int minIngroupSampleNumWithNonMissingDataToInclude) {
		super(locusIterator, locusFilter, locusTypeNamePredicateMap, locusTypeNamesProducer, selectedOutgroupSampleIndexPopLevelMap,
				minOutgroupPopLevels, ingroupSampleIndices, minIngroupSampleNumWithNonMissingDataToInclude,
				true,//boolean onlyIncludeBiallelicLoci
				true, //toStoreCalculatedSingleLocusData
				true//boolean toIncludeSeparateSFSForAncestralAlleleOfPresenceAndAbsenceType
				);
		// TODO Auto-generated constructor stub
	}
	
	
	@Override
	protected Pair<Boolean, Boolean> checkAlleleTypeOfGenotypeCalledOutgroupSamples(T locus) {
		return locus.checkAlleleTypeOfGenotypeCalledOutgroupSamples(this.getSelectedOutgroupSampleIndexPopLevelMap().keySet());
	}
	
	@Override
	protected List<Integer> getSampleIndicesWithNonMissingData(T locus, Collection<Integer> targetSampleIndices) {
		return locus.getSampleIndicesWithNonMissingData(targetSampleIndices);
	}
	
	@Override
	protected boolean locusIsVariant(T locus) {
		return locus.isVariant(this.getAllSampleIndices());
	}
	
	/**
	 * currently all variant locus of {@link SingleLocusBase} subtype are biallelic
	 * 		sv, mei, rgene
	 * ==========================
	 * this may not be the case in future version?
	 */
	@Override
	protected boolean locusIsBiallelic(T locus) {
		return locus.isBiallelic(this.getAllSampleIndices());
	}
	
	@Override
	protected int getPresenceTypeAlleleCount(T locus, Collection<Integer> targetSampleIndices) {
		return locus.getTotalPresenceAlleleCount(targetSampleIndices);
	}
	
	@Override
	protected int getAbsenceTypeAlleleCount(T locus, Collection<Integer> targetSampleIndices) {
		return locus.getTotalAbsenceAlleleCount(targetSampleIndices);
	}


}
