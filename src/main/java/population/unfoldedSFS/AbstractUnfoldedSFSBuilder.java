package population.unfoldedSFS;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;
import java.util.function.Predicate;

import basic.Pair;
import htsjdk.variant.variantcontext.VariantContext;
import population.utils.SingleLocusBase;

/**
 * abstract class to build unfolded SFS for a set of variant loci
 * 
 * the main task is to perform the following steps for each locus:
 * 		1. check if the locus is variant or not (continue if variant, skip if not)
 * 			a locus is non-variant only if there all outgroup samples and ingroup samples with called genotype are homogzygous of the same allele;
 * 			otherwise, the locus is variant;
 * 		
 * 		2. infer the ancestral allele state of the variant locus
 * 			the strategy is constrained by the {@link #minOutgroupPopLevels} and {@link #minIngroupSampleNumWithNonMissingDataToInclude}
 * 		
 * 		3. calculate derived allele frequency and proportion for the locus
 * 			the derived allele frequency is the number of derived alleles for the locus, which is useless to construct unfolded SFS since some samples may have missing data (genotype ./.)
 * 			thus, only the proportion of derived allele of a locus in the ingroup sample is useful to construct the SFS;
 * 		
 * 		4. build bins for unfolded SFS using all loci with inferred ancestral allele state and calculated derived allele proportion
 * 			note that for each locus type, a separate unfolded SFS will be constructed;
 * 			
 * ==================================================
 * all unfolded SFS calculation should be subclass of this class for the following variation types:
 * 		1. SNP and indel (including snpEff annotated SNPs)
 * 			use the {@link VariantContext} for variant loci
 * 		2. SV
 * 		3. MEI
 * 		4. R gene presence/absence
 * 
 * ==================================
 * derived allele frequency: the number of derived allele of a specific locus
 * derived allele proportion: the proportion of derived alleles of all alleles for a specific locus
 * 
 * @author tanxu
 * 
 * @param <T> type of a locus for a specific genomic feature in a population, could be a {@link VariantContext} of htsjdk or a subtype of {@link SingleLocusBase} of this package
 */
public abstract class AbstractUnfoldedSFSBuilder<T> {
	/**
	 * iterator for all variant loci to calculate the unfolded SFS
	 * each loci should contain genotype for both ingroup samples and outgroup samples selected for ancient allele state inference
	 */
	private final Iterator<T> locusIterator;
	
	/**
	 * 
	 */
	private final Predicate<T> locusFilter;
	
	/**
	 * map from type name to the Predicate of each locus;
	 * the {@link Predicate}s are not required to form a complete covering of all loci; 
	 * locus not identified to be any of the types will be excluded from the constructed SFS
	 * 
	 * 
	 * it is allowed for a locus to be of multiple types
	 * 		for example, for a R gene insert locus with homologous reference R gene of LRR-PKinase domain composition, the locus can be of the following types simutaneously
	 * 			1. exact LRR-Pkinase
	 * 			2. LRR-containing
	 * 			3. Pkinase-containing
	 * 
	 * 
	 * either {@link #locusTypeNamePredicateMap} or {@link #locusTypeNamesProducer} should be provided;
	 * the difference is that types of locus given by this function are known beforehand while while the types are not known until runtime  if {@link #locusTypeNamesProducer} is provided
	 */
	protected final Map<String, Predicate<T>> locusTypeNamePredicateMap;
	
	/**
	 * function on how the types of a locus can be produced from the locus;
	 * 
	 * locus not identified to be any of the types will be excluded from the constructed SFS
	 * 
	 * this should be provided if types of variations are not known beforehand, for example annotation types of SNP sites by snpEff;
	 * 
	 * either {@link #locusTypeNamePredicateMap} or {@link #locusTypeNamesProducer} should be non-null;
	 * 
	 * the difference is that types of locus given by this function is not known until runtime while the types are known beforehand if {@link #locusTypeNamePredicateMap} is provided
	 */
	protected final Function<T, List<String>> locusTypeNamesProducer;
	
	/**
	 * map from sample index to the population level of a set of selected outgroup samples to infer the ancient allele state
	 * the sample index must be consistent with the map between sample index and the genotype
	 * 		for {@link VariantContext}, the {@link VariantContext#getGenotype(int)} has first sample's genotype with index 0
	 * 		for {@link SingleLocusBase#getSampleIndexGenotypeMap()} has first sample's index as 1  
	 * 
	 */
	private final Map<Integer, Integer> selectedOutgroupSampleIndexPopLevelMap;
	
	
	/**
	 * minimal population levels of outgroup required to infer the ancient allele state;
	 * must be positive integer
	 * 
	 * for example, 
	 * 		if 1, only outgroup sample(s) of a single population level with called homozygous genotype is enough to infer the ancient allele state;
	 * 		if 2, at least outgroup samples of two or more population levels with called homozygous genotype are required to infer the ancient allele state;
	 * 
	 */
	private final int minOutgroupPopLevels;
	
	/**
	 * 
	 */
	private final List<Integer> ingroupSampleIndices;
	
	/**
	 * only when the number of ingroup sample of a locus with non-missing data equal to or more than this value, the locus will be included in the SFS construction
	 * must be positive integer
	 */
	private final int minIngroupSampleNumWithNonMissingDataToInclude;
	
	/**
	 * whether or not only consider biallelic loci;
	 * 
	 */
	private final boolean onlyIncludeBiallelicLoci;
	
	/**
	 * whether or not to store the inferred ancient allele state and the derived allele freq(proportion) of each locus to
	 * {@link #locusDerivedAlleleProportionMap} and {@link #locusPresenceAsDerivedAlleleMap}
	 * 
	 * 
	 */
	private final boolean toStoreCalculatedSingleLocusData;
	
	/**
	 * whether or not to build unfolded SFS separately for loci with inferred ancestral state of PRESENCE and ABSENCE of the genomic feature
	 * if true, for loci of a type with name "T", three unfoded SFS will be built:
	 * 		1. SFS using all loci of that type
	 * 		2. SFS using loci of that type with inferred ancestral state being "PRESENCE" of the genomic feature with type name as 'T-Presence'
	 * 		3. SFS using loci of that type with inferred ancestral state being "ABSENCE" of the genomic feature with type name as 'T-Absence'
	 * if false, only the first SFS will be built;
	 * =============================
	 * this field is primarily for genomic features of a functional structure such as R gene or mobile element to facilitate investigating difference in loss and gain of the feature in the population
	 * it makes little sense to build such SFSs for features such as SV and SNP!!!!!!!!!!
	 */
	private final boolean toIncludeSeparateSFSForAncestralAlleleOfPresenceAndAbsenceType;
	
	/////////////////////
//	/**
//	 * with of bin,equal to 1/(2 * number of ingroup sample)
//	 */
//	private double binWidth;
	
	/**
	 * indices of all ingroup samples and selected outgroup samples
	 */
	private List<Integer> allSampleIndices;
	/**
	 * a 
	 */
	protected List<String> orderedLocusTypeNames;
	
	/**
	 * map from the bin index to Bin;
	 * 
	 * total number of bin is 2n, where n is the number of ingroup samples
	 * the first bin index is 1, the last bin index is 2n;
	 * the ith bin has derived allele proportion in range ((i-1)/2n, i/2n]
	 * 
	 * for more details, see {@link CompositeDerivedAlleleProportionBin}
	 */
	private Map<Integer, CompositeDerivedAlleleProportionBin> binIndexMap;
	
	
	/**
	 * map from the locus type name to the total number of locus in the ingroup that are included in the SFS construction;
	 */
	private Map<String, Integer> typeNameTotalLocusNumberIncludedInSFSMap;
	
	
	/**
	 * map from the type of the variant to the number of locus included in SFS and ancestral state is inferred as PRESENCE of the variant
	 * for MEI and R gene insertion, the PRESENCE type variant means the presence of the ME/R gene
	 */
	private Map<String, Integer> typeNameTotalNumOfLocusWithPresenceTypeAncestralAlleleStateIncludedInSFSMap;

	/**
	 * map from the type of the variant to the number of locus included in SFS and ancestral state is inferred as ABSENCE of the variant
	 * for MEI and R gene insertion, the ABSENCE type variant means the absence of the ME/R gene
	 */
	private Map<String, Integer> typeNameTotalNumOfLocusWithAbsenceTypeAncestralAlleleStateIncludedInSFSMap;
	
	
	/**
	 * map from the locus type to whether 
	 * 1. the presence type allele of the {@link SingleLocusBase} is detected as derived allele and the absence type allele of the {@link SingleLocusBase} is detected as ancient allele
	 * 2. the absence type allele of the {@link SingleLocusBase} is detected as derived allele and the presence type allele of the {@link SingleLocusBase} is detected as ancient allele
	 * 
	 * the map value is true if the first scenario is met; false otherwise;
	 */
	private Map<T, Boolean> locusPresenceAsDerivedAlleleMap;
	
	
	/**
	 * map from the {@link SingleLocusBase} to the calculated proportion of derived allele
	 */
	private Map<T, Double> locusDerivedAlleleProportionMap;
	
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
	 * @param toStoreCalculatedSingleLocusData
	 */
	public AbstractUnfoldedSFSBuilder(
			Iterator<T> locusIterator,
			Predicate<T> locusFilter,
			Map<String, Predicate<T>> locusTypeNamePredicateMap, 
			Function<T, List<String>> locusTypeNamesProducer,
			Map<Integer, Integer> selectedOutgroupSampleIndexPopLevelMap,
			int minOutgroupPopLevels,
			List<Integer> ingroupSampleIndices,
			int minIngroupSampleNumWithNonMissingDataToInclude,
			boolean onlyIncludeBiallelicLoci,
			boolean toStoreCalculatedSingleLocusData,
			boolean toBuildSeparateSFSForAncestralAlleleOfPresenceAndAbsenceType
			) {
		super();
		if(minOutgroupPopLevels<=0) {
			throw new IllegalArgumentException("given minOutgroupPopLevels must be positive integer!");
		}
		if(minIngroupSampleNumWithNonMissingDataToInclude<=0) {
			throw new IllegalArgumentException("given minIngroupSampleNumWithNonMissingDataToInclude must be positive integer!");
		}
		if(locusTypeNamePredicateMap==null&&locusTypeNamesProducer==null) {
			throw new IllegalArgumentException("given locusTypeNamePredicateMap and locusTypeNameProducer cannot both be null!");
		}
		if(locusTypeNamePredicateMap!=null&&locusTypeNamesProducer!=null) {
			throw new IllegalArgumentException("given locusTypeNamePredicateMap and locusTypeNameProducer cannot both be non-null!");
		}
		
		////////////
		this.locusIterator=locusIterator;
		this.locusFilter=locusFilter;
		this.locusTypeNamePredicateMap=locusTypeNamePredicateMap;
		this.locusTypeNamesProducer=locusTypeNamesProducer;
		this.selectedOutgroupSampleIndexPopLevelMap=selectedOutgroupSampleIndexPopLevelMap;
		this.minOutgroupPopLevels=minOutgroupPopLevels;
		this.ingroupSampleIndices=ingroupSampleIndices;
		this.minIngroupSampleNumWithNonMissingDataToInclude = minIngroupSampleNumWithNonMissingDataToInclude;
		this.onlyIncludeBiallelicLoci=onlyIncludeBiallelicLoci;
		this.toStoreCalculatedSingleLocusData=toStoreCalculatedSingleLocusData;
		this.toIncludeSeparateSFSForAncestralAlleleOfPresenceAndAbsenceType=toBuildSeparateSFSForAncestralAlleleOfPresenceAndAbsenceType;
		/////////////////////////////
		this.prepare();
		this.buildBins();
		this.run();
	}
	
	private void prepare() {
		this.typeNameTotalLocusNumberIncludedInSFSMap=new HashMap<>();
		this.typeNameTotalNumOfLocusWithAbsenceTypeAncestralAlleleStateIncludedInSFSMap=new HashMap<>();
		this.typeNameTotalNumOfLocusWithPresenceTypeAncestralAlleleStateIncludedInSFSMap=new HashMap<>();
		
		this.orderedLocusTypeNames=new ArrayList<>();
		
		////////////////////
		if(this.toStoreCalculatedSingleLocusData) {
			this.locusPresenceAsDerivedAlleleMap=new HashMap<>();
			this.locusDerivedAlleleProportionMap=new HashMap<>();
		}
		
	}
	
	private void buildBins() {
		int totalAlleleCount=(2*this.ingroupSampleIndices.size());
		
//		this.binWidth=(double)1/totalAlleleCount;
		
		this.binIndexMap=new LinkedHashMap<>();
		
		for(int i=0;i<this.ingroupSampleIndices.size()*2;i++) {
//			double start = this.binWidth*i; // this will lead to incorreact assignment; for example, for 46/48 derived allele, the proportion calculated is '0.9583333333333334', but the corresponding bin range calculated by this method is (0.9375, 0.9583333333333333], thus the locus can not be correctly assigned!
//			double end=this.binWidth*(i+1);
			double start = (double)i/totalAlleleCount; //this is consistent with how derived allele proportion is calculated so that it can be correctly assigned to the bin!!!!!!!!!!
			double end=(double)(i+1)/totalAlleleCount;//this is consistent with how derived allele proportion is calculated so that it can be correctly assigned to the bin!!!!!!!!!!
			
			CompositeDerivedAlleleProportionBin bin=new CompositeDerivedAlleleProportionBin(i+1, start, end);
			
			this.binIndexMap.put(i, bin);
		}
	}
	
	
	/**
	 * return the number of allele representing the presence of the genomic feature (SV, SNP, MEI, R gene insert, etc) for the given set of samples in the given locus
	 * 
	 * @param locus
	 * @param targetSampleIndices
	 * @return
	 */
	protected abstract int getPresenceTypeAlleleCount(T locus, Collection<Integer> targetSampleIndices);
	
	/**
	 * return the number of allele representing the absence of the genomic feature (SV, SNP, MEI, etc) for the given set of samples in the given locus
	 * 
	 * @param locus
	 * @param targetSampleIndices
	 * @return
	 */
	protected abstract int getAbsenceTypeAlleleCount(T locus, Collection<Integer> targetSampleIndices);
	
	/**
	 * check and return whether the given locus is a variant locus or not (in terms of the selected outgroup samples and ingroup samples)
	 * 		note that a locus may be variant for a set of ingroup and outgroup samples while be non-variant for a different set of ingroup and outgroup samples
	 * 
	 * this method is needed because only variant locus is included in unfolded SFS construction;
	 * 		for non-variant locus, the derived allele frequency and proportion is always 0
	 * 
	 * ===============================
	 * implementation:
	 * 1. if all outgroup samples are of missing genotype, return false;
	 * 2. if all ingrooup samples are of missing genotype, return false;
	 * 3. if all outgroup and ingroup samples with called genotype are homozygous of the same allele, return false;
	 * 4. otherwise, return true;
	 * 
	 * @param locus
	 * @return
	 */
	protected abstract boolean locusIsVariant(T locus);
	
	/**
	 * check and return whether the given locus is a biallelic locus or not (in terms of the selected outgroup samples and ingroup samples)
	 * 		note that a locus may be biallelic for a set of ingroup and outgroup samples while be multi-allelic or single-allelic for a different set of ingroup and outgroup samples
	 * 
	 * 
	 * @param locus
	 * @return
	 */
	protected abstract boolean locusIsBiallelic(T locus);
	
	/**
	 * check and return the state of alleles of outgroup samples of the given locus;
	 * 		the result is used to infer the ancient allele state of the given variant locus
	 * 		
	 * @param locus
	 * @return first boolean value indicates whether all outgroup samples with called non-missing genotype are homozygous absence of variant in the given locus; second boolean value indicates whether all outgroup samples with called non-missing genotype are homozygous presence of variant in the given locus
	 * 
	 */
	protected abstract Pair<Boolean, Boolean> checkAlleleTypeOfGenotypeCalledOutgroupSamples(T locus);
	
	/**
	 * find out and return the list of sample indices with non-missing data/called genotype of the given locus and given sample indices
	 * 
	 * @param locus
	 * @return
	 */
	protected abstract List<Integer> getSampleIndicesWithNonMissingData(T locus, Collection<Integer> targetSampleIndices);
	
	
	private int counter=0;
//	private int i=0;
	/**
	 * iterate through all variant locus and calculate the unfolded SFS
	 */
	private void run() {
		T locus;
		while(this.locusIterator.hasNext()&&(locus=this.locusIterator.next())!=null) {
			counter++;
			if(counter%1000000==0) {
				int num=counter/1000000;
				System.out.println(num+" M records found and processed!");
			}
			
			if(!this.locusFilter.test(locus)) {
				continue;
			}
			
			if(!this.locusIsVariant(locus)) {//skip locus that is not variant
				continue;
			}
			
			if(this.onlyIncludeBiallelicLoci&&!this.locusIsBiallelic(locus)) {
				continue;//skip non-biallelic sites if only include biallelic ones
			}
			
			//check outgroup sample pop level
			List<Integer> outgroupSampleIndicesWithNonMissingData = 
					this.getSampleIndicesWithNonMissingData(locus, this.selectedOutgroupSampleIndexPopLevelMap.keySet());
			
			Set<Integer> outgroupPopLevelWithNonMissingData=new HashSet<>();
			for(int index:outgroupSampleIndicesWithNonMissingData) {
				outgroupPopLevelWithNonMissingData.add(this.selectedOutgroupSampleIndexPopLevelMap.get(index));
			}
			if(outgroupPopLevelWithNonMissingData.size()<this.minOutgroupPopLevels) {//min pop level is not met
				continue;
			}
			
			//check if outgroup sample with called genotype are all homozygous of the same type;
			Pair<Boolean, Boolean> outgroupGenotypes = this.checkAlleleTypeOfGenotypeCalledOutgroupSamples(locus);
			boolean allOutgroupSampleWithAbsenceGenotype = outgroupGenotypes.getFirst();
			boolean allOutgroupSampleWithPresenceGenotype = outgroupGenotypes.getSecond();
			
			//
			if(!allOutgroupSampleWithAbsenceGenotype && !allOutgroupSampleWithPresenceGenotype) {//skip the locus if no unambiguous ancestral allele state can be inferred
				continue;
			}
			
//			i++;
//			System.out.println(i);
			
			///check ingroup sample with non-missing data;
			List<Integer> ingroupSampleIndicesWithNonMissingData= this.getSampleIndicesWithNonMissingData(locus, this.ingroupSampleIndices);
			if(ingroupSampleIndicesWithNonMissingData.size()<this.minIngroupSampleNumWithNonMissingDataToInclude) {
				continue;
			}
			
			
			//calculate derived allele frequency(proportion)
//			int alleleCountWithPresenceInIngroupSample = this.getPresenceTypeAlleleCount(locus, this.ingroupSampleIndices);
//			int alleleCountWithAbsenceInIngroupSample = this.getAbsenceTypeAlleleCount(locus, this.ingroupSampleIndices);
//			Double derivedAlleleProportion=null;
			Integer derivedAlleleCount=null;
			if(allOutgroupSampleWithPresenceGenotype) {//derived allele is absence type
//				derivedAlleleProportion = 1-(double)alleleCountWithPresenceInIngroupSample/(ingroupSampleIndicesWithNonMissingData.size()*2);
				derivedAlleleCount=this.getAbsenceTypeAlleleCount(locus, this.ingroupSampleIndices);
			}else if(allOutgroupSampleWithAbsenceGenotype) {//derived allele is the presence type
//				derivedAlleleProportion = (double)alleleCountWithPresenceInIngroupSample/(ingroupSampleIndicesWithNonMissingData.size()*2);
				derivedAlleleCount=this.getPresenceTypeAlleleCount(locus, this.ingroupSampleIndices);
			}
			
			//find out the index of the bin
			Integer binIndexByAlleleCount=this.getBinIndex(derivedAlleleCount, ingroupSampleIndicesWithNonMissingData.size());
			
			
			///
			//first look up the type of the locus
			boolean ancestralAlleleStateIsPresenceOfTheFeature=allOutgroupSampleWithPresenceGenotype;
			
			List<String> types = this.getTypes(locus);
			if(types.isEmpty()) { //skip the locus if it is not identified to be any of the required types
				continue;
			}
			
			for(String type:types) {
				this.binIndexMap.get(binIndexByAlleleCount).addOneToType(type, ancestralAlleleStateIsPresenceOfTheFeature);
				
				//////////////////////
				if(this.toStoreCalculatedSingleLocusData) {
					this.locusPresenceAsDerivedAlleleMap.put(locus, !ancestralAlleleStateIsPresenceOfTheFeature);
					this.locusDerivedAlleleProportionMap.put(locus, (double)derivedAlleleCount/(2*ingroupSampleIndicesWithNonMissingData.size()));
				}
				
				//update the total number of locus included in the calculation of SFS of the corresponding type
				if(!this.typeNameTotalLocusNumberIncludedInSFSMap.containsKey(type)) {
					this.typeNameTotalLocusNumberIncludedInSFSMap.put(type,0);
					this.typeNameTotalNumOfLocusWithAbsenceTypeAncestralAlleleStateIncludedInSFSMap.put(type, 0);
					this.typeNameTotalNumOfLocusWithPresenceTypeAncestralAlleleStateIncludedInSFSMap.put(type, 0);
				}
				
				this.typeNameTotalLocusNumberIncludedInSFSMap.put(type,this.typeNameTotalLocusNumberIncludedInSFSMap.get(type)+1);
				if(ancestralAlleleStateIsPresenceOfTheFeature) {
					this.typeNameTotalNumOfLocusWithPresenceTypeAncestralAlleleStateIncludedInSFSMap.put(
							type,this.typeNameTotalNumOfLocusWithPresenceTypeAncestralAlleleStateIncludedInSFSMap.get(type)+1);
				}else {
					this.typeNameTotalNumOfLocusWithAbsenceTypeAncestralAlleleStateIncludedInSFSMap.put(
							type,this.typeNameTotalNumOfLocusWithAbsenceTypeAncestralAlleleStateIncludedInSFSMap.get(type)+1);
				}
			}	
			
		}
	}
	
	/**
	 * retrieve the bin index with the given derived allele count and ingroup sample num with non-missing data
	 * @param derivedAlleleCount
	 * @return
	 */
	private int getBinIndex(int derivedAlleleCount, int ingroupSampleNumWithNonMissingData) {
		if(ingroupSampleNumWithNonMissingData==this.ingroupSampleIndices.size()) {//all ingroup sample have non-missing data
			for(int index:this.binIndexMap.keySet()) {
				if(this.binIndexMap.get(index).getDerivedAlleleCountWhenAllIngroupSamplesHaveNonMissingData()==derivedAlleleCount) {
					return index;
				}
			}
		}else {
			double derivedAlleleProportion=(double)derivedAlleleCount/(ingroupSampleNumWithNonMissingData*2);
			for(int index:this.binIndexMap.keySet()) {
				if(this.binIndexMap.get(index).isCovering(derivedAlleleProportion)) {
					return index;
				}
			}
		}
		
		throw new IllegalArgumentException("given derivedAlleleCount is outside of boundary!");
	}
	
	/**
	 * find out and return the types of the given locus
	 * 
	 * @param locus
	 * @return
	 */
	protected List<String> getTypes(T locus) {
		List<String> ret=new ArrayList<>();
		
		if(this.locusTypeNamePredicateMap!=null) {
			for(String type:this.locusTypeNamePredicateMap.keySet()) {
				if(this.locusTypeNamePredicateMap.get(type).test(locus)) {
					ret.add(type);
				}
			}
		}else {
			ret = this.locusTypeNamesProducer.apply(locus);
		}
		
		/////////////////
		if(!ret.isEmpty()) {
			for(String type:ret) {
				if(!this.orderedLocusTypeNames.contains(type))
					this.orderedLocusTypeNames.add(type);
			}
		}
		
		return ret;
	}
	
	/**
	 * write the calculated SFS to files
	 * @throws IOException
	 */
	public void writeToFiles(Path outDir, String outFileInfoPrefix) {
		/////////////
		Path summaryStatisticsFile = Path.of(outDir.toString(),outFileInfoPrefix.concat(".summary.statistics.txt"));
		this.checkOutputFile(summaryStatisticsFile);
		
		Path outputSFSTableFile=Path.of(outDir.toString(),outFileInfoPrefix.concat(".unfolded.SFS.table.txt"));
		this.checkOutputFile(outputSFSTableFile);
		
		Collections.sort(this.orderedLocusTypeNames);
		
		try {
			///////////////output summary statistics
			BufferedWriter bw = new BufferedWriter(new FileWriter(summaryStatisticsFile.toString()));
//			bw.append("#variant locus num in ingroup samples:");
//			bw.newLine();
//			for(String type:this.typeNameLocusNumInIngroupMap.keySet()) {
//				String line=type.concat("\t").concat(Integer.toString(this.typeNameLocusNumInIngroupMap.get(type)));
//				
//				bw.append(line);
//				bw.newLine();
//			}
			bw.append("#locus num included in SFS:");
			bw.newLine();
			for(String type:this.typeNameTotalLocusNumberIncludedInSFSMap.keySet()) {
				String line=type.concat("\t").concat(Integer.toString(this.typeNameTotalLocusNumberIncludedInSFSMap.get(type)));
				
				bw.append(line);
				bw.newLine();
			}
			bw.flush();
			bw.close();
			
			////////////////////////////////
			//write to unfolded SFS table file
			BufferedWriter sfsWriter = new BufferedWriter(new FileWriter(outputSFSTableFile.toString()));
			sfsWriter.append(this.getSFSTableFileHeaderLine());
			sfsWriter.newLine();
			
			for(int binIndex:this.binIndexMap.keySet()) {
				CompositeDerivedAlleleProportionBin bin=this.binIndexMap.get(binIndex);
				
				for(String type:this.orderedLocusTypeNames) {
					StringBuilder sb = new StringBuilder();
					if(this.typeNameTotalLocusNumberIncludedInSFSMap.get(type)>0) {
						//output the bin of SFS for all loci
						double proportion=(double)bin.getLocusNum(type)/this.typeNameTotalLocusNumberIncludedInSFSMap.get(type);
						
						sb.append(binIndex).append("\t")
						.append(bin.getStart()).append("\t")
						.append(bin.getEnd()).append("\t")
						.append(type).append("\t")
						.append("all").append("\t")
						.append(proportion);
						
						sfsWriter.append(sb.toString());
						sfsWriter.newLine();
						
						if(this.toIncludeSeparateSFSForAncestralAlleleOfPresenceAndAbsenceType) {
							//loci with Presence type ancestral allele 
							sb=new StringBuilder();
							proportion=
									(double)bin.getLocusNumWithPresenceTypeAncestralAllele(type)/this.typeNameTotalNumOfLocusWithPresenceTypeAncestralAlleleStateIncludedInSFSMap.get(type);
							
							sb.append(binIndex).append("\t")
							.append(bin.getStart()).append("\t")
							.append(bin.getEnd()).append("\t")
							.append(type).append("\t")
							.append("Presence").append("\t")
							.append(proportion);
							
							sfsWriter.append(sb.toString());
							sfsWriter.newLine();
							
							//loci with Absence type ancestral allele
							sb=new StringBuilder();
							proportion=
									(double)bin.getLocusNumWithAbsenceTypeAncestralAllele(type)/this.typeNameTotalNumOfLocusWithAbsenceTypeAncestralAlleleStateIncludedInSFSMap.get(type);
							
							sb.append(binIndex).append("\t")
							.append(bin.getStart()).append("\t")
							.append(bin.getEnd()).append("\t")
							.append(type).append("\t")
							.append("Absence").append("\t")
							.append(proportion);
							
							sfsWriter.append(sb.toString());
							sfsWriter.newLine();
						}
					}
				}
			}
			
			sfsWriter.flush();
			sfsWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	protected void checkOutputFile(Path outFile) {
		if(outFile.toFile().exists()) {
			System.out.println("output file "+outFile.toString()+" already exists, delete it..." );
			outFile.toFile().delete();
		}
	}
	
	/**
	 * 
	 * @return
	 */
	protected String getSFSTableFileHeaderLine() {
		StringBuilder sb = new StringBuilder();
		
		sb.append("bin_index").append("\t")
		.append("bin_start").append("\t")
		.append("bin_end").append("\t")
		.append("type").append("\t")
		.append("ancestral_state").append("\t")
		.append("proportion");
		
		return sb.toString();
	}
	
	////////////////////////////
	protected List<Integer> getAllSampleIndices(){
		if(this.allSampleIndices==null) {
			this.allSampleIndices=new ArrayList<>();
			this.allSampleIndices.addAll(this.getIngroupSampleIndices());
			this.allSampleIndices.addAll(this.getSelectedOutgroupSampleIndexPopLevelMap().keySet());
		}
		return this.allSampleIndices;
	}
	/**
	 * @return the selectedOutgroupSampleIndexPopLevelMap
	 */
	protected final Map<Integer, Integer> getSelectedOutgroupSampleIndexPopLevelMap() {
		return selectedOutgroupSampleIndexPopLevelMap;
	}

	/**
	 * @return the ingroupSampleIndices
	 */
	protected final List<Integer> getIngroupSampleIndices() {
		return ingroupSampleIndices;
	}
	
	/////////////////////////////////////
	/**
	 * @return the locusPresenceAsDerivedAlleleMap
	 */
	public final Map<T, Boolean> getLocusPresenceAsDerivedAlleleMap() {
		return locusPresenceAsDerivedAlleleMap;
	}
	
	/**
	 * @return the locusDerivedAlleleProportionMap
	 */
	public final Map<T, Double> getLocusDerivedAlleleProportionMap() {
		return locusDerivedAlleleProportionMap;
	}
	

	/**
	 * @return the typeNameTotalLocusNumberIncludedInSFSMap
	 */
	public Map<String, Integer> getTypeNameTotalLocusNumberIncludedInSFSMap() {
		return typeNameTotalLocusNumberIncludedInSFSMap;
	}

	///////////////////////////

	/**
	 * @return the typeNameTotalNumOfLocusWithPresenceTypeAncestralAlleleStateIncludedInSFSMap
	 */
	public Map<String, Integer> getTypeNameTotalNumOfLocusWithPresenceTypeAncestralAlleleStateIncludedInSFSMap() {
		return typeNameTotalNumOfLocusWithPresenceTypeAncestralAlleleStateIncludedInSFSMap;
	}

	/**
	 * @return the typeNameTotalNumOfLocusWithAbsenceTypeAncestralAlleleStateIncludedInSFSMap
	 */
	public Map<String, Integer> getTypeNameTotalNumOfLocusWithAbsenceTypeAncestralAlleleStateIncludedInSFSMap() {
		return typeNameTotalNumOfLocusWithAbsenceTypeAncestralAlleleStateIncludedInSFSMap;
	}


}
