package population.sv.writer;

import java.util.ArrayList;
import java.util.List;

import population.popHierarchy.PopulationStructureFileReader;
import population.sv.utils.SimpleSVLocus;
import population.utils.Genotype;


/**
 * extract {@link SimpleSVLocus}s from a given group based on the allele types in ingroup samples
 * 
 * see {@link #mustHaveTwoAlleleTypes} for details
 * 
 * @author tanxu
 *
 */
public class IngroupSVLocusExtractor {
	
	private final List<SimpleSVLocus> svLocusList;
	
	private final PopulationStructureFileReader dellyRunPopulationStructureFileReader;
	
	/**
	 * if true, only SimpleSVLocus with two alleles types found in ingroup samples will be extracted
	 * 		1. presence of SV - non ref variant 
	 * 		2. absence of SV - reference type)
	 * if false, any SimpleSVLocus with one or more alleles of type presence of SV(non-ref variant) will be extracted
	 */
	private final boolean mustHaveTwoAlleleTypes;
	
	
	////////////////////////
	private List<SimpleSVLocus> extractedSVLocusList;
	
	
	public IngroupSVLocusExtractor(List<SimpleSVLocus> svLocusList,
			PopulationStructureFileReader dellyRunPopulationStructureFileReader, boolean mustHaveTwoAlleleTypes) {
		super();
		this.svLocusList = svLocusList;
		this.dellyRunPopulationStructureFileReader = dellyRunPopulationStructureFileReader;
		this.mustHaveTwoAlleleTypes = mustHaveTwoAlleleTypes;
		
		this.run();
	}





	void run() {
		this.extractedSVLocusList=new ArrayList<>();
		
		for(SimpleSVLocus locus:this.svLocusList) {
			boolean svPresenceAlleleFound=false;
			boolean svAbsenceAlleleFound=false;
			
			for(int ingroupSampleIndex:this.dellyRunPopulationStructureFileReader.getAllIngroupSampleIndices()) {
				if(locus.getSampleIndexGenotypeMap().get(ingroupSampleIndex).equals(Genotype.MISSING)){
					//
				}else if(locus.getSampleIndexGenotypeMap().get(ingroupSampleIndex).equals(Genotype.ABSENCE)) {
					svAbsenceAlleleFound=true;
				}else if(locus.getSampleIndexGenotypeMap().get(ingroupSampleIndex).equals(Genotype.PRESENCE)) {
					svPresenceAlleleFound=true;
				}else if(locus.getSampleIndexGenotypeMap().get(ingroupSampleIndex).equals(Genotype.HETER)) {
					svPresenceAlleleFound=true;
					svAbsenceAlleleFound=true;
				}
			}
			
			if(this.mustHaveTwoAlleleTypes) {
				if(svPresenceAlleleFound && svAbsenceAlleleFound) {
					this.extractedSVLocusList.add(locus);
				}
			}else {
				if(svPresenceAlleleFound) {
					this.extractedSVLocusList.add(locus);
				}
			}
		}
	}





	/**
	 * @return the extractedSVLocusList
	 */
	public List<SimpleSVLocus> getExtractedSVLocusList() {
		return extractedSVLocusList;
	}
	
	
}
