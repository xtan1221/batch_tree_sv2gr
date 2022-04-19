package population.snp.snpEff;

import java.nio.file.Path;
import java.util.List;

import population.popHierarchy.PopulationStructureFileReader;

public class SnpEffAnnotationTypeStatSummaryTest {
	
	public static void main(String[] args) {
//		Path vcfFile=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\snp annotation\\snpEff\\sorghum\\output\\sorghum.27.samples.filtered.snp.sites.chrom1-10.snpEff.annotation.first200k.lines.vcf");
		Path vcfFile=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\snp annotation\\snpEff\\sorghum\\output\\sorghum.27.samples.filtered.snp.sites.snpEff.annotation.first10000lines.vcf");
		Path populationHierachyTableFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\snp annotation\\snpEff\\sorghum\\output\\sorghum.27.sample.SnpEff.sample_index.population_level.table");		
		PopulationStructureFileReader populationStructureFileReader=new PopulationStructureFileReader(populationHierachyTableFile,true);
		//
		List<Integer> ingroupSamples = populationStructureFileReader.getOrderedSampleIndices();
		System.out.println(ingroupSamples);
		Path outputDir = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\snp annotation\\snpEff\\sorghum\\output\\test");
		
		SnpEffAnnotationTypeStatSummary summary=new SnpEffAnnotationTypeStatSummary(
				vcfFile,
				ingroupSamples,
				outputDir
				);
	}
}
