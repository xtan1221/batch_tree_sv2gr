package population.popGeneticsStatistics.nucdiv;

import java.nio.file.Path;

public class PiMeanAndSDCalculatorTest {

	
	public static void main(String[] args) {
		Path file = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\population_genetic_statistics\\nucleotide_diversity\\testing_data\\12.indica.all.sites.filtered.QUAL30.GQ10.NonMissing3.coverage4-10X.non.variant.and.SNP.sites.windowed.pi");
		
//		Path file = Path.of("/scratch/tanxu/reseq/rice/gatk_joint_snp_genotype_calling/result/12.indica.all.sites.from.scratch/filtered/pi/vcftools/12.indica.all.sites.filtered.QUAL30.GQ10.NonMissing3.coverage3-10X.non.variant.and.SNP.sites.sites.pi");
		
		
		boolean onlyIncludeNon0Values = false;
		int targetColIndex = 4;
		
		
		PiMeanAndSDCalculator ca=new PiMeanAndSDCalculator(file, onlyIncludeNon0Values, targetColIndex);
	}
}
