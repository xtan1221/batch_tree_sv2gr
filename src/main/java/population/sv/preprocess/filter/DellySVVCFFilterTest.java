package population.sv.preprocess.filter;

import java.nio.file.Path;

public class DellySVVCFFilterTest {
	static Path inputVCFFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\maize\\output\\delly\\maize.46.samples.merged.genotype.sv.filtered.vcf");
	static int maxSVLen = 20000;
	static double maxBoundaryUncertainty=0.5;
	static Path outputVCFFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\test\\filtered.maize.delly.sv.vcf");
	
	
	public static void main(String[] args) {
		DellySVVCFFilter filter= new DellySVVCFFilter(
				inputVCFFile,
				maxSVLen, maxBoundaryUncertainty, 
				outputVCFFile
				);
	}
}
