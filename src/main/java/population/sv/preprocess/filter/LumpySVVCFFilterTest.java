package population.sv.preprocess.filter;

import java.nio.file.Path;

public class LumpySVVCFFilterTest {
	static Path inputVCFFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\maize\\output\\lumpy\\46.maize.samples.joint.sv.lumpy.breakpoint.type.sv.filtered.svtyper.genotyped.vcf");
	static int maxSVLen = 20000;
	static double minQual = 20;
	static double maxBoundaryUncertainty = 0.5;
	static Path outputVCFFile=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\test\\filtered.sv.vcf");
	
	public static void main(String[] args) {
		LumpySVVCFFilter filter= new LumpySVVCFFilter(
				inputVCFFile,
				maxSVLen, minQual, maxBoundaryUncertainty, 
				outputVCFFile
				);
	}
}
