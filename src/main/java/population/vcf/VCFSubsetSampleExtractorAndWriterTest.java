package population.vcf;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Predicate;

import htsjdk.variant.variantcontext.VariantContext;

public class VCFSubsetSampleExtractorAndWriterTest {
	
	
	
	public static void main(String[] args) {
		Path inVCF=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\sorghum\\output\\lumpy\\27.sorghum.samples.joint.sv.lumpy.breakpoint.type.sv.filtered.vcf");
		
		Path outVCF=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\sorghum\\output\\test.sv.vcf");
		
		Predicate<VariantContext> locusFilter = e->{return true;};
		
		List<String> samples=new ArrayList<>();
		samples.add("3237");
		samples.add("11282");
		
		VCFSubsetSampleExtractorAndWriter ew=new VCFSubsetSampleExtractorAndWriter(
				inVCF,//Path inVCFFile, 
				samples,//List<String> targetSamples,
				locusFilter,//Predicate<VariantContext> locusFilter, 
				outVCF//Path outputVCFFile
				);
	}
}
