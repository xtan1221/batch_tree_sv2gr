package population.vcf;

import java.nio.file.Path;

/**
 * 
 * @author tanxu
 * 
 */
public class LumpyRawJointCalledVCFFilePostprocessorTest {
	
	
	
	public static void main(String[] args) {
		//sorghum
//		Path inVCF = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\sorghum\\output\\lumpy\\27.sorghum.samples.joint.sv.lumpy.vcf");
//		Path outVCF = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\sorghum\\output\\lumpy\\27.sorghum.samples.joint.sv.lumpy.sample.reordered.to.delly.breakpoint.filtered.vcf");
		
		
		//rice
		Path inVCF = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\rice\\output\\lumpy\\40.rice.samples.joint.sv.lumpy.vcf");
		Path outVCF = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\rice\\output\\lumpy\\40.rice.samples.joint.sv.lumpy.sample.reordered.to.delly.breakpoint.filtered.vcf");
		
		//maize
//		Path inVCF = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\maize\\output\\lumpy\\46.maize.samples.joint.sv.lumpy.vcf");
//		Path outVCF = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\maize\\output\\lumpy\\46.maize.samples.joint.sv.lumpy.sample.reordered.to.delly.breakpoint.filtered.vcf");
				
		
		LumpyRawJointCalledVCFFilePostprocessor r=new LumpyRawJointCalledVCFFilePostprocessor(inVCF, outVCF);
		
		
		
	}
}
