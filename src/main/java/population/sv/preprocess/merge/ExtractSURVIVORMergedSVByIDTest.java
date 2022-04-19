package population.sv.preprocess.merge;

import java.nio.file.Path;

public class ExtractSURVIVORMergedSVByIDTest {
	static Path mergedSvIDVcfFile=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\test\\27.sorghum.delly.lumpy.sv.union.merged.50.bp.min.2.bp.vcf");
	static Path dellyCalledSvVcfFile=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\sorghum\\output\\delly\\27.sorghum.samples.merged.genotype.sv.filtered.all.filter.applied.vcf");
	static Path lumpyCalledSvVcfFile=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\sorghum\\output\\lumpy\\27.sorghum.samples.joint.sv.lumpy.svtyper.genotyped.all.filter.applied.sample.reordered.vcf");
	static Path outVcfFile=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\test\\27.sorghum.delly.lumpy.sv.merged.vcf");
	
	public static void main(String[] args) {
		ExtractSURVIVORMergedSVByID merge = new ExtractSURVIVORMergedSVByID(
				mergedSvIDVcfFile,
				dellyCalledSvVcfFile,
				lumpyCalledSvVcfFile,
				outVcfFile
				);
		
	}
}
