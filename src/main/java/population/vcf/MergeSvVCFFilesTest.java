package population.vcf;

import java.nio.file.Path;

public class MergeSvVCFFilesTest {
	
	static Path inputSVVcfFilePathListFile=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\test\\input.vcf.file.path.list.txt");
	static Path outputMergedVCFFile=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\test\\merged.vcf");
	
//	static Path inputSVVcfFilePathListFile=Path.of("/scratch/tanxu/reseq/maize/sv/lumpy/bam_with_RefGen_v5/joint_called_by_lumpy/file.path.list.txt");
//	static Path outputMergedVCFFile=Path.of("/scratch/tanxu/reseq/maize/sv/lumpy/bam_with_RefGen_v5/joint_called_by_lumpy/46.maize.samples.joint.sv.lumpy.breakpoint.type.sv.filtered.genotyped.vcf");
	
	
	
	
	public static void main(String[] args) {
		MergeSvVCFFiles m = new MergeSvVCFFiles(inputSVVcfFilePathListFile, outputMergedVCFFile);
	}
	
}
