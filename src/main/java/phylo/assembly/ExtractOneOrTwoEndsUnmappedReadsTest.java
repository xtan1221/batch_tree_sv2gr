package phylo.assembly;

import java.nio.file.Path;

public class ExtractOneOrTwoEndsUnmappedReadsTest {

	
	static void test1() {
		Path inputBamFile = Path.of("/scratch/tanxu/reseq/sb/single_sample_hc_pipeline_all_32/result/bam_files/SRR562030/SRR562030.trimmed.coord.sorted.dup.marked.rg.added.bam");
		
		Path outputBamFile = Path.of("/scratch/tanxu/reseq/test/assembly/reads_extraction/unmapped/SRR562030/SRR562030.trimmed.coord.sorted.dup.marked.rg.added.not.both.ends.mapped.q10.bam");
		
		int minMapq = 10;
		
		ExtractOneOrTwoEndsUnmappedReads extractor = new ExtractOneOrTwoEndsUnmappedReads(inputBamFile, outputBamFile, minMapq);
		
		extractor.run();
	}
	
	static void test2() {
		Path inputBamFile = 
				Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\data processing and analysis\\SV detection\\result\\analysis\\SRR486615.trimmed.coord.sorted.dup.marked.rg.added.chr01.1-10000000.out.bam");
		Path outputBamFile = 
				Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\assembly\\reads extraction\\data\\test.out.bam");
		
		int minMapq = 10;
		
		ExtractOneOrTwoEndsUnmappedReads extractor = new ExtractOneOrTwoEndsUnmappedReads(inputBamFile, outputBamFile, minMapq);
		
		extractor.run();
	}
	
	
	public static void main(String[] args) {
		test1();
//		test2();
	}
}
