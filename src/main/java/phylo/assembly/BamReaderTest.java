package phylo.assembly;

import java.nio.file.Path;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamReaderTest {
	static String readname="SRR486615.108695557";
	static Path inputSortedBamFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\data processing and analysis\\SV detection\\result\\analysis\\SRR486615.trimmed.coord.sorted.dup.marked.rg.added.chr01.1-1000000.bam");
	
	
	public static void test() {
		
		SamReader reader = SamReaderFactory.makeDefault().open(inputSortedBamFile.toFile());
		SamReader reader2 = SamReaderFactory.makeDefault().open(inputSortedBamFile.toFile());
		
		reader.iterator().forEachRemaining(r->{
			if(r.getSupplementaryAlignmentFlag()) {
				
				SAMRecord mate = reader2.queryMate(r);
				
				if(mate!=null) {
					SAMRecord mate2 = SAMRecordQueryUtils.queryNonSupplementaryMate(reader2, mate); //reader2.queryMate(mate); 
					if(mate2!=null) {
						System.out.println("=======================");
						System.out.println(r.getSAMString().trim());
						System.out.println(mate.getSAMString().trim());
						System.out.println(mate2.getSAMString().trim());
//						if(mate2.getSupplementaryAlignmentFlag()) {
//							System.out.println("=======================");
//							System.out.println(r.getSAMString().trim());
//							System.out.println(mate.getSAMString().trim());
//							System.out.println(mate2.getSAMString().trim());
//						}
					}
				}
			}
		});
		
	}
	
	
	public static void main(String[] args) {
		test();
	}
}
