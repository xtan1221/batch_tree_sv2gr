package genomics.geneAnnotation;

import java.nio.file.Path;

public class TranscriptGenomicRegionBedReaderTest {

	
	public static void main(String[] args) {
		
		Path transcriptGenomicRegionBedFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\data processing and analysis\\input data\\sorghum\\reference\\Sbicolor\\v3.1.1\\annotation\\transcript.genomic.region.bed");
		
		TranscriptGenomicRegionBedReader reader = new TranscriptGenomicRegionBedReader(transcriptGenomicRegionBedFile);
		
		System.out.println(reader.getTranscriptIDGenomicRegionMap().size());
	}
}
