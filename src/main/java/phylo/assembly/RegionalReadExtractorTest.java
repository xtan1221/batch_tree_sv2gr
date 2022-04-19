package phylo.assembly;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import phylo.ref.Region;

public class RegionalReadExtractorTest {
	
	public static void main(String[] args) {
		Path inputSortedBamFile = 
				Path.of("/scratch/tanxu/reseq/sb/single_sample_hc_pipeline_all_32/result/bam_files/SRR486615/SRR486615.trimmed.coord.sorted.dup.marked.rg.added.bam");
		List<Region> regionList = new ArrayList<>(); 
		regionList.add(new Region("Chr01", 1972000, 1983000));
		Path outputBamFile = Path.of("/scratch/tanxu/reseq/test/assembly/reads_extraction/SRR486615.Chr01.1972000.1983000.bam");
		int minMapq = 10;
		
		RegionalReadExtractor extractor = new RegionalReadExtractor(inputSortedBamFile, regionList, outputBamFile, minMapq);
		
		extractor.run();
	}
}
