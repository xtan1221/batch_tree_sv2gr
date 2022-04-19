package phylo.ref;

import java.io.IOException;
import java.nio.file.Path;

public class RegionFileWriterTest {
	public static void main(String[] args) throws IOException {
//		Path chromLenFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv2gr_java\\data\\Sbicolor_454_v3.0.1.chrom.length.txt");
//		Path chromLenFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\phylogeny\\rice\\build_full_genome_chrom_trees\\region_files\\Osativa_323_v7.0.chrom.length.txt");
		Path chromLenFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\phylogeny\\maize\\build_batch_local_tree\\region_files\\chrom1-10.length.txt");
		
		Path outputRegionFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\phylogeny\\maize\\build_batch_local_tree\\region_files\\zm.RefGen_v5.chrom1-10.windowsize.500k.txt");
		
		////////////////
		RefPartition rp = new RefPartition(
				chromLenFile, //Path chromLengthFile, 
				10000000,//int minChromLen, 
				500000,//int windowSize, 
				true//boolean toKeepTrailingShortWindow
				);
		
		rp.getWindowList().forEach(r->{
			System.out.println(r);
		});
		
		///////////////////////////
		RegionFileWriter writer = new RegionFileWriter(outputRegionFile);
		
		writer.run(rp.getWindowList());
	}
}
