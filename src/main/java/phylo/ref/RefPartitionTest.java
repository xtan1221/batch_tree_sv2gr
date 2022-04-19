package phylo.ref;

import java.nio.file.Path;

import phylo.batch.BcfUtils;

public class RefPartitionTest {
	public static void main(String[] args) {
		Path chromLenFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv2gr_java\\data\\Sbicolor_454_v3.0.1.chrom.length.txt");
		
		RefPartition rp = new RefPartition(
				chromLenFile, //Path chromLengthFile, 
				10000000,//int minChromLen, 
				100000000,//int windowSize, 
				true//boolean toKeepTrailingShortWindow
				);
		
		rp.getWindowList().forEach(r->{
			System.out.println(r);
		});
		
		System.out.println(BcfUtils.buildBcftoolsRegionString(rp.getWindowList()));
	}
}
