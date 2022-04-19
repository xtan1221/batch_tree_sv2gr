package phylo.tree.fullgenome;

import java.nio.file.Path;

public class FullGenomeAndChromTreeFromBcfFileManagerTest {
	
	public static void main(String[] args) {
		Path bcfFile = Path.of("/scratch/tanxu/reseq/sb/gatk_joint_snp_genotype_calling/result/sorghum.1.prop.all.24.bicolor.all.sites/filtered/sorghum.1.prop.all.24.bicolor.all.sites.filtered.by.coverage.4-4X.filtered.by.coverage.bcf");
		Path regionBedFile = Path.of("/scratch/tanxu/reseq/sb/reference/v3.1.1/assembly/Sbicolor_454_v3.0.1.1-10.chrom.bed");
		Path rootOutDir = Path.of("/scratch/tanxu/reseq/test/sv2gr/full_genome_chrom_tree/jc_nj_tree");
		String outgroupName = "3237";
		int regionLen = 1000000;
		int minQUALForVariantSites = 30; 
		int threadNum = 20;
		
		String[] ssss = new String[7];
		ssss[0]=bcfFile.toString();
		ssss[1]=regionBedFile.toString();
		ssss[2]=rootOutDir.toString();
		ssss[3]=outgroupName;
		ssss[4]=Integer.toString(regionLen);
		ssss[5]=Integer.toString(minQUALForVariantSites);
		ssss[6]=Integer.toString(threadNum);
//		FullGenomeAndChromTreeFromBcfFileManager manager = 
//				new FullGenomeAndChromTreeFromBcfFileManager(
//						bcfFile, regionBedFile, rootOutDir, outgroupName, 
//						regionLen, minQUALForVariantSites, threadNum);
//		manager.run();
		
		FullGenomeAndChromTreeFromBcfFileManager.main(ssss);
	}
	
	
}
