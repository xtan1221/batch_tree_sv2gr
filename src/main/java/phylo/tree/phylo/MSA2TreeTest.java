package phylo.tree.phylo;

import java.io.File;
import java.nio.file.Path;

import htsjdk.variant.vcf.VCFFileReader;
import phylo.alignment.MultipleAlignment;
import phylo.exception.NewickTreeFileNotCreatedException;
import phylo.tree.phylo.PhylipBashRunner.DistanceModel;
import phylo.tree.reader.NewickFileFormatType;
import phylo.vcf.Vcf2AlignmentSeqBuilder;
import population.vcf.utils.GenotypeRecoderFactory;
import population.vcf.utils.VariantContextFilterFactory;

public class MSA2TreeTest {

	
	
	public static void testPhylipTreeBuilder(MultipleAlignment msa) {
		
		String dataName = "test";
		Path tmpDir = Path.of("/scratch/tanxu/reseq/test/sv2gr/");
		DistanceModel distanceModel = DistanceModel.Jukes_Cantor;
		
		
		PhylipBasedTreeBuilder builder = new PhylipBasedTreeBuilder(distanceModel);
		
		try {
			builder.build(dataName, msa, tmpDir);
			System.out.println(builder.getTree().toFullNewickString(NewickFileFormatType.SIMPLE_NEWICK_2));
		} catch (NewickTreeFileNotCreatedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}
	
	public static void testMegaXTreeBuilder(MultipleAlignment msa) {
		
		String dataName = "test";
		Path tmpDir = Path.of("/scratch/tanxu/reseq/test/sv2gr/megax/");
		Path maoFile = Path.of("/home/tanxu/phylogeny/megaX/infer_NJ_nucleotide_pairwise_deletion.mao");
		
		MegaXBasedTreeBuilder builder = new MegaXBasedTreeBuilder(maoFile);
		
		try {
			builder.build(dataName, msa, tmpDir);
			System.out.println(builder.getTree().toFullNewickString(NewickFileFormatType.SIMPLE_NEWICK_2));
		} catch (NewickTreeFileNotCreatedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	
	public static void main(String[] args) {
		File vcfFile=new File("/scratch/tanxu/reseq/sb/gatk_joint_snp_genotype_calling/result/sorghum.1.prop.all.24.bicolor.all.sites/sorghum.1.prop.all.24.bicolor.all.sites.jointly.called.raw.first.1M.vcf");
//		File vcfFile=new File("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv2gr_java\\vcf\\sorghum.1.prop.all.24.bicolor.all.sites.jointly.called.raw.first.200k.vcf");
//		File vcfFile=new File("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv2gr_java\\vcf\\toy.vcf");
//		
		VCFFileReader reader = new VCFFileReader(vcfFile, false);
		
		//
		Vcf2AlignmentSeqBuilder.firstN=1000;
		
		int maxMissingGenotype=0;
		int minQUALForVariantSites = 30;
		Vcf2AlignmentSeqBuilder builder = new Vcf2AlignmentSeqBuilder(
//				VariantContextFilterFactory.onlySNPSites().and(VariantContextFilterFactory.allIndividualsAreGenericHomozygous()).and(VariantContextFilterFactory.maxMissingCount(maxMissingGenotype)), 
				VariantContextFilterFactory.nonIndelSite().and(VariantContextFilterFactory.filterVariantSitesByQual(minQUALForVariantSites)).and(VariantContextFilterFactory.allIndividualsAreGenericHomozygous()).and(VariantContextFilterFactory.maxMissingCount(maxMissingGenotype)),
//				GenotypeRecoderFactory.biBaseSeqRecoder()
				GenotypeRecoderFactory.singleBaseSeqRecoder()
				);
		
		builder.run(reader);
		
		
//		testPhylipTreeBuilder(builder.getMultipleAlignment());
		testMegaXTreeBuilder(builder.getMultipleAlignment());
	}
}
