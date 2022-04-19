package phylo.vcf;

import java.io.File;

import htsjdk.variant.vcf.VCFFileReader;
import population.vcf.utils.GenotypeRecoderFactory;
import population.vcf.utils.VariantContextFilterFactory;

public class Vcf2AlignmentSeqBuilderTest {
	
	public static void main(String[] args) {
		
		File vcfFile=new File("/scratch/tanxu/reseq/test/sv2gr/batch/sb.chrom1-10.merged.no.missing.nj.bs.tree/tmp/1/1.vcf");
//		File vcfFile=new File("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv2gr_java\\vcf\\sorghum.1.prop.all.24.bicolor.all.sites.jointly.called.raw.first.200k.vcf");
//		File vcfFile=new File("C:\\Users\\tanxu\\Desktop\\scratch\\tmp_vcf\\1.vcf");
//		
		VCFFileReader reader = new VCFFileReader(vcfFile, false);
		
		//
		Vcf2AlignmentSeqBuilder.firstN=2000000;
		
		int maxMissingGenotype=0;
		int minQUALForVariantSites = 30;
		Vcf2AlignmentSeqBuilder builder = new Vcf2AlignmentSeqBuilder(
//				VariantContextFilterFactory.onlySNPSites().and(VariantContextFilterFactory.allIndividualsAreGenericHomozygous()).and(VariantContextFilterFactory.maxMissingCount(maxMissingGenotype)), 
				VariantContextFilterFactory.nonIndelSite().and(VariantContextFilterFactory.nonMixedTypeSite()).and(VariantContextFilterFactory.filterVariantSitesByQual(minQUALForVariantSites)).and(VariantContextFilterFactory.allIndividualsAreGenericHomozygous()).and(VariantContextFilterFactory.maxMissingCount(maxMissingGenotype)),
//				GenotypeRecoderFactory.biBaseSeqRecoder()
				GenotypeRecoderFactory.singleBaseSeqRecoder()
				);
		
		builder.run(reader);
		builder.print();
		
//		System.out.println(builder.getMultipleAlignment().getAlignmentLen());
	}
}
