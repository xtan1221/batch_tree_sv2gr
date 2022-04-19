package phylo.vcf;

import java.io.File;
import java.nio.file.Path;

import htsjdk.variant.vcf.VCFFileReader;
import phylo.alignment.AligmentFileWriter;
import phylo.alignment.FastaAlignmentFileWriter;
import phylo.alignment.MultipleAlignment;
import population.vcf.utils.GenotypeRecoderFactory;
import population.vcf.utils.VariantContextFilterFactory;

public class VCF2FastaTest {

	
	public static void main(String[] args) {
		File vcfFile=new File("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv2gr_java\\vcf\\sorghum.1.prop.all.24.bicolor.all.sites.jointly.called.raw.first.200k.vcf");
//		File vcfFile=new File("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv2gr_java\\vcf\\toy.vcf");
//		
		VCFFileReader reader = new VCFFileReader(vcfFile, false);
		
		//
		Vcf2AlignmentSeqBuilder.firstN=10000;
//		
		int maxMissingGenotype=0;
		
		Vcf2AlignmentSeqBuilder builder = new Vcf2AlignmentSeqBuilder(
//				VariantContextFilterFactory.onlySNPSites().and(VariantContextFilterFactory.allIndividualsAreGenericHomozygous()).and(VariantContextFilterFactory.maxMissingCount(maxMissingGenotype)), 
				VariantContextFilterFactory.nonIndelSite().and(VariantContextFilterFactory.allIndividualsAreGenericHomozygous()).and(VariantContextFilterFactory.maxMissingCount(maxMissingGenotype)),
//				GenotypeRecoderFactory.biBaseSeqRecoder()
				GenotypeRecoderFactory.singleBaseSeqRecoder()
				);
		
		builder.run(reader);
		builder.print();
		
		Path outputFilePath = Path.of("C:\\Users\\tanxu\\Desktop\\scratch\\tmp_output");
		String outFileBaseName = "sorghum.1.prop.all.24.bicolor.all.sites.jointly.called.raw.first.10000";
		AligmentFileWriter writer = new FastaAlignmentFileWriter(
				new MultipleAlignment(builder.getSampleNameSequenceStringMap()), 
				outputFilePath, 
				outFileBaseName);
		
		writer.write();
		
	}
}
