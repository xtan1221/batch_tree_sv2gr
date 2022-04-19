package phylo.tree.dist.diff;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import htsjdk.variant.vcf.VCFFileReader;
import population.vcf.utils.GenotypeRecoder;
import population.vcf.utils.GenotypeRecoderFactory;
import population.vcf.utils.VariantContextFilter;
import population.vcf.utils.VariantContextFilterFactory;

public class PairwiseDiffMatrixFromVcfTest {

	
	public static void main(String[] args) throws Exception {
		//////////////////
//		Path bcfFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv2gr_java\\vcf\\sorghum.1.prop.all.24.bicolor.all.sites.jointly.called.raw.first.200k.vcf");
		Path bcfFile = Path.of("C:\\Users\\tanxu\\Desktop\\scratch\\tmp_vcf\\sorghum.1.prop.all.24.bicolor.all.sites.filtered.by.coverage.4-4X.filtered.by.coverage.first100k.vcf");
		/////////////////
		List<String> orderedSampleNameList = new ArrayList<>();
		VCFFileReader reader = new VCFFileReader(bcfFile, false);
		reader.getFileHeader().getSampleNamesInOrder().forEach(n->{
			orderedSampleNameList.add(n);
		});
		reader.close();
		
		////////////////////////////
		String regionListString = "";
		String uniqueRegionIdentifier = "1"; 
		
		VariantContextFilter variantFilter = VariantContextFilterFactory.nonIndelSite().and(VariantContextFilterFactory.nonMixedTypeSite()).and(VariantContextFilterFactory.filterVariantSitesByQual(30)).and(VariantContextFilterFactory.allIndividualsAreGenericHomozygous()).and(VariantContextFilterFactory.maxMissingCount(0)); 
		GenotypeRecoder genotypeRecoder = GenotypeRecoderFactory.singleBaseSeqRecoder();//GenotypeRecoder genotypeRecoder
		Path regionDiffMatrixOutDir = Path.of("C:\\Users\\tanxu\\Desktop\\scratch\\tmp_dist_out\\");
		Path tmpOutDir=Path.of("C:\\Users\\tanxu\\Desktop\\scratch\\tmp_dist_out\\");
		
		
		PairwiseDiffMatrixFromVcf2 collector = 
				new PairwiseDiffMatrixFromVcf2(
						orderedSampleNameList, bcfFile, regionListString, uniqueRegionIdentifier, 
						variantFilter, genotypeRecoder,
						regionDiffMatrixOutDir,  tmpOutDir);
		
		System.out.println(collector.call());
	}
}
