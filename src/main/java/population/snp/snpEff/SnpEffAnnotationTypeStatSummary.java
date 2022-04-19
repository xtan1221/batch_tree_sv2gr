package population.snp.snpEff;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import population.vcf.utils.VariantContextFilterFactory;

/**
 * find out the number of sites with at least one target sample have non-reference allele for each type of SnpEff annotation type;
 * 
 * @author tanxu
 *
 */
public class SnpEffAnnotationTypeStatSummary {
	private final Path vcfFile;
	
	/**
	 * only sites with at least one sample in this given list will be included in the summary;
	 */
	private final List<Integer> targetSampleIndexList;
	
	private final Path outputStatSummaryFileDir;
	/////////////////////////////////
	
	private Path outputStatSummaryFile;
	/**
	 * 
	 */
	private Map<String, Integer> annotatedTypeSiteNumMap;
	
	public SnpEffAnnotationTypeStatSummary(
			Path vcfFile, 
			List<Integer> targetSampleIndexList, 
			Path outputStatSummaryFileDir) {
		super();
		this.vcfFile = vcfFile;
		this.targetSampleIndexList = targetSampleIndexList;
		this.outputStatSummaryFileDir = outputStatSummaryFileDir;
		
		//////////////////////////////
		this.prepare();
		this.run();
		try {
			this.writeToFile();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	void prepare() {
//		//////////////////////////
		this.outputStatSummaryFile = Path.of(this.outputStatSummaryFileDir.toString(),"summary.statistics.txt");
		this.checkOutputFile(outputStatSummaryFile);
		
	}
	
	private void checkOutputFile(Path outFile) {
		if(outFile.toFile().exists()) {
			System.out.println("output file "+outFile.toString()+" already exists, delete it..." );
			outFile.toFile().delete();
		}
	}
	
	
	private int counter=0;
	void run() {
		this.annotatedTypeSiteNumMap=new HashMap<>();
		
		VCFFileReader reader = new VCFFileReader(this.vcfFile, false);
		Stream<VariantContext> stream=reader.iterator().stream();
		
		stream.forEach(vc->{
			if( //must be biallelic SNP sites with good quality
					VariantContextFilterFactory.nonIndelSite()//not indel, can be skipped since only SNP sites are included in the vcf file annotated by SnpEff
					.and(VariantContextFilterFactory.nonMixedTypeSite())//not mixed type
					.and(VariantContextFilterFactory.notLowQual()) //FILTER column does not contain 'LowQual'
					.and(VariantContextFilterFactory.biallelicSite()) //must be biallelic!!!!!!!!!!!!!!
					.test(vc)) {
				//passes all filters
			}else {
				return;//skip
			}
			
			//debug
//			if(vc.getCommonInfo().getAttribute("ANN").toString().contains("synonymous_variant")) {
//				if(!SnpEffVariantContextUtils.getSNPAnnotatedTypes(vc).contains("synonymous_variant")) {
//					System.out.println(vc.toString());
//				}
//			}
//			System.out.println(vc.toString());
			
			counter++;
			if(counter%1000000==0) {
				int num=counter/1000000;
				System.out.println(num+" M records found and processed!");
			}
			
			//check if at least one target sample has non-reference allele
			boolean targetSampleWithNonRefAlleleFound=false;
			for(int index:this.targetSampleIndexList) {
//				System.out.println(index);
				Genotype gt=vc.getGenotype(index);
				if(!gt.isNoCall()) {//
					if(!gt.isHomRef()) {//not homo reference
						targetSampleWithNonRefAlleleFound=true;
						break;
					}
				}
			}
			
			//
			if(targetSampleWithNonRefAlleleFound) {
				for(String annotatedType: SnpEffVariantContextUtils.getSNPAnnotatedTypes(vc)) {
					if(annotatedType.isEmpty()) {//skip empty annotated types
						continue;
					}
					
					if(!this.annotatedTypeSiteNumMap.containsKey(annotatedType)) {
						this.annotatedTypeSiteNumMap.put(annotatedType, 0);
					}
					//
					this.annotatedTypeSiteNumMap.put(annotatedType, this.annotatedTypeSiteNumMap.get(annotatedType)+1);
				}
			}
		});
		
		reader.close();
	}
	
	
	void writeToFile() throws IOException {
		////////////
		//write to summary file
		BufferedWriter summaryWriter = new BufferedWriter(new FileWriter(this.outputStatSummaryFile.toString()));
		
		for(String type: this.annotatedTypeSiteNumMap.keySet()) {
			summaryWriter.append(type+"\t"+this.annotatedTypeSiteNumMap.get(type));
			summaryWriter.newLine();
		}
		
		summaryWriter.flush();
		summaryWriter.close();
	}
}
