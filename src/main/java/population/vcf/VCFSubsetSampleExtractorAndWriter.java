package population.vcf;

import java.nio.file.Path;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Stream;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import population.vcf.utils.VariantContextUtils;


/**
 * create a new VCF file containing only a subset of samples from a given VCF file;
 * 
 * only loci with at least one of the target samples with called non-ref allele will be included in the output vcf file;
 * 
 * @author tanxu
 * 
 */
public class VCFSubsetSampleExtractorAndWriter {
	/**
	 * 
	 */
	private final Path inVCFFile;
	/**
	 * 
	 */
	private final List<String> targetSamples;
	/**
	 * additional filter for each original variant locus
	 */
	private final Predicate<VariantContext> locusFilter;
	/**
	 * 
	 */
	private final Path outputVCFFile;
	
	////////////////////////////////
	private VCFFileReader reader;
	private VariantContextWriter writer;
	
	/**
	 * 
	 * @param inVCFFile
	 * @param targetSamples
	 * @param locusFilter
	 * @param outputVCFFile
	 */
	public VCFSubsetSampleExtractorAndWriter(
			Path inVCFFile, List<String> targetSamples,
			Predicate<VariantContext> locusFilter, Path outputVCFFile) {
		super();
		this.inVCFFile = inVCFFile;
		this.targetSamples = targetSamples;
		this.locusFilter = locusFilter;
		this.outputVCFFile = outputVCFFile;
		
		/////////////////////
		if(this.outputVCFFile.toFile().exists()) {
			System.out.println("given outputVCFFile already exists, delete it...");
			this.outputVCFFile.toFile().delete();
		}
		
		this.prepare();
		this.run();
	}
	
	
	void prepare() {
		this.reader = new VCFFileReader(this.inVCFFile, false);
		//////////
		VariantContextWriterBuilder builder = 
				new VariantContextWriterBuilder().setOutputFile(this.outputVCFFile.toFile())
				.setReferenceDictionary(this.reader.getFileHeader().getSequenceDictionary()).setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
		this.writer =builder.build();
		/////////////////modify the header (line containing the sample names)
		VCFHeader header = this.reader.getHeader();
		VCFHeader newHeader=new VCFHeader(header.getMetaDataInInputOrder(), this.targetSamples); //new VCFHeader with the target samples only
		this.writer.writeHeader(newHeader); //write the header section of the vcf file
	}
	
	void run() {
		
		
		////////////////process each locus
		Stream<VariantContext> stream=reader.iterator().stream();
		stream.forEach(vc->{
			if(!this.locusFilter.test(vc)) {
				return;
			}
			
			///////////////
			VariantContext extractedVC=VariantContextUtils.subsetGenotypes(vc, targetSamples);
			
			if(extractedVC!=null) {
				this.writer.add(extractedVC);
			}
		});
		
		this.reader.close();
		this.writer.close();
	}
	
}
