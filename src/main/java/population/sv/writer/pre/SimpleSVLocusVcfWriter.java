package population.sv.writer.pre;

import java.nio.file.Path;
import java.util.Collections;
import java.util.List;

import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import population.sv.utils.SimpleSVLocus;

/**
 * write a set of {@link SimpleSVLocus} with non-null {@link VariantContext} to a vcf file
 * 
 * 
 * @author tanxu
 *
 */
public class SimpleSVLocusVcfWriter {
	private final List<SimpleSVLocus> svLocusList;
	private final VCFFileReader VCFFileReader;
	private final Path outVCFFile;
	
	////////////////////////////////
	private VariantContextWriter writer;
	
	
	public SimpleSVLocusVcfWriter(List<SimpleSVLocus> svLocusList, VCFFileReader VCFFileReader, Path outVCFFile) {
		super();
		this.svLocusList = svLocusList;
		this.VCFFileReader = VCFFileReader;
		this.outVCFFile = outVCFFile;
		
		if(this.outVCFFile.toFile().exists()) {
			System.out.println("given outVCFFile exists, delete it...");
			this.outVCFFile.toFile().delete();
		}
		
		
		///////////////////////
		this.sortSimpleSVLocus();
		this.run();
	}


	void sortSimpleSVLocus() {
		Collections.sort(this.svLocusList, (a,b)->{
			if(a.getChrom().equals(b.getChrom())) {
				return a.getStart()-b.getStart();
			}else {
				return a.getChrom().compareTo(b.getChrom());
			}
		});
		
		
		
	}

	void run() {
		VariantContextWriterBuilder builder = 
				new VariantContextWriterBuilder().setOutputFile(this.outVCFFile.toFile())
				.setReferenceDictionary(this.VCFFileReader.getFileHeader().getSequenceDictionary()).setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);

		this.writer =builder.build();
		
		this.writer.writeHeader(this.VCFFileReader.getHeader());
		
		for(SimpleSVLocus locus:this.svLocusList) {
			if(locus.getVariantContext()==null)
				throw new IllegalArgumentException("no VariantContext found in SimpleSVLocus!");
			writer.add(locus.getVariantContext());
		}
		
		writer.close();
	}
	
	
	
}
