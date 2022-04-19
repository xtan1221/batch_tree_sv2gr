package population.sv.preprocess.merge;

import java.nio.file.Path;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Stream;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import population.sv.preprocess.filter.DellySVVCFFilter;
import population.sv.preprocess.filter.LumpySVVCFFilter;
import population.vcf.ReorderSamples;


/**
 * the SURVIVOR output vcf file only contains the SV IDs from input VCF files without the full genotype information of all samples
 * 		the IDs from input SV vcf files should be all unique, otherwise, error will occur!!!!
 * 
 * thus, the IDs should be retrieved to extract the full SV records from the input VCF files to build a VCF file containing the merged SVs
 * 
 * 
 * @author tanxu
 *
 */
public class ExtractSURVIVORMergedSVByID {
	/**
	 * merged SVs from {@link #dellyCalledSvVcfFile} and {@link #lumpyCalledSvVcfFile} by SURVIVOR;
	 * 
	 * the 'ID' column should be used to identify whether the SV is from the Delly caller or LUMPY caller and extract the corresponding SV to output to {@link #outVcfFile}
	 */
	private final Path mergedSvIDVcfFile;
	
	/**
	 * the original SV vcf file called by Delly that is used as input for SURVIVOR;
	 * 
	 * more specifically, the SV vcf file that is filtered with {@link DellySVVCFFilter}
	 */
	private final Path dellyCalledSvVcfFile;
	
	/**
	 * the original SV vcf file called by LUMPY that is used as input for SURVIVOR
	 * 
	 * specifically, the VCF file should be processed to 
	 * 1. filter out 'BND' type SVs and 
	 * 2. genotyped by SVTyper, 
	 * 3. then filtered by {@link LumpySVVCFFilter}
	 * 4. samples reordered by {@link ReorderSamples} to be consistent with the sample orders in {@link #dellyCalledSvVcfFile}
	 */
	private final Path lumpyCalledSvVcfFile;

	/**
	 * output VCF file contains the merged SVs by SURVIVOR extracted from {@link #dellyCalledSvVcfFile} and {@link #lumpyCalledSvVcfFile} using the 'ID' column from {@link #mergedSvIDVcfFile}
	 */
	private final Path outVcfFile;
	
	
	//////////////////////
	private VCFFileReader reader;
	private VariantContextWriter writer;
	private Set<String> mergedSvIDs;
	private Set<String> mergedSVIDs2;
	
	public ExtractSURVIVORMergedSVByID(
			Path mergedSvIDVcfFile, 
			Path dellyCalledSvVcfFile, Path lumpyCalledSvVcfFile,
			Path outVcfFile) {
		super();
		this.mergedSvIDVcfFile = mergedSvIDVcfFile;
		this.dellyCalledSvVcfFile = dellyCalledSvVcfFile;
		this.lumpyCalledSvVcfFile = lumpyCalledSvVcfFile;
		this.outVcfFile = outVcfFile;
		
		if(this.outVcfFile.toFile().exists()) {
			System.out.println("outVcfFile already exists, delete it...");
			this.outVcfFile.toFile().delete();
		}
		
		//////////////////////////////////
		this.readMergedSVIDs();
		this.prepareOutVCF();
		this.readDellySVVCFFile();
		this.readLumpySVVCFFile();
		this.close();
		this.report();
	}

	void readMergedSVIDs() {
		this.mergedSvIDs=new HashSet<>();
		this.mergedSVIDs2=new HashSet<>();
		
		this.reader = new VCFFileReader(this.mergedSvIDVcfFile, false);
		
		////////////////process each locus
		Stream<VariantContext> stream=this.reader.iterator().stream();
		stream.forEach(vc->{
			String id=vc.getID();
			
			if(this.mergedSvIDs.contains(id))
				throw new IllegalArgumentException("duplicate SV ID found:"+id);
			
			this.mergedSvIDs.add(id);
			this.mergedSVIDs2.add(id);
		});
		
		this.reader.close();
	}
	
	/**
	 * write the header section to output vcf file with the header section of {@link #dellyCalledSvVcfFile}
	 */
	void prepareOutVCF() {
		System.out.println("write out the header section of output vcf file...");
		this.reader = new VCFFileReader(this.dellyCalledSvVcfFile, false);
		//////////
		VariantContextWriterBuilder builder = 
				new VariantContextWriterBuilder().setOutputFile(this.outVcfFile.toFile())
				.setReferenceDictionary(this.reader.getFileHeader().getSequenceDictionary()).setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
		
		this.writer =builder.build();
		
		/////////////////modify the header (line containing the sample names)
		
		VCFHeader header = this.reader.getHeader();
		this.writer.writeHeader(header); //write the header section of the vcf file
		
		this.reader.close();
	}
	
	/**
	 * must be invoked after {@link #prepareOutVCF()}
	 */
	void readDellySVVCFFile() {
		this.reader = new VCFFileReader(this.dellyCalledSvVcfFile, false);
		
		////////////////process each locus
		Stream<VariantContext> stream=this.reader.iterator().stream();
		stream.forEach(vc->{
			String id=vc.getID();
			
			if(this.mergedSvIDs.contains(id)) {
				this.writer.add(vc);
				this.mergedSVIDs2.remove(id);
			}
				
			
		});
		
		this.reader.close();
	}
	
	
	
	void readLumpySVVCFFile() {
		this.reader = new VCFFileReader(this.lumpyCalledSvVcfFile, false);
		
		////////////////process each locus
		Stream<VariantContext> stream=this.reader.iterator().stream();
		stream.forEach(vc->{
			String id=vc.getID();
			
			if(this.mergedSvIDs.contains(id)) {
				this.writer.add(vc);
				this.mergedSVIDs2.remove(id);
			}
				
		});
		
		this.reader.close();
	}
	
	
	/**
	 * close {@link #writer}
	 */
	void close() {
		this.writer.close();
	}
	
	void report() {
		System.out.println("unidentified SVs:"+this.mergedSVIDs2.size());
		this.mergedSVIDs2.forEach(i->{
			System.out.println(i);
		});
	}
}
