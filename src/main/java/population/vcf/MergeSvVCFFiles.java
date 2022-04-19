package population.vcf;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Stream;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/**
 * merge a list of SV vcf files into a single one by removing any duplicate SVs
 * 
 * 1. same chrom
 * 2. same SV type
 * 3. same start and end coordinates
 * 
 * 
 * the VCF header of the first input vcf file will be used as the header of the output merged vcf file
 * 
 * @author tanxu
 *
 */
public class MergeSvVCFFiles {
	/**
	 * each line contains the full path to one input vcf file of SVs
	 * the first one must contain the reference dictionary lines as following
	 * ##reference=/scratch/tanxu/reseq/sb/reference/v3.1.1/assembly/Sbicolor_454_v3.0.1.fa
		##contig=<ID=Chr01,length=80884392>
		##contig=<ID=Chr02,length=77742459>
		##contig=<ID=Chr03,length=74386277>
		##contig=<ID=Chr04,length=68658214>
		##contig=<ID=Chr05,length=71854669>
		##contig=<ID=Chr06,length=61277060>
		... ...
	 */
	private final Path inputSVVcfFilePathListFile;
	
	/**
	 * 
	 */
	private final Path outputMergedVCFFile;
	
	////////////////////////////////
	private List<Path> inputSVVcfFiles;
	private Set<SVRecord> records;
	
	////////////////////////////////
	private VCFFileReader inputVCFReader;
	private VariantContextWriter writer;

	
	public MergeSvVCFFiles(Path inputSVVcfFilePathListFile, Path outputMergedVCFFile) {
		super();
		this.inputSVVcfFilePathListFile = inputSVVcfFilePathListFile;
		this.outputMergedVCFFile = outputMergedVCFFile;
		
		if(this.outputMergedVCFFile.toFile().exists()) {
			System.out.println("outputMergedVCFFile already exists, delete it...");
			this.outputMergedVCFFile.toFile().delete();
		}
		
		/////////////////////
		this.readInInputVCFFilePaths();
		this.prepareOutVCF();
		this.run();
	}

	void readInInputVCFFilePaths() {
		System.out.println("read in input vcf file paths...");
		this.inputSVVcfFiles = new ArrayList<>();
		
		try{////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(this.inputSVVcfFilePathListFile.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty()) { //first line is id	family	order
					continue;
				}
				
				this.inputSVVcfFiles.add(Path.of(line));
			}
	
			lineReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/**
	 * write the header section to output vcf file with the header section of the first input vcf file
	 */
	void prepareOutVCF() {
		System.out.println("write out the header section of output vcf file...");
		this.inputVCFReader = new VCFFileReader(this.inputSVVcfFiles.get(0), false);
		//////////
		VariantContextWriterBuilder builder = 
				new VariantContextWriterBuilder().setOutputFile(this.outputMergedVCFFile.toFile())
				.setReferenceDictionary(this.inputVCFReader.getFileHeader().getSequenceDictionary()).setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
		this.writer =builder.build();
		
		/////////////////modify the header (line containing the sample names)
		
		VCFHeader header = this.inputVCFReader.getHeader();
		this.writer.writeHeader(header); //write the header section of the vcf file
		
		this.inputVCFReader.close();
	}
	
	
	/**
	 * read the input vcf files one by one and filter out duplicate SVs
	 */
	void run() {
		System.out.println("scan each input vcf file and filter out duplicate SVs...");
		this.records=new HashSet<>();
		
		for(Path inVCF:this.inputSVVcfFiles) {
			this.inputVCFReader = new VCFFileReader(inVCF, false);
			////////////////process each locus
			Stream<VariantContext> stream=this.inputVCFReader.iterator().stream();
			stream.forEach(vc->{
				String chrom=vc.getContig();
				int start = vc.getStart();
				Integer end=vc.getAttribute("END")==null?null:Integer.parseInt((String)vc.getAttribute("END")); 
				String type=(String)vc.getAttribute("SVTYPE"); 
				
				if(end==null) {
					throw new IllegalArgumentException("null SV end found!");
				}
				
				SVRecord record=new SVRecord(chrom, type, start, end);
				
				if(this.records.contains(record)) {
					System.out.println("duplicate SV found:"+record.toString());
				}else {
					this.records.add(record);
					this.writer.add(vc);
				}
			});
			
			this.inputVCFReader.close();
		}
		
		
		this.writer.close();
	}
	
	
	static class SVRecord{
		private final String chrom;
		private final String type;
		private final int start;
		private final int end;
		public SVRecord(String chrom, String type, int start, int end) {
			super();
			this.chrom = chrom;
			this.type = type;
			this.start = start;
			this.end = end;
		}
		
		
		@Override
		public String toString() {
			return "SVRecord [chrom=" + chrom + ", type=" + type + ", start=" + start + ", end=" + end + "]";
		}


		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((chrom == null) ? 0 : chrom.hashCode());
			result = prime * result + end;
			result = prime * result + start;
			result = prime * result + ((type == null) ? 0 : type.hashCode());
			return result;
		}
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (!(obj instanceof SVRecord))
				return false;
			SVRecord other = (SVRecord) obj;
			if (chrom == null) {
				if (other.chrom != null)
					return false;
			} else if (!chrom.equals(other.chrom))
				return false;
			if (end != other.end)
				return false;
			if (start != other.start)
				return false;
			if (type == null) {
				if (other.type != null)
					return false;
			} else if (!type.equals(other.type))
				return false;
			return true;
		}
	}
}
