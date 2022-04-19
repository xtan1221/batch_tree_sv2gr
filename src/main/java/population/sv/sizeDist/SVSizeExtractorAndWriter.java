package population.sv.sizeDist;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.List;

import population.sv.utils.SimpleSVLocus;
import population.unfoldedSFS.SingleLocusUnfoldedSFSBuilder;


/**
 * 
 * @author tanxu
 *
 */
public class SVSizeExtractorAndWriter {
	private final List<SimpleSVLocus> SVLoci;
	private final SingleLocusUnfoldedSFSBuilder<SimpleSVLocus> unfoldedSFSBuilder;
	private final Path outDir;
	
	///////////////////////////
	private Path outFile;
	
	public SVSizeExtractorAndWriter(List<SimpleSVLocus> sVLoci, SingleLocusUnfoldedSFSBuilder<SimpleSVLocus> unfoldedSFSBuilder, Path outDir) {
		super();
		SVLoci = sVLoci;
		this.unfoldedSFSBuilder=unfoldedSFSBuilder;
		this.outDir = outDir;
		
		/////////////////////////////
		this.prepare();
		this.run();
	}

	
	void prepare() {
		this.outFile=Path.of(this.outDir.toString(), "sv.type.size.tsv");
		
		if(this.outFile.toFile().exists()) {
			System.out.println("given out file already exists, delete it...");
			this.outFile.toFile().delete();
		}
	}
	
	
	void run() {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(this.outFile.toFile()));
			writer.append(this.getOutFileHeaderLine());
			writer.newLine();
			
			for(SimpleSVLocus sv:this.SVLoci) {
				StringBuilder sb=new StringBuilder();
				
				sb.append(sv.getChrom()).append("\t")
				.append(sv.getStart()).append("\t")
				.append(sv.getEnd()).append("\t")
				.append(sv.getSize()).append("\t")
				.append(sv.getType().toString()).append("\t")
				.append(!this.unfoldedSFSBuilder.getLocusDerivedAlleleProportionMap().containsKey(sv)?"N/A":this.unfoldedSFSBuilder.getLocusDerivedAlleleProportionMap().get(sv));
				
				writer.append(sb.toString());
				writer.newLine();
			}
			
			writer.flush();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	
	private String getOutFileHeaderLine() {
		//chrom	start	end	size	type
		StringBuilder sb=new StringBuilder();
		
		sb.append("chrom").append("\t")
		.append("start").append("\t")
		.append("end").append("\t")
		.append("size").append("\t")
		.append("type").append("\t")
		.append("derived_allele_prop");
		
		return sb.toString();
	}
}
