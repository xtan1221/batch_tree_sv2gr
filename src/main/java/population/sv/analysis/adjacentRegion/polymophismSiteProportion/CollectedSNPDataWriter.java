package population.sv.analysis.adjacentRegion.polymophismSiteProportion;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;

import population.sv.utils.SimpleSVLocus;

/**
 * post process of the windows collected by {@link AllSiteVCFFileWindowedSNPCollector} to generate output table files that can be directly used to perform statistical anaylsis with R
 * 
 * 
 * @author tanxu
 *
 */
public class CollectedSNPDataWriter {
	private final SVAdjacentRegionWindowBuilder builder;
	private final AllSiteVCFFileWindowedSNPCollector collector;
	private final Path outDir;
	
	//////////////////////////////////////
	private Path outFullGenomeWindowInforTableFile;
	private Path outSVLocusWindowInfoTableFile;
	
	public CollectedSNPDataWriter(
			SVAdjacentRegionWindowBuilder builder,
			AllSiteVCFFileWindowedSNPCollector collector, 
			Path outDir) {
		super();
		this.builder = builder;
		this.collector = collector;
		this.outDir = outDir;
		
		
		//
		this.prepare();
		this.writeToFullGenomeWindowInfoTableFile();
		this.writeToSVLocusAdjacentWindowTableFile();
	}


	void prepare() {
//		this.outSummaryTableFile = Path.of(this.outDir.toString(), "summary.txt");
		this.outFullGenomeWindowInforTableFile=Path.of(this.outDir.toString(), 
				"full_genome_windows_".concat(this.buildOutFileInfor()).concat(".txt"));
		this.outSVLocusWindowInfoTableFile= Path.of(this.outDir.toString(), 
				"sv_adjacent_windows_".concat(this.buildOutFileInfor()).concat(".txt"));
		
//		if(this.outSummaryTableFile.toFile().exists()) {
//			System.out.println("outSummaryTableFile already exists, delete it...");
//			this.outSummaryTableFile.toFile().delete();
//		}
		if(this.outFullGenomeWindowInforTableFile.toFile().exists()) {
			System.out.println("outFullGenomeWindowInforTableFile already exists, delete it...");
			this.outFullGenomeWindowInforTableFile.toFile().delete();
		}
		if(this.outSVLocusWindowInfoTableFile.toFile().exists()) {
			System.out.println("outSVLocusWindowInfoTableFile already exists, delete it...");
			this.outSVLocusWindowInfoTableFile.toFile().delete();
		}
	}
	
	private String buildOutFileInfor() {
		StringBuilder sb = new StringBuilder();
		
		sb.append("num_").append(this.builder.getWindowNum()).append("_size_").append(this.builder.getWindowSize());
		
		return sb.toString();
	}
	
	
	/**
	 * 
	 */
	void writeToFullGenomeWindowInfoTableFile() {
		
		 try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(this.outFullGenomeWindowInforTableFile.toString()));
			writer.append(this.getFullGenomeWindowInfoTableFileHeaderLine());
			writer.newLine();
			
			for(Window w:this.collector.getFullGenomeWindows()) {
				StringBuilder sb=new StringBuilder();
				sb.append(w.getChrom()).append("\t")
				.append(w.getStart()).append("\t")
				.append(w.getEnd()).append("\t")
				.append(w.getTotalSiteNum()).append("\t")
				.append(w.getTotalPolymophSiteNum()).append("\t")
				.append(w.getTotalSiteNum()>0?(double)w.getTotalPolymophSiteNum()/w.getTotalSiteNum():"N/A");
				
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
	
	private String getFullGenomeWindowInfoTableFileHeaderLine() {
		StringBuilder sb=new StringBuilder();
		sb.append("chrom").append("\t")
		.append("start").append("\t")
		.append("end").append("\t")
		.append("total_sites").append("\t")
		.append("total_snp_sites").append("\t")
		.append("propertion");
		
		return sb.toString();
	}
	
	/**
	 * 
	 */
	void writeToSVLocusAdjacentWindowTableFile() {
	    try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(this.outSVLocusWindowInfoTableFile.toString()));
			writer.append(this.getSVLocusAdjacentWindowTableFileHeaderLine());
			writer.newLine();
			
			for(SimpleSVLocus locus:this.builder.getSvLeftSideWindowsMap().keySet()) {
				StringBuilder sb = new StringBuilder();
				sb.append(locus.getType().toString()).append("\t")
				.append(locus.getChrom()).append("\t")
				.append(locus.getStart()).append("\t")
				.append(locus.getEnd()).append("\t")
				.append(locus.getEnd()-locus.getStart()+1);
				
				//add left windows
				for(int i=0;i<this.builder.getWindowNum();i++) {
					Window w = this.builder.getSvLeftSideWindowsMap().get(locus).get(i);
					
					if(w==null||w.getTotalSiteNum()==0) {
						sb.append("\t").append("N/A");
					}else {
						sb.append("\t").append((double)w.getTotalPolymophSiteNum()/w.getTotalSiteNum());
					}
				}
				
				//add right windows
				for(int i=0;i<this.builder.getWindowNum();i++) {
					Window w = this.builder.getSvRightSideWindowsMap().get(locus).get(i);
					
					if(w==null||w.getTotalSiteNum()==0) {
						sb.append("\t").append("N/A");
					}else {
						sb.append("\t").append((double)w.getTotalPolymophSiteNum()/w.getTotalSiteNum());
					}
				}
				
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
	
	
	
	private String getSVLocusAdjacentWindowTableFileHeaderLine() {
		StringBuilder sb=new StringBuilder();
		sb.append("sv_type").append("\t")
		.append("chrom").append("\t")
		.append("start").append("\t")
		.append("end").append("\t")
		.append("size");
		
		//add left window columns
		for(int i=0;i<this.builder.getWindowNum();i++) {
			sb.append("\t").append("left_").append(i+1);
		}
		
		//add right window columns
		for(int i=0;i<this.builder.getWindowNum();i++) {
			sb.append("\t").append("right_").append(i+1);
		}
		
		return sb.toString();
	}
	
}
