package population.sv.analysis.adjacentRegion.piAndThetaW;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;

import population.sv.utils.SimpleSVLocus;
import population.unfoldedSFS.SingleLocusUnfoldedSFSBuilder;

/**
 * post process of the windows collected by {@link WindowPiAndThetaWCollector} to generate output table files that can be directly used to perform statistical anaylsis with R
 * 
 * 
 * @author tanxu
 *
 */
public class CollectedPiAndThetaWDataWriter {
	private final SVAdjacentRegionWindowBuilder builder;
	private final WindowPiAndThetaWCollector collector;
	private final SingleLocusUnfoldedSFSBuilder<SimpleSVLocus> SVUnfoldedSFSBuilder;
	private final Path outDir;
	
	//////////////////////////////////////
	private Path outFullGenomeWindowInforTableFile;
	private Path outSVLocusWindowMeanPiTableFile;
	private Path outSVLocusWindowMeanThetaWTableFile;
	
	public CollectedPiAndThetaWDataWriter(
			SVAdjacentRegionWindowBuilder builder,
			WindowPiAndThetaWCollector collector, 
			SingleLocusUnfoldedSFSBuilder<SimpleSVLocus> SVUnfoldedSFSBuilder,
			Path outDir) {
		super();
		this.builder = builder;
		this.collector = collector;
		this.SVUnfoldedSFSBuilder=SVUnfoldedSFSBuilder;
		this.outDir = outDir;
		
		//
		this.prepare();
		this.writeToFullGenomeWindowInfoTableFile();
		this.writeToSVLocusAdjacentWindowMeanPiTableFile();
		this.writeToSVLocusAdjacentWindowMeanThetaWTableFile();
	}
	

	void prepare() {
		this.outFullGenomeWindowInforTableFile=Path.of(this.outDir.toString(), 
				"full_genome_windows_".concat(this.buildOutFileInfor()).concat(".txt"));
		this.outSVLocusWindowMeanPiTableFile= Path.of(this.outDir.toString(), 
				"sv_adjacent_windows_mean_pi_".concat(this.buildOutFileInfor()).concat(".txt"));
		this.outSVLocusWindowMeanThetaWTableFile= Path.of(this.outDir.toString(), 
				"sv_adjacent_windows_mean_thetaW_".concat(this.buildOutFileInfor()).concat(".txt"));
		
		if(this.outFullGenomeWindowInforTableFile.toFile().exists()) {
			System.out.println("outFullGenomeWindowInforTableFile already exists, delete it...");
			this.outFullGenomeWindowInforTableFile.toFile().delete();
		}
		if(this.outSVLocusWindowMeanPiTableFile.toFile().exists()) {
			System.out.println("outSVLocusWindowInfoTableFile already exists, delete it...");
			this.outSVLocusWindowMeanPiTableFile.toFile().delete();
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
			
			for( PiAndThetaWWindow w:this.collector.getFullGenomeWindows()) {
				StringBuilder sb=new StringBuilder();
				sb.append(w.getChrom()).append("\t")
				.append(w.getStart()).append("\t")
				.append(w.getEnd()).append("\t")
				.append(w.getTotalSiteNum()).append("\t")
				.append(w.getTotalPi()).append("\t")
				.append(w.getTotalThetaW()).append("\t")
				.append(w.getTotalSiteNum()>0?(double)w.getTotalPi()/w.getTotalSiteNum():"N/A").append("\t")
				.append(w.getTotalSiteNum()>0?(double)w.getTotalThetaW()/w.getTotalSiteNum():"N/A");
				
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
		.append("total_pi").append("\t")
		.append("total_thetaw").append("\t")
		.append("mean_pi").append("\t")
		.append("mean_thetaw");
		
		return sb.toString();
	}
	
	/**
	 * 
	 */
	void writeToSVLocusAdjacentWindowMeanPiTableFile() {
	    try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(this.outSVLocusWindowMeanPiTableFile.toString()));
			writer.append(this.getSVLocusAdjacentWindowTableFileHeaderLine());
			writer.newLine();
			
			for(SimpleSVLocus locus:this.builder.getSvLeftSideWindowsMap().keySet()) {
				StringBuilder sb = new StringBuilder();
				sb.append(locus.getType().toString()).append("\t")
				.append(locus.getChrom()).append("\t")
				.append(locus.getStart()).append("\t")
				.append(locus.getEnd()).append("\t")
				.append(locus.getEnd()-locus.getStart()+1).append("\t")
				.append(this.SVUnfoldedSFSBuilder.getLocusDerivedAlleleProportionMap().containsKey(locus)?
						this.SVUnfoldedSFSBuilder.getLocusDerivedAlleleProportionMap().get(locus):"N/A");
				
				//add left windows
				for(int i=0;i<this.builder.getWindowNum();i++) {
					PiAndThetaWWindow w = this.builder.getSvLeftSideWindowsMap().get(locus).get(i);
					
					if(w==null||w.getTotalSiteNum()==0) {
						sb.append("\t").append("N/A");
					}else {
						sb.append("\t").append((double)w.getTotalPi()/w.getTotalSiteNum());
					}
				}
				
				//add right windows
				for(int i=0;i<this.builder.getWindowNum();i++) {
					PiAndThetaWWindow w = this.builder.getSvRightSideWindowsMap().get(locus).get(i);
					
					if(w==null||w.getTotalSiteNum()==0) {
						sb.append("\t").append("N/A");
					}else {
						sb.append("\t").append((double)w.getTotalPi()/w.getTotalSiteNum());
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
	
	/**
	 * 
	 */
	void writeToSVLocusAdjacentWindowMeanThetaWTableFile() {
	    try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(this.outSVLocusWindowMeanThetaWTableFile.toString()));
			writer.append(this.getSVLocusAdjacentWindowTableFileHeaderLine());
			writer.newLine();
			
			for(SimpleSVLocus locus:this.builder.getSvLeftSideWindowsMap().keySet()) {
				StringBuilder sb = new StringBuilder();
				sb.append(locus.getType().toString()).append("\t")
				.append(locus.getChrom()).append("\t")
				.append(locus.getStart()).append("\t")
				.append(locus.getEnd()).append("\t")
				.append(locus.getEnd()-locus.getStart()+1).append("\t")
				.append(this.SVUnfoldedSFSBuilder.getLocusDerivedAlleleProportionMap().containsKey(locus)?
						this.SVUnfoldedSFSBuilder.getLocusDerivedAlleleProportionMap().get(locus):"N/A");
				
				//add left windows
				for(int i=0;i<this.builder.getWindowNum();i++) {
					PiAndThetaWWindow w = this.builder.getSvLeftSideWindowsMap().get(locus).get(i);
					
					if(w==null||w.getTotalSiteNum()==0) {
						sb.append("\t").append("N/A");
					}else {
						sb.append("\t").append((double)w.getTotalThetaW()/w.getTotalSiteNum());
					}
				}
				
				//add right windows
				for(int i=0;i<this.builder.getWindowNum();i++) {
					PiAndThetaWWindow w = this.builder.getSvRightSideWindowsMap().get(locus).get(i);
					
					if(w==null||w.getTotalSiteNum()==0) {
						sb.append("\t").append("N/A");
					}else {
						sb.append("\t").append((double)w.getTotalThetaW()/w.getTotalSiteNum());
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
		.append("size").append("\t")
		.append("derived_allele_prop");
		
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
