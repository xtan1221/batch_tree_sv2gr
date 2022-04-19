package population.sv.analysis.assignSiteToClosestSV.postprocess;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Map;

import basic.NumericRange;
import basic.Pair;
import population.sv.analysis.assignSiteToClosestSV.postprocess.utils.BackgroundNucSitesCollectorAndWriter;
import population.sv.analysis.assignSiteToClosestSV.postprocess.utils.WindowedNucSitesCollectorAndWriter;
import population.sv.utils.SimpleSVType;
import statistics.StatTestUtils;

public class GroupSummaryFileWriter {
	/**
	 * output dir for summary file 
	 */
	private final Path outDir;
	
	/**
	 * 
	 */
	private final BackgroundNucSitesCollectorAndWriter backgroundNucSitesCollectorAndWriter;
	/**
	 * 
	 */
	private final WindowedNucSitesCollectorAndWriterBatchBuilder groupProcessorAndWriterBatchBuilder;
	
	
	//////////////////////////
	/**
	 * 
	 */
	private Path outBackgroundNucSitesSummaryFile;
	/**
	 * output summary with each contains the data for a {@link WindowedNucSitesCollectorAndWriter}
	 * 
	 * the columns are described in {@link #getOutFileHeaderLine()}
	 */
	private Path outTargetWindowsSummaryTSVFile;
	
	
	
	/**
	 * map from sv type to 
	 * 		map from sv size range to 
	 * 			map from closest sv derived allele proportion range to
	 * 				map from dist to closest sv range to the pair of 
	 * 					{@link WindowedNucSitesCollectorAndWriter} of nuc sites at the left side of closest sv
	 * 					{@link WindowedNucSitesCollectorAndWriter} of nuc sites at the right side of closest sv
	 */
	private Map<SimpleSVType, Map<NumericRange,Map<NumericRange, Map<NumericRange, Pair<WindowedNucSitesCollectorAndWriter,WindowedNucSitesCollectorAndWriter>>>>> groupMap;

	public GroupSummaryFileWriter(
			Path outDir, 
			BackgroundNucSitesCollectorAndWriter backgroundNucSitesCollectorAndWriter,
			WindowedNucSitesCollectorAndWriterBatchBuilder groupProcessorAndWriterBatchBuilder) {
		super();
		this.outDir = outDir;
		this.backgroundNucSitesCollectorAndWriter = backgroundNucSitesCollectorAndWriter;
		this.groupProcessorAndWriterBatchBuilder = groupProcessorAndWriterBatchBuilder;
		
		
		///////////////////////////
		this.prepare();
		this.writeToBackgroundNucSitesSummaryFile();
		this.writeToWindowSummaryTSVFile();
	}


	void prepare() {
		this.outBackgroundNucSitesSummaryFile=Path.of(this.outDir.toString(),"background.summary.txt");
		if(this.outBackgroundNucSitesSummaryFile.toFile().exists()) {
			System.out.println("outBackgroundNucSitesSummaryFile exists, delete it...");
			this.outBackgroundNucSitesSummaryFile.toFile().delete();
		}
		
		//////////////////
		this.outTargetWindowsSummaryTSVFile=Path.of(this.outDir.toString(),"windows.summary.tsv");
		if(this.outTargetWindowsSummaryTSVFile.toFile().exists()) {
			System.out.println("outSummaryTSVFile exists, delete it...");
			this.outTargetWindowsSummaryTSVFile.toFile().delete();
		}
	}
	
	
	void writeToBackgroundNucSitesSummaryFile() {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(this.outBackgroundNucSitesSummaryFile.toString()));
			///////////////////
			StringBuilder sb=new StringBuilder();
			sb.append("total sites:").append(this.backgroundNucSitesCollectorAndWriter.getTotalSites()).append("\n")
			.append("mean pi:").append(this.backgroundNucSitesCollectorAndWriter.getMeanOfPi()).append("\n")
			.append("sd of pi:").append(this.backgroundNucSitesCollectorAndWriter.getSumStatOfPi().getStandardDeviation()).append("\n")
			.append("mean theta_w:").append(this.backgroundNucSitesCollectorAndWriter.getMeanOfThetaW()).append("\n")
			.append("sd of theta_w:").append(this.backgroundNucSitesCollectorAndWriter.getSumStatOfThetaW().getStandardDeviation())
			;
			writer.append(sb.toString());
			
			///////////////////
			writer.flush();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	void writeToWindowSummaryTSVFile() {
		this.groupMap=this.groupProcessorAndWriterBatchBuilder.getGroupMap();
		
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(this.outTargetWindowsSummaryTSVFile.toString()));
			writer.append(this.getOutFileHeaderLine());
			writer.newLine();

			for(SimpleSVType svType:this.groupMap.keySet()) {
				for(NumericRange svSizeRange:this.groupMap.get(svType).keySet()){
					for(NumericRange closestSVDap:this.groupMap.get(svType).get(svSizeRange).keySet()) {
						for(NumericRange distToClosestSV:this.groupMap.get(svType).get(svSizeRange).get(closestSVDap).keySet()) {
							WindowedNucSitesCollectorAndWriter left = this.groupMap.get(svType).get(svSizeRange).get(closestSVDap).get(distToClosestSV).getFirst();
							WindowedNucSitesCollectorAndWriter right = this.groupMap.get(svType).get(svSizeRange).get(closestSVDap).get(distToClosestSV).getSecond();
							
							writer.append(this.getDataLineString(left));
							writer.newLine();
							writer.append(this.getDataLineString(right));
							writer.newLine();
						}
					}
				}
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
	 * @param group
	 * @return
	 */
	String getDataLineString(WindowedNucSitesCollectorAndWriter group) {
		StringBuilder sb = new StringBuilder();
		
		sb.append(group.getTargetSVType().toString()).append("\t")
		.append(group.getSvSizeRange().getMin()).append("\t")
		.append(group.getSvSizeRange().getMax()).append("\t")
		.append(group.getDistToClosestSVRange().getMin()).append("\t")
		.append(group.getDistToClosestSVRange().getMax()).append("\t")
		.append(group.getSvDAPRange().getMin()).append("\t")
		.append(group.getSvDAPRange().getMax()).append("\t")
		//========================
		.append(group.getTotalSites()).append("\t")
		
		.append(group.getTotalSites()==0?"N/A":group.getMeanOfPi()).append("\t")
		.append(group.getTotalSites()==0?"N/A":group.getVarianceOfPi()).append("\t")
		.append(group.getTotalSites()==0?"N/A":StatTestUtils.calculateTwoSampleTTestStatistics(group.getSumStatOfPi(), this.backgroundNucSitesCollectorAndWriter.getSumStatOfPi())).append("\t")
		
		.append(group.getTotalSites()==0?"N/A":group.getMeanOfThetaW()).append("\t")
		.append(group.getTotalSites()==0?"N/A":group.getVarianceOfThetaW()).append("\t")
		.append(group.getTotalSites()==0?"N/A":StatTestUtils.calculateTwoSampleTTestStatistics(group.getSumStatOfThetaW(), this.backgroundNucSitesCollectorAndWriter.getSumStatOfThetaW())).append("\t")
		//=========================
		.append(group.getWindowIndex()).append("\t")
		.append(group.buildDescriptiveString()).append("\t")  //TODO
		.append(group.isToWriteNucRecordToFile()?group.getNucSiteRecordOutFile().getFileName():"N/A");
		
		return sb.toString();
	}
	
	/**
	 * columns:
	 * 
	 * 1: sv_type	
	 * 2: min_size  //inclusive
	 * 3: max_size  //inclusive
	 * 4: min_dist
	 * 5: max_dist
	 * 6: min_DAP  //derived allele prop
	 * 7: max_DAP  //derived allele prop
	 * ===========================
	 * 8: total_sites  //total nuc sites with calculated pi and theta_w
	 * 9: mean_pi
	 * 10: variance_of_pi
	 * 11: z_statistics_pi		//z-statistics of the two-sample z test for means of pi of the nuc sites in the current window and the background; can be used to calculate the p-value;
	 * 							//for detailed information about two-sample z test for mean, see https://cran.r-project.org/web/packages/distributions3/vignettes/two-sample-z-test.html
	 * 12: mean_theta_w
	 * 13: variance_of_theta_w
	 * 14: z_statistics_theta_w //z-statistics of the two-sample z test for means of pi of the nuc sites in the current window and the background; can be used to calculate the p-value;
	 * 							//for detailed information about two-sample z test for mean, see https://cran.r-project.org/web/packages/distributions3/vignettes/two-sample-z-test.html
	 * ====================
	 * 15: window_index //  for figure plot x axis legend
	 * 16: group_name //for figure plot legend
	 * 17: nuc_site_file_name //not needed
	 * 
	 * 
	 * @return
	 */
	String getOutFileHeaderLine() {
		StringBuilder sb = new StringBuilder();
		
		sb.append("sv_type").append("\t")
		.append("min_size").append("\t")
		.append("max_size").append("\t")
		.append("min_dist").append("\t")
		.append("max_dist").append("\t")
		.append("min_DAP").append("\t")
		.append("max_DAP").append("\t")
		//========================
		.append("total_sites").append("\t")
		
		.append("mean_pi").append("\t")
		.append("variance_of_pi").append("\t")
		.append("z_statistics_pi").append("\t") //z-statistics of the two-sample z test for means of pi of the nuc sites in the current window and the background; can be used to calculate the p-value
		
		.append("mean_theta_w").append("\t")
		.append("variance_of_theta_w").append("\t")
		.append("z_statistics_theta_w").append("\t") //z-statistics of the two-sample z test for means of theta_w of the nuc sites in the current window and the background; can be used to calculate the p-value
		
		//=========================
		.append("window_index").append("\t")
		.append("group_name").append("\t")
		.append("nuc_site_file_name");
		
		return sb.toString();
	}
}
