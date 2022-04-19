package population.sv.analysis.assignSiteToClosestSV.postprocess;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import population.sv.analysis.assignSiteToClosestSV.NucSiteToClosestSVAssignerAndWriter;
import population.sv.analysis.assignSiteToClosestSV.NucSiteToClosestSVAssignerAndWriter.NucSiteRecord;
import population.sv.analysis.assignSiteToClosestSV.postprocess.utils.BackgroundNucSitesCollectorAndWriter;
import population.sv.analysis.assignSiteToClosestSV.postprocess.utils.WindowedNucSitesCollectorAndWriter;


/**
 * classify all nuc sites with calculated pi and theta w and assigned to the closest SVs into groups defined by given {@link WindowedNucSitesCollectorAndWriter}s
 * 
 * @author tanxu
 * 
 */
public class NucSiteAssignedToClosestSVExtractor {
	public static int counter=0;
	/**
	 * file built by {@link NucSiteToClosestSVAssignerAndWriter}
	 * 
	 * contains nucleotide sites with calculated pi and theta w and assigned to the closest SVs
	 */
	private final Path nucWithPiAndThetaWAndAssignedToSVFile;
	
	/**
	 * the nuc sites that are not included in any {@link #targetWindowedNucSitesCollectorAndWriters} and will be used as background to perform statistical testing
	 * 
	 */
	private final BackgroundNucSitesCollectorAndWriter backgroundNucSitesCollectorAndWriter;
	
	/**
	 * contains the built target {@link WindowedNucSitesCollectorAndWriter} for nuc sites in a specific window
	 */
	private final WindowedNucSitesCollectorAndWriterBatchBuilder groupProcessorAndWriterBatchBuilder;
	
	
	/////////////////////////////
	/**
	 * 
	 */
	private List<WindowedNucSitesCollectorAndWriter> targetWindowedNucSitesCollectorAndWriters;
	
	public NucSiteAssignedToClosestSVExtractor(
			Path nucWithPiAndThetaWAndAssignedToSVFile,
			BackgroundNucSitesCollectorAndWriter backgroundNucSitesCollectorAndWriter,
			WindowedNucSitesCollectorAndWriterBatchBuilder groupProcessorAndWriterBatchBuilder) throws IOException {
		super();
		this.nucWithPiAndThetaWAndAssignedToSVFile = nucWithPiAndThetaWAndAssignedToSVFile;
		this.backgroundNucSitesCollectorAndWriter=backgroundNucSitesCollectorAndWriter;
		this.groupProcessorAndWriterBatchBuilder = groupProcessorAndWriterBatchBuilder;
		
		//////////////////////////////
		this.prepare();
		this.run();
		this.closeAllWriters();
	}
	
	void prepare() {
		this.targetWindowedNucSitesCollectorAndWriters=new ArrayList<>();
		
		this.targetWindowedNucSitesCollectorAndWriters.addAll(this.groupProcessorAndWriterBatchBuilder.getAllGroups());
	}
	
	
	//////////////////////////////////////
	private int count=0;
	private int reportLength=1000000;
	void run() {
		try {
			BufferedReader lineReader = new BufferedReader(new FileReader(this.nucWithPiAndThetaWAndAssignedToSVFile.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty()||line.startsWith("#")) { //first line is id	family	order
					continue;
				}
				
				///////////////////////
				count++;
				if(count % reportLength == 0) {
					System.out.println(count/reportLength+"M sites processed...");
				}
				
				NucSiteRecord record=NucSiteRecord.fromDataLineString(line);
				if(record==null) {
					continue;
				}
				
				
				int addedWindNum=0;
				List<WindowedNucSitesCollectorAndWriter> addedGroups=new ArrayList<>();
				for(WindowedNucSitesCollectorAndWriter group:this.targetWindowedNucSitesCollectorAndWriters) {
					if(group.process(record)) {
						addedWindNum++;
						addedGroups.add(group);
					}
				}
				
//				if(addedWindNum>=2) {
//					System.out.println("debug");
//				}
				
				if(addedWindNum==0) {//the NucSiteRecord is not belonging to any WindowedNucSitesCollectorAndWriter, add it to backgroundNucSitesCollectorAndWriter 
					this.backgroundNucSitesCollectorAndWriter.process(record);
				}
				
//				if(counter>1000000) {
//					System.out.println("reached111111");
//				}
			}
			
			lineReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	/**
	 * 
	 * @throws IOException
	 */
	void closeAllWriters() throws IOException {
		this.backgroundNucSitesCollectorAndWriter.closeWriter();
		
		for(WindowedNucSitesCollectorAndWriter w:this.targetWindowedNucSitesCollectorAndWriters) {
			w.closeWriter();
		}
	}
}
