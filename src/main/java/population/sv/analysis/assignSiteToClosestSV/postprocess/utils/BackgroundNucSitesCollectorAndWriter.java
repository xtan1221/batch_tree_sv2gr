package population.sv.analysis.assignSiteToClosestSV.postprocess.utils;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import population.sv.analysis.assignSiteToClosestSV.NucSiteToClosestSVAssignerAndWriter.NucSiteRecord;
import population.sv.analysis.assignSiteToClosestSV.postprocess.NucSiteAssignedToClosestSVExtractor;

/**
 * collector and writer for all nucleotide sites that does not meet the conditions of any {@link WindowedNucSitesCollectorAndWriter} in a run of {@link NucSiteAssignedToClosestSVExtractor}
 * 		see {@link #test(NucSiteRecord)} method
 * 
 * @author tanxu
 *
 */
public class BackgroundNucSitesCollectorAndWriter extends NucSitesCollectorAndWriterBase{
	
	
	/**
	 * 
	 * @param outDir
	 * @throws IOException
	 */
	public BackgroundNucSitesCollectorAndWriter(boolean toWriteNucRecordToFile, Path outDir) throws IOException {
		super(toWriteNucRecordToFile, outDir);
		// TODO Auto-generated constructor stub
		
		//////////////////////////////
		this.buildOutputNucSiteRecordTSVFile();
	}
	
	@Override
	protected void buildOutputNucSiteRecordTSVFile() throws IOException {
		if(this.toWriteNucRecordToFile) {
			this.nucSiteRecordOutFile = Path.of(this.outDir.toString(), "background.nuc.records.tsv");
			this.writer = new BufferedWriter(new FileWriter(this.nucSiteRecordOutFile.toString()));
		}
	}
	
	@Override
	public boolean process(NucSiteRecord record) throws IOException {
		if(this.toWriteNucRecordToFile) {
			//write the nuc site data line to the output file
			this.writer.append(record.toDataLine());
			this.writer.newLine();
		}
		
		//update summary statistics
		NucSiteAssignedToClosestSVExtractor.counter++;
		this.totalSites++;
		this.sumStatOfPi.addValue(record.getPi());
		this.sumStatOfThetaW.addValue(record.getThetaW());
		
		return true;
	}
	
}