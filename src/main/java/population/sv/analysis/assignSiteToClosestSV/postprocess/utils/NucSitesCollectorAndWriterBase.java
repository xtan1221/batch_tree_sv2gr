package population.sv.analysis.assignSiteToClosestSV.postprocess.utils;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Path;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import population.sv.analysis.assignSiteToClosestSV.NucSiteToClosestSVAssignerAndWriter.NucSiteRecord;

/**
 * base class 
 * @author tanxu
 * 
 */
public abstract class NucSitesCollectorAndWriterBase{
	/**
	 * whether or not to write nuc sites to output file
	 * true only for debug and testing since the test statistic is directly calculated from the mean and variance of nuc sites of each group;
	 * 		two sample Z test for mean
	 */
	protected final boolean toWriteNucRecordToFile;
	/**
	 * output dir for the file to write the detailed data for each nuc site of this group
	 */
	protected final Path outDir;
	
	/////////////////////////////////////
	/**
	 * file to write the detailed data for {@link NucSiteRecord}s that pass {@link #test(NucSiteRecord)} of this {@link NucSitesCollectorAndWriterBase}
	 * each line contains the {@link NucSiteRecord#toDataLine()} for a {@link NucSiteRecord}
	 */
	protected Path nucSiteRecordOutFile;
	
	protected BufferedWriter writer;
	
	/////////////////////////
	/**
	 * summary information for all {@link NucSiteRecord}s that pass {@link #test(NucSiteRecord)} 
	 */
	protected int totalSites;
	/**
	 * 
	 */
	protected SummaryStatistics sumStatOfPi;
	/**
	 * 
	 */
	protected SummaryStatistics sumStatOfThetaW;
	
	public NucSitesCollectorAndWriterBase(
			boolean toWriteNucRecordToFile,
			Path outDir
			) throws IOException {
		super();
		this.toWriteNucRecordToFile=toWriteNucRecordToFile;
		this.outDir = outDir;
		
		///////////////
		this.totalSites=0;
		this.sumStatOfPi=new SummaryStatistics();
		this.sumStatOfThetaW=new SummaryStatistics();
		
	}
	
	
	/**
	 * build the {@link #nucSiteRecordOutFile}
	 * @throws IOException
	 */
	protected abstract void buildOutputNucSiteRecordTSVFile() throws IOException;
	
	
	/**
	 * check the given {@link NucSiteRecord}, if it should be included in this {@link NucSitesCollectorAndWriterBase}, 
	 * 1. write it to the {@link #nucSiteRecordOutFile}
	 * 
	 * 2. update the {@link #totalSites}, {@link #sumStatOfPi} and {@link #sumStatOfThetaW}
	 * @param record
	 * @throws IOException
	 */
	public abstract boolean process(NucSiteRecord record) throws IOException;
	
	
	/**
	 * 
	 * @throws IOException
	 */
	public void closeWriter() throws IOException {
		if(this.toWriteNucRecordToFile) {
			this.writer.flush();
			this.writer.close();
		}
	}
	
	
	/**
	 * @return the totalSites
	 */
	public final int getTotalSites() {
		return totalSites;
	}
	
	/**
	 * @return the summedPi
	 */
	public final double getSummedPi() {
		return this.sumStatOfPi.getSum();
	}

	/**
	 * @return the meanOfPi
	 */
	public final double getMeanOfPi() {
		return this.sumStatOfPi.getMean();
	}

	/**
	 * @return the sDOfPi
	 */
	public final double getVarianceOfPi() {
		return this.sumStatOfPi.getVariance();
	}

	/**
	 * @return the summedThetaW
	 */
	public final double getSummedThetaW() {
		return this.sumStatOfThetaW.getSum();
	}

	/**
	 * @return the meanOfThetaW
	 */
	public final double getMeanOfThetaW() {
		return this.sumStatOfThetaW.getMean();
	}
	
	/**
	 * @return the sDOfPi
	 */
	public final double getVarianceOfThetaW() {
		return this.sumStatOfThetaW.getVariance();
	}
	
	/**
	 * @return the nucSiteRecordOutFile
	 */
	public final Path getNucSiteRecordOutFile() {
		return nucSiteRecordOutFile;
	}

	/**
	 * @return the sumStatOfPi
	 */
	public final SummaryStatistics getSumStatOfPi() {
		return sumStatOfPi;
	}

	/**
	 * @return the sumStatOfThetaW
	 */
	public final SummaryStatistics getSumStatOfThetaW() {
		return sumStatOfThetaW;
	}


	/**
	 * @return the toWriteNucRecordToFile
	 */
	public boolean isToWriteNucRecordToFile() {
		return toWriteNucRecordToFile;
	}
	
	
}