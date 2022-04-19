package population.sv.analysis.assignSiteToClosestSV.postprocess.utils;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import basic.NumericRange;
import population.sv.analysis.assignSiteToClosestSV.NucSiteToClosestSVAssignerAndWriter.NucSiteRecord;
import population.sv.analysis.assignSiteToClosestSV.postprocess.NucSiteAssignedToClosestSVExtractor;
import population.sv.utils.SimpleSVType;

/**
 * collector and writer for all nucleotide sites that meet the conditions specified in this {@link WindowedNucSitesCollectorAndWriter};
 * 		see {@link #test(NucSiteRecord)} method
 * 
 * ===========================
 * in detail, perform the following two tasks
 * 
 * 1. for all nucleotide sites that pass {@link #test(NucSiteRecord)}, write them to the {@link #outgroupSamplesOutFile} with each line contains the data string as {@link NucSiteRecord#toDataLine()}
 * 
 * 2. calculate the summary information for calculated pi and theta w for the qualified {@link NucSiteRecord}s
 * 		1. total number of sites
 * 			{@link #totalSites}
 * 		2. summed pi of these sites
 * 			{@link #sumStatOfPi}
 * 		3. mean pi of these sites
 * 			{@link #sumStatOfPi}
 * 		4. sd of pi of these sites
 * 			{@link #sumStatOfPi}
 * 		5. summed theta w of these sites
 * 			{@link #sumStatOfThetaW}
 * 		6. mean theta w of these sites
 * 			{@link #sumStatOfThetaW}
 * 		7. sd of theta w of these sites
 * 			{@link #sumStatOfThetaW}
 * 
 * @author tanxu
 *
 */
public class WindowedNucSitesCollectorAndWriter extends NucSitesCollectorAndWriterBase{
	
	/**
	 * target SV type 
	 */
	private final SimpleSVType targetSVType;
	
	/**
	 * range for sv size 
	 */
	private final NumericRange svSizeRange;
	
	/**
	 * range for sv derived allele proportion 
	 */
	private final NumericRange svDAPRange;
	
	////////////////////////////////////////key infor of this window
	/**
	 * range for the distance of the site to the closest SV's bounary
	 */
	private final NumericRange distToClosestSVRange;
	
	/**
	 * if true, only nuc sites at the left side of its closest SV will be included in this group;
	 * if false, only nuc sites at the right side of its closest SV will be included in this group
	 */
	private final boolean atLeftSideOfClosestSV;
	
	/**
	 * the absolute value indicate the relative position of {@link #distToClosestSVRange} among all the ranges
	 * 		for example, for dist to closest SV ranges (0,100], (100,200], (200,300], (300,400], (500,600], the indices are 
	 * 				0.5, 1.5, 2.5, 3.5
	 * 
	 * if {@link #atLeftSideOfClosestSV} is true, this value is negative
	 * if {@link #atLeftSideOfClosestSV} is false, this value is positive
	 * ==============================
	 * the window index is designed in the way so that it can be directly used to create plot with R
	 * 
	 */
	private final double windowIndex;
	
	
	
	
	public WindowedNucSitesCollectorAndWriter(
			boolean toWriteNucRecordToFile, Path outDir, 
			SimpleSVType targetSVType, NumericRange svSizeRange, NumericRange svDAPRange, 
			NumericRange distToClosestSVRange, double windowIndex, boolean atLeftSideOfClosestSV) throws IOException {
		super(toWriteNucRecordToFile, outDir);
		this.targetSVType = targetSVType;
		this.svSizeRange = svSizeRange;
		this.svDAPRange = svDAPRange;
		this.distToClosestSVRange = distToClosestSVRange;
		this.windowIndex=windowIndex;
		this.atLeftSideOfClosestSV = atLeftSideOfClosestSV;
		
		///////////////
		this.totalSites=0;
		this.sumStatOfPi=new SummaryStatistics();
		this.sumStatOfThetaW=new SummaryStatistics();
		
		//////////////////////////////
		this.buildOutputNucSiteRecordTSVFile();

	}
	
	@Override
	protected void buildOutputNucSiteRecordTSVFile() throws IOException {
		if(this.toWriteNucRecordToFile) {
			this.nucSiteRecordOutFile = Path.of(this.outDir.toString(), this.buildDescriptiveString().concat(".nuc.records.tsv"));
			this.writer = new BufferedWriter(new FileWriter(this.nucSiteRecordOutFile.toString()));
		}
	}
	
	
	/**
	 * build and return the descriptive string of this {@link WindowedNucSitesCollectorAndWriter}
	 * @return
	 */
	public String buildDescriptiveString() {
		StringBuilder sb=new StringBuilder();
		sb.append(this.targetSVType.toString())
		.append(".size.").append(this.svSizeRange.getMin()).append("-").append(this.svSizeRange.getMax())
		.append(".dap.").append(this.svDAPRange.getMin()).append("-").append(this.svDAPRange.getMax())
		.append(".dist.").append(this.distToClosestSVRange.getMin()).append("-").append(this.distToClosestSVRange.getMax())
		.append(".side.").append(this.atLeftSideOfClosestSV?"left":"right")
		;
		
		return sb.toString();
	}
	
	@Override
	public boolean process(NucSiteRecord record) throws IOException {
		boolean ret = test(record);
		if( ret) {
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
		}
		
		return ret;
	}
	
	/**
	 * 
	 * @param record
	 * @return
	 */
	private boolean test(NucSiteRecord record) {
//		boolean test1=this.atLeftSideOfClosestSV&&record.isAtLeftSide();
//		boolean test2=!this.atLeftSideOfClosestSV&&!record.isAtLeftSide();
//		
//		boolean ret=this.targetSVType.equals(record.getClosestSVType()) &&
//				this.svSizeRange.covers(record.getClosestSVSize()) &&
//				this.svDAPRange.covers(record.getClosestSVDerivedAlleleProportion()) &&
//				this.distToClosestSVRange.covers(record.getDistToClosestSVBoundary()) && 
//				((this.atLeftSideOfClosestSV&&record.isAtLeftSide()) || (!this.atLeftSideOfClosestSV&&!record.isAtLeftSide()))
//				;
		
		return this.targetSVType.equals(record.getClosestSVType()) &&
				this.svSizeRange.covers(record.getClosestSVSize()) &&
				this.svDAPRange.covers(record.getClosestSVDerivedAlleleProportion()) &&
				this.distToClosestSVRange.covers(record.getDistToClosestSVBoundary()) && 
				((this.atLeftSideOfClosestSV&&record.isAtLeftSide()) || (!this.atLeftSideOfClosestSV&&!record.isAtLeftSide()))
				;
	}

	/**
	 * @return the targetSVType
	 */
	public final SimpleSVType getTargetSVType() {
		return targetSVType;
	}

	/**
	 * @return the windowIndex
	 */
	public final double getWindowIndex() {
		return windowIndex;
	}

	/**
	 * @return the svSizeRange
	 */
	public final NumericRange getSvSizeRange() {
		return svSizeRange;
	}
	
	/**
	 * @return the svDAPRange
	 */
	public final NumericRange getSvDAPRange() {
		return svDAPRange;
	}

	/**
	 * @return the distToClosestSVRange
	 */
	public final NumericRange getDistToClosestSVRange() {
		return distToClosestSVRange;
	}
	
}