package population.sv.analysis.adjacentRegion.postprocess;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

import basic.NumericRange;
import population.sv.analysis.adjacentRegion.postprocess.CollectedWindowedDataFileReader.SVRecord;
import population.sv.utils.SimpleSVType;

public class GroupStatSummary {
	private final SimpleSVType svtype;
	private final NumericRange sizeRange;
	private final NumericRange derivedAllelePropRange;
	
	private final List<SVRecord> svRecords;
	private final int windowNum;
	
	
	///////////////////////////
	private List<Double> leftWindowMeans;
	private List<Double> rightWindowMeans;
	private List<Double> leftWindowSDs;
	private List<Double> rightWindowSDs;
	
	
	public GroupStatSummary(
			SimpleSVType svtype,
			NumericRange sizeRange,
			NumericRange derivedAllelePropRange,
			List<SVRecord> svRecords, int windowNum) {
		super();
		this.svtype=svtype;
		this.sizeRange=sizeRange;
		this.derivedAllelePropRange=derivedAllelePropRange;
		this.svRecords = svRecords;
		this.windowNum = windowNum;
		
		this.calculate();
	}
	
	
	void calculate() {
		this.leftWindowMeans=new ArrayList<>();
		this.rightWindowMeans=new ArrayList<>();
		this.leftWindowSDs=new ArrayList<>();
		this.rightWindowSDs=new ArrayList<>();
		
		
		//left windows
		for(int i=0;i<this.windowNum;i++) {
			SummaryStatistics stat = new SummaryStatistics();
			
			for(SVRecord sv:this.svRecords) {
				if(sv.getLeftSideWindowProportion().get(i)!=null)
					stat.addValue(sv.getLeftSideWindowProportion().get(i));
			}
			
			this.leftWindowMeans.add(stat.getMean());
			this.leftWindowSDs.add(stat.getStandardDeviation());
		}
		
		//right windows
		for(int i=0;i<this.windowNum;i++) {
			SummaryStatistics stat = new SummaryStatistics();
			
			for(SVRecord sv:this.svRecords) {
				if(sv.getRightSideWindowProportion().get(i)!=null)
					stat.addValue(sv.getRightSideWindowProportion().get(i));
			}
			this.rightWindowMeans.add(stat.getMean());
			this.rightWindowSDs.add(stat.getStandardDeviation());
		}
	}


	/**
	 * @return the leftWindowMeans
	 */
	public List<Double> getLeftWindowMeans() {
		return leftWindowMeans;
	}


	/**
	 * @return the rightWindowMeans
	 */
	public List<Double> getRightWindowMeans() {
		return rightWindowMeans;
	}


	/**
	 * @return the leftWindowSDs
	 */
	public List<Double> getLeftWindowSDs() {
		return leftWindowSDs;
	}


	/**
	 * @return the rightWindowSDs
	 */
	public List<Double> getRightWindowSDs() {
		return rightWindowSDs;
	}
	
	
}
