package population.sv.analysis.adjacentRegion.piAndThetaW;

import genomics.utils.SimpleGenomicRegion;


/**
 * 
 * @author tanxu
 *
 */
public class PiAndThetaWWindow extends SimpleGenomicRegion{
	/**
	 * 
	 */
	private int totalSiteNum;
	private double totalPi;
	private double totalThetaW;
	
	/**
	 * 
	 * @param chrom
	 * @param start
	 * @param end
	 */
	public PiAndThetaWWindow(String chrom, int start, int end) {
		super(chrom, start, end);
		// TODO Auto-generated constructor stub
		

		this.totalPi=0d;
		this.totalThetaW=0d;
		this.totalSiteNum=0;
	}
	

	public void addOneToTotalSite() {
		this.totalSiteNum++;
	}
	
	public void addPiAndThetaW(double pi, double thetaW) {
		this.totalPi+=pi;
		this.totalThetaW+=thetaW;
	}
	
	/**
	 * @return the totalSiteNum
	 */
	public int getTotalSiteNum() {
		return totalSiteNum;
	}
	
	/**
	 * @return the totalPi
	 */
	public final double getTotalPi() {
		return totalPi;
	}

	/**
	 * @return the totalThetaW
	 */
	public final double getTotalThetaW() {
		return totalThetaW;
	}

	
}
