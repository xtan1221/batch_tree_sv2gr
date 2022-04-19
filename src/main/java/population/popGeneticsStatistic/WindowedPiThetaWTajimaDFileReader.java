package population.popGeneticsStatistic;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import genomics.utils.SimpleGenomicRegion;

/**
 * 
 * @author tanxu
 *
 */
public class WindowedPiThetaWTajimaDFileReader {
	/**
	 * output file of ANGSD for slide window analysis of tajimaD
	 * each line contains the calculated theta w, pi and tajimaD for each window
	 * 
	 * col1: chrom
	 * col2: window center position
	 * col3: summed thetaW of all sites in the window with calculated thetaW
	 * col4: summed pi of all sites in the window with calculated pi
	 * col5: tajimaD
	 * col6: number of sites with calculated tajimaD (and pi, thetaW)
	 */
	private final Path angsdSlideWindowOutFile;
	private final int windSize;
	
	
	/////////////////////////////
	private List<WindowedPiThetaWTajimaDRecord> windowedPiThetaWTajimaDRecords;
	
	public WindowedPiThetaWTajimaDFileReader(Path angsdSlideWindowOutFile, int windSize) {
		super();
		this.angsdSlideWindowOutFile = angsdSlideWindowOutFile;
		this.windSize = windSize;
		
		////////////////////////
		this.run();
	}



	void run() {
		this.windowedPiThetaWTajimaDRecords=new ArrayList<>();
		try {
			////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(this.angsdSlideWindowOutFile.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty()) { //first line is id	family	order
					continue;
				}
				
				//Chr	WinCenter	tW	tP	Tajima	nSites
				//Chr01	30000	68.528524	56.635117	-0.632988	20000
				//Chr01	50000	111.148012	100.674711	-0.345030	19992
				String[] splits=line.split("\\s+");
				
				if(splits[0].equals("Chr")) {//header line
					continue;
				}
				
				String chrom=splits[0];
				int windCenter=Integer.parseInt(splits[1]);
				int windStart=windCenter-this.windSize/2+1;
				int windEnd=windStart+this.windSize-1;
				double totalThetaW=Double.parseDouble(splits[2]);
				double totalPi=Double.parseDouble(splits[3]);
				double tajimaD=Double.parseDouble(splits[4]);
				int totalSites=Integer.parseInt(splits[5]);
				
				if(totalSites==0) {
					continue;
				}
				
				WindowedPiThetaWTajimaDRecord record = 
						new WindowedPiThetaWTajimaDRecord(
								chrom, windStart, windEnd, 
								windCenter, 
								totalThetaW, totalPi, tajimaD, 
								totalSites);
				
				this.windowedPiThetaWTajimaDRecords.add(record);
			}
			
			lineReader.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	
	/**
	 * @return the windowedPiThetaWTajimaDRecords
	 */
	public final List<WindowedPiThetaWTajimaDRecord> getWindowedPiThetaWTajimaDRecords() {
		return windowedPiThetaWTajimaDRecords;
	}


	
	/**
	 * 
	 * @author tanxu
	 *
	 */
	public static class WindowedPiThetaWTajimaDRecord extends SimpleGenomicRegion{
		private final int centerPos;
		private final double summedThetaW;
		private final double summedPi;
		private final double tajimaD;
		private final int siteNum;
		
		
		public WindowedPiThetaWTajimaDRecord(
				String chrom, int start, int end, 
				int centerPos, 
				double summedThetaW, double summedPi, double tajimaD, 
				int siteNum) {
			super(chrom, start, end);
			this.centerPos = centerPos;
			this.summedThetaW = summedThetaW;
			this.summedPi = summedPi;
			this.tajimaD = tajimaD;
			this.siteNum = siteNum;
		}


		/**
		 * @return the centerPos
		 */
		public final int getCenterPos() {
			return centerPos;
		}


		/**
		 * @return the windSize
		 */
		public final int getWindSize() {
			return this.getLen();
		}

		
		/**
		 * @return the summedThetaW
		 */
		public final double getSummedThetaW() {
			return summedThetaW;
		}


		/**
		 * @return the summedPi
		 */
		public final double getSummedPi() {
			return summedPi;
		}


		/**
		 * @return the tajimaD
		 */
		public final double getTajimaD() {
			return tajimaD;
		}


		/**
		 * @return the siteNum
		 */
		public final int getSiteNum() {
			return siteNum;
		}
		
	}
}
