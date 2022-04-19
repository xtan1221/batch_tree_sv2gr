package population.sv.analysis.adjacentRegion.piAndThetaW;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import genomics.chrom.ChromLenReader;
import genomics.utils.GenomicRegionUtils;

/**
 * 
 * 1.  for each given window, collect the following information from a file containing all the sites with calculated nucleotide diversity (pi) and Watterson estimator (theta W) for genetic diversity
 * 		
 * 		1. number of sites with calculated pi and theta W
 * 		
 * 		2. summed value of pi and theta W of all sites in the window
 * 
 * 2. for given window size, calculate the same infor for every window on target genome
 * 		the result is used to estimate the mean and sd of the full genome
 * 
 * ========================
 * output file columns
 * col1 = chrom
 * col2 = start
 * col3 = end
 * col4 = total site num with non-missing data (see {@link #minSampleNumWithNonMissingDataToInclude})
 * col5 = total site num with SNP
 * 
 * 
 * @author tanxu
 * 
 */
public class WindowPiAndThetaWCollector {
	/**
	 * maximal allowed full genome windows to estimate full genome mean and sd of the proportion of polymorphism sites;
	 */
	static int MAX_FULL_GENOME_WINDOW_NUM = 100000;
	
	/////////////////////////
	/**
	 * file containing with calculated pi (column 'tp') and theta W (column 'tw') for each window of size 1 and step size 1
	 * 
	 * specifically, the per site theta file produced by ANGSD
	 * 		see http://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests#Step_1:_Finding_a_.27global_estimate.27_of_the_SFS
	 * 
	 * column 1 is the chrom name
	 * column 2 is the position
	 * column 3 is the watterson W
	 * column 4 is the pi
	 * column 5 is the tajimaD
	 * column 6 is the number of sites with calculated tw and tp
	 * 		since the window size is 1, thus if this column is 0, there is no tw and tp calculated for this site
	 * =================
	 * note that column 3 and 4 may be '-inf', for which the pi/theta W should be 0
	 */
	private final Path perSiteThetaFile;
	
	
	/**
	 * a list of target windows 
	 */
	private final Map<String,List<PiAndThetaWWindow>> chromTargetWindowsMap;
	
	/**
	 * 
	 */
	private final int windowSize;
	
	/**
	 * 
	 */
	private final ChromLenReader chromLenReader;
	
	
	/////////////
	private int totalSitesWithCalculatedPiAndThetaW=0;

	
	private double totalPi=0;
	private double totalThetaW=0;
	
	/**
	 * the set of {@link Window}s built based on window size from full genome;
	 * this is built to facilitate calculate the mean and sd of propertion of sites with SNP
	 */
	private List<PiAndThetaWWindow> fullGenomeWindows;
	
	/**
	 * the map from the chrom name to all windows on the chrom including both these from {@link #chromTargetWindowsMap} and {@link #fullGenomeWindows}
	 */
	private Map<String,List<PiAndThetaWWindow>> chromAllWindowsMap;
	
	public WindowPiAndThetaWCollector(
			Path perSiteThetaFile,
			Map<String, List<PiAndThetaWWindow>> chromBuiltWindowsMap, 
			int windowSize, ChromLenReader chromLenReader			
			) {
		super();
		this.perSiteThetaFile=perSiteThetaFile;
		this.chromTargetWindowsMap = chromBuiltWindowsMap;
		this.windowSize=windowSize;
		this.chromLenReader=chromLenReader;
		
		/////////////////////////////
		this.buildFullGenomeWindows();
		this.prepareWindowsForCalculation();
		this.run();
	}
	
	
	void buildFullGenomeWindows() {
		System.out.println("build full genome windows...");
		
		this.chromAllWindowsMap=new HashMap<>();
		for(String chrom:this.chromTargetWindowsMap.keySet()) {
			this.chromAllWindowsMap.put(chrom, new ArrayList<>());
			this.chromAllWindowsMap.get(chrom).addAll(this.chromTargetWindowsMap.get(chrom));
		}
		
		////////////
		for(String chrom:this.chromLenReader.getChromNameLengthMap().keySet()) {
			if(!this.chromAllWindowsMap.containsKey(chrom)) {
				this.chromAllWindowsMap.put(chrom, new ArrayList<>());
			}
		}
		
		//////////////////////////////
		this.fullGenomeWindows = new ArrayList<>();
		
		int totalChromLen=0;
		for(String chrom:this.chromLenReader.getChromNameLengthMap().keySet()) {
			totalChromLen+=this.chromLenReader.getChromNameLengthMap().get(chrom);
		}
		
		int winNum=totalChromLen/this.windowSize;
		if(winNum>MAX_FULL_GENOME_WINDOW_NUM) {//
			System.out.println("randomly selected window number for mean and sd calculation:"+MAX_FULL_GENOME_WINDOW_NUM);
			Random rand=new Random();
			List<String> chromNames = new ArrayList<>();
			chromNames.addAll(this.chromLenReader.getChromNameLengthMap().keySet());
			
			for(int i=0;i<MAX_FULL_GENOME_WINDOW_NUM;i++) {
				String chrom=chromNames.get(rand.nextInt(chromNames.size()));
				
				int start = rand.nextInt(this.chromLenReader.getChromNameLengthMap().get(chrom));
				PiAndThetaWWindow w=new PiAndThetaWWindow(chrom, start, start+this.windowSize-1);
				
				if(!this.chromAllWindowsMap.get(chrom).contains(w))
					this.chromAllWindowsMap.get(chrom).add(w);
				
				this.fullGenomeWindows.add(w);
			}
		}else {
			System.out.println("use all windows to calculate mean and sd:"+winNum);
			for(String chrom:this.chromLenReader.getChromNameLengthMap().keySet()) {
				for(int i=1;i<this.chromLenReader.getChromNameLengthMap().get(chrom);i+=this.windowSize) {
					PiAndThetaWWindow w=new PiAndThetaWWindow(chrom, i, i+this.windowSize);
					if(!this.chromAllWindowsMap.get(chrom).contains(w))
						this.chromAllWindowsMap.get(chrom).add(w);
					
					this.fullGenomeWindows.add(w);
				}
			}
		}
	}
	
	void prepareWindowsForCalculation() {
		for(String chrom:this.chromAllWindowsMap.keySet()) {
			Collections.sort(
					this.chromAllWindowsMap.get(chrom),
					GenomicRegionUtils.sorterByChromAndStartPos()
					);
		}
	}
	
	
	////////////////////
	private int count=0;
	private int reportLength=1000000;
	void run() {
		System.out.println("start reading per site Pi and ThetaW file ....");
		currentCoveringWindows = new ArrayList<>();
		
		try {
			////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(this.perSiteThetaFile.toFile()));
			String line = null;
			
			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty()) { //first line is id	family	order
					continue;
				}
				
				//Chr     WinCenter       tW      tP      Tajima  nSites
				//Chr01   20      0.004730        0.003807        -0.041221       1
				String[] splits = line.split("\\s+");
				if(splits[0].equals("Chr")) {//header line
					continue;
				}
				///////////////////////
				count++;
				if(count % reportLength == 0) {
					System.out.println(count/reportLength+"M sites processed...");
				}
				//////////////////////
				int nSites=Integer.parseInt(splits[5]);
				if(nSites==0) {//no pi and theta w calculated for the current site, skip
					continue;
				}
				
				//
				totalSitesWithCalculatedPiAndThetaW++;
				
				String chrom=splits[0];
				int pos=Integer.parseInt(splits[1]);
				
//				double thetaW=splits[2].equals("-inf")?0:Double.parseDouble(splits[2]);
//				double pi=splits[3].equals("-inf")?0:Double.parseDouble(splits[3]);
				
				double thetaW=Double.parseDouble(splits[2]);
				double pi=Double.parseDouble(splits[3]);
				
				this.totalPi+=pi;
				this.totalThetaW+=thetaW;
				
				////////////
				this.updateWindow(chrom, pos);
				
				for(PiAndThetaWWindow w:currentCoveringWindows) {
					w.addOneToTotalSite();
					w.addPiAndThetaW(pi, thetaW);
				}
			}
			
			lineReader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	///////////////////////////
	private List<PiAndThetaWWindow> currentCoveringWindows;
	private String currentChrom=null;
	private int currentWindowIndex=-1; //the index of the last window added to the currentCoveringWindows
	
	/**
	 * check if any {@link PiAndThetaWWindow} should be added to and/or removed from the {@link #currentCoveringWindows}
	 * 
	 * if yes, do it
	 * @param chrom
	 * @param pos
	 */
	void updateWindow(String chrom, int pos) {
//		System.out.println(chrom+"\t"+pos);
		
		if(currentChrom==null || !currentChrom.equals(chrom)) {
			//initialize or reset
			currentChrom=chrom;
			currentWindowIndex=-1;
			currentCoveringWindows.clear();
		}
		
		//the chrom is not included for any window calculation, skip it
		if(!this.chromAllWindowsMap.containsKey(chrom)) {
			return;
		}
		
		////check if any current windows need to be removed
		List<PiAndThetaWWindow> windowsToBeRemoved=new ArrayList<>();
		for(PiAndThetaWWindow w:this.currentCoveringWindows) {
			if(!w.covers(pos)) {
				windowsToBeRemoved.add(w);
//				System.out.println(w.toString());
			}
		}
		
		this.currentCoveringWindows.removeAll(windowsToBeRemoved);
		
		////////////////
		//check if any new windows need to be added
		List<PiAndThetaWWindow> windowsToBeAdded=new ArrayList<>();
		int nextWindowIndex=this.currentWindowIndex+1;
		int newWindowAdded=0;
		if(this.chromAllWindowsMap.get(currentChrom).size()>nextWindowIndex) {
			PiAndThetaWWindow nextWindow=this.chromAllWindowsMap.get(currentChrom).get(nextWindowIndex);
			while(nextWindow.covers(pos)) {
				windowsToBeAdded.add(nextWindow);
				nextWindowIndex++;
				newWindowAdded++;
				if(this.chromAllWindowsMap.get(currentChrom).size()>nextWindowIndex) {//there is at least one more window
					nextWindow=this.chromAllWindowsMap.get(currentChrom).get(nextWindowIndex);
				}else {
					break;
				}
			}
		}
		
		this.currentWindowIndex+=newWindowAdded;
		this.currentCoveringWindows.addAll(windowsToBeAdded);
	}
	
	///////////////////
	/**
	 * @return the chromBuiltWindowsMap
	 */
	public Map<String, List<PiAndThetaWWindow>> getChromBuiltWindowsMap() {
		return chromTargetWindowsMap;
	}
	
	
	/**
	 * @return the fullGenomeWindows
	 */
	public List<PiAndThetaWWindow> getFullGenomeWindows() {
		return fullGenomeWindows;
	}
	/**
	 * @return the totalSitesWithCalculatedPiAndThetaW
	 */
	public final int getTotalSitesWithCalculatedPiAndThetaW() {
		return totalSitesWithCalculatedPiAndThetaW;
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
