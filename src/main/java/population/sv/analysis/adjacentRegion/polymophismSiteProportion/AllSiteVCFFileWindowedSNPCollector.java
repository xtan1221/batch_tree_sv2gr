package population.sv.analysis.adjacentRegion.polymophismSiteProportion;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.stream.Stream;

import genomics.chrom.ChromLenReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import population.vcf.utils.VariantContextFilterFactory;

/**
 * 1.  for each given window, collect the following information from a vcf file containing all sites (snp+indel+non-variant sites)
 * 
 * 		1. number of polymorphism sites
 * 
 * 		2. number of site with non-missing data
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
public class AllSiteVCFFileWindowedSNPCollector {
	/**
	 * maximal allowed full genome windows to estimate full genome mean and sd of the proportion of polymorphism sites;
	 */
	static int MAX_FULL_GENOME_WINDOW_NUM = 100000;
	
	/////////////////////////
	/**
	 * 
	 */
	private final Path allSitesVCFFile;
	
	/**
	 * the list of index of samples in the {@link #snpSitesOnlyVCFFile} to be included in the calculation;
	 * 
	 * first sample has index 1;
	 */
	private final List<Integer> includedSampleIndexList;
	
	/**
	 * a list of windows 
	 */
	private final Map<String,List<Window>> chromTargetWindowsMap;
	
	/**
	 * 
	 */
	private final int windowSize;
	
	/**
	 * 
	 */
	private final ChromLenReader chromLenReader;
	
	
	/**
	 * only site with number of samples with non-missing data no less than this value is counted as non-missing site
	 * all other sites are treated as missing data (thus not included in total sites and total )
	 */
	private final int minSampleNumWithNonMissingDataToInclude;
	
	
	////////////////
	private VCFFileReader reader;
	
	
	/////////////
	private int totalSiteNumWithNonMissingData=0;
	private int totalPolymorphismSiteNum=0;
	
	/**
	 * the set of {@link Window}s built based on window size from full genome;
	 * this is built to facilitate calculate the mean and sd of proportion of sites with SNP
	 */
	private List<Window> fullGenomeWindows;
	
	/**
	 * the map from the chrom name to all windows on the chrom including both these from {@link #chromTargetWindowsMap} and {@link #fullGenomeWindows}
	 */
	private Map<String,List<Window>> chromAllWindowsMap;
	
	public AllSiteVCFFileWindowedSNPCollector(
			Path allSitesVCFFile, List<Integer> includedSampleIndexList,
			Map<String, List<Window>> chromBuiltWindowsMap, 
			int windowSize, ChromLenReader chromLenReader,
			int minSampleNumWithNonMissingDataToInclude
			) {
		super();
		this.allSitesVCFFile = allSitesVCFFile;
		this.includedSampleIndexList = includedSampleIndexList;
		this.chromTargetWindowsMap = chromBuiltWindowsMap;
		this.windowSize=windowSize;
		this.chromLenReader=chromLenReader;
		this.minSampleNumWithNonMissingDataToInclude = minSampleNumWithNonMissingDataToInclude;
		
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
				Window w=new Window(chrom, start, start+this.windowSize);
				
				if(!this.chromAllWindowsMap.get(chrom).contains(w))
					this.chromAllWindowsMap.get(chrom).add(w);
				
				this.fullGenomeWindows.add(w);
			}
		}else {
			System.out.println("use all windows to calculate mean and sd:"+winNum);
			for(String chrom:this.chromLenReader.getChromNameLengthMap().keySet()) {
				for(int i=1;i<this.chromLenReader.getChromNameLengthMap().get(chrom);i+=this.windowSize) {
					Window w=new Window(chrom, i, i+this.windowSize);
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
					(a,b)->{
						return a.getStart()-b.getStart();
						}
					);
		}
	}
	
	////////////////////
	private int count=0;
	private int reportLength=1000000;
	void run() {
		System.out.println("start reading all site vcf file ....");
		this.reader = new VCFFileReader(this.allSitesVCFFile, false);
		Stream<VariantContext> stream=reader.iterator().stream();
		currentCoveringWindows = new ArrayList<>();
		
		stream.forEach(vc->{
			count++;
			if(count % reportLength == 0) {
				System.out.println(count/reportLength+"M sites processed...");
			}
			
			///update windows
			this.updateWindow(vc.getContig(), vc.getStart());
//	
//			if(vc.getStart()==154479) {
//				System.out.println();
//			}
			
			if(VariantContextFilterFactory.nonIndelSite()//not indel
					.and(VariantContextFilterFactory.nonMixedTypeSite())//not mixed type
					.and(VariantContextFilterFactory.notLowQual()) //FILTER column does not contain 'LowQual'
					.test(vc)) {
				//passes all filters
			}else {
				return;//skip
			}
			
			////////////////////check each sample
			int numberOfSampleWithNonMissingData=0;
			Set<String> calledAlleles = new HashSet<>(); //all non-missing allele types in the selected sample
			
			for(int index:this.includedSampleIndexList) {
				Genotype gt = vc.getGenotype(index-1); //the given sample index starts from 1, thus need to adjust to 0-based here
				if(gt.isNoCall()) {
					//missing data
				}else {
					numberOfSampleWithNonMissingData++;
					for(Allele allele:gt.getAlleles()) {
						if(allele.isCalled()) {
							calledAlleles.add(allele.getDisplayString());
						}
					}
				}
			}
			
			if(numberOfSampleWithNonMissingData>=this.minSampleNumWithNonMissingDataToInclude) {
				totalSiteNumWithNonMissingData++;
				
				for(Window w:currentCoveringWindows) {
					totalPolymorphismSiteNum++;
					w.addOneToTotalSite();
					if(calledAlleles.size()>=2) {
						w.addOneToPolymorphSite();
					}
				}
			}
			
		});
		
	}
	
	///////////////////////////
	private List<Window> currentCoveringWindows;
	private String currentChrom=null;
	private int currentWindowIndex=-1; //the index of the last window added to the currentCoveringWindows
	
	///check if any window need to be write to file
	///check if any window need to be added
	
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
		List<Window> windowsToBeRemoved=new ArrayList<>();
		for(Window w:this.currentCoveringWindows) {
			if(w.getEnd()<pos) {
				windowsToBeRemoved.add(w);
//				System.out.println(w.toString());
			}
		}
		
		this.currentCoveringWindows.removeAll(windowsToBeRemoved);
		
		////////////////
		//check if any new windows need to be added
		List<Window> windowsToBeAdded=new ArrayList<>();
		int nextWindowIndex=this.currentWindowIndex+1;
		int newWindowAdded=0;
		if(this.chromAllWindowsMap.get(currentChrom).size()>nextWindowIndex) {
			Window nextWindow=this.chromAllWindowsMap.get(currentChrom).get(nextWindowIndex);
			while(nextWindow.getStart()<=pos) {
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
	public Map<String, List<Window>> getChromBuiltWindowsMap() {
		return chromTargetWindowsMap;
	}
	

	/**
	 * @return the totalSiteNumWithNonMissingData
	 */
	public int getTotalSiteNumWithNonMissingData() {
		return totalSiteNumWithNonMissingData;
	}


	/**
	 * @return the totalPolymorphismSiteNum
	 */
	public int getTotalPolymorphismSiteNum() {
		return totalPolymorphismSiteNum;
	}
	
	

	/**
	 * @return the fullGenomeWindows
	 */
	public List<Window> getFullGenomeWindows() {
		return fullGenomeWindows;
	}

}
