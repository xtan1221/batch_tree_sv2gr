package population.sv.analysis.assignSiteToClosestSV;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import genomics.utils.GenomicRegionUtils;
import genomics.utils.Position;
import population.sv.utils.SimpleSVLocus;
import population.sv.utils.SimpleSVType;
import population.unfoldedSFS.SingleLocusUnfoldedSFSBuilder;


/**
 * assign each site in a per site data file with calculated Pi and theta w to the closest SVs and output all information to a file
 * 
 * @author tanxu
 *
 */
public class NucSiteToClosestSVAssignerAndWriter2 {
	/**
	 * 
	 */
	private final List<SimpleSVLocus> targetSVs;
	/**
	 * 
	 */
	private final SingleLocusUnfoldedSFSBuilder<SimpleSVLocus> SVUnfoldedSFSBuilder;
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
	 * 
	 */
	private final Path outDir;
	
	//////////////////////////
	/**
	 * the output tsv file contains the extracted closest SV information for each site in the given {@link #perSiteThetaFile} with calculated pi and/or theta w
	 * the columns are as in {@link #getOutFileHeaderLineString()}
	 */
	private Path outTSVFile;
	private Map<String, List<SimpleSVLocus>> chromSVsMap;
	
	
	public NucSiteToClosestSVAssignerAndWriter2(
			List<SimpleSVLocus> targetSVs, 
			SingleLocusUnfoldedSFSBuilder<SimpleSVLocus> SVUnfoldedSFSBuilder,
			Path perSiteThetaFile, Path outDir) {
		super();
		this.targetSVs = targetSVs;
		this.SVUnfoldedSFSBuilder=SVUnfoldedSFSBuilder;
		this.perSiteThetaFile = perSiteThetaFile;
		this.outDir = outDir;
		
		//////////////////////////////
		this.prepare();
		this.run();
	}
	
	
	private void prepare() {
		this.outTSVFile=Path.of(this.outDir.toString(), "nuc.site.assigned.to.closest.SV.pi.thetaW.tsv");

		if(this.outTSVFile.toFile().exists()) {
			System.out.println("given outTSVFile already exists, delete it...");
			this.outTSVFile.toFile().delete();
		}
		
		
		////////////////////
		this.chromSVsMap = new HashMap<>();
		
		this.targetSVs.forEach(sv->{
			if(!this.chromSVsMap.containsKey(sv.getChrom())) {
				this.chromSVsMap.put(sv.getChrom(), new ArrayList<>());
			}
			this.chromSVsMap.get(sv.getChrom()).add(sv);
		});
		
		for(String chrom:this.chromSVsMap.keySet()) {
			Collections.sort(
					this.chromSVsMap.get(chrom),
					GenomicRegionUtils.sorterByChromAndStartPos()
					);
		}
	}
	
	
	/////////////////////////////////
	private int count=0;
	private int reportLength=1000000;
	
	void run() {
		try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(this.outTSVFile.toString()));
			writer.append(this.getOutFileHeaderLineString());
			writer.newLine();
			
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
				
				String chrom=splits[0];
				int pos=Integer.parseInt(splits[1]);
				
				double thetaW=Double.parseDouble(splits[2]);
				double pi=Double.parseDouble(splits[3]);
				
				
				NucSiteRecord record = this.query(chrom, pos, pi, thetaW);
				if(record!=null) {
					writer.append(record.toDataLine());
					writer.newLine();
				}
			}
			
			lineReader.close();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	
	private NucSiteRecord query(String chrom, int pos, double pi, double thetaW) {
		if(!this.chromSVsMap.containsKey(chrom)) {
			return null;
		}
		
		/////////////////
		SimpleSVLocus closestSV=null;
		Integer distToClosestSVBoundary=null;
		Boolean atLeftSideOfClosestSV=null;
		
		
		//termination conditions
		//1. there is no SV left on the chrom (trivial)
		//2. a SV is found to be located at the right side of the given position with dist > distToClosestSVBoundary and atLeftSideOfClosestSV is true
		for(SimpleSVLocus sv:this.chromSVsMap.get(chrom)) {
			if(sv.covers(pos)) {//skip SV that covers the given position
				continue; 
			}
			//terminate if a SV is found to be located at the right side of the given position with dist > distToClosestSVBoundary and atLeftSideOfClosestSV is true
			if(sv.getStart()>pos && distToClosestSVBoundary!=null && atLeftSideOfClosestSV && sv.getStart()-pos>distToClosestSVBoundary) {
				break;
			}
			
			////current SV is either at the left side or right side of the given position
			/////////update if needed
			if(sv.getEnd()<pos) {//SV is at the left side of the target pos, note that all previous checked SVs must be fully at the left side of the target position or covering the target position
				int dist = pos-sv.getEnd();
				if(closestSV==null||dist<distToClosestSVBoundary) {//first one or a closer SV is found
					closestSV=sv;
					distToClosestSVBoundary=dist;
					atLeftSideOfClosestSV=false;
				}
			}else {//SV is at the right side of the target pos
				int dist = sv.getStart()-pos;
				if(closestSV==null||dist<distToClosestSVBoundary) {//first one or a closer SV is found
					closestSV=sv;
					distToClosestSVBoundary=dist;
					atLeftSideOfClosestSV=true;
				}
			}
		}
		
		/////////////////////
		if(closestSV==null) {
			return null;
		}else {
			return new NucSiteRecord(
					chrom, pos, closestSV.getType(), closestSV.getSize(),
					this.SVUnfoldedSFSBuilder.getLocusDerivedAlleleProportionMap().containsKey(closestSV)?this.SVUnfoldedSFSBuilder.getLocusDerivedAlleleProportionMap().get(closestSV):null,
					atLeftSideOfClosestSV,
					atLeftSideOfClosestSV?closestSV.getStart():closestSV.getEnd(), 
					distToClosestSVBoundary, 
					pi, thetaW
					);
		}
	}
	
	/**
	 * 
	 * @return
	 */
	private String getOutFileHeaderLineString() {
		StringBuilder sb=new StringBuilder();
		sb.append("#chrom").append("\t")
		.append("pos").append("\t")
		.append("closestSVType").append("\t")
		.append("closestSVSize").append("\t")
		.append("closestSVDerivedAlleleProp").append("\t")
		.append("closestSVBoundary").append("\t")
		.append("isAtLeftOfClosestSV").append("\t")
		.append("distToClosestSVBound").append("\t")
		.append("pi").append("\t")
		.append("thetaW");
		
		return sb.toString();
	}
	
	
	/////////////////////////////////////
	/**
	 * a nucleotide site and the set of information to its closest SV
	 * @author tanxu
	 *
	 */
	public static class NucSiteRecord extends Position{
		private final SimpleSVType closestSVType;
		private final int closestSVSize;
		private final Double closestSVDerivedAlleleProportion;
		private final boolean atLeftSideOfClosestSV;
		private final int closestSVBoundary;
		private final int distToClosestSVBoundary;
		private final double pi;
		private final double thetaW;
		
		public NucSiteRecord(
				String chrom, int coordinate, SimpleSVType closestSVType, int closestSVSize,
				Double closestSVDerivedAlleleProportion, boolean atLeftSideOfClosestSV, int closestSVBoundary,
				int distToClosestSVBoundary, double pi, double thetaW) {
			super(chrom, coordinate);
			this.closestSVType = closestSVType;
			this.closestSVSize = closestSVSize;
			this.closestSVDerivedAlleleProportion = closestSVDerivedAlleleProportion;
			this.atLeftSideOfClosestSV = atLeftSideOfClosestSV;
			this.closestSVBoundary = closestSVBoundary;
			this.distToClosestSVBoundary = distToClosestSVBoundary;
			this.pi = pi;
			this.thetaW = thetaW;
		}
		
		/**
		 * build and return a {@link NucSiteRecord} from the given data line string from a TSV file with the columns as described in {@link #getOutFileHeaderLineString()}
		 * @param dataLine
		 * @return
		 */
		public static NucSiteRecord fromDataLineString(String dataLine) {
			String[] splits = dataLine.split("\t");
			if(splits.length!=10) {
				System.out.println("skip an abnormal data line:"+dataLine);
				return null;
			}
			String chrom = splits[0];
			int coordinate = Integer.parseInt(splits[1]);
			SimpleSVType closestSVType = SimpleSVType.valueOf(splits[2]);
			int closestSVSize = Integer.parseInt(splits[3]);
			Double closestSVDerivedAlleleProportion = splits[4].equals("N/A")?null:Double.parseDouble(splits[4]);
			int closestSVBoundary = Integer.parseInt(splits[5]);
			boolean atLeftSideOfClosestSV = Boolean.parseBoolean(splits[6]);
			int distToClosestSVBoundary = Integer.parseInt(splits[7]);
			Double pi = splits[8].equals("N/A")?null:Double.parseDouble(splits[8]);
			Double thetaW = splits[9].equals("N/A")?null:Double.parseDouble(splits[9]);
			
			return new NucSiteRecord(
					chrom, coordinate, closestSVType, closestSVSize,
					closestSVDerivedAlleleProportion, atLeftSideOfClosestSV, closestSVBoundary,
					distToClosestSVBoundary, pi, thetaW);
		}
		
		
		/**
		 * build a tab-delimited data line string that is ready to output to a tsv file;
		 * 
		 * must be consistent with columns in {@link #getOutFileHeaderLineString()}
		 * 
		 * @return
		 */
		public String toDataLine() {
			StringBuilder sb=new StringBuilder();
			sb.append(this.getChrom()).append("\t")
			.append(this.getStart()).append("\t")
			.append(this.closestSVType.toString()).append("\t")
			.append(this.closestSVSize).append("\t")
			.append(this.closestSVDerivedAlleleProportion==null?"N/A":this.closestSVDerivedAlleleProportion).append("\t")
			.append(this.closestSVBoundary).append("\t")
			.append(this.atLeftSideOfClosestSV).append("\t")
			.append(this.distToClosestSVBoundary).append("\t")
			.append(this.pi).append("\t")
			.append(this.thetaW);
			
			return sb.toString();
		}
		

		/**
		 * @return the closestSVType
		 */
		public final SimpleSVType getClosestSVType() {
			return closestSVType;
		}

		/**
		 * @return the closestSVSize
		 */
		public final int getClosestSVSize() {
			return closestSVSize;
		}

		/**
		 * @return the closestSVDerivedAlleleProportion
		 */
		public final Double getClosestSVDerivedAlleleProportion() {
			return closestSVDerivedAlleleProportion;
		}

		/**
		 * @return the closetSVBoundary
		 */
		public final int getClosetSVBoundary() {
			return closestSVBoundary;
		}

		/**
		 * @return the atLeftSide
		 */
		public final boolean isAtLeftSide() {
			return atLeftSideOfClosestSV;
		}

		/**
		 * @return the distToClosestSVBoundary
		 */
		public final int getDistToClosestSVBoundary() {
			return distToClosestSVBoundary;
		}

		/**
		 * @return the pi
		 */
		public final double getPi() {
			return pi;
		}

		/**
		 * @return the thetaW
		 */
		public final double getThetaW() {
			return thetaW;
		}
	}
}
