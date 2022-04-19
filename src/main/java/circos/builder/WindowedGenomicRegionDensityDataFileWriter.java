package circos.builder;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import circos.WindowedNumericValueDataFileWriter;
import circos.WindowedNumericValueRecord;
import genomics.utils.AbstractGenomicRegion;
import genomics.utils.GenomicRegionUtils;

/**
 * 
 * @author tanxu
 *
 */
public class WindowedGenomicRegionDensityDataFileWriter<T extends AbstractGenomicRegion> {
	private final List<T> genomicRegions;
	/**
	 * target chroms
	 */
	private final Map<String, Integer> targetChromNameLengthMap;
	private final int windSize;
	private final String genomicRegionType;
	private final Path outDir;
	
	////////////////////////////////
	private Path outFile;
	private int totalGenomicRegionNum;
	private Map<String, List<T>> chromNameSortedGenomicRegionListMap;
	private Map<String, Map<Integer, Integer>> chromWindowIndexGenomicRegionNumMapMap;
	
	private WindowedNumericValueDataFileWriter windowedNumericValueDataFileWriter;
	
	
	public WindowedGenomicRegionDensityDataFileWriter(
			List<T> genomicRegions,
			Map<String, Integer> targetChromNameLengthMap, 
			int windSize, 
			String genomicRegionType,
			Path outDir
			) {
		super();
		this.genomicRegions = genomicRegions;
		this.targetChromNameLengthMap = targetChromNameLengthMap;
		this.windSize = windSize;
		this.genomicRegionType=genomicRegionType;
		this.outDir = outDir;
		
		/////////////////////////////
		this.prepare();
		this.run();
	}
	
	void prepare() {
		this.outFile=Path.of(this.outDir.toString(), this.genomicRegionType.concat(".wind.size.").concat(Integer.toString(this.windSize)).concat(".density.tsv"));
//		this.outFile=Path.of(this.outDir.toString(), this.genomicRegionType.concat(".density.tsv"));
		
		if(this.outFile.toFile().exists()) {
			System.out.println("outfile exists, delete it...");
			this.outFile.toFile().delete();
		}
		///////////////
		this.chromNameSortedGenomicRegionListMap=new HashMap<>();
		for(T region:this.genomicRegions) {
			if(!this.chromNameSortedGenomicRegionListMap.containsKey(region.getChrom())) {
				this.chromNameSortedGenomicRegionListMap.put(region.getChrom(), new ArrayList<>());
			}
			this.chromNameSortedGenomicRegionListMap.get(region.getChrom()).add(region);
		}
		
		/////////////////////////
		this.totalGenomicRegionNum=this.genomicRegions.size();
		
		////////////
		this.chromWindowIndexGenomicRegionNumMapMap=new HashMap<>();
		for(String chrom:this.targetChromNameLengthMap.keySet()) {
			this.chromWindowIndexGenomicRegionNumMapMap.put(chrom, new HashMap<>());
			
			int windIndex=0;
			int windEnd=this.windSize;
			while(windEnd<this.targetChromNameLengthMap.get(chrom)) {
				this.chromWindowIndexGenomicRegionNumMapMap.get(chrom).put(windIndex, 0);
				
				windIndex++;
				windEnd+=this.windSize;
			}
			this.chromWindowIndexGenomicRegionNumMapMap.get(chrom).put(windIndex, 0);
		}
		
		//////////////
		for(String chrom:this.targetChromNameLengthMap.keySet()) {
			List<T> regions=this.chromNameSortedGenomicRegionListMap.get(chrom);
			Collections.sort(regions, GenomicRegionUtils.sorterByChromAndStartPos());
			
			for(T region:regions) {
				int windIndex=region.getCenter()/this.windSize;
				
				this.chromWindowIndexGenomicRegionNumMapMap.get(chrom).put(windIndex, this.chromWindowIndexGenomicRegionNumMapMap.get(chrom).get(windIndex)+1);
			}
		}
	}
	
	/**
	 * write to file
	 */
	void run() {
		
		//////////////
		try {
			this.windowedNumericValueDataFileWriter=new WindowedNumericValueDataFileWriter(this.outFile);
			
			for(String chrom:this.chromWindowIndexGenomicRegionNumMapMap.keySet()) {
				for(int windIndex:this.chromWindowIndexGenomicRegionNumMapMap.get(chrom).keySet()) {
					int windStart = windIndex*this.windSize+1;
					int windEnd=(windIndex+1)*this.windSize>this.targetChromNameLengthMap.get(chrom)?
							this.targetChromNameLengthMap.get(chrom):(windIndex+1)*this.windSize;
					int regionNum=this.chromWindowIndexGenomicRegionNumMapMap.get(chrom).get(windIndex);
					WindowedNumericValueRecord record=
							new WindowedNumericValueRecord(
									chrom, windStart, windEnd,
									(double)regionNum/this.totalGenomicRegionNum //region density
									);
					
					this.windowedNumericValueDataFileWriter.add(record);
				}
			}
			
			this.windowedNumericValueDataFileWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
}
