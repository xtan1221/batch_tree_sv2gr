package population.sv.analysis.adjacentRegion.postprocess;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Map;

import basic.NumericRange;
import population.sv.utils.SimpleSVType;

public class GroupStatSummaryWriter {
	private final Path outTableFile;
	private final int windowNum;
	/**
	 * map from the sv type to the 
	 * 		map from sv size range to the 
	 * 			map from derived allele proportion range to the {@link GroupStatSummary}
	 */
	private Map<SimpleSVType, Map<NumericRange, Map<NumericRange,GroupStatSummary>>> svTypeSizeRangeDAPRangeGroupStatSummaryMapMapMap;
	
	
	
	public GroupStatSummaryWriter(
			Path outTableFile, int windowNum,
			Map<SimpleSVType, Map<NumericRange, Map<NumericRange,GroupStatSummary>>> svTypeSizeRangeDAPRangeGroupStatSummaryMapMapMap) {
		super();
		this.outTableFile = outTableFile;
		this.windowNum = windowNum;
		this.svTypeSizeRangeDAPRangeGroupStatSummaryMapMapMap = svTypeSizeRangeDAPRangeGroupStatSummaryMapMapMap;
		if(this.outTableFile.toFile().exists()) {
			System.out.println("output table file exits, delete it...");
			this.outTableFile.toFile().delete();
		}
		
		this.write();
	}

	void write() {
		 try {
			BufferedWriter writer = new BufferedWriter(new FileWriter(this.outTableFile.toString()));
			writer.append(this.getHeaderLine());
			writer.newLine();
			
			for(SimpleSVType svType:this.svTypeSizeRangeDAPRangeGroupStatSummaryMapMapMap.keySet()) {
				for(NumericRange sizeRange:this.svTypeSizeRangeDAPRangeGroupStatSummaryMapMapMap.get(svType).keySet()) {
					for(NumericRange dapRange:this.svTypeSizeRangeDAPRangeGroupStatSummaryMapMapMap.get(svType).get(sizeRange).keySet()) {
						GroupStatSummary group=this.svTypeSizeRangeDAPRangeGroupStatSummaryMapMapMap.get(svType).get(sizeRange).get(dapRange);
						
						//left windows
						for(int i=0;i<this.windowNum;i++) {
							StringBuilder sb=new StringBuilder();
							sb.append(svType.toString()).append("\t")
							.append(sizeRange.getMin()).append("\t")
							.append(sizeRange.getMax()).append("\t")
							.append(dapRange.getMin()).append("\t")
							.append(dapRange.getMax()).append("\t")
							.append(this.buildGroupName(svType, sizeRange, dapRange)).append("\t")
							.append(this.buildLeftWindowIndex(i)).append("\t")
							.append(group.getLeftWindowMeans().get(i)).append("\t")
							.append(group.getLeftWindowSDs().get(i));
							
							writer.append(sb.toString());
							writer.newLine();
						}
						
						//right windows
						for(int i=0;i<this.windowNum;i++) {
							StringBuilder sb=new StringBuilder();
							sb.append(svType.toString()).append("\t")
							.append(sizeRange.getMin()).append("\t")
							.append(sizeRange.getMax()).append("\t")
							.append(dapRange.getMin()).append("\t")
							.append(dapRange.getMax()).append("\t")
							.append(this.buildGroupName(svType, sizeRange, dapRange)).append("\t")
							.append(this.buildRightWindowIndex(i)).append("\t")
							.append(group.getRightWindowMeans().get(i)).append("\t")
							.append(group.getRightWindowSDs().get(i));
							
							writer.append(sb.toString());
							writer.newLine();
						}
					}

				}
			}
			
			writer.flush();
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	private String buildGroupName(SimpleSVType svType, NumericRange sizeRange, NumericRange dapRange) {
		return svType.toString().concat("_")
				.concat(Double.toString(sizeRange.getMin())).concat("_")
				.concat(Double.toString(sizeRange.getMax())).concat("_")
				.concat(Double.toString(dapRange.getMin())).concat("_")
				.concat(Double.toString(dapRange.getMax()))
				;
	}
	
	/**
	 * given index = 0 indicate the right-most window
	 * @param index
	 * @return
	 */
	private String buildLeftWindowIndex(int index) {
		return Double.toString((double)-index-0.5);
	}
	
	/**
	 * given index = 0 indicate the left-most window
	 * @param index
	 * @return
	 */
	private String buildRightWindowIndex(int index) {
		return Double.toString((double)index+0.5);
	}
	
	private String getHeaderLine() {
		StringBuilder sb=new StringBuilder();
		
		sb.append("sv_type").append("\t")
		.append("min_size").append("\t")
		.append("max_size").append("\t")
		.append("min_derived_allele_prop").append("\t")
		.append("max_derived_allele_prop").append("\t")
		.append("group_name").append("\t")
		.append("window_index").append("\t")
		.append("mean").append("\t")
		.append("sd");
		
		return sb.toString();
	}
	
}
