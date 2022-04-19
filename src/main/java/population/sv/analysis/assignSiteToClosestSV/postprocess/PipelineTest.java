package population.sv.analysis.assignSiteToClosestSV.postprocess;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import basic.NumericRange;
import population.sv.analysis.assignSiteToClosestSV.postprocess.utils.BackgroundNucSitesCollectorAndWriter;
import population.sv.utils.SimpleSVType;


/**
 * testing pipeline 
 * @author tanxu
 *
 */
public class PipelineTest {
	static boolean toWriteNucRecordToFile = false;
	static Path outDir=Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\sorghum\\sv_nuc_site_assign_to_closest_sv\\test\\output");
	//////////////////////
	static List<SimpleSVType> svTypes;
	static List<NumericRange> svSizeRanges;
	static List<NumericRange> closestSVDerivedAllelePropRanges;
	static List<NumericRange> orderedDistToClosestSVRanges;
	static WindowedNucSitesCollectorAndWriterBatchBuilder groupProcessorAndWriterBatchBuilder;
	
	////////////////
	static Path nucWithPiAndThetaWAndAssignedToSVFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\sorghum\\sv_nuc_site_assign_to_closest_sv\\test\\nuc.site.assigned.to.closest.SV.pi.thetaW.first1M.lines.tsv");
//	static Path nucWithPiAndThetaWAndAssignedToSVFile = Path.of("C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\sv\\sorghum\\sv_nuc_site_assign_to_closest_sv\\test\\test.tsv");

	static BackgroundNucSitesCollectorAndWriter backgroundNucSitesCollectorAndWriter;
	static NucSiteAssignedToClosestSVExtractor nucSiteAssignedToClosestSVExtractor;
	
	////////////////
	static GroupSummaryFileWriter groupSummaryFileWriter;
	
	
	/////////////////////////////////
	/**
	 * 
	 */
	static void prepareGroupProcessorAndWriterBatchBuilder() {
		svTypes=new ArrayList<>();
//		svTypes.add(SimpleSVType.DEL);
		svTypes.addAll(SVAdjacentNucSiteNumericRangeUtils.getAllSVTypes());
		
		/////////
		svSizeRanges=new ArrayList<>();
		svSizeRanges.add(SVAdjacentNucSiteNumericRangeUtils.getSVFullSizeRanges());
//		svSizeRanges.addAll(SVAdjacentNucSiteNumericRangeUtils.getSVSizeRanges1To10000WithSize1000());
		svSizeRanges.addAll(SVAdjacentNucSiteNumericRangeUtils.getSVSizeRanges1To1000WithSize100());
		
		//////////
		closestSVDerivedAllelePropRanges=new ArrayList<>();
		closestSVDerivedAllelePropRanges.add(SVAdjacentNucSiteNumericRangeUtils.getSVFullDerivedAllelePropRange());
		closestSVDerivedAllelePropRanges.addAll(SVAdjacentNucSiteNumericRangeUtils.getSVDerivedAllelePropRanges());
		
		////////
		orderedDistToClosestSVRanges=new ArrayList<>();
		orderedDistToClosestSVRanges.addAll(SVAdjacentNucSiteNumericRangeUtils.getOrderedDistToClosestSVRanges(100,100));
		
	}
	
	
	
	static void buildGroupProcessorAndWriterBatchBuilder() throws IOException {
		groupProcessorAndWriterBatchBuilder = 
				new WindowedNucSitesCollectorAndWriterBatchBuilder(
						toWriteNucRecordToFile,
						outDir,//Path outDir, 
						svTypes,//List<SimpleSVType> svTypes, 
						svSizeRanges,//List<NumericRange> svSizeRanges,
						closestSVDerivedAllelePropRanges,//List<NumericRange> closestSVDerivedAllelePropRanges, 
						orderedDistToClosestSVRanges//List<NumericRange> orderedDistToClosestSVRanges
						);
	}
	
	
	///////////////
	static void buildBackgroundNucSitesCollectorAndWriter() throws IOException {
		backgroundNucSitesCollectorAndWriter=
				new BackgroundNucSitesCollectorAndWriter(toWriteNucRecordToFile, outDir);
	}
	
	static void buildNucSiteAssignedToClosestSVExtractor() throws IOException {
		nucSiteAssignedToClosestSVExtractor = 
				new NucSiteAssignedToClosestSVExtractor(
						nucWithPiAndThetaWAndAssignedToSVFile,//Path nucWithPiAndThetaWAndAssignedToSVFile,
						backgroundNucSitesCollectorAndWriter,//BackgroundNucSitesCollectorAndWriter backgroundNucSitesCollectorAndWriter,
						groupProcessorAndWriterBatchBuilder//GroupProcessorAndWriterBatchBuilder groupProcessorAndWriterBatchBuilder
						);
	}
	
	static void buildGroupSummaryFileWriter() {
		groupSummaryFileWriter = 
				new GroupSummaryFileWriter(
						outDir,//Path outDir, 
						backgroundNucSitesCollectorAndWriter, //BackgroundNucSitesCollectorAndWriter backgroundNucSitesCollectorAndWriter,
						groupProcessorAndWriterBatchBuilder//GroupProcessorAndWriterBatchBuilder groupProcessorAndWriterBatchBuilder
						);
	}
	
	static void run() throws IOException {
		prepareGroupProcessorAndWriterBatchBuilder();
		buildGroupProcessorAndWriterBatchBuilder();
		buildBackgroundNucSitesCollectorAndWriter();
		buildNucSiteAssignedToClosestSVExtractor();
		buildGroupSummaryFileWriter();
	}
	
	
	public static void main(String[] args) {
		try {
			run();
			System.out.println(NucSiteAssignedToClosestSVExtractor.counter);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}

