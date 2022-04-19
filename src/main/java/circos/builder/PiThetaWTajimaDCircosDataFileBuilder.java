package circos.builder;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Set;

import circos.WindowedNumericValueDataFileWriter;
import circos.WindowedNumericValueRecord;
import population.popGeneticsStatistic.WindowedPiThetaWTajimaDFileReader;
import population.popGeneticsStatistic.WindowedPiThetaWTajimaDFileReader.WindowedPiThetaWTajimaDRecord;


/**
 * 
 * @author tanxu
 *
 */
public class PiThetaWTajimaDCircosDataFileBuilder {
	private final WindowedPiThetaWTajimaDFileReader windowedPiThetaWTajimaDFileReader;
	
	private final int windSize;
	
	private final Path outDir;
	/**
	 * target chroms
	 */
	private final Set<String> targetChroms;
	
	//////////////////////////
	private Path windowedMeanPiDataFile;
	private Path windowedMeanThetaWDataFile;
	private Path windowedMeanTajimaDDataFile;
	
	private WindowedNumericValueDataFileWriter windowedMeanPiDataFileWriter;
	private WindowedNumericValueDataFileWriter windowedMeanThetaWDataFileWriter;
	private WindowedNumericValueDataFileWriter windowedMeanTajimaDDataFileWriter;
	
	public PiThetaWTajimaDCircosDataFileBuilder(
			WindowedPiThetaWTajimaDFileReader windowedPiThetaWTajimaDFileReader, int windSize, Path outDir, Set<String> targetChroms) {
		super();
		this.windowedPiThetaWTajimaDFileReader = windowedPiThetaWTajimaDFileReader;
		this.windSize = windSize;
		this.outDir = outDir;
		this.targetChroms=targetChroms;
		
		//////////////////////////
		this.prepare();
		this.run();
	}

	void prepare() {
		this.windowedMeanPiDataFile=Path.of(this.outDir.toString(), "wind.size.".concat(Integer.toString(this.windSize)).concat(".mean.pi.tsv"));
		this.windowedMeanThetaWDataFile=Path.of(this.outDir.toString(), "wind.size.".concat(Integer.toString(this.windSize)).concat(".mean.thetaW.tsv"));
		this.windowedMeanTajimaDDataFile=Path.of(this.outDir.toString(), "wind.size.".concat(Integer.toString(this.windSize)).concat(".mean.tajimaD.tsv"));
		
		if(this.windowedMeanPiDataFile.toFile().exists()) {
			System.out.println("windowedMeanPiDataFile already exists, delete it...");
			this.windowedMeanPiDataFile.toFile().delete();
		}
		if(this.windowedMeanThetaWDataFile.toFile().exists()) {
			System.out.println("windowedMeanThetaWDataFile already exists, delete it...");
			this.windowedMeanThetaWDataFile.toFile().delete();
		}
		if(this.windowedMeanTajimaDDataFile.toFile().exists()) {
			System.out.println("windowedMeanTajimaDDataFile already exists, delete it...");
			this.windowedMeanTajimaDDataFile.toFile().delete();
		}
		
		
	}
	
	void run() {
		
		try {
			this.windowedMeanPiDataFileWriter = new WindowedNumericValueDataFileWriter(this.windowedMeanPiDataFile);
			this.windowedMeanThetaWDataFileWriter=new WindowedNumericValueDataFileWriter(this.windowedMeanThetaWDataFile);
			this.windowedMeanTajimaDDataFileWriter=new WindowedNumericValueDataFileWriter(this.windowedMeanTajimaDDataFile);
			
			
			for(WindowedPiThetaWTajimaDRecord record:this.windowedPiThetaWTajimaDFileReader.getWindowedPiThetaWTajimaDRecords()) {
				String chrom=record.getChrom();
				if(!this.targetChroms.contains(chrom)) {
					continue;
				}
				
				WindowedNumericValueRecord piRecord=
						new WindowedNumericValueRecord(
								chrom, record.getStart(), record.getEnd(), 
								record.getSummedPi()/record.getSiteNum());
				this.windowedMeanPiDataFileWriter.add(piRecord);
				
				WindowedNumericValueRecord thetaWRecord=
						new WindowedNumericValueRecord(
								chrom, record.getStart(), record.getEnd(), 
								record.getSummedThetaW()/record.getSiteNum());
				this.windowedMeanThetaWDataFileWriter.add(thetaWRecord);
				
				WindowedNumericValueRecord tajimaDRecord=
						new WindowedNumericValueRecord(
								chrom, record.getStart(), record.getEnd(), 
								record.getTajimaD());
				this.windowedMeanTajimaDDataFileWriter.add(tajimaDRecord);
			}
			
			///////////////
			this.windowedMeanPiDataFileWriter.close();
			this.windowedMeanThetaWDataFileWriter.close();
			this.windowedMeanTajimaDDataFileWriter.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
}
