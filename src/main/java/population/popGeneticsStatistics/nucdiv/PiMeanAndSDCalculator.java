package population.popGeneticsStatistics.nucdiv;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Path;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;



public class PiMeanAndSDCalculator {
	private final Path tableFile;
	private final boolean onlyIncludeNon0Values;
	private final int targetColIndex;
	
	private SummaryStatistics sumStat;
	 
	public PiMeanAndSDCalculator(Path tableFile, boolean onlyIncludeNon0Values, int targetColIndex) {
		super();
		this.tableFile = tableFile;
		this.onlyIncludeNon0Values=onlyIncludeNon0Values;
		this.targetColIndex=targetColIndex;
		
		this.run();
	}


	void run() {
		this.sumStat=new SummaryStatistics();
		
		try {
			////////////////////
			BufferedReader lineReader = new BufferedReader(new FileReader(this.tableFile.toFile()));
			String line = null;

			while ((line = lineReader.readLine()) != null) {
				if(line.isEmpty()||line.startsWith("CHROM")) {
					continue;
				}
				//CHROM   POS     PI
				//Chr1    1001    0
				//Chr1    1002    0.2
				String[] splits=line.split("\\s+");
				double value=Double.parseDouble(splits[targetColIndex]);
				
				if(value==0&&this.onlyIncludeNon0Values) {//skip sites with pi=0 (non-variant)
					continue;
				}
				
				this.sumStat.addValue(value);
			}
			
			lineReader.close();
		} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
		}
		
		
		System.out.println(this.sumStat.getMean());
		System.out.println(this.sumStat.getStandardDeviation());
	}
}
