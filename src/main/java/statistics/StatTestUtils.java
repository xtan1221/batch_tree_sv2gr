package statistics;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;

public class StatTestUtils {
	
	
	/**
	 * calculate and return the t statistics for two sample t-test for equality of means
	 * 
	 * for details, see https://stattrek.com/statistics/dictionary.aspx?definition=two-sample%20t-test
	 * @param sample1
	 * @param sample2
	 * @return
	 */
	public static double calculateTwoSampleTTestStatistics(SummaryStatistics sample1, SummaryStatistics sample2) {
		return (sample1.getMean()-sample2.getMean())/Math.sqrt(sample1.getVariance()/sample1.getN()+sample2.getVariance()/sample2.getN());
	}
	
	public static void main(String[] args) {
		SummaryStatistics sample1 = new SummaryStatistics();
		//3,7,11,0,7,0,4,5,6,2,4,7,2,9
		sample1.addValue(3);
		sample1.addValue(7);
		sample1.addValue(11);
		sample1.addValue(0);
		sample1.addValue(7);
		sample1.addValue(0);
		sample1.addValue(4);
		sample1.addValue(5);
		sample1.addValue(6);
		sample1.addValue(2);
		sample1.addValue(4);
		sample1.addValue(7);
		sample1.addValue(2);
		sample1.addValue(9);
		
		SummaryStatistics sample2 = new SummaryStatistics();
		//5,5,4,5,4,5,7,2,6,2,2,7,2,6,4,2,5,2
		sample2.addValue(5);
		sample2.addValue(5);
		sample2.addValue(4);
		sample2.addValue(5);
		sample2.addValue(4);
		sample2.addValue(5);
		sample2.addValue(7);
		sample2.addValue(2);
		sample2.addValue(6);
		sample2.addValue(2);
		sample2.addValue(2);
		sample2.addValue(7);
		sample2.addValue(2);
		sample2.addValue(6);
		sample2.addValue(4);
		sample2.addValue(2);
		sample2.addValue(5);
		sample2.addValue(2);
		
		System.out.println(calculateTwoSampleTTestStatistics(sample1, sample2));
	}
}
