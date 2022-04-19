package phylo.batch;

import java.nio.file.Path;
import java.util.Collections;
import java.util.List;

import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProcessExecutor;
import phylo.ref.Region;


public class BcfUtils {
	public static final Log log=Log.getInstance(BcfUtils.class);
	/**
	 * 
	 * 
	 * bcftools view -r Chr01:1-2000,Chr01:10000-20000 -O v .bcf > test1.vcf
	 */
	public static void queryRegionIntoVcfFileFromIndexedBcfFile(String bcfFilePathString, String regionListString, String outputVcfFilePathString) {
		
		String commandLineString="bcftools view -r ".concat(regionListString).concat(" -O v ").concat(bcfFilePathString).concat(" -o ").concat(outputVcfFilePathString);
		log.info("run command to query regional vcf: '".concat(commandLineString).concat("'"));
		
		ProcessExecutor.execute(commandLineString); //this will do exactly what is done below plus blocking
		
//		try {
//			Process process=Runtime.getRuntime().exec(commandLineString);
//			
//			BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
//		    String line = "";
//		    while ((line = reader.readLine()) != null) {
//		        System.out.println(line);
//		    }
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
		
//		int exitValue;
//		try {
//			exitValue = rt.exec(commandLineString).exitValue();
//			return exitValue;
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		
//		return 1;
	}
	
	/**
	 * 
	 * @param bcfFile
	 * @return
	 */
	public static void indexBcfFile(Path bcfFile) {
		String commandLineString="bcftools index ".concat(bcfFile.toString());
		log.info("run command to index bcf file: '".concat(commandLineString).concat("'"));
		
		ProcessExecutor.execute(commandLineString);
//		try {
//			Process process=Runtime.getRuntime().exec(commandLineString);
//			
//			BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
//		    String line = "";
//		    while ((line = reader.readLine()) != null) {
//		        System.out.println(line);
//		    }
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
	}
	
	/**
	 * extract the header section of the bcf file and output to the given vcf file
	 * @param bcfFilePathString
	 * @param headerVcfFilePathString
	 */
	public static void queryHeaderRegionFromBcfFile(String bcfFilePathString, String headerVcfFilePathString) {
		String commandLineString="bcftools view -h ".concat(" -O v ").concat(bcfFilePathString).concat(" -o ").concat(headerVcfFilePathString);
		log.info("run command to query header section from bcf to output vcf: '".concat(commandLineString).concat("'"));
		
		ProcessExecutor.execute(commandLineString); //this will do exactly what is done below plus blocking
	}
	
	/**
	 * build and return Comma-separated list of regions
	 * 		chr|chr:pos|chr:beg-end|chr:beg-[,…]
	 * 
	 * @param regionList
	 * @return
	 */
	public static String buildBcftoolsRegionString(List<Region> regionList) {
		if(regionList==null||regionList.isEmpty())
			throw new IllegalArgumentException("given regionList cannot be null or empty!");
		
		Collections.sort(regionList);
		
		StringBuilder sb = new StringBuilder();
		
		boolean firstAdded = false;
		for(Region r:regionList) {
			if(firstAdded) {
				sb.append(",");
			}else {
				firstAdded=true;
			}
			sb.append(r.toString());
		}
		
		return sb.toString();
	}
	
}
