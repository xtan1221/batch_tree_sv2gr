package phylo.batch;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import phylo.ref.Region;

public class BcfUtilsTest {
	
	/**
	 * java this bcfFile
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
//		String bcfFilePathString=args[0];
//		String outputVcfFilePathString=args[1];
//		
//		List<Region> regionList = new ArrayList<>();
//		regionList.add(new Region("Chr01", 1, 10000));
//		regionList.add(new Region("Chr01", 111111, 1000110));
//		
//		BcfUtils.queryRegionIntoVcfFile(bcfFilePathString, regionList, outputVcfFilePathString);
		System.out.println("=================");
		String commandLine="bcftools view -r Chr01:1-10000,Chr01:111111-1000110 -O v /scratch/tanxu/reseq/test/batch_local_tree/sorghum.1.versi.1.timo.1.prop.all.24.bicolor.jointly.called.raw.first1000000lines.bcf -o /scratch/tanxu/reseq/test/batch_local_tree/out.vcf";
		
		Process process=Runtime.getRuntime().exec(commandLine);
		
		BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
	    String line = "";
	    while ((line = reader.readLine()) != null) {
	        System.out.println(line);
	    }
		
	    System.out.println("=================");
		commandLine="mkdir test";
		process=Runtime.getRuntime().exec(commandLine);
		
		reader = new BufferedReader(new InputStreamReader(process.getInputStream()));
	    while ((line = reader.readLine()) != null) {
	        System.out.println(line);
	    }
	}
	
	
}
