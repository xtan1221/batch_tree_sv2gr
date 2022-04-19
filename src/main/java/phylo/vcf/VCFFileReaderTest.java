package phylo.vcf;

import java.io.File;

import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;


/**
 * test class for VCFFileReader class in htsjdk lib
 * @author tanxu
 *
 */
public class VCFFileReaderTest {
	private String bcfFilePath="C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\snp filter\\test_data\\sorghum.1.versi.1.timo.1.prop.all.24.bicolor.jointly.called.raw.first1000000lines.bcf";
	private String bcfIndex="C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\snp filter\\test_data\\sorghum.1.versi.1.timo.1.prop.all.24.bicolor.jointly.called.raw.first1000000lines.bcf.csi";
	private File bcfFile=new File(bcfFilePath);
	private File bcfIndexFile=new File(bcfIndex);
	
	static String vcfFilePath="C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\snp filter\\test_data\\sorghum.1.versi.1.timo.1.prop.all.24.bicolor.jointly.called.raw.first1000000lines.vcf";
	static File vcfFile=new File(vcfFilePath);
	
	static String vcfgzFilePath="C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\snp filter\\test_data\\sorghum.1.versi.1.timo.1.prop.all.24.bicolor.jointly.called.raw.first1000000lines.vcf.gz";
	static File vcfgzFile=new File(vcfgzFilePath);
	static String vcfgzIndex="C:\\Users\\tanxu\\Desktop\\reseq-structrual variation\\2021\\sap2test pipeline codes\\snp filter\\test_data\\sorghum.1.versi.1.timo.1.prop.all.24.bicolor.jointly.called.raw.first1000000lines.vcf.gz.csi";
	static File vcfgzIndexFile=new File(vcfgzIndex);
	
	
	static VCFFileReader reader;
	
	
	VCFFileReaderTest(){
//		this.reader = new VCFFileReader(bcfFile,bcfIndexFile); //bcf2codec of htsjdk is not updated to the current vcf specification (bcf2.2?), thus cannot parse the new bcf files coded by the specification
		this.reader = new VCFFileReader(vcfFile, false); //vcf file does not have index!
//		this.reader = new VCFFileReader(vcfgzFile, vcfgzIndexFile); //vcf.gz file is not supported?
	}
	
	
	
	
	void testGetFileHeader() {
		VCFHeader VCFHeader = this.reader.getHeader();
		VCFHeader.getGenotypeSamples().forEach(e->{
			System.out.println(e);
		});
	}
	
	void testQuery() {
		System.out.println(reader.isQueryable()); //apparently vcf file is not queryable
	}
	
	static int counter=0;
	void testIterator() {
		reader.iterator().forEachRemaining(v->{
			if(v.isBiallelic()&&v.isSNP()) {
				System.out.println(v.toString());
				counter++;
				if(counter>10)
					System.exit(0);;
			}
		});
	}
	
	
	public static void main(String[] args) {
		VCFFileReaderTest test = new VCFFileReaderTest();
		System.out.println();
//		test.testGetFileHeader();
		
//		test.testQuery();
		test.testIterator();
	}
}
