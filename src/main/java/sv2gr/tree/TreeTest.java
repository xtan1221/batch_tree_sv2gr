package sv2gr.tree;

import phylo.tree.reader.NewickFileFormatType;

public class TreeTest {
	
	public static void main(String[] args) {
		String newickString="(3237:0.00994728,((((((13123:7.681E-4,(13158:6.3E-6,(Gp0009948:0.0,Gp0018195:0.0)1.0:7.1E-7)1.0:5.939E-5)0.94:4.785E-5,13153:7.7132E-4)0.71:6.97E-5,(((((((11657:1.8126E-4,Gp0018187:1.904E-4)1.0:2.1885E-4,Gp0018197:1.9498E-4)1.0:7.484E-5,(Gp0018183:7.936E-5,Gp0018186:3.75E-5)1.0:3.6364E-4)1.0:3.5426E-4,11757:4.656E-4)1.0:1.8781E-4,Gp0018192:2.8049E-4)0.8:2.312E-5,Gp0018209:6.0678E-4)1.0:4.9755E-4,(13068:2.0196E-4,13148:4.2692E-4)1.0:4.1038E-4)0.62:5.867E-5)0.53:7.493E-5,11647:4.4317E-4)0.79:1.9181E-4,((((11652:8.1422E-4,Gp0018189:6.1966E-4)1.0:1.6319E-4,(13138:6.6993E-4,Gp0009978:5.5325E-4)1.0:9.3673E-4)0.99:1.8875E-4,13128:0.00101802)0.96:2.1528E-4,13133:9.5886E-4)0.77:1.2088E-4)0.76:2.0615E-4,(13143:0.00108893,Gp0018212:3.8241E-4)1.0:8.0569E-4):0.00994728);";
		
		Tree tree=Tree.fromNewickString(newickString, NewickFileFormatType.SIMPLE_NEWICK_2);
		
		System.out.println();
	}
	
	
	
}
