package phylo.tree.distance_ditched;

public class Test {
	public static void main(String[] args) {
		int t = 4;
		
		for(int i=0;i<t-1;i++) {
			for(int j=i+1;j<t;j++) {
				System.out.println(i+ " " + j);
			}
		}
	}
}
