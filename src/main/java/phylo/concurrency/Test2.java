package phylo.concurrency;

public class Test2 {
	
	
	public static void main(String[] args) {
		int cores = Runtime.getRuntime().availableProcessors();
		System.out.println(cores);
		
		int a = 10;
		int b = 15;
		
		System.out.println((double)a/b);
	}
}
