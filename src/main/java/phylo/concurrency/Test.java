package phylo.concurrency;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class Test {
	
	static class SimpleRun implements Callable<Integer> {
		private final int id;
		SimpleRun(int id){
			this.id=id;
		}
		@Override
		public Integer call() {
			try {
				Thread.sleep(3000); //suppose each task needs 3 seconds to finish
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			System.out.println("hello from:"+id);
			return 1;
		}
	}
	
	static void test() throws InterruptedException {
		ExecutorService executorService = Executors.newFixedThreadPool(10000); 
		
		List<Future<Integer>> futureList=new ArrayList<>();
		for(int i=0;i<100;i++) {
			SimpleRun sr=new SimpleRun(i);
			futureList.add(executorService.submit(sr));
		}
		executorService.shutdown();
		
		boolean allDone=false;
		while(!allDone) {
			Thread.sleep(1000);
			System.out.println("still running...");
			allDone=true;
			for(Future<Integer> future:futureList) {
				allDone=allDone&&future.isDone();
			}
		}
		
		System.out.println("all done!");
	}
	
	
	public static void main(String[] args) {
		System.out.println("start");
		long startTime=System.nanoTime();
		
		try {
			test();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		long elapsedNanos = System.nanoTime() - startTime;
		
		System.out.println("elpased time:"+TimeUtils.getReadableTimeFromNanoTime(elapsedNanos));
	}
}
