package varbin;

import java.io.InputStream;
import java.io.OutputStream;

//public class Pipe implements java.lang.Runnable {
public class Pipe {
	private InputStream stdout;
	private OutputStream stdin;
	
	public Pipe(InputStream stdout, OutputStream stdin){
		this.stdout = stdout;
		this.stdin = stdin;
	}
	
	//@Override
	public void run() {
		System.out.println("in pipe run");
		try {
			byte[] bArray = new byte[512];
			int readLength = 1;
			// keep reading till no more bites (readLength = -1 means eof)
			while (readLength > -1) {
				readLength = this.stdout.read(bArray, 0, bArray.length);
				if (readLength > -1) {
					this.stdin.write(bArray, 0, readLength);
				}
				System.out.println("in pipe loop");
			}
		} catch (Exception e) {
			System.out.println("broken pipe");
			throw new RuntimeException("Broken pipe", e);
		} finally {
			System.out.println("pipe closing");
			try{
				this.stdout.close();
			} catch (Exception e){
			}
			try{
				this.stdin.close();
				System.out.println("pipe closing");
			} catch (Exception e){
			}
		}
	}
	
	public static InputStream multiPipe(Process... processes) throws java.lang.InterruptedException{
		Process p1;
		Process p2;
		for (int i = 0; i < processes.length - 1; i++){
			p1 = processes[i];
			p2 = processes[i+1];
			//new Thread(new Pipe(p1.getInputStream(), p2.getOutputStream())).start();
		}
		Process pLast = processes[processes.length - 1];
		pLast.waitFor();
		return pLast.getInputStream();
		
	}
}
