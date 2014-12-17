package varbin;


import java.io.*;

public class StdIOExample2 {
	  public StdIOExample2(String argv) {
		    try {
		      String line;
		      String lineRead;
		      //OutputStream stdin = null;
		      //InputStream stderr = null;
		      //InputStream stdout = null;

		      // launch the command and grab stdin/stdout and stderr
		      
		      
		      ProcessBuilder pb = new ProcessBuilder("sed", "s/hello/goodbye/");
		      Process process = pb.start();
		      //Process process = Runtime.getRuntime().exec(" sed s/hello/goodbye/ ");
		      BufferedReader stdout = new BufferedReader(new  InputStreamReader(process.getInputStream()));
		      BufferedReader stderr = new BufferedReader(new InputStreamReader(process.getErrorStream()));
		      BufferedWriter stdin = new BufferedWriter(new OutputStreamWriter(process.getOutputStream()));

		      line = "bla bla bla hello bla" + "\n";   
		      System.out.println("[Stdin] " + line);
		      //while ((lineRead = stdout.readLine()) != null) {
		    	  System.out.println("[Stdout] " + line);
		      //}
		      stdin.write(line);
		      stdin.flush();


		      line = "yadda yadda hello yadda" + "\n";
		      System.out.println("[Stdin] " + line);
		      //while ((lineRead = stdout.readLine()) != null) {
		    	  System.out.println("[Stdout] " + line);
		      //}
		      stdin.write(line);
		      stdin.flush();


		      stdin.close();
		      
		      System.out.println("Clean up if any output in stdout.");
		      while ((line = stdout.readLine()) != null) {
		        System.out.println ("[Stdout] " + line);
		      }
		      stdout.close();

		      System.out.println("Clean up if any output in stderr.");
		      while ((line = stderr.readLine()) != null) {
		        System.out.println ("[Stderr] " + line);
		      }
		      stderr.close();
		      

		      System.out.println("Exit value: " + process.exitValue());

		      process.destroy();
		      System.out.println("stdout: " + stdout.toString());
		    } catch (Exception err) {
		      err.printStackTrace();
		    }
		  }
}
