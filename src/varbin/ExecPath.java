package varbin;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class ExecPath {
	private String nameStr;
	private String pathStr;

	public ExecPath(String nameStr) {
		this.nameStr = nameStr;
		this.pathStr = null;

		// TODO Auto-generated constructor stub
	}

	public String findPath() {
		String sysPath = System.getenv("PATH");

		List<String> dirs = Arrays.asList(sysPath.split(File.pathSeparator));
		System.out.println(sysPath);
		System.out.println(dirs.toString());

		for (String dir : dirs) {
			File file = new File(dir, this.nameStr);
			if (file.isFile()) {
				this.pathStr = file.getAbsolutePath();
				break;
			}
		}
		
		if (this.pathStr == null) {
			System.err.println("ERROR: path for executable " + this.nameStr + " not found");
			System.exit(1);
		}
		return pathStr;
	}
}