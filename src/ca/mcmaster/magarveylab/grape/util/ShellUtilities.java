package ca.mcmaster.magarveylab.grape.util;

import java.io.IOException;

public class ShellUtilities {
	public static void runCommand(String command) {
		Runtime rt = Runtime.getRuntime();
		Process pr = null;
		try {
			pr = rt.exec(new String[]{"/bin/bash","-c", command});
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		try {
			pr.waitFor();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}
}
