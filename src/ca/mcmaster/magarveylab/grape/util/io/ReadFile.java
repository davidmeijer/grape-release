package ca.mcmaster.magarveylab.grape.util.io;

import java.io.File;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Scanner;

public class ReadFile {

	@SuppressWarnings("resource")
	public static Map<String, String> read(File file) {
		Map<String, String> inputData = new LinkedHashMap<String, String>();
		Scanner nrpScanner = null;
		try {	
			nrpScanner = new Scanner(file);
		} catch (Exception e) {
			e.printStackTrace();
		}
		while(nrpScanner.hasNextLine()) {
			String line = nrpScanner.nextLine();
			
			String[] vals = line.split("\t");
			if(vals.length < 2) {
				continue;
			}
			if(vals[0].startsWith("#")) {
				continue;
			}
			try {
				inputData.put(vals[0], vals[1]);
			} catch (Exception e) {
				System.err.println("Bad SMILES on line " + line);				
			}
		}
		return inputData;
	}
}
