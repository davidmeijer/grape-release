package ca.mcmaster.magarveylab.grape.util.io;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.grape.enums.DomainEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.AminoAcidEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.TailoringDomainEnums;
import ca.mcmaster.magarveylab.grape.mol.chem.matcher.SubstrateMatcher;

public class ReadAminoAcidFile {
	
	/**
	 * Read Amino Acid input. 
	 */
	public static void read(String aminoAcidPath, SubstrateMatcher sm) {
		// Generate array of amino acids and NRP IAtomContainers from the file.
		
		// Read amino acid file
		Scanner aaScanner = null;
		try {	
			aaScanner = new Scanner(new File(aminoAcidPath));
		} catch (Exception e) {
			e.printStackTrace();
		}
		while(aaScanner.hasNextLine()) {
			String line = aaScanner.nextLine();
			String[] vals = line.split("\t");
			if(vals[0].startsWith("#")) {
				continue;
			}
			if(vals.length < 4) {
				System.err.println("Found line with fewer than four values: " + line);
				continue;
			}
			try {
				String name = vals[0];
				IAtomContainer mol = SmilesIO.readSmilesTemplates(vals[2]);
				AminoAcidEnums aminoAcidEnum = DomainEnums.getAminoAcidEnumFromAbbreviation(vals[3]);
				List<TailoringDomainEnums> tailoringDomains = new ArrayList<TailoringDomainEnums>();
				for(int i = 4; i < vals.length; i++) {
					// These are tailoring domains
					TailoringDomainEnums tailoringDomain = DomainEnums.getTailoringDomainFromAbbreviation(vals[i]);
					if(tailoringDomain == null) {
						System.err.println("Warning: In amino acid file, found unknown tailoring domain " + vals[i] + " - skipping");
					}
					else {
						tailoringDomains.add(tailoringDomain);
					}
				}
				sm.addAminoAcid(mol, name, aminoAcidEnum, tailoringDomains);
			} catch (Exception e) {
				System.err.println("Found bad line: " + line);
				e.printStackTrace();
			}	
		}
	}
}
