package ca.mcmaster.magarveylab.grape;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

import ca.mcmaster.magarveylab.grape.enums.AcylAdenylatingSubstrates;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.SmallPKunits;
import ca.mcmaster.magarveylab.grape.util.io.SmilesIO;

public class NamesMain {

	static ArrayList<String> inputMoleculeNames = new ArrayList<String>();
	static ArrayList<String> inputMoleculeSmiles = new ArrayList<String>();
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
//		readInput("allAntibiotics_22jul.smi");
//		StringBuilder sb = new StringBuilder();
//		Set<String> allNames = new HashSet<String>();
//		for(int i = 0; inputMoleculeNames.size() > i; i++){
//			String cleanName = SmilesIO.getCleanFileName(inputMoleculeNames.get(i));
//			
//			String cleanNameFixed = cleanName;
//			int counter = 0;
//			while(allNames.contains(cleanNameFixed)){
//				counter ++;
//				cleanNameFixed = "dup_" + counter + "_" +cleanName;
//			}
//			allNames.add(cleanNameFixed);
//			sb.append(inputMoleculeNames.get(i) + "\t" +
//					inputMoleculeSmiles.get(i) + "\t" +
//					cleanNameFixed + "\n");
//		}
		BufferedWriter bw = new BufferedWriter(new FileWriter("monomers_pk_acyl.tsv"));
		bw.append("pks\n");
		for(SmallPKunits pk : SmallPKunits.values()){
			bw.append(pk.name() + "\t");
			for(int i = 0; i < pk.pkDomains().length; i++){
				if(i > 0){
					bw.append(" -- ");
				}
				bw.append(pk.pkDomains()[i].name() + " (" + pk.pkSubstrates()[i].name() + ")");
			}
			bw.append("\t" + pk.smiles() + "\n");
		}
		bw.append("acyl\n");
		for(AcylAdenylatingSubstrates acyl : AcylAdenylatingSubstrates.values()){
			bw.append(acyl.fullName() + "\t" + acyl.abbreviation() + "\t");
			for(int i = 0; i < acyl.smiles().length; i++){
				if(i > 0){
					bw.append(" -- ");
				}
				bw.append(acyl.smiles()[i]);
			}
			bw.append("\n");
		}
		
		bw.close();
	}
	
	public static void readInput(String filePath){
		// Read NRP file
		Scanner nrpScanner = null;
		try {	
			nrpScanner = new Scanner(new File(filePath));
		} catch (Exception e) {
			e.printStackTrace();
		}
		while(nrpScanner.hasNextLine()) {
			String line = nrpScanner.nextLine();
			
			String[] vals = line.split("\t");
			if(vals.length < 2) {
				System.out.println("More than 2 columns");
				continue;
			}
			inputMoleculeNames.add(vals[0]);
			inputMoleculeSmiles.add(vals[1]);
		}
	}
}
