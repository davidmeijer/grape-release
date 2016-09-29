package ca.mcmaster.magarveylab.grape.util.io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


import ca.mcmaster.magarveylab.grape.enums.AcylAdenylatingSubstrates;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.AminoAcidEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.SugarModificationsEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.TailoringDomainEnums;
import ca.mcmaster.magarveylab.grape.nrp.chem.ChemicalAbstraction;
import ca.mcmaster.magarveylab.grape.nrp.chem.Fragment;
import ca.mcmaster.magarveylab.grape.nrp.chem.Fragment.FragmentType;
import ca.mcmaster.magarveylab.grape.util.ChemicalUtilities;

public class TextOutput {

	/**
	 * Given a String of text and a String representing a file path it will write the text to the file
	 * @param text the text to be written to the file
	 * @return filename the name of the file to be written to
	 * @throws IOException
	 */
	public static void writeFile(String text, String filename) throws IOException {
		// make output folder if it does not already exist
		new File("text/").mkdirs();
		
		File file = new File("text/" + filename);
		if(!file.exists()) {
			try {
				file.createNewFile();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		BufferedWriter bw = new BufferedWriter(new FileWriter(file.getAbsoluteFile()));
		bw.write(text);
		bw.close();
	}

	/**
	 * Create an output text string
	 * @param monomerFragments
	 * @throws IOException 
	 */
	public static String createOutputFile(ChemicalAbstraction chemicalAbstraction) {
		if(chemicalAbstraction.getErrorMessage() != null){
			return chemicalAbstraction.getErrorMessage();
		}
		List<Fragment> monomerFragments = chemicalAbstraction.getMonomerFragments();
		List<ArrayList<Fragment>> orderedMonomerFragmentSequences = ChemicalUtilities.getOrderedMonomerFragmentSequences(monomerFragments);
				
		// Useful debugging code for displaying Fragment connectivity
		/*
		for(int i = 0; i < chemicalAbstraction.getMonomerFragments().size(); i++) {
			Fragment frag = chemicalAbstraction.getMonomerFragments().get(i);
			System.out.print(i + " ");
			if(frag.getFragmentType() == FragmentType.AMINO_ACID) {
				System.out.print(frag.getAminoAcidDomains());
			} else {
				System.out.print(frag.getPkDomains());
			}
			int fragIndexAfterC = -1;
			int fragIndexAfterN = -1;
			for(int j = 0; j < chemicalAbstraction.getMonomerFragments().size(); j++) {
				Fragment frag2 = chemicalAbstraction.getMonomerFragments().get(j);
				if(frag2.getMolecule().contains(frag.getAtomAfterCTerminus())) {
					fragIndexAfterC = j;
				}
				if(frag2.getMolecule().contains(frag.getAtomAfterNTerminus())) {
					fragIndexAfterN = j;
				}
			}
			if(fragIndexAfterC > -1)
				System.out.print("After C: " + fragIndexAfterC + chemicalAbstraction.getMonomerFragments().get(fragIndexAfterC).getAminoAcidDomains());
			if(fragIndexAfterN > -1)
				System.out.print("After N: " + fragIndexAfterN + chemicalAbstraction.getMonomerFragments().get(fragIndexAfterN).getAminoAcidDomains());
			System.out.println();
		}
		*/
		
		StringBuilder outputText = new StringBuilder();
		
		outputText.append("Chemical_type\t" + chemicalAbstraction.getChemicalType() + "\n" + "Chemical_subtype\t" + chemicalAbstraction.getChemicalSubType() + "\n");
		
		
		if(chemicalAbstraction.getChemicalScaffoldName() != null){
			outputText.append("Chemical_scaffold\t" + chemicalAbstraction.getChemicalScaffoldName() + "\n");
		}
		
		for(List<Fragment> orderedSequence : orderedMonomerFragmentSequences) {
			outputText.append("Ordered sequence: \n");
			for(Fragment m : orderedSequence) {
				if(m.getFragmentType() == FragmentType.AMINO_ACID) {
					outputText.append("\t" + m.getFragmentType() + "\t");
					if(m.getAminoAcidDomains().size() < 1) {
						outputText.append("Unk\t");
					}
					else {
						outputText.append(m.getAminoAcidDomains().get(0).getAbbreviation() + "\t");
					}
					if(!m.hasKnownStart()){
						outputText.append("cyclic\t");
					}
					for(TailoringDomainEnums t : m.getTailoringDomains()) {
						outputText.append(t.getAbbreviation() + "\t");
					}
				}
				else if(m.getFragmentType() == FragmentType.SUGAR) {
					outputText.append("\t" + m.getFragmentType() + "\t");
					if(m.getSugarNames().size() > 0) {
						for(String sugarName : m.getSugarNames()){
							outputText.append(sugarName + "\t");
						}
					}
					for(SugarModificationsEnums s : m.getSugarModifications()) {
						outputText.append(s.getAbbreviation() + "\t");
					}
				}
				else if(m.getFragmentType() == FragmentType.MULTIPLE_AMINO_ACID_PIECE) {
					if(m.getAminoAcidDomains().size() < 1){
						outputText.append("\t" + m.getFragmentType() + "\t" + m.getNumNitrogens());
					}else{
						for(int i = 0 ; m.getAminoAcidDomains().size() > i ; i++){
							AminoAcidEnums aa = m.getAminoAcidDomains().get(i);
							outputText.append("\t" + FragmentType.AMINO_ACID + "\t" + aa.getAbbreviation() + "\t");
							if(m.getAminoAcidDomains().size() - 1 > i) outputText.append("\n");
						}
					}
					if(!m.hasKnownStart()){
						outputText.append("cyclic\t");
					}
					for(AcylAdenylatingSubstrates starter : m.getStarters()){
						outputText.append(starter.fullName() + "\t");
					}
					for(TailoringDomainEnums t : m.getTailoringDomains()) {
						outputText.append(t.getAbbreviation() + "\t");
					}
				}
				else if(m.getFragmentType() == FragmentType.POLYKETIDE) {
					for(int i = 0; i < m.getPkDomains().size(); i++) {
						outputText.append("\t" + m.getFragmentType() + "\t" + m.getPkDomains().get(i) + "-");
						outputText.append(m.getLoadingUnits().get(i).getAbbreviation() + "\t");
						if(i != m.getPkDomains().size() - 1) {
							outputText.append("\n");
						}
					}
					for(TailoringDomainEnums t : m.getTailoringDomains()) {
						outputText.append(t.getAbbreviation() + "\t");
					}
					for(AcylAdenylatingSubstrates starter : m.getStarters()){
						outputText.append(starter.fullName() + "\t");
					}
				}
				else if(m.getFragmentType() == FragmentType.FA_OR_PK) {
					for(int i = 0; i < m.getPkDomains().size(); i++) {
						outputText.append("\t" + m.getFragmentType() + "\t" + m.getPkDomains().get(i) + "-");
						outputText.append(m.getLoadingUnits().get(i).getAbbreviation() + "\t");
						if(i != m.getPkDomains().size() - 1) {
							outputText.append("\n");
						}
					}
					for(TailoringDomainEnums t : m.getTailoringDomains()) {
						outputText.append(t.getAbbreviation() + "\t");
					}
					for(AcylAdenylatingSubstrates starter : m.getStarters()){
						outputText.append(starter.fullName() + "\t");
					}
				}
				else {
					outputText.append("\t" + m.getFragmentType() + "\t");
					for(TailoringDomainEnums t : m.getTailoringDomains()) {
						outputText.append(t.getAbbreviation() + "\t");
					}
					for(AcylAdenylatingSubstrates starter : m.getStarters()){
						outputText.append(starter.fullName() + "\t");
					}
				}
	
				outputText.append("\n");
			}
		}
		return outputText.toString();
		
	}

}
