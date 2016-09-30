package ca.mcmaster.magarveylab.grape.util.io;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.grape.enums.DomainEnums.PolyKetideDomainEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.SugarModificationsEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.TailoringDomainEnums;
import ca.mcmaster.magarveylab.grape.mol.chem.ChemicalAbstraction;
import ca.mcmaster.magarveylab.grape.mol.chem.Fragment;
import ca.mcmaster.magarveylab.grape.mol.chem.Fragment.FragmentType;
import ca.mcmaster.magarveylab.grape.util.ShellUtilities;

public class ImageOutput {

	/**
	 * Create a graphical image report for NRPs
	 * @param nrp
	 * @param monomerFragments
	 * @param index
	 * @throws IOException
	 * @throws CDKException 
	 */
	public static void createReport(String outputPath, String name, String smiles, String id, ChemicalAbstraction chemicalAbstraction) throws IOException, CDKException {
		outputPath = outputPath + "image/";
		new File(outputPath).mkdirs();
		IAtomContainer currentMolecule = null;
		if(name == null){
			name = id;
		}
		String filename = SmilesIO.getCleanFileName(name);
		if(chemicalAbstraction.getErrorMessage() != null){
			File output = new File(outputPath + filename + ".txt");
			PrintWriter writer = new PrintWriter(output, "UTF-8");
			writer.println(chemicalAbstraction.getErrorMessage());
			writer.close();
			return;
		}
		
		List<Fragment> monomerFragments = chemicalAbstraction.getMonomerFragments();
		
		currentMolecule = SmilesIO.readSmilesTemplates(smiles);
		SmilesIO.drawMolecule(outputPath, currentMolecule, "tempCurrentNrp");
		// Make images
		int nrows, ncols;
		int numImages = monomerFragments.size()+1;
		if(numImages <= 6) {
			nrows = 3;
			ncols = 2;
		}
		else if(numImages <= 12) {
			nrows = 4;
			ncols = 3;
		}
		else if(numImages <= 24) {
			nrows = 6;
			ncols = 4;
		} 
		else if (numImages <= 64) {
			nrows = 8;
			ncols = 8;
		}
		else {
			nrows = 12;
			ncols = 12;
		}
		
		for(int i = 0; i < monomerFragments.size(); i++) {
			SmilesIO.drawMolecule(outputPath, monomerFragments.get(i).getAtomContainer(), "tempTempGeneratedFragment_" + i);
			
			String label = "";
			if(monomerFragments.get(i).getFragmentType() == FragmentType.POLYKETIDE) {
				for(int k = 0; k < monomerFragments.get(i).getPkDomains().size(); k++) {
					List<PolyKetideDomainEnums> p = monomerFragments.get(i).getPkDomains().get(k);
					label += p + ": " + monomerFragments.get(i).getLoadingUnits().get(k).getAbbreviation() + "\n";
				}
				for(TailoringDomainEnums tailoringDomain:monomerFragments.get(i).getTailoringDomains()) {
					label += " " + tailoringDomain.toString();
				}
			}
			else if(monomerFragments.get(i).getFragmentType() == FragmentType.AMINO_ACID) {
				monomerFragments.get(i).getFragmentType().toString();
				if(monomerFragments.get(i).getAminoAcidDomains() != null) {
					label += " " + monomerFragments.get(i).getAminoAcidDomains();
				}
				else {
					label += " " + "Unknown"; 
				}
				for(TailoringDomainEnums tailoringDomain:monomerFragments.get(i).getTailoringDomains()) {
					label += " " + tailoringDomain.toString();
				}
				if(monomerFragments.get(i).getTanimotoScore() > 0) {
					//label += " " + monomerFragments.get(i).getTanimotoScore();
				}
			}
			else if(monomerFragments.get(i).getFragmentType() == FragmentType.SUGAR) {
				if(monomerFragments.get(i).getSugarNames().size() > 0) {
					label += " " + monomerFragments.get(i).getSugarNames().get(0);
				}
				for(SugarModificationsEnums s : monomerFragments.get(i).getSugarModifications()) {
					label += " " + s.getAbbreviation();
				}
			}
			else if(monomerFragments.get(i).getFragmentType() == FragmentType.MULTIPLE_AMINO_ACID_PIECE) {
				
				if(monomerFragments.get(i).getAminoAcidDomains().size() > 0) {
					label += monomerFragments.get(i).getAminoAcidDomains();
				}
			}
			else {
				label += monomerFragments.get(i).getFragmentType().toString();
				for(TailoringDomainEnums tailoringDomain:monomerFragments.get(i).getTailoringDomains()) {
					label += " " + tailoringDomain.toString();
				}
			}
			/*
			if(monomerFragments.get(i).getMonomerType() == FragmentType.SUGAR) {
				label += " " + monomerFragments.get(i).getSugarName();
				label += monomerFragments.get(i).getTanimotoScore();
			}
			*/
			label = label.replaceAll("[\"']", "_");
			
			String command = "convert -background white -fill black -pointsize 72 label:'" + label + "' " + outputPath + "tempLabel_"+i+".png; ";
			
			//TODO: add monomertype label. See:
			//convert Geodiamolide_A.png -gravity South   -background White   -splice 0x90  -pointsize 72          -annotate +0+2 'Faerie Dragon'   out.png
			
			command += "convert " + outputPath + "tempTempGeneratedFragment_" + i + ".png " + outputPath + "tempLabel_" + i + ".png -append " + outputPath + "tempGeneratedFragment_"+ i+".png;";
			
			
			ShellUtilities.runCommand(command);
		}
			
		String command = "";
		command += "montage " + outputPath + "tempCurrentNrp* " + outputPath + "tempGeneratedFragment_* -mode Concatenate -tile " + ncols + "x" + nrows + " " + outputPath + "temp_output_nolabel.png; ";
		command += "convert " + outputPath + "temp_output_nolabel.png -gravity South  -background White   -splice 0x180  -pointsize 72     -annotate +0+2 'Name: " + name + "\nChemical type: " + chemicalAbstraction.getChemicalType() + "'   " + outputPath + "" +filename + ".png; ";
		command += "rm " + outputPath + "temp*;";
		
		ShellUtilities.runCommand(command);
	}

}
