package ca.mcmaster.magarveylab.grape.util.io;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.interfaces.IMolecule;

import ca.mcmaster.magarveylab.grape.enums.DomainEnums.PolyKetideDomainEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.SugarModificationsEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.TailoringDomainEnums;
import ca.mcmaster.magarveylab.grape.nrp.chem.ChemicalAbstraction;
import ca.mcmaster.magarveylab.grape.nrp.chem.Fragment;
import ca.mcmaster.magarveylab.grape.nrp.chem.Fragment.FragmentType;
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
	public static void createReport(String name, String smiles, String id, ChemicalAbstraction chemicalAbstraction) throws IOException, CDKException {
		new File("image/").mkdirs();
		IMolecule currentMolecule = null;
		if(name == null){
			name = id;
		}
		String filename = SmilesIO.getCleanFileName(name);
		if(chemicalAbstraction.getErrorMessage() != null){
			File output = new File("image/" +filename + ".txt");
			PrintWriter writer = new PrintWriter(output, "UTF-8");
			writer.println(chemicalAbstraction.getErrorMessage());
			writer.close();
			return;
		}
		
		List<Fragment> monomerFragments = chemicalAbstraction.getMonomerFragments();
		
		currentMolecule = SmilesIO.readSmiles(smiles);
		SmilesIO.drawMolecule(currentMolecule, "tempCurrentNrp");
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
		else {
			nrows = 6;
			ncols = 4;
		}
		
		for(int i = 0; i < monomerFragments.size(); i++) {
			SmilesIO.drawMolecule(monomerFragments.get(i).getMolecule(), "tempTempGeneratedFragment_" + i);
			
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
			
			String command = "convert -background white -fill black -pointsize 72 label:'" + label + "' image/tempLabel_"+i+".png; ";
			
			//TODO: add monomertype label. See:
			//convert Geodiamolide_A.png -gravity South   -background White   -splice 0x90  -pointsize 72          -annotate +0+2 'Faerie Dragon'   out.png
			
			command += "convert image/tempTempGeneratedFragment_" + i + ".png image/tempLabel_" + i + ".png -append image/tempGeneratedFragment_"+ i+".png;";
			
			
			ShellUtilities.runCommand(command);
		}
			
		String command = "";
		command += "montage image/tempCurrentNrp* image/tempGeneratedFragment_* -mode Concatenate -tile " + ncols + "x" + nrows + " image/temp_output_nolabel.png; ";
		command += "convert image/temp_output_nolabel.png -gravity South  -background White   -splice 0x180  -pointsize 72     -annotate +0+2 'Name: " + name + "\nChemical type: " + chemicalAbstraction.getChemicalType() + "'   image/" +filename + ".png; ";
		command += "rm image/temp*;";
		
		ShellUtilities.runCommand(command);
	}

}
