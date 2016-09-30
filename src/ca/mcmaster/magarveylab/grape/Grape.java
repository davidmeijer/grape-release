package ca.mcmaster.magarveylab.grape;

import java.io.IOException;
import java.util.Map;
import java.util.Map.Entry;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.grape.nrp.chem.ChemicalAbstraction;
import ca.mcmaster.magarveylab.grape.nrp.chem.MoleculePredictor;
import ca.mcmaster.magarveylab.grape.util.ChemicalUtilities;
import ca.mcmaster.magarveylab.grape.util.io.ImageOutput;
import ca.mcmaster.magarveylab.grape.util.io.JsonOutput;
import ca.mcmaster.magarveylab.grape.util.io.ReadSmilesFile;
import ca.mcmaster.magarveylab.grape.util.io.TextOutput;

public class Grape {
	
	MoleculePredictor predictor;
	private GrapeConfig gc;

	public Grape(GrapeConfig gc){
		this.gc = gc;
		predictor =  new MoleculePredictor(gc.getAminoAcidsPath());
	}
	
	/**
	 * Run the analysis
	 * @throws InvalidSmilesException 
	 */
	public void run() throws IOException, InvalidSmilesException{
		
		if(!gc.getSmiles().isEmpty()){ //Single compound to run
			ChemicalAbstraction chemicalAbstraction = ChemicalUtilities.getChemicalAbstractionFromSmiles(gc.getSmiles(), gc.getName(), predictor);
			if(gc.txt()){
				String textOutput = TextOutput.createOutputFile(chemicalAbstraction);
				if(!gc.getID().isEmpty()){
					TextOutput.writeFile(textOutput, gc.getOutputPath(), gc.getID() + ".txt");
				}else{
					TextOutput.writeFile(textOutput, gc.getOutputPath(), gc.getName() + ".txt");
				}	
			}
			if(gc.json()){
				JsonOutput.writeJSON(gc.getOutputPath(), gc.getName(), gc.getSmiles(), gc.getID(), chemicalAbstraction);
			}
			if(gc.image()){
				try {
					ImageOutput.createReport(gc.getOutputPath(), gc.getName(), gc.getSmiles(), gc.getID(), chemicalAbstraction);
				} catch (CDKException e) {
					System.err.println("Could not draw molecule");
				}
			}							
		}else if(gc.getFile() != null){ //many compounds to run
			Map<String, String> inputData = ReadSmilesFile.read(gc.getFile());
			for(Entry<String, String> single : inputData.entrySet()){
				String name = single.getKey();
				String smiles = single.getValue().replace(" ", "");
				System.out.println("Working on: " + name);
				ChemicalAbstraction chemicalAbstraction = ChemicalUtilities.getChemicalAbstractionFromSmiles(smiles, name, predictor);
				if(chemicalAbstraction.errored()){
					System.err.println(chemicalAbstraction.getErrorMessage());
					continue;
				}
				if(gc.txt()){
					String textOutput = TextOutput.createOutputFile(chemicalAbstraction);
					TextOutput.writeFile(textOutput, gc.getOutputPath(), name + ".txt");
				}
				if(gc.json()){
					JsonOutput.writeJSON(gc.getOutputPath(), name, smiles, name, chemicalAbstraction);
				}
				if(gc.image()){
					try{
						ImageOutput.createReport(gc.getOutputPath(), name, smiles, name, chemicalAbstraction);
					}catch(Exception e){
						System.err.println("Could not draw: " + name);
					}
				}
			}
		}
	}
}
