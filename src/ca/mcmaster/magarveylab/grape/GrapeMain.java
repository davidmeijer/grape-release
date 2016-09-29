package ca.mcmaster.magarveylab.grape;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;

import ca.mcmaster.magarveylab.grape.classify.InputEntry;
import ca.mcmaster.magarveylab.grape.nrp.chem.ChemicalAbstraction;
import ca.mcmaster.magarveylab.grape.nrp.chem.NRPPredictor;
import ca.mcmaster.magarveylab.grape.util.ChemicalUtilities;
import ca.mcmaster.magarveylab.grape.util.io.ImageOutput;
import ca.mcmaster.magarveylab.grape.util.io.JsonOutput;
import ca.mcmaster.magarveylab.grape.util.io.ReadFile;
import ca.mcmaster.magarveylab.grape.util.io.SmilesIO;
import ca.mcmaster.magarveylab.grape.util.io.TextOutput;

/* TODO imports not used
import ca.mcmaster.magarveylab.grape.pk.modules.Modifications;
import com.ibm.icu.text.SimpleDateFormat;
import ca.mcmaster.magarveylab.grape.enums.MoleculeClasses.ChemicalSubType;
*/

/**
 * The main class for GRAPE
 * @author gmchen
 */
@Deprecated
public class GrapeMain {

	static Map<String, String> data = new HashMap<String, String>();
	public static String currentName = "";
	/**
	 * Main method.
	 * @param args first argument should be a file path to a file with tab delimited data: name tab smiles. With one line per entry.
	 * @throws IOException 
	 * @throws CDKException 
	 * @throws ExceptionSmilesNotConnected 
	 * @throws ExceptionSmilesTooLarge 
	 */
	public static void main(String[] args) throws IOException, CDKException{
		
		if((args.length) == 0) {
			System.out.println("File name not provided - exiting.");
			System.exit(1);
		}
		
		String inputFileName = args[0];
		data = ReadFile.read(new File(inputFileName));
		NRPPredictor predictor = new NRPPredictor("data/monomer_types/amino_acids_with_domains.txt");
		runAllMolecules(predictor);
	}
	/**
	 * Run all molecules through reverse synthesis
	 * @throws IOException
	 * @throws CDKException
	 */
	private static void runAllMolecules(NRPPredictor predictor) throws IOException, CDKException{
		
		int x = 0;
		for(Entry<String, String> set : data.entrySet()) {
			String name = set.getKey();
			currentName = name;
			System.out.println("Name: " + name);
			
			ChemicalAbstraction chemicalAbstraction = ChemicalUtilities.getChemicalAbstractionFromSmiles(set.getValue(), predictor);
			InputEntry inputEntry = new InputEntry(set.getValue(), SmilesIO.getCleanFileName(name), 1);
			if(chemicalAbstraction != null) {
		
				String text = TextOutput.createOutputFile(chemicalAbstraction);
				
				//TextOutput.writeFile(text, SmilesIO.getCleanFileName(name) + ".txt");

				JsonOutput.writeJSON(".",inputEntry.getName(), inputEntry.getSmiles(), String.valueOf(x), chemicalAbstraction);

				ImageOutput.createReport(".",inputEntry.getName(), inputEntry.getSmiles(), String.valueOf(x), chemicalAbstraction);
			}
			x ++;
		}
	}

}
