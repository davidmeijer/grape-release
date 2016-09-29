package ca.mcmaster.magarveylab.grape;

import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.openscience.cdk.exception.InvalidSmilesException;

/**
 * Main method for grape for command line interface.
 * @author prees cdejong
 * 
 */
//TODO: add output directory
 public class CliMain {
	/**
	 * Execute a Grape search.
	 * Usage: java -jar grape.jar -file sequence.fa [options]
	 * @param args	command line arguments
	 */
	public static void main(String[] args) {
		Options options = createOptions();
		String header = "GRAPE - Genetic Retro Assembly Prediction Engine, is a tool designed to deconstruct natural products into fragments that represent their biosynthetic components.";
		String footer = "Magarvey Lab 2015. Written by Chris Dejong and Greg Chen.";
		
		try {
			CommandLineParser parser = new GnuParser();
			CommandLine line = parser.parse(options, args, true);

			
			if (line.hasOption("h") || line.getOptions().length == 0) {
				HelpFormatter formatter = new HelpFormatter();
				formatter.printHelp("grape", header, options, footer, true);
			} else {
				GrapeConfig gc = GrapeConfig.getConfig();
			
				if (line.hasOption("v")) { 
					System.out.println(header);
					System.out.println("Grape: version " + gc.getVersion());
					System.out.println(footer);
				} else {
					line = parser.parse(options, args);
					parseCommandLine(line);
					if(gc.getSmiles().isEmpty() && gc.getFile() == null){
						throw new IllegalArgumentException();
					}
					//run grape here
					Grape grape = new Grape(gc);
					grape.run();
				}
			} 
		} catch (ParseException e) {
			System.err.println("Error: parsing command line arguments");
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("grape", header, options, footer, true);			
		} catch (IllegalArgumentException e) {
			System.err.println("Error: must specify input sequence file or a SMILES!");
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("grape", header, options, footer, true);
		} catch (InvalidSmilesException e) {
			System.err.println("Error: bad smiles when trying to Draw");
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("grape", header, options, footer, true);
		} catch (IOException e) {
			System.err.println("Error: cannot write files");
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("grape", header, options, footer, true);
		} finally {
			System.exit(0);
		}
	}

	/**
	 * Create command line options. 
	 * @return	newly created options
	 */
	@SuppressWarnings("static-access")
	public static Options createOptions() {
		Options options = new Options();
		
		// construct options with values
		Option file = OptionBuilder.withLongOpt("input file").hasArg().withArgName("FILE")
				.withDescription("Set the path to a data input file.").create("f");
		
		Option smiles = OptionBuilder.withLongOpt("smiles").hasArg().withArgName("SMILES")
				.withDescription("Pass in a molecule in SMILES format. Note: SMILES generally needs to be quoted depending on your Shell").create("s");
		
		Option id = OptionBuilder.withLongOpt("id").hasArg().withArgName("ID")
				.withDescription("Give the molelecule an ID number.").create("i");
		
		Option name = OptionBuilder.withLongOpt("name").hasArg().withArgName("NAME")
				.withDescription("Give the molecule a name.").create("n");

		Option aminoAcid = OptionBuilder.withLongOpt("aminoAcidPath").hasArg().withArgName("OUTPUT")
				.withDescription("Set amino acid file path.").create("a");
		
		Option output = OptionBuilder.withLongOpt("output").hasArg().withArgName("OUTPUT")
				.withDescription("Set an output folder path.").create("o");
		
		//construct boolean options
		Option images = new Option("img", "Images", false, "Create graphical images of output");
		Option json = new Option("json", "JSON", false, "Create JSON file output");
		Option text = new Option("txt", "text", false, "Create text file output");
		
		Option help = new Option("h", "help", false, "Print this message");
		Option version = new Option("v", "version", false, "Print the current version and exit");
	
		options.addOption(file);
		options.addOption(smiles);
		options.addOption(id);
		options.addOption(name);
		options.addOption(aminoAcid);
		options.addOption(output);

		options.addOption(images);
		options.addOption(json);
		options.addOption(text);

		options.addOption(help);
		options.addOption(version);

		return options;
	}
	
	/**
	 * Parse the command line input. 
	 * @param line	line to parse
	 * @throws ParseException
	 */
	public static void parseCommandLine(CommandLine line) throws ParseException {
		
		GrapeConfig gc = GrapeConfig.getConfig();
		//parse all other options here  and set them to GrapeConfig object.
		if (line.hasOption("f")) {
			gc.setFile(line.getOptionValue("f"));
		}
		if(line.hasOption("s")){
			gc.setSmiles(line.getOptionValue("s"));
		}
		if(line.hasOption("i")){
			gc.setID(line.getOptionValue("i"));
		}
		if(line.hasOption("n")){
			gc.setName(line.getOptionValue("n"));
		}
		if(line.hasOption("o")){
			gc.setOutputPath(line.getOptionValue("o"));
		}
		if(line.hasOption("a")){
			gc.setAminoAcidsPath(line.getOptionValue("a"));
		}
		if(line.hasOption("img")){
			gc.setImageTrue();
		}
		if(line.hasOption("json")){
			gc.setJsonTrue();
		}
		if(line.hasOption("txt")){
			gc.setTxtTrue();
		}
	}
}
