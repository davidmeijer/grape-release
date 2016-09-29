package ca.mcmaster.magarveylab.grape;

import java.io.File;

public class GrapeConfig {
	
	private static GrapeConfig config = new GrapeConfig();
	
	private String version = "2.9.1";
	private File file = null;
	private String smiles = new String();
	private String name = new String();
	private String id = new String();
	private boolean image = false;
	private boolean json = false;
	private boolean txt = false;
	private String outputPath = "./";
	private String aminoAcidsPath = "data/amino_acids_with_domains.txt";
	
	private GrapeConfig() {}
	
	//Getters
	public static GrapeConfig getConfig() {
		return config;
	}
	public String getVersion() {
		return version;
	}
	public File getFile() {
		return file;	
	}
	public String getSmiles() {
		return smiles;
	}
	public String getName(){
		return name;
	}
	public String getID() {
		return id;
	}
	public boolean image(){
		return image;
	}
	public boolean json(){
		return json;
	}
	public boolean txt(){
		return txt;
	}
	public String getOutputPath(){
		return outputPath;
	}
	
	public String getAminoAcidsPath() {
		return aminoAcidsPath;
	}
	
	//Setters
	public void setFile(String filePath) {
		this.file = new File(filePath);
	}
	public void setSmiles(String smiles) {
		this.smiles = smiles;
	}
	public void setName(String name){
		this.name = name;
	}
	public void setID(String id) {
		this.id = id;
	}
	public void setImageTrue(){
		image = true;
	}
	public void setJsonTrue(){
		json = true;
	}
	public void setTxtTrue(){
		txt = true;
	}
	public void setOutputPath(String outputPath){
		this.outputPath = outputPath + "/";
	}
	public void setAminoAcidsPath(String aminoAcidsPath) {
		this.aminoAcidsPath = aminoAcidsPath;
	}
}
