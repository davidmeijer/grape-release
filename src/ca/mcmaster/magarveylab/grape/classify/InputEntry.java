package ca.mcmaster.magarveylab.grape.classify;


public class InputEntry {
	private final String smiles;
	private final String name;
	private final int id;
	private boolean badSmiles = false;
	
	public InputEntry(String smiles, String name, int id) {
		this.smiles = smiles;
		this.name = name;
		this.id = id;
	}
	
	public InputEntry(String smiles, String name, String id) {
		this.smiles = smiles;
		this.name = name;
		this.id = Integer.parseInt(id);
	}
	
	public InputEntry(String smiles, String id) {
		this.smiles = smiles;
		this.name = null;
		this.id = Integer.parseInt(id);
	}
	
	public String getSmiles() {
		return this.smiles;
	}
	
	public String getName() {
		return this.name;
	}
	
	public int getId() {
		return this.id;
	}
	
	public void badSmiles() {
		badSmiles = true;
	}
	
	public boolean isBadSmiles() {
		return badSmiles;
	}
}
