package ca.mcmaster.magarveylab.grape.util;

import java.util.ArrayList;

public class MonomerMatch {
	
	private ArrayList<String> monomerNames;
	private double tanimotoScore;
	public MonomerMatch(ArrayList<String> monomerNames, double tanimotoScore) {
		this.monomerNames = monomerNames;
		this.tanimotoScore = tanimotoScore;
	}
	public ArrayList<String> getMonomerNames() {
		return monomerNames;
	}
	public void setMonomerNames(ArrayList<String> monomerNames) {
		this.monomerNames = monomerNames;
	}
	public double getTanimotoScore() {
		return tanimotoScore;
	}
	public void setTanimotoScore(double tanimotoScore) {
		this.tanimotoScore = tanimotoScore;
	}
}
