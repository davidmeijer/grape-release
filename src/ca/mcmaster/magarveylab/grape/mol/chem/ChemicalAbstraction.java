package ca.mcmaster.magarveylab.grape.mol.chem;

import java.util.ArrayList;
import java.util.List;

import ca.mcmaster.magarveylab.grape.enums.MoleculeClasses.ChemicalSubType;
import ca.mcmaster.magarveylab.grape.enums.MoleculeClasses.ChemicalType;

public class ChemicalAbstraction {
	private List<Fragment> monomerFragments  = new ArrayList<Fragment>();
	private ChemicalType chemicalType;
	private String name;
	private ChemicalSubType chemicalSubType = ChemicalSubType.NONE;
	private String chemicalScaffoldName = null;
	private String errorMessage;
	private boolean error = false;
	
	public ChemicalAbstraction(String name) {
		this.name = name;
	}

	public String getName() {
		return name;
	}

	/**
	 * @return the monomerFragments
	 */
	public List<Fragment> getMonomerFragments() {
		return monomerFragments;
	}

	/**
	 * @param monomerFragments the monomerFragments to set
	 */
	public void setMonomerFragments(List<Fragment> monomerFragments) {
		if (monomerFragments != null) {
			this.monomerFragments = monomerFragments;
		}
	}


	/**
	 * @return the chemicalType
	 */
	public ChemicalType getChemicalType() {
		return chemicalType;
	}

	/**
	 * @param chemicalType the chemicalType to set
	 */
	public void setChemicalType(ChemicalType chemicalType) {
		this.chemicalType = chemicalType;
	}

	/**
	 * @return the chemicalSubType
	 */
	public ChemicalSubType getChemicalSubType() {
		return chemicalSubType;
	}

	/**
	 * @param chemicalType the chemicalSubType to set
	 */
	public void setChemicalSubType(ChemicalSubType chem) {
		this.chemicalSubType = chem;
	}
	
	public void setChemicalScaffoldName(String chemicalScaffoldName){
		this.chemicalScaffoldName = chemicalScaffoldName;
	}
	
	public String getChemicalScaffoldName(){
		return chemicalScaffoldName;
	}

	public void setErrorMessage(String errorMessage) {
		this.errorMessage = errorMessage;	
		error = true;
	}
	public String getErrorMessage(){
		return errorMessage;
	}
	public boolean errored(){
		return error;
	}
}
