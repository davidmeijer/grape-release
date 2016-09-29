package ca.mcmaster.magarveylab.grape.nrp.chem;

import java.util.ArrayList;
import java.util.List;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IMoleculeSet;

import ca.mcmaster.magarveylab.grape.enums.AcylAdenylatingSubstrates;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.AminoAcidEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.PolyKetideDomainEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.SugarModificationsEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.TailoringDomainEnums;
import ca.mcmaster.magarveylab.grape.enums.KnownOtherEnums;
import ca.mcmaster.magarveylab.grape.pk.modules.PKsubstrate;
import ca.mcmaster.magarveylab.grape.pk.modules.Modifications;
import ca.mcmaster.magarveylab.grape.util.ShellUtilities;
import ca.mcmaster.magarveylab.grape.util.io.SmilesIO;

//import ca.mcmaster.magarveylab.grape.util.DomainEnums.AcyltransferaseSubstratesEnum;
//import ca.mcmaster.magarveylab.grape.util.DomainEnums.AdenylationSubstratesEnum;

/**
 * Fragment class specifies objects containing information from NRP monomers: molecule object,
 * most likely identity, connectivity, pointers to important atoms, and pointers to connected fragments.
 * @author gmchen
 */
public class Fragment 
{	
	public enum FragmentType {
	    AMINO_ACID, MULTIPLE_AMINO_ACID_PIECE, FATTY_ACID, SUGAR, POLYKETIDE, FA_OR_PK, CYCLOHEXANE, UNKNOWN_OTHER, KNOWN_OTHER, ACYL_ADENYLATING
	}
	// Note: FA_OR_PK is the generic type of identified PK units that may not be PKs - most likely fatty acids, possibly modified amino acids or others
	
	private FragmentType fragmenType;
	
	private IMolecule molecule = null;
	private IAtom thiazoleS = null;
	private IAtom betaLactamS = null;
	private IAtom oxazoleO = null;
	private List<IAtom> aminoNs = new ArrayList<IAtom>();
	private List<IAtom> aminoCs = new ArrayList<IAtom>();
	private IAtom lactoneCarboxylC = null;
	private IAtom lactoneHydroxylO = null;
	private IAtom lactoneHydroxylC = null;
	// References to atoms in other monomer fragments, used to retain graph structure of the fragments.
	private IAtom atomAfterCTerminus = null;
	private IAtom atomAfterNTerminus = null;
	private IAtom atomInAttachedSugar = null;
	private IAtom atomOppositeSugarBond = null; // If this is a sugar, this atom is in the piece to which the sugar is attached
	private IAtom atomAfterLactoneCarboxyl = null;
	private IAtom atomAfterLactoneHydroxyl = null;
	private List<String> sugarNames = new ArrayList<String>();
	private double tanimotoScore = -1;
	// A boolean to track whether connectivity of this monomer fragment is reliable.
	private boolean ignoreConnectivity = false;
	// A list of domains
	private List<AminoAcidEnums> aminoAcidDomains = new ArrayList<AminoAcidEnums>();
	private KnownOtherEnums knownOther = null;
	private List<TailoringDomainEnums> tailoringDomains = new ArrayList<TailoringDomainEnums>();
	private List<List<PolyKetideDomainEnums>> pkDomains = new ArrayList<List<PolyKetideDomainEnums>>();
	private List<SugarModificationsEnums> sugarModifications = new ArrayList<SugarModificationsEnums>();
	private List<AcylAdenylatingSubstrates> starters = new ArrayList<AcylAdenylatingSubstrates>();
	private Modifications pkModifications;
	private List<PKsubstrate> loadingUnits;
	// For multiple amino acid piece, track number of nitrogens
	private int numNitrogens = -1;
	//private List<AcyltransferaseSubstratesEnum> acyltranferaseDomains = new ArrayList<AcyltransferaseSubstratesEnum>();
	//Set to false when this fragment is part of a larger cyclic list of ordered monomer fragments. This means the start is unknown
	private boolean knownStart = true;
	
	/**
	 * Consstruct a monomer fragment from a molecule
	 * @param molecule
	 */
	public Fragment(IMolecule molecule) {
		this.molecule = molecule;
	}
	
	/**
	 * Partition this Fragment into monomer fragments based on the connected components of the molecule. Contents of domains
	 * are not preserved, and therefore will need to be repopulated.
	 * @return An ArrayList of monomer fragments based on connected components.
	 */
	public ArrayList<Fragment> partitionIntoMonomerFragments() {
		ArrayList<Fragment> monomerFragments = new ArrayList<Fragment>();
		
		IMoleculeSet partitions = ConnectivityChecker.partitionIntoMolecules(molecule);
		
		if(partitions.getMoleculeCount() == 1) {
			monomerFragments.add(this);
			return(monomerFragments);
		}
		
		for(int i = 0; i < partitions.getMoleculeCount(); i++) {
			IMolecule currentMolecule = partitions.getMolecule(i);
			Fragment newFragment = new Fragment(currentMolecule);
			if(currentMolecule.contains(thiazoleS)) {
				newFragment.setThiazoleS(thiazoleS);
			}
			if(currentMolecule.contains(oxazoleO)) {
				newFragment.setOxazoleO(oxazoleO);
			}
			for(IAtom aminoN : aminoNs){
				if(currentMolecule.contains(aminoN)) {
					newFragment.addAminoN(aminoN);
					newFragment.setAtomAfterNTerminus(atomAfterNTerminus);
				}	
			}
			for(IAtom aminoC : aminoCs){
				if(currentMolecule.contains(aminoC)) {
					newFragment.addAminoC(aminoC);
					newFragment.setAtomAfterCTerminus(atomAfterCTerminus);
				}
			}
			if(currentMolecule.contains(lactoneCarboxylC)) {
				newFragment.setLactoneCarboxylC(lactoneCarboxylC);
				newFragment.setAtomAfterLactoneCarboxyl(atomAfterLactoneCarboxyl);
			}
			if(currentMolecule.contains(lactoneHydroxylO)) {
				newFragment.setLactoneHydroxylO(lactoneHydroxylO);
				newFragment.setAtomAfterLactoneHydroxyl(atomAfterLactoneHydroxyl);
				newFragment.setLactoneHydroxylC(lactoneHydroxylC);
			}
			if(currentMolecule.contains(atomInAttachedSugar)) {
				newFragment.setAtomInAttachedSugar(atomInAttachedSugar);
			}
			if(currentMolecule.contains(atomOppositeSugarBond)) {
				newFragment.setAtomOppositeSugarBond(atomOppositeSugarBond);
			}
			if(currentMolecule.contains(betaLactamS)) {
				newFragment.setBetaLactamS(betaLactamS);
			}
			monomerFragments.add(newFragment);
		}
		return monomerFragments;
	}
	
	
	
	/**
	 * Given a list of MonomerFragments, identify the one that contains the atom after the N terminus
	 * @param monomerFragmentList
	 * @return
	 */
	public Fragment getFragmentAfterNTerminus(List<Fragment> monomerFragmentList) {
		if(atomAfterNTerminus == null) {
			return null;
		}
		for(Fragment m : monomerFragmentList) {
			if(m.getMolecule().contains(atomAfterNTerminus)) {
				return(m);
			}
		}
		return null;
	}
	
	/**
	 * Given a list of MonomerFragments, identify the one that contains the atom after the C terminus
	 * @param monomerFragmentList
	 * @return
	 */
	public Fragment getFragmentAfterCTerminus(List<Fragment> monomerFragmentList) {
		if(atomAfterCTerminus == null) {
			return null;
		}
		for(Fragment m : monomerFragmentList) {
			if(m.getMolecule().contains(atomAfterCTerminus)) {
				return(m);
			}
		}
		return null;
	}
	
	public IMolecule getMolecule() {
		return this.molecule;
	}
	public void setMolecule(IMolecule newMolecule) {
		this.molecule = newMolecule;
	}
	public List<IAtom> getAminoNs() {
		return aminoNs;
	}

	public void addAminoN(IAtom aminoN) {
		if(aminoNs.size() > 0) {
			this.fragmenType = FragmentType.MULTIPLE_AMINO_ACID_PIECE;
		}
		aminoNs.add(aminoN);
	}

	public List<IAtom> getAminoCs() {
		return(aminoCs);
	}

	public void addAminoC(IAtom aminoC) {
		if(aminoCs.size() > 0) {
			this.fragmenType = FragmentType.MULTIPLE_AMINO_ACID_PIECE;
		}
		aminoCs.add(aminoC);
	}

	public IAtom getLactoneCarboxylC() {
		return lactoneCarboxylC;
	}

	public void setLactoneCarboxylC(IAtom lactoneCarboxylC) {
		this.lactoneCarboxylC = lactoneCarboxylC;
	}

	public IAtom getLactoneHydroxylO() {
		return lactoneHydroxylO;
	}

	public void setLactoneHydroxylO(IAtom lactoneHydroxylO) {
		this.lactoneHydroxylO = lactoneHydroxylO;
	}
	
	/**
	 * @return the carbonToSugar
	 */
	public IAtom getAtomInAttachedSugar() {
		return atomInAttachedSugar;
	}

	/**
	 * @param carbonToSugar the carbonToSugar to set
	 */
	public void setAtomInAttachedSugar(IAtom atomAfterSugarBond) {
		this.atomInAttachedSugar = atomAfterSugarBond;
	}

	/**
	 * @return the thiazoleS
	 */
	public IAtom getThiazoleS() {
		return thiazoleS;
	}

	/**
	 * @param thiazoleS the thiazoleS to set
	 */
	public void setThiazoleS(IAtom thiazoleS) {
		this.thiazoleS = thiazoleS;
	}

	/**
	 * @return the oxazoleO
	 */
	public IAtom getOxazoleO() {
		return oxazoleO;
	}

	/**
	 * @param oxazoleO the oxazoleO to set
	 */
	public void setOxazoleO(IAtom oxazoleO) {
		this.oxazoleO = oxazoleO;
	}

	/**
	 * @return the atomAfterCTerminus
	 */
	public IAtom getAtomAfterCTerminus() {
		return atomAfterCTerminus;
	}

	/**
	 * @param atomAfterCTerminus the atomAfterCTerminus to set
	 */
	public void setAtomAfterCTerminus(IAtom atomAfterCTerminus) {
		this.atomAfterCTerminus = atomAfterCTerminus;
	}

	/**
	 * @return the atomAfterNTerminus
	 */
	public IAtom getAtomAfterNTerminus() {
		return atomAfterNTerminus;
	}

	/**
	 * @param atomAfterNTerminus the atomAfterNTerminus to set
	 */
	public void setAtomAfterNTerminus(IAtom atomAfterNTerminus) {
		this.atomAfterNTerminus = atomAfterNTerminus;
	}

	/**
	 * @return the atomAfterLactoneCarboxyl
	 */
	public IAtom getAtomAfterLactoneCarboxyl() {
		return atomAfterLactoneCarboxyl;
	}

	/**
	 * @param atomAfterLactoneCarboxyl the atomAfterLactoneCarboxyl to set
	 */
	public void setAtomAfterLactoneCarboxyl(IAtom atomAfterLactoneCarboxyl) {
		this.atomAfterLactoneCarboxyl = atomAfterLactoneCarboxyl;
	}

	/**
	 * @return the atomAfterLactoneHydroxyl
	 */
	public IAtom getAtomAfterLactoneHydroxyl() {
		return atomAfterLactoneHydroxyl;
	}

	/**
	 * @param atomAfterLactoneHydroxyl the atomAfterLactoneHydroxyl to set
	 */
	public void setAtomAfterLactoneHydroxyl(IAtom atomAfterLactoneHydroxyl) {
		this.atomAfterLactoneHydroxyl = atomAfterLactoneHydroxyl;
	}

	/**
	 * @return the adenylationDomains
	 *//*
	public List<AdenylationSubstratesEnum> getAdenylationDomains() {
		return adenylationDomains;
	}

	/**
	 * @param adenylationDomains the adenylationDomains to set
	 *//*
	public void setAdenylationDomains(
			List<AdenylationSubstratesEnum> adenylationDomains) {
		this.adenylationDomains = adenylationDomains;
	}

	/**
	 * @return the tailoringDomains
	 */
	public List<TailoringDomainEnums> getTailoringDomains() {
		return tailoringDomains;
	}

	/**
	 * @param tailoringDomains the tailoringDomains to set
	 */
	public void setTailoringDomains(List<TailoringDomainEnums> tailoringDomains) {
		this.tailoringDomains = tailoringDomains;
	}
	
	public void setKnownOther(KnownOtherEnums knownOther){
		this.knownOther = knownOther;
	}
	
	public KnownOtherEnums getKnownOther(){
		return knownOther;
	}

	/**
	 * @return the acyltranferaseDomains
	 *//*
	public List<AcyltransferaseSubstratesEnum> getAcyltranferaseDomains() {
		return acyltranferaseDomains;
	}

	/**
	 * @param acyltranferaseDomains the acyltranferaseDomains to set
	 *//*
	public void setAcyltranferaseDomains(
			List<AcyltransferaseSubstratesEnum> acyltranferaseDomains) {
		this.acyltranferaseDomains = acyltranferaseDomains;
	}*/

	/**
	 * Set the amino acid domain
	 * @param newDomain
	 */
	public void addAminoAcidDomain(AminoAcidEnums newDomain) {
		aminoAcidDomains.add(newDomain);
	}
	
	public void addAminoAcidDomains(AminoAcidEnums[] newDomains) {
		for(AminoAcidEnums newDomain : newDomains){
			aminoAcidDomains.add(newDomain);
		}
	}
	
	/**
	 * Get the amino acid domain
	 * @return The amino acid domain
	 */
	public List<AminoAcidEnums> getAminoAcidDomains() {
		return aminoAcidDomains;
	}
	/**
	 * Add a tailoring domain, if it has not already been added.
	 * @param newDomain
	 */
	public void addTailoringDomain(TailoringDomainEnums newDomain) {
		if(!tailoringDomains.contains(newDomain)) {
			tailoringDomains.add(newDomain);
		}
	}
	
	public void addStarter(AcylAdenylatingSubstrates starter){
		starters.add(starter);
	}
	
	public List<AcylAdenylatingSubstrates> getStarters(){
		return starters;
	}
	/**
	 * Add a sugar modification, if it has not already been added.
	 * @param sugarModification
	 */
	public void addSugarModification(SugarModificationsEnums sugarModification) {
		if(!sugarModifications.contains(sugarModification)) {
			sugarModifications.add(sugarModification);
		}
	}
	/**
	 * Add an acyltransferase domain, if it ha not already been added.
	 * @param newDomain
	 *//*
	public void addAcyltransferaseDomain(AcyltransferaseSubstratesEnum newDomain) {
		if(!acyltranferaseDomains.contains(newDomain)) {
			acyltranferaseDomains.add(newDomain);
		}
	}*/

	/**
	 * @return the sugarModifications
	 */
	public List<SugarModificationsEnums> getSugarModifications() {
		return sugarModifications;
	}

	/**
	 * @param sugarModifications the sugarModifications to set
	 */
	public void setSugarModifications(List<SugarModificationsEnums> sugarModifications) {
		this.sugarModifications = sugarModifications;
	}

	/**
	 * @return the fragmenType
	 */
	public FragmentType getFragmentType() {
		return fragmenType;
	}

	/**
	 * @param fragmenType the fragmenType to set
	 */
	public void setFragmentType(FragmentType fragmentType) {
		this.fragmenType = fragmentType;
	}

	/**
	 * @return the tanimotoScore
	 */
	public double getTanimotoScore() {
		return tanimotoScore;
	}

	/**
	 * @param tanimotoScore the tanimotoScore to set
	 */
	public void setTanimotoScore(double tanimotoScore) {
		this.tanimotoScore = tanimotoScore;
	}

	/**
	 * @return the pkDomains
	 */
	public List<List<PolyKetideDomainEnums>> getPkDomains() {
		return pkDomains;
	}

	/**
	 * Set the PK domains. Remove any END domains
	 * @param pkDomains the pkDomains to set
	 */
	public void setPkDomains(List<List<PolyKetideDomainEnums>> polyketideDomains) {
		this.pkDomains = polyketideDomains;
	}

	/**
	 * @return the numNitrogens
	 */
	public int getNumNitrogens() {
		return numNitrogens;
	}

	/**
	 * @param numNitrogens the numNitrogens to set
	 */
	public void setNumNitrogens(int numNitrogens) {
		this.numNitrogens = numNitrogens;
	}

	/**
	 * @return the sugarName
	 */
	public List<String> getSugarNames() {
		return sugarNames;
	}

	/**
	 * @param sugarName the sugarName to set
	 */
	public void addSugarName(String sugarName) {
		this.sugarNames.add(sugarName);
	}
	
	/**
	 * @return the lactoneHydroxylC
	 */
	public IAtom getLactoneHydroxylC() {
		return lactoneHydroxylC;
	}

	/**
	 * @param lactoneHydroxylC the lactoneHydroxylC to set
	 */
	public void setLactoneHydroxylC(IAtom lactoneHydroxylC) {
		this.lactoneHydroxylC = lactoneHydroxylC;
	}

	/**
	 * @return the pkModifications
	 */
	public Modifications getPkModifications() {
		return pkModifications;
	}

	/**
	 * @param pkModifications the pkModifications to set
	 */
	public void setPkModifications(Modifications pkModifications) {
		this.pkModifications = pkModifications;
	}

	/**
	 * @return the loadingUnits
	 */
	public List<PKsubstrate> getLoadingUnits() {
		return loadingUnits;
	}

	/**
	 * @param loadingUnits the loadingUnits to set
	 */
	public void setLoadingUnits(List<PKsubstrate> loadingUnits) {
		this.loadingUnits = loadingUnits;
	}

	/**
	 * @return the atomOppositeSugarBond
	 */
	public IAtom getAtomOppositeSugarBond() {
		return atomOppositeSugarBond;
	}

	/**
	 * @param atomOppositeSugarBond the atomOppositeSugarBond to set
	 */
	public void setAtomOppositeSugarBond(IAtom atomOppositeSugarBond) {
		this.atomOppositeSugarBond = atomOppositeSugarBond;
	}
	
	/**
	 * Set to false when this fragment is part of a larger cyclic list of ordered monomer fragments. This means the start is unknown
	 * @param b
	 */
	public void setKnownStart(boolean b) {
		knownStart = b;
	}
	
	/**
	 * false when this fragment is part of a larger cyclic list of ordered monomer fragments. This means the start is unknown.
	 * @return
	 */
	public Boolean hasKnownStart(){
		return knownStart;
	}
	
	/**
	 * If this monomer fragment has an attached sugar, return that sugar from a list
	 * @param monomerFragments
	 * @return
	 */
	public Fragment getAttachedSugar(List<Fragment> monomerFragments) {
		if(this.getAtomInAttachedSugar() == null) {
			return null;
		}
		for(Fragment m : monomerFragments) {
			if(m.getMolecule().contains(this.atomInAttachedSugar)) {
				return m;
			}
		}
		return null;
	}
	
	/**
	 * If this monomer fragment is a sugar attached to another fragment, return that other fragment from a list
	 */
	public Fragment getFragmentOppositeSugarBond(List<Fragment> monomerFragments) {
		if(this.getAtomOppositeSugarBond() == null) {
			return null;
		}
		for(Fragment m : monomerFragments) {
			if(m.getMolecule().contains(this.atomOppositeSugarBond)) {
				return m;
			}
		}
		return null;
	}

	/**
	 * @return the betaLactamS
	 */
	public IAtom getBetaLactamS() {
		return betaLactamS;
	}

	/**
	 * @param betaLactamS the betaLactamS to set
	 */
	public void setBetaLactamS(IAtom betaLactamS) {
		this.betaLactamS = betaLactamS;
	}

	public static void drawMonomerFragments(
			List<Fragment> monomerFragments, String name) {
		name = SmilesIO.getCleanFileName(name);
		int numImages = monomerFragments.size();
		int nrows, ncols;
	
		if (numImages <= 12) {
			nrows = 4;
			ncols = 3;
		} else {
			nrows = 6;
			ncols = 4;
		}
		for (int i = 0; i < monomerFragments.size(); i++) {
			SmilesIO.drawMolecule(monomerFragments.get(i).getMolecule(),
					"temp_step_" + i);
		}
		String command = "montage images/temp_step_* -mode concatenate -tile "
				+ ncols + "x" + nrows + " images/temp_step_montage.png; ";
		command += "convert -background white -fill black -pointsize 72 label:'"
				+ name + "' images/temp_step_label.png; ";
		command += "convert -colorspace RGB images/temp_step_label.png images/temp_step_montage.png -append images/"
				+ name + ".png; ";
		command += "rm images/temp_step_*;";
	
		ShellUtilities.runCommand(command);
	}
}
