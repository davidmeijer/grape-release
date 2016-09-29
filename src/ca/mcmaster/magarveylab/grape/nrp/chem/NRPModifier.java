package ca.mcmaster.magarveylab.grape.nrp.chem;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;

import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.isomorphism.mcss.RMap;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;

import ca.mcmaster.magarveylab.grape.GrapeMain;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.PolyKetideDomainEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.TailoringDomainEnums;
import ca.mcmaster.magarveylab.grape.nrp.chem.Fragment.FragmentType;
import ca.mcmaster.magarveylab.grape.pk.modules.PKsubstrate;
import ca.mcmaster.magarveylab.grape.util.ChemicalUtilities;
import ca.mcmaster.magarveylab.grape.util.io.SmilesIO;

/**
 * Perform modifications to NRPs, breaking them into constitudent monomer units
 * 
 * @author gmchen cDejong
 */
public class NRPModifier {

	private final static boolean imageDump = false;
	private final static boolean showBonds = false;
	private final static boolean printSteps = false;

	private boolean fungal = false; // set true if fungal chemistry should be done
	private List<Fragment> monomerFragments;

	// This global list is necessary for determining the presence and
	// directionality of lactone bonds in the case that there is one lactone
	// whose adjacent monomers are separated in the peptide cleavage stage.
	private List<IAtom> lactoneCarboxylCList;
	private List<IAtom> thiazoleSList;
	private List<IAtom> betaLactamSList;
	private List<IAtom> oxazoleOList;
	private List<IAtom> sugarOxygens;
	private ChemicalAbstraction chemicalAbstraction;

	/**
	 * Constructor for this NRPModifier. Typical usage: construct an NRPModifier
	 * from an original NRP IMolecule, and get the output from method
	 * getMonomerFragments.
	 * 
	 * @param originalMolecule
	 */
	public NRPModifier(IMolecule originalMolecule, ChemicalAbstraction chemicalAbstraction) {
		lactoneCarboxylCList = new ArrayList<IAtom>();
		thiazoleSList = new ArrayList<IAtom>();
		oxazoleOList = new ArrayList<IAtom>();
		sugarOxygens = new ArrayList<IAtom>();
		monomerFragments = new ArrayList<Fragment>();
		monomerFragments.add(new Fragment(originalMolecule));
		this.chemicalAbstraction = chemicalAbstraction; 
		performAllNrpModifications();
	}
	
	public NRPModifier(IMolecule originalMolecule, ChemicalAbstraction chemicalAbstraction, boolean fungal) {
		lactoneCarboxylCList = new ArrayList<IAtom>();
		thiazoleSList = new ArrayList<IAtom>();
		oxazoleOList = new ArrayList<IAtom>();
		sugarOxygens = new ArrayList<IAtom>();
		monomerFragments = new ArrayList<Fragment>();
		monomerFragments.add(new Fragment(originalMolecule));
		this.chemicalAbstraction = chemicalAbstraction;
		this.fungal = fungal;
		performAllNrpModifications();
	}
	/**
	 * Sets the modifier to also do fungal breakdowns
	 */
	public void fungal(){
		fungal = true;
	}

	/**
	 * Get the monomer fragments from the modifications performed by the
	 * NRPModifier
	 * 
	 * @return
	 */
	public List<Fragment> getMonomerFragments() {
		return monomerFragments;
	}

	/**
	 * Perform all NRP modifications in the intended order
	 */
	private void performAllNrpModifications() { //TODO start splitting these up into different classes
		// if(false)
		if (!imageDump) {
			if(fungal){
				breakFungalSubstructures();
			}
			breakDisulfideBridges(); //add fungal stuff to this
			breakAromaticEthers();
			breakAdjoinedAromaticRings();
			processBetaLactamLike();
			processMonobactams();
			breakLinearEthers();
			processUniqueSubstructures();
			processThiazoles();
			processOxazoles();
			modifyImine();
			breakThioesters();
			openLactoneRings();
			processUreido();
			breakEsterLinkages();
			breakPeptideBonds();
			breakSugarGroups();
			breakSulfateGroups();
			processAAMethylations();
			processAACHydroxyliations();
			processChlorinations();
			findEpoxides();
			processPeptideEpoxiKetones(); 
		}

		if (imageDump) {
			Fragment.drawMonomerFragments(monomerFragments, GrapeMain.currentName
					+ "_00-Original");
			breakDisulfideBridges();
			breakAromaticEthers();
			breakAdjoinedAromaticRings();
			Fragment.drawMonomerFragments(monomerFragments, GrapeMain.currentName
					+ "_01-Chemical_Bridging");
			processBetaLactamLike();
			processThiazoles();
			processOxazoles();
			Fragment.drawMonomerFragments(monomerFragments, GrapeMain.currentName
					+ "_02-Heterocycles");
			breakThioesters();
			openLactoneRings();
			processUreido();
			breakEsterLinkages();
			breakPeptideBonds();
			Fragment.drawMonomerFragments(monomerFragments, GrapeMain.currentName
					+ "_03-Core_Bonds");
			breakSugarGroups();
			breakSulfateGroups();
			processAAMethylations();
			processChlorinations();
			findEpoxides();
			Fragment.drawMonomerFragments(monomerFragments, GrapeMain.currentName
					+ "_04-Tailoring");
		}

		// Update thiazole and oxazole fields
		for (Fragment m : monomerFragments) {
			for (IAtom t : thiazoleSList) {
				if (m.getMolecule().contains(t)) {
					m.setThiazoleS(t);
					m.addTailoringDomain(TailoringDomainEnums.THIAZOLE);
				}
			}
			for (IAtom o : oxazoleOList) {
				if (m.getMolecule().contains(o)) {
					m.setOxazoleO(o);
					m.addTailoringDomain(TailoringDomainEnums.OXAZOLE);
				}
			}
			for (IAtom o : sugarOxygens) {
				if (m.getMolecule().contains(o)) {
					m.setFragmentType(FragmentType.SUGAR);
				}
			}
			if(m.getBetaLactamS() != null) {
				m.addTailoringDomain(TailoringDomainEnums.SULFUR_BETA_LACTAM);
			}
		}
	}

	/**
	 * Breaks fungal specific substructures
	 */
	private void breakFungalSubstructures() {
		processAminoDerivedTricyclicRing();
		processPrenylations();
		processDielsAlder();
		process2Pyridone();
		processBenzoquinone();
		processDecalin();
		processTetramicAcid();
	}
	private void processPrenylations() {
		
		SMARTSQueryTool querytool = null;
		try {
			querytool = new SMARTSQueryTool("C[C;D2]=[C;D3]([C;D1])[C;D1]");
		} catch (CDKException e) {
			e.printStackTrace();
		}
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		// go through each fragment.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IMolecule mol = fragment.getMolecule();
			try {
				while(querytool.matches(mol)){
					List<Integer> matches = querytool.getMatchingAtoms().get(0);
					fragment.addTailoringDomain(TailoringDomainEnums.PRENYLATION);
										
					//get the atoms that were matched
					List<IAtom> matchedAtoms = new ArrayList<IAtom>();
					for(Integer match : matches){
						matchedAtoms.add(mol.getAtom(match));
					}
					
					//add hydroxyl to the connected carbon to the group
					for(IAtom atom : matchedAtoms){
						for(IAtom connectedAtom : mol.getConnectedAtomsList(atom)){
							if(!matchedAtoms.contains(connectedAtom)){ //Atom that is not in the list must be the one that had a hydroxyl before
								Atom oxygen = new Atom("O");
								oxygen.setAtomTypeName("O.sp3");
								mol.addAtom(oxygen);
								mol.addBond(
										fragment.getMolecule().getAtomNumber(connectedAtom),
										fragment.getMolecule().getAtomNumber(oxygen),
										IBond.Order.SINGLE);
							}
						}
					}
					
					//remove the matched atoms
					for(IAtom atom : matchedAtoms){
						mol.removeAtomAndConnectedElectronContainers(atom);
					} // add hydroxyl
				}
			} catch (CDKException e) {
				e.printStackTrace();
			}
		}
	
	}
	private void processTetramicAcid() { // CC(=O)C1=CCNC1=O & CC(=O)C1CCNC1=O // check if dealing with epoxide first //C\C(O)=C1\CCNC1=O
		IMolecule[] templates = new IMolecule[3]; 
		try {
			templates[0] = SmilesIO.readSmiles("CC(=O)C1=CCNC1=O");
			templates[1] = SmilesIO.readSmiles("CC(=O)C1CCNC1=O"); 
			templates[2] = SmilesIO.readSmiles("C\\C(O)=C1\\CCNC1=O");
		} catch(Exception e) {
			e.printStackTrace();
		}
		IBond[] templateBondToBreak = new IBond[]{
				templates[0].getBond(templates[0].getAtom(3),templates[0].getAtom(4)),
				templates[1].getBond(templates[1].getAtom(3),templates[1].getAtom(4)),
				templates[2].getBond(templates[2].getAtom(3),templates[2].getAtom(4))
		};
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IMolecule mol = fragment.getMolecule();
			List<IBond> templateBondMatches = new ArrayList<IBond>();
			int templateIndex = 0; //template index to check
			boolean match = false; //change to true if a match comes up in the next statement
			while(!match && templates.length > templateIndex){
				templateBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(templates[templateIndex], templateBondToBreak[templateIndex], mol);
				if(templateBondMatches.size() == 0){
					templateIndex ++;
				}else{
					match = true;
				}
			}
			if(!match){ //if no matches breaks the loop
				break;
			}
			for(IBond bond : templateBondMatches){
				for(IAtom atom2 : mol.getConnectedAtomsList(bond.getAtom(0))){ // remove epoxide oxygen
					if(atom2.getAtomicNumber() == 8){
						mol.removeAtomAndConnectedElectronContainers(atom2);
						break;
					}
				}
				SmilesIO.drawMoleculeHighlightingBonds(mol, "tets", bond);
				mol.removeBond(bond); //break bond that matched template
				for(IAtom atom : bond.atoms()){
					if(mol.getConnectedAtomsCount(atom) == 1){
						Atom dub = new Atom("O");
						dub.setAtomTypeName("O.sp2");
						Atom single = new Atom("O");
						single.setAtomTypeName("O.sp3");
						mol.addAtom(dub);
						mol.addBond(new Bond(atom, dub, IBond.Order.DOUBLE));
						mol.addAtom(single);
						mol.addBond(new Bond(atom, single, IBond.Order.SINGLE));
						fragment.addAminoC(atom);
					}
				}
			}
		}
	}
	private void processDecalin() {
		IMolecule[] templates = new IMolecule[3];
		
		try {
			templates[0] = SmilesIO.readSmiles("CC1C=CC2CCCCC2C1C(O)=C");
			templates[1] = SmilesIO.readSmiles("CC1C=CC2CCCCC2\\C1=C(/C)O");
			templates[2] = SmilesIO.readSmiles("CC1C=CC2CCCCC2C1C(C)=O");
		} catch(Exception e) {
			e.printStackTrace();
		}
		IBond[] templateBondToBreak1 = new IBond[]{
				templates[0].getBond(templates[0].getAtom(1),templates[0].getAtom(10)),
				templates[1].getBond(templates[1].getAtom(1),templates[1].getAtom(10)),
				templates[2].getBond(templates[2].getAtom(1),templates[2].getAtom(10))
		};
		
		IBond[] templateBondToBreak2 = new IBond[]{
				templates[0].getBond(templates[0].getAtom(4),templates[0].getAtom(9)),
				templates[1].getBond(templates[1].getAtom(4),templates[1].getAtom(9)),
				templates[2].getBond(templates[2].getAtom(4),templates[2].getAtom(9))
		};
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IMolecule mol = fragment.getMolecule();
			List<IBond> templateBondMatches1 = new ArrayList<IBond>();
			int templateIndex = 0; //template index to check
			boolean match = false; //change to true if a match comes up in the next statement
			while(!match && templates.length > templateIndex){
				templateBondMatches1 = ChemicalUtilities.findMatchingBondsFromTemplate(templates[templateIndex], templateBondToBreak1[templateIndex], mol);
				if(templateBondMatches1.size() == 0){
					templateIndex ++;
				}else{
					match = true;
				}
			}
			if(!match){ //if no matches breaks the loop
				break;
			}
			List<IBond> templateBondMatches2 = ChemicalUtilities.findMatchingBondsFromTemplate(templates[templateIndex], templateBondToBreak2[templateIndex], mol);
			for(IBond bond : templateBondMatches1){
				mol.removeBond(bond);
			}
			for(IBond bond : templateBondMatches2){
				mol.removeBond(bond);
			}
		}
	}
	private void processBenzoquinone() {
		IMolecule benzoquinoneTemplate = null;
		
		try {
			benzoquinoneTemplate = SmilesIO.readSmiles("CC1=C(O)C(=O)C(C)=C(O)C1=O");
		} catch(Exception e) {
			e.printStackTrace();
		}
		IBond templateBondToBreak1 = benzoquinoneTemplate.getBond(benzoquinoneTemplate.getAtom(1),benzoquinoneTemplate.getAtom(2));
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IMolecule mol = fragment.getMolecule();
			List<IBond> templateBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(benzoquinoneTemplate, templateBondToBreak1, mol);
			for(IBond bond : templateBondMatches){
				mol.removeBond(bond);
				for(IAtom atom : bond.atoms()){
					for(IAtom connectedAtom : mol.getConnectedAtomsList(atom)){ // add oxygen to carboyxlic carbon
						if(connectedAtom.getAtomicNumber() == 8){
							Atom oxygen = new Atom("O");
							oxygen.setAtomTypeName("O.sp2");
							mol.addAtom(oxygen);
							mol.addBond(
									fragment.getMolecule().getAtomNumber(atom),
									fragment.getMolecule().getAtomNumber(oxygen),
									IBond.Order.DOUBLE);
							break;
						}
					}
				}
			}
			if(!ConnectivityChecker.isConnected(fragment.getMolecule())){
				IMoleculeSet splitFragments = ConnectivityChecker.partitionIntoMolecules(fragment.getMolecule());
				int indexToReplace = monomerFragments.indexOf(fragment);
				monomerFragments.remove(indexToReplace);
				Fragment firstFrag = new Fragment(splitFragments.getMolecule(0));
				monomerFragments.add(indexToReplace, firstFrag);
				for(int i = 1; splitFragments.getMoleculeCount() > i; i++){
					Fragment frag = new Fragment(splitFragments.getMolecule(i));
					monomerFragments.add(indexToReplace + 1, frag);
				}
			}
		}
	}
	private void processAminoDerivedTricyclicRing() { //Need to add an aditional breakdown here potentially
		IMolecule template1 = null;
		IMolecule template2 = null;
		try {
			template1 = SmilesIO.readSmiles("O=C1CN2C(=O)C3=C(C=CC=C3)N=C2CN1");
			template2 = SmilesIO.readSmiles("O=C1CN2C(=O)C3=CC=CC=C3N=C2CN1"); //Two templates are use because aromaticity is not functioning properly for matching
		} catch(Exception e) {
			e.printStackTrace();
		}
		IBond templateBondToBreak1 = template1.getBond(template1.getAtom(12),template1.getAtom(13));
		IBond templateBondToBreak2 = template2.getBond(template2.getAtom(12),template2.getAtom(13));
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IMolecule mol = fragment.getMolecule();
			List<IBond> templateBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template1, templateBondToBreak1, mol);
			if(templateBondMatches.size() == 0){
				templateBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template2, templateBondToBreak2, mol);
			}
			for(IBond bond : templateBondMatches){
				mol.removeBond(bond);
				for(IAtom atom : bond.atoms()){
					if(atom.getAtomicNumber() == 6){
						Atom oxygen = new Atom("O");
						oxygen.setAtomTypeName("O.sp2");
						mol.addAtom(oxygen);
						mol.addBond(
								fragment.getMolecule().getAtomNumber(atom),
								fragment.getMolecule().getAtomNumber(oxygen),
								IBond.Order.DOUBLE);
						
					}
				}
			}
			if(!ConnectivityChecker.isConnected(fragment.getMolecule())){
				IMoleculeSet splitFragments = ConnectivityChecker.partitionIntoMolecules(fragment.getMolecule());
				int indexToReplace = monomerFragments.indexOf(fragment);
				monomerFragments.remove(indexToReplace);
				Fragment firstFrag = new Fragment(splitFragments.getMolecule(0));
				monomerFragments.add(indexToReplace, firstFrag);
				for(int i = 1; splitFragments.getMoleculeCount() > i; i++){
					Fragment frag = new Fragment(splitFragments.getMolecule(i));
					monomerFragments.add(indexToReplace + 1, frag);
				}
			}
		}
	}
	private void process2Pyridone() {
		IMolecule template = null;
		try {
			template = SmilesIO.readSmiles("CC1=C(O)C(=CNC1=O)C1=CC=C(O)C=C1");
			template = SmilesIO.readSmiles(SmilesIO.generateSmiles(template));
		} catch(Exception e) {
			e.printStackTrace();
		}
		IBond templateBondToBreak = template.getBond(template.getAtom(12),template.getAtom(14));
		IBond templateBondToRemoveCarbon = template.getBond(template.getAtom(2),template.getAtom(3));
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IMolecule mol = fragment.getMolecule();
			List<IBond> templateBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateBondToBreak, mol);
			if(templateBondMatches.size() == 0){
				break;
			}
			List<IBond> templateBondMatches2 = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateBondToRemoveCarbon, mol);
			for(IBond bond : templateBondMatches){
				mol.removeBond(bond);
				for(IAtom atom : bond.atoms()){
					boolean carboxylicCarbon = false;
					for(IAtom atom2 : mol.getConnectedAtomsList(atom)){
						if(atom2.getAtomicNumber() == 6){
							carboxylicCarbon = true;
						}
					}
					if(carboxylicCarbon){
						Atom oxygen = new Atom("O");
						oxygen.setAtomTypeName("O.sp2");
						mol.addAtom(oxygen);
						mol.addBond(
								fragment.getMolecule().getAtomNumber(atom),
								fragment.getMolecule().getAtomNumber(oxygen),
								IBond.Order.DOUBLE);
						fragment.addAminoC(atom);
					}
				}
			}
			for(IBond bond : templateBondMatches2){ //Remove the carbon and reattach the nitrogen
				IAtom carbonToRemove = null;
				IAtom nitrogen = null;
				IAtom carbon = null; //carbon to attach to nitrogen
				for(IAtom atom : bond.atoms()){
					if(atom.getAtomicNumber() == 6){
						carbonToRemove = atom;
					}else{
						nitrogen = atom;
					}
				}
				for(IAtom atom : mol.getConnectedAtomsList(carbonToRemove)){
					if(atom.getAtomicNumber() == 6){
						carbon = atom;
						 break;
					}
				}
				mol.removeAtomAndConnectedElectronContainers(carbonToRemove);
				mol.addBond(
						fragment.getMolecule().getAtomNumber(carbon),
						fragment.getMolecule().getAtomNumber(nitrogen),
						IBond.Order.SINGLE);
				
			}
		}
		
	}
	private void processDielsAlder() {
		IMolecule dielsAlderTemplate1 = null;
		IMolecule dielsAlderTemplate2 = null;
		try {
			dielsAlderTemplate1 = SmilesIO.readSmiles("O=C1NCC2CC=CCC12");
			dielsAlderTemplate2 = SmilesIO.readSmiles("O=C1NCC2CC3OC3CC12");
		} catch(Exception e) {
			e.printStackTrace();
		}
		
		//TODO clean these up into arrays
		//no epoxide
		IBond templateBondToBreak1a = dielsAlderTemplate1.getBond(dielsAlderTemplate1.getAtom(9),dielsAlderTemplate1.getAtom(4));
		IBond templateBondToBreak1b = dielsAlderTemplate1.getBond(dielsAlderTemplate1.getAtom(4),dielsAlderTemplate1.getAtom(5));
		IBond templateBondToBreak1c = dielsAlderTemplate1.getBond(dielsAlderTemplate1.getAtom(5),dielsAlderTemplate1.getAtom(6));
		IBond templateBondToBreak1d = dielsAlderTemplate1.getBond(dielsAlderTemplate1.getAtom(6),dielsAlderTemplate1.getAtom(7));
		IBond templateBondToBreak1e = dielsAlderTemplate1.getBond(dielsAlderTemplate1.getAtom(7),dielsAlderTemplate1.getAtom(8));
		IBond templateBondToBreak1f = dielsAlderTemplate1.getBond(dielsAlderTemplate1.getAtom(8),dielsAlderTemplate1.getAtom(9));
		
		//epoxide
		IBond templateBondToBreak2a = dielsAlderTemplate2.getBond(dielsAlderTemplate2.getAtom(10),dielsAlderTemplate2.getAtom(4));
		IBond templateBondToBreak2b = dielsAlderTemplate2.getBond(dielsAlderTemplate2.getAtom(4),dielsAlderTemplate2.getAtom(5));
		IBond templateBondToBreak2c = dielsAlderTemplate2.getBond(dielsAlderTemplate2.getAtom(5),dielsAlderTemplate2.getAtom(6));
		IBond templateBondToBreak2d = dielsAlderTemplate2.getBond(dielsAlderTemplate2.getAtom(6),dielsAlderTemplate2.getAtom(8));
		IBond templateBondToBreak2e = dielsAlderTemplate2.getBond(dielsAlderTemplate2.getAtom(8),dielsAlderTemplate2.getAtom(9));
		IBond templateBondToBreak2f = dielsAlderTemplate2.getBond(dielsAlderTemplate2.getAtom(9),dielsAlderTemplate2.getAtom(10));
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IMolecule mol = fragment.getMolecule();
			List<IBond> templateBondMatches1 = ChemicalUtilities.findMatchingBondsFromTemplate(dielsAlderTemplate1, templateBondToBreak1a, mol);
			List<IBond> templateBondMatches2 = null;
			List<IBond> templateBondMatches3 = null;
			List<IBond> templateBondMatches4 = null;
			List<IBond> templateBondMatches5 = null;
			List<IBond> templateBondMatches6 = null;
			if(templateBondMatches1.size() == 0){
				templateBondMatches1 = ChemicalUtilities.findMatchingBondsFromTemplate(dielsAlderTemplate2, templateBondToBreak2a, mol);
				if(templateBondMatches1.size() == 0){
					break;
				}else{
					templateBondMatches2 = ChemicalUtilities.findMatchingBondsFromTemplate(dielsAlderTemplate2, templateBondToBreak2b, mol);
					templateBondMatches3 = ChemicalUtilities.findMatchingBondsFromTemplate(dielsAlderTemplate2, templateBondToBreak2c, mol);
					templateBondMatches4 = ChemicalUtilities.findMatchingBondsFromTemplate(dielsAlderTemplate2, templateBondToBreak2d, mol);
					templateBondMatches5 = ChemicalUtilities.findMatchingBondsFromTemplate(dielsAlderTemplate2, templateBondToBreak2e, mol);
					templateBondMatches6 = ChemicalUtilities.findMatchingBondsFromTemplate(dielsAlderTemplate2, templateBondToBreak2f, mol);
				}
			}else{
				templateBondMatches2 = ChemicalUtilities.findMatchingBondsFromTemplate(dielsAlderTemplate1, templateBondToBreak1b, mol);
				templateBondMatches3 = ChemicalUtilities.findMatchingBondsFromTemplate(dielsAlderTemplate1, templateBondToBreak1c, mol);
				templateBondMatches4 = ChemicalUtilities.findMatchingBondsFromTemplate(dielsAlderTemplate1, templateBondToBreak1d, mol);
				templateBondMatches5 = ChemicalUtilities.findMatchingBondsFromTemplate(dielsAlderTemplate1, templateBondToBreak1e, mol);
				templateBondMatches6 = ChemicalUtilities.findMatchingBondsFromTemplate(dielsAlderTemplate1, templateBondToBreak1f, mol);
			}
			for(IBond bond : templateBondMatches1){ //break
				mol.removeBond(bond);
			}
			for(IBond bond : templateBondMatches2){ //break
				mol.removeBond(bond);
			}
			for(IBond bond : templateBondMatches3){ //double
				bond.setOrder(Order.DOUBLE);
			}
			for(IBond bond : templateBondMatches4){ //single, remove epoxide
				bond.setOrder(Order.SINGLE);
				for(IAtom atom2 : mol.getConnectedAtomsList(bond.getAtom(0))){ // remove epoxide oxygen
					if(atom2.getAtomicNumber() == 8){
						mol.removeAtomAndConnectedElectronContainers(atom2);
						break;
					}
				}
			}
			for(IBond bond : templateBondMatches5){ //double
				bond.setOrder(Order.DOUBLE);
			}
			for(IBond bond : templateBondMatches6){ //break
				mol.removeBond(bond);
			}
			if(!ConnectivityChecker.isConnected(fragment.getMolecule())){
				IMoleculeSet splitFragments = ConnectivityChecker.partitionIntoMolecules(fragment.getMolecule());
				int indexToReplace = monomerFragments.indexOf(fragment);
				monomerFragments.remove(indexToReplace);
				Fragment firstFrag = new Fragment(splitFragments.getMolecule(0));
				monomerFragments.add(indexToReplace, firstFrag);
				for(int i = 1; splitFragments.getMoleculeCount() > i; i++){
					Fragment frag = new Fragment(splitFragments.getMolecule(i));
					monomerFragments.add(indexToReplace + 1, frag);
				}
			}
		}
	}
	
	private void processUniqueSubstructures() {
		processKendoLikeSubstructures();
		processPieriLikeSubstructures();
		processAnthramycinLikeSubstructures();
		breakHybridDiCystines();
	}

	private void processPeptideEpoxiKetones() {
		IMolecule peptideEpoxiKetoneTemplate = null;
				
		try {
			peptideEpoxiKetoneTemplate = SmilesIO.readSmiles("CC1(CO1)C(=O)CN");
		} catch(Exception e) {
			e.printStackTrace();
		}
		IBond templateBondToBreak = peptideEpoxiKetoneTemplate.getBond(peptideEpoxiKetoneTemplate.getAtom(1),peptideEpoxiKetoneTemplate.getAtom(4));
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IMolecule mol = fragment.getMolecule();
			List<IBond> templateBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(peptideEpoxiKetoneTemplate, templateBondToBreak, mol);
			for(IBond bond : templateBondMatches){
				mol.removeBond(bond);
				for(IAtom atom : bond.atoms()){
					if(mol.getConnectedAtomsCount(atom) == 2){
						
						Atom oxygen = new Atom("O");
						oxygen.setAtomTypeName("O.sp3");
						mol.addAtom(oxygen);
						mol.addBond(
								fragment.getMolecule().getAtomNumber(atom),
								fragment.getMolecule().getAtomNumber(oxygen),
								IBond.Order.SINGLE);
					}
				}
			}

			IMoleculeSet molFrags = ConnectivityChecker.partitionIntoMolecules(mol);
			if (molFrags.getAtomContainerCount() > 1){
				int indexToReplace = monomerFragments.indexOf(fragment);
				monomerFragments.remove(indexToReplace);
				IMolecule[] molPiece = new IMolecule[2];
				try {
					molPiece[0] = SmilesIO.readSmiles("CC1CO1");
					molPiece[1] = SmilesIO.readSmiles("OCC1CO1");
				} catch(Exception e) {
					e.printStackTrace();
				}
				 int addCount = 0;
				for(int i = 0; molFrags.getAtomContainerCount() > i; i++){
					IMolecule molFrag = molFrags.getMolecule(i);
					Fragment fragmentToAdd = null;
					if(ChemicalUtilities.getTanimotoScore(molFrag, molPiece[0]) >= 0.8 || ChemicalUtilities.getTanimotoScore(molFrag, molPiece[0]) >= 0.8){
						fragmentToAdd = new Fragment(molFrag);
						fragmentToAdd.addTailoringDomain(TailoringDomainEnums.C_METHYLTRANSFERASE);
						fragmentToAdd.setFragmentType(FragmentType.POLYKETIDE);
						List<PolyKetideDomainEnums> singlePolyketideDomains = new ArrayList<PolyKetideDomainEnums>();
						List<List<PolyKetideDomainEnums>> polyketideDomains = new ArrayList<List<PolyKetideDomainEnums>>();
						singlePolyketideDomains.add(PolyKetideDomainEnums.KETONE);
						polyketideDomains.add(singlePolyketideDomains);
						fragmentToAdd.setPkDomains(polyketideDomains);
						List<PKsubstrate> loadingUnits = new ArrayList<PKsubstrate>();
						loadingUnits.add(PKsubstrate.MALONYL);
						fragmentToAdd.setLoadingUnits(loadingUnits);
					}else{
						fragmentToAdd = new Fragment(molFrag);
					}
					monomerFragments.add(indexToReplace + addCount, fragmentToAdd);
					addCount += 1;
				}
			}
		}
	}

	private void processMonobactams() {
		
		IMolecule monobactamTemplate = null;
		
		try {
			monobactamTemplate = SmilesIO.readSmiles("O=C1CCN1");
		} catch(Exception e) {
			e.printStackTrace();
		}
		IBond templateBondToBreak = monobactamTemplate.getBond(monobactamTemplate.getAtom(3),monobactamTemplate.getAtom(4));
		IBond templateBondToBreakSulfate = monobactamTemplate.getBond(monobactamTemplate.getAtom(4),monobactamTemplate.getAtom(1));
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IMolecule mol = fragment.getMolecule();
			List<IBond> templateBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(monobactamTemplate, templateBondToBreak, mol);
			List<IBond> templateBondMatchesSulfate = ChemicalUtilities.findMatchingBondsFromTemplate(monobactamTemplate, templateBondToBreakSulfate, mol);
			if(templateBondMatches.size() == 0 && templateBondMatchesSulfate.size() == 0){
				continue;
			}
			boolean containsSulfate = false;
			for(IBond bond : templateBondMatchesSulfate){ //check if there is a sulfate
				for(IAtom atom : bond.atoms()){
					if(atom.getAtomicNumber() == 7){
						for(IAtom connectedAtom : mol.getConnectedAtomsList(atom)){
							if(connectedAtom.getAtomicNumber() == 16){
								containsSulfate = true;
							}
						}
					}
				}
			}
			if(containsSulfate){
				for(IBond bond : templateBondMatchesSulfate){
					for(IAtom atom : bond.atoms()){
						if(atom.getAtomicNumber() == 6){
							mol.removeBond(bond);
							Atom oxygen = new Atom("O");
							oxygen.setAtomTypeName("O.sp3");
							mol.addAtom(oxygen);
							mol.addBond(
									fragment.getMolecule().getAtomNumber(atom),
									fragment.getMolecule().getAtomNumber(oxygen),
									IBond.Order.SINGLE);
							fragment.addAminoC(atom);
						}else {
							fragment.addAminoN(atom);
						}
					}
				}
				continue;
			}
			for(IBond bond : templateBondMatches){
				mol.removeBond(bond);
				for(IAtom atom : bond.atoms()){
					if(atom.getAtomicNumber() == 6){
						Atom oxygen = new Atom("O");
						oxygen.setAtomTypeName("O.sp3");
						mol.addAtom(oxygen);
						mol.addBond(
								fragment.getMolecule().getAtomNumber(atom),
								fragment.getMolecule().getAtomNumber(oxygen),
								IBond.Order.SINGLE);
					}
				}
			}
		}
	}

	private void breakLinearEthers() {
		
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IMolecule mol = fragment.getMolecule();
			List<IAtom> cyclicAtoms = ChemicalUtilities.getAtomsInRings(mol);
			List<IAtom> linearEtherOxygens = new ArrayList<IAtom>(); 
			for(IAtom oxiCandidate : mol.atoms()){
				if(oxiCandidate.getAtomicNumber() != 8 
						|| mol.getConnectedAtomsCount(oxiCandidate) != 2 
						|| cyclicAtoms.contains(oxiCandidate)
						|| mol.getBondOrderSum(oxiCandidate) != 2
						) continue;
				boolean properEther = true;
				for(IAtom connectedC : mol.getConnectedAtomsList(oxiCandidate)){
					int numOxygens = 0;
					for(IAtom connectedAtom : mol.getConnectedAtomsList(connectedC)){
						if(connectedAtom.getAtomicNumber() == 8) numOxygens ++;
					}
					if(numOxygens > 1){
						properEther = false; 
						break; //Stop checking the atoms
					}
					if(connectedC.getAtomicNumber() != 6 //Must be connected to two carbonsw
					 || mol.getConnectedAtomsCount(connectedC) < 2){
						properEther = false; //Carbons must not be terminal
						break; //Stop checking the atoms
					}
				}
				if(properEther) linearEtherOxygens.add(oxiCandidate);
			}
			for(IAtom etherOxygen : linearEtherOxygens){
				List<IAtom>connectedCs = new ArrayList<IAtom>(); 
				for(IAtom connectedC : mol.getConnectedAtomsList(etherOxygen)){
					connectedCs.add(connectedC);
				}
				mol.removeAtomAndConnectedElectronContainers(etherOxygen);
				for(IAtom connectedC : connectedCs){
					Atom oxygen = new Atom("O");
					oxygen.setAtomTypeName("O.sp3");
					mol.addAtom(oxygen);
					mol.addBond(
							fragment.getMolecule().getAtomNumber(connectedC),
							fragment.getMolecule().getAtomNumber(oxygen),
							IBond.Order.SINGLE);
				}
			}
			if(!ConnectivityChecker.isConnected(fragment.getMolecule())){
				IMoleculeSet splitFragments = ConnectivityChecker.partitionIntoMolecules(fragment.getMolecule());
				int indexToReplace = monomerFragments.indexOf(fragment);
				monomerFragments.remove(indexToReplace);
				Fragment firstFrag = new Fragment(splitFragments.getMolecule(0));
				monomerFragments.add(indexToReplace, firstFrag);
				for(int i = 1; splitFragments.getMoleculeCount() > i; i++){
					Fragment frag = new Fragment(splitFragments.getMolecule(i));
					monomerFragments.add(indexToReplace + 1, frag);
				}
			}
		}
	}

	private void breakHybridDiCystines() {
		IMolecule template = null;
		try {template = SmilesIO.readSmiles("CC1=NC(=CS1)C1=NC(C)=CS1");} catch (Exception e) {}
		IBond templateToBreakBond = template.getBond(template.getAtom(8),template.getAtom(9));
		IMolecule template2 = null;
		try {template2 = SmilesIO.readSmiles("CC1CSC(N1)C1CSC(C)=N1");} catch (Exception e) {}
		IBond templateToBreakBond2 = template2.getBond(template2.getAtom(1),template2.getAtom(0));
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {

			Fragment fragment = q.poll();
			List<IBond> templateDoubleBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateToBreakBond, fragment.getMolecule());
			templateDoubleBondMatches.addAll(ChemicalUtilities.findMatchingBondsFromTemplate(template2, templateToBreakBond2, fragment.getMolecule()));
			if(templateDoubleBondMatches.size() < 1) continue;
			IMolecule mol = fragment.getMolecule();
			for(IBond bond : templateDoubleBondMatches){
				mol.removeBond(bond);
				for(IAtom atom : bond.atoms()){
					IMolecule carboxylicAcid = null;
					try {carboxylicAcid = SmilesIO.readSmiles("C(O)=O");} catch (Exception e) {}
					IAtom carboxilicC = carboxylicAcid.getAtom(0);
					mol.add(carboxylicAcid);
					mol.addBond(
							mol.getAtomNumber(atom),
							mol.getAtomNumber(carboxilicC),
							IBond.Order.SINGLE);
				}
			}
			if(!ConnectivityChecker.isConnected(fragment.getMolecule())){
				IMoleculeSet splitFragments = ConnectivityChecker.partitionIntoMolecules(fragment.getMolecule());
				int indexToReplace = monomerFragments.indexOf(fragment);
				// there must be two fragments
				Fragment firstFrag = new Fragment(splitFragments.getMolecule(0));
				Fragment secondFrag = new Fragment(splitFragments.getMolecule(1));
				monomerFragments.remove(indexToReplace);
				monomerFragments.add(indexToReplace, firstFrag);
				monomerFragments.add(indexToReplace + 1, secondFrag);
			}
		}
	}
	
	private void processAnthramycinLikeSubstructures() {
		IMolecule template = null;
		try {template = SmilesIO.readSmiles("O=C1NCC=NC=C1");} catch (Exception e) {}
		IBond templateToBreakBond = template.getBond(template.getAtom(4),template.getAtom(5));
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IMolecule mol = fragment.getMolecule();
			
			List<IBond> templateDoubleBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateToBreakBond, mol);
			
			for(IBond bond : templateDoubleBondMatches){
				bond.setOrder(Order.SINGLE);
				IAtom carbon = null;
				if(bond.getAtom(0).getAtomicNumber() == 6){
					carbon = bond.getAtom(0);
				}else{
					carbon = bond.getAtom(1);
				}
				Atom oxygen = new Atom("O");
				oxygen.setAtomTypeName("O.sp2");
				mol.addAtom(oxygen);
				mol.addBond(new Bond(carbon, oxygen, IBond.Order.DOUBLE));
			}
		}
	}

	private void processPieriLikeSubstructures() {
		IMolecule template = null;
		try {template = SmilesIO.readSmiles("CC1=C(C)C(=O)C(O)=C(O)N1");} catch (Exception e) {}
		IBond templateToBreakBond = template.getBond(template.getAtom(0),template.getAtom(1));
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {

			Fragment fragment = q.poll();
			IMolecule molClone = null;
			try {molClone = fragment.getMolecule().clone();}catch (CloneNotSupportedException e) {}
			List<IBond> templateDoubleBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateToBreakBond, molClone);
			
			if(templateDoubleBondMatches.size() < 1) continue;
			
			for(IBond bond : templateDoubleBondMatches){
				molClone.removeBond(bond);
				
				if(!ConnectivityChecker.isConnected(molClone)){
					IMoleculeSet splitFragments = ConnectivityChecker.partitionIntoMolecules(molClone);
					if(splitFragments.getMoleculeCount() == 2){
						if(splitFragments.getMolecule(0).getAtomCount() > splitFragments.getMolecule(1).getAtomCount()){
							molClone.addBond(bond);
							continue;
						}
						if(splitFragments.getMolecule(0).getAtomCount() > splitFragments.getMolecule(1).getAtomCount()){
							molClone = splitFragments.getMolecule(0);
						}else{
							molClone = splitFragments.getMolecule(1);
						}
						IAtom terminalCarbon = null;
						if(molClone.contains(bond.getAtom(0))){
							terminalCarbon = bond.getAtom(0);
						}else if(molClone.contains(bond.getAtom(1))){
							terminalCarbon = bond.getAtom(1);
						}else{
							System.out.println("Something went wrong with method processPieriLikeSubstructures -- This should never happen, this method was ignored for the final breakdown of this compound");
							break;
						}
						IMolecule pieceToAdd = null;
						try {pieceToAdd = SmilesIO.readSmiles("CC(C=O)C(=O)CC(O)=O");} catch (Exception e) {}
						IAtom atomToConnect = pieceToAdd.getAtom(2);
						molClone.add(pieceToAdd);
						molClone.addBond(molClone.getAtomNumber(terminalCarbon), molClone.getAtomNumber(atomToConnect), IBond.Order.SINGLE);
						fragment.setMolecule(molClone);
					}else{
						molClone.addBond(bond);
						continue;
					}					
				}else{ //Shouldn't be connected
					molClone.addBond(bond);
					continue;
				}	
			}
			
		}
		
	}

	private void processKendoLikeSubstructures() {

		IMolecule template = null;
		try {template = SmilesIO.readSmiles("CC1=C2OC(C)(O)C=C2C(C)=C(O)C1=O");} catch (Exception e) {}
		IBond templateDoubleBond = template.getBond(template.getAtom(7),template.getAtom(8));
		IBond templateEtherBond = template.getBond(template.getAtom(3),template.getAtom(4));
		IBond templateHydroxylBond = template.getBond(template.getAtom(4),template.getAtom(6));
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		// go through each nrp fragment. In current implementation there should
		// be only one.
		while (!q.isEmpty()) {

			Fragment fragment = q.poll();
			List<IBond> templateDoubleBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateDoubleBond, fragment.getMolecule());
			
			if(templateDoubleBondMatches.size() < 1) continue;
			
			List<IBond> templateEtherBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateEtherBond, fragment.getMolecule());
			List<IBond> templateHydroxylBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateHydroxylBond, fragment.getMolecule());
			
			for(IBond bond : templateDoubleBondMatches){
				fragment.getMolecule().removeBond(bond);
				for(IAtom atom : bond.atoms()){
					if(fragment.getMolecule().getConnectedAtomsCount(atom) == 2){ //Add ketone to carbon in the starter unit
						Atom ketone = new Atom("O");
						ketone.setAtomTypeName("O.sp2");
						fragment.getMolecule().addAtom(ketone);
						fragment.getMolecule().addBond(
								fragment.getMolecule().getAtomNumber(ketone),
								fragment.getMolecule().getAtomNumber(atom),
								IBond.Order.DOUBLE);
					}else if(fragment.getMolecule().getConnectedAtomsCount(atom) == 1){ //Add carboxcylic acid to carbon not in starter unit
						IMolecule carboxylicAcid = null;
						try {carboxylicAcid = SmilesIO.readSmiles("C(O)=O");} catch (Exception e) {}
						IAtom carboxilicC = carboxylicAcid.getAtom(0);
						fragment.getMolecule().add(carboxylicAcid);
						fragment.getMolecule().addBond(
								fragment.getMolecule().getAtomNumber(atom),
								fragment.getMolecule().getAtomNumber(carboxilicC),
								IBond.Order.SINGLE);
						
					}else{
						System.out.println("Something went wrong with method processKendoLikeSubstructures (Double bond match doesn't have appropriate carbons): reattached broken bond"); //TODO make exception
						fragment.getMolecule().addBond(bond);
					}				
				}		
			}
			for(IBond bond : templateEtherBondMatches){
				fragment.getMolecule().removeBond(bond);
				if(!ConnectivityChecker.isConnected(fragment.getMolecule())){
					fragment.getMolecule().addBond(bond);
					System.out.println("Something went wrong with method processKendoLikeSubstructures (Breaking ether makes molecule not connected): reattached broken bond");
				}
			}
			for(IBond bond : templateHydroxylBondMatches){
				bond.setOrder(IBond.Order.DOUBLE);
			}
		}
	}

	private void modifyImine() {
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		// go through each nrp fragment. In current implementation there should
		// be only one.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();

			for (int j = 0; j < fragment.getMolecule().getAtomCount(); j++) {
				if (!(fragment.getMolecule().getAtom(j).getAtomTypeName()
						.equals("N.sp2"))) {
					continue;
				}
				IAtom imineN = fragment.getMolecule().getAtom(j);
				List<Set<IAtom>> rings = ChemicalUtilities.getSmallestRings(fragment.getMolecule());
				boolean inSmallRing = false;
				for(Set<IAtom> ringAtoms : rings){
					if(ringAtoms.contains(imineN)){
						inSmallRing = true;
						break;
					}
				}
				if(inSmallRing) continue;
				
				IAtom imineC = null;
				for (IAtom c : fragment.getMolecule().getConnectedAtomsList(imineN)) {
					if (c.getAtomTypeName().equals("C.sp2")) {

						if (fragment.getMolecule()
										.getBond(imineN, c)
										.getOrder() == IBond.Order.DOUBLE
										&& fragment.getMolecule().getConnectedAtomsList(c).size() == 2){
							imineC = c;
							
						}
					}
				}
				if(imineC == null) continue;
				// Then check if breaking the imineN and imineC leaves
				// one piece.
				IBond removedBond = fragment.getMolecule().removeBond(
						imineN, imineC);

				if (ConnectivityChecker.isConnected(fragment.getMolecule())) {
					//Do the removal of double bond and create single
					fragment.getMolecule().addBond(
							fragment.getMolecule().getAtomNumber(imineN),
							fragment.getMolecule().getAtomNumber(imineC),
							IBond.Order.SINGLE);
					Atom ketone = new Atom("O");
					ketone.setAtomTypeName("O.sp2");
					fragment.getMolecule().addAtom(ketone);
					fragment.getMolecule().addBond(
							fragment.getMolecule().getAtomNumber(ketone),
							fragment.getMolecule().getAtomNumber(imineC),
							IBond.Order.DOUBLE);
					
				}else{
					//re-add the original bond
					fragment.getMolecule().addBond(removedBond);
				}
			}
		}
		
	}

	/**
	 * Break peptide bonds. This method modifies monomerFragments and should be
	 * run after thiazole, oxazole, and lactone processing.
	 * 
	 * @throws CDKException
	 */
	private void breakPeptideBonds() {
		// The resulting arraylist is ordered C to N
		// For each nrp fragment, identify all peptide bonds/beta
		// carbons/nitrogens. Then add these to the queue.

		int numStandardPeptideBonds = 0;
		int numNonStandardPeptideBonds = 0;
		ArrayList<IBond> peptideBonds = new ArrayList<IBond>();

		List<Set<IAtom>> smallRings = new ArrayList<Set<IAtom>>();
		List<IAtom> atomsInRings = new ArrayList<IAtom>();
		
		for (Fragment frag : monomerFragments) {
			IMolecule m = frag.getMolecule();
			
			// Keep track of rings with six or fewer atoms. If the C-N or C=N bond participates in this ring, break it but do not
			// change the connectivity annotation
			smallRings.addAll(ChemicalUtilities.getSmallestRings(frag.getMolecule()));
			
			// Keep track of atoms in rings of any size
			atomsInRings.addAll(ChemicalUtilities.getAtomsInRings(m));
			
			// Prepare peptide template 1 (amino acids with C terminus) and
			// peptide template 2 (amino acids with a nitrogen attached to C
			// terminus)
			IMolecule[] peptideTemplates = new IMolecule[5];
			IBond[] peptideTemplateBonds = new IBond[5];
			try {
				peptideTemplates[0] = SmilesIO.readSmiles("CCNC(C)=O");
				peptideTemplates[1] = SmilesIO.readSmiles("CC(=O)NC=C");
				peptideTemplates[2] = SmilesIO.readSmiles("CC\\N=C(/C)O");
				peptideTemplates[3] = SmilesIO.readSmiles("CCNC(C)OC");
				peptideTemplates[4] = SmilesIO.readSmiles("C\\C(O)=N/C=C");
				//peptideTemplates[5] = SmilesIO.readSmiles("CC(=O)NC=C");
				// peptideTemplates.add(SmilesIO.readSmiles("CC(=O)NCC=O"));
				// peptideTemplates.add(SmilesIO.readSmiles("CC(O)NCC=O"));
				// peptideTemplates.add(SmilesIO.readSmiles("CC(=O)NCCO"));
				// peptideTemplates.add(SmilesIO.readSmiles("CC(O)NCCO"));
				// peptideTemplates.add(SmilesIO.readSmiles("C\\C(O)=N\\CCO"));
				// peptideTemplates.add(SmilesIO.readSmiles("C\\C(O)=N\\CC=O"));
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			// Find the template bond
			for (int i = 0; i < peptideTemplates.length; i++) {
				IMolecule peptideTemplate = peptideTemplates[i];
				for (int j = 0; j < peptideTemplate.getBondCount(); j++) {
					if (peptideTemplate.getBond(j).getAtom(0).getAtomicNumber() == 6
							&& peptideTemplate.getBond(j).getAtom(1)
									.getAtomicNumber() == 7
							|| peptideTemplate.getBond(j).getAtom(0)
									.getAtomicNumber() == 7
							&& peptideTemplate.getBond(j).getAtom(1)
									.getAtomicNumber() == 6) {
						IAtom currentCarbon = null;
						IAtom currentNitrogen = null;
						if (peptideTemplate.getBond(j).getAtom(0)
								.getAtomicNumber() == 6) {
							currentCarbon = peptideTemplate.getBond(j).getAtom(
									0);
							currentNitrogen = peptideTemplate.getBond(j)
									.getAtom(1);
						} else {
							currentCarbon = peptideTemplate.getBond(j).getAtom(
									1);
							currentNitrogen = peptideTemplate.getBond(j)
									.getAtom(0);
						}
						boolean carbonNextToOxygen = false;
						for (IAtom a : peptideTemplate
								.getConnectedAtomsList(currentCarbon)) {
							if (a.getAtomicNumber() == 8) {
								carbonNextToOxygen = true;
							}
						}
						if (!carbonNextToOxygen) {
							continue;
						}
						peptideTemplateBonds[i] = (peptideTemplate.getBond(j));
					}
				}
			}
			
			for (int x = 0; x < peptideTemplates.length; x++) {
				IMolecule peptideTemplate = peptideTemplates[x];
				IBond peptideTemplateBond = peptideTemplateBonds[x];
				List<List<RMap>> templateMatchMap = null;
				try {
					templateMatchMap = UniversalIsomorphismTester
							.getSubgraphMaps(m, peptideTemplate);
				} catch (CDKException e) {
					e.printStackTrace();
				}
				for (int i = 0; i < templateMatchMap.size(); i++) {
					// None of the carbons attached to another carbon can be
					// terminal.
					boolean foundWrongTerminalCarbon = false;
					for (int j = 0; j < templateMatchMap.get(i).size(); j++) {
						IBond currentTemplateMatchBond = peptideTemplate
								.getBond(templateMatchMap.get(i).get(j)
										.getId2());
						for (int k = 0; k < 2; k++) {
							// For each bond in the molecule structure
							IAtom a = m.getBond(
									templateMatchMap.get(i).get(j).getId1())
									.getAtom(k);
							if (m.getConnectedAtomsCount(a) == 1
									&& a.getAtomicNumber() == 6) {
								for (IAtom connectedAtom : m
										.getConnectedAtomsList(a)) {
									if (connectedAtom.getAtomicNumber() == 6) {
										foundWrongTerminalCarbon = true;
									}
								}
							}
						}
					}
					if (foundWrongTerminalCarbon) {
						continue;
					}
					for (int j = 0; j < templateMatchMap.get(i).size(); j++) {
						IBond currentTemplateMatchBond = peptideTemplate
								.getBond(templateMatchMap.get(i).get(j)
										.getId2());
						if (currentTemplateMatchBond == peptideTemplateBond) {
							IBond candidatePeptideBond = m
									.getBond(templateMatchMap.get(i).get(j)
											.getId1());

							IAtom currentCarbon = null;
							IAtom currentNitrogen = null;
							if (candidatePeptideBond.getAtom(0)
									.getAtomicNumber() == 6) {
								currentCarbon = candidatePeptideBond.getAtom(0);
								currentNitrogen = candidatePeptideBond
										.getAtom(1);
							} else {
								currentCarbon = candidatePeptideBond.getAtom(1);
								currentNitrogen = candidatePeptideBond
										.getAtom(0);
							}
							// Check that this nitrogen has not been assigned to
							// another peptide bond already

							boolean nitrogenAlreadyPresent = false;
							for (IBond prevBond : peptideBonds) {
								if (prevBond.contains(currentNitrogen)) {
									nitrogenAlreadyPresent = true;
									break;
								}
							}
							if (nitrogenAlreadyPresent && frag.getMolecule().getConnectedAtomsCount(currentNitrogen) < 3) {
								break;
							}

							if (x == 0) {
								numStandardPeptideBonds++;
							} else {
								numNonStandardPeptideBonds++;
							}
							
							peptideBonds.add(candidatePeptideBond);
						}
					}
				}
			}

			if (showBonds) {
				if (peptideBonds.isEmpty()) {
					SmilesIO.drawMolecule(m, GrapeMain.currentName
							+ "peptide_bonds");
				} else {
					SmilesIO.drawMoleculeHighlightingBonds(
							m,
							SmilesIO.getCleanFileName(GrapeMain.currentName
									+ "peptide_bonds"), peptideBonds);
				}
			}
		}
		
		// Swap the order of peptide bonds so that bonds with the carbon in a cylcle are processed first.
		// This takes precedence over bonds involved in amino acid side chains 
		
		int lastUnvisitedIndex = peptideBonds.size()-1; // All elements from lastUnvisitedIndex onward have a cyclic carbon
		for(int i = 0; i < lastUnvisitedIndex; i++) {
			IAtom peptideC = null;
			if(peptideBonds.get(i).getAtom(0).getAtomicNumber() == 6) {
				peptideC = peptideBonds.get(i).getAtom(0);
			} else {
				peptideC = peptideBonds.get(i).getAtom(1);
			}
			if(!atomsInRings.contains(peptideC)) {
				Collections.swap(peptideBonds, i, lastUnvisitedIndex);
				lastUnvisitedIndex--;
				i--;
			}
		}
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		while (!q.isEmpty()) {
			Fragment subfragment = q.poll();
			// get peptide bond index or move on if there are no peptide bonds
			int peptideBondIndex = -1;
			int bondCount = subfragment.getMolecule().getBondCount();
			for (int i = 0; i < bondCount; i++) {
				if (peptideBonds.contains(subfragment.getMolecule().getBond(i))) {
					peptideBondIndex = i;
					break;
				}
			}
			if (peptideBondIndex == -1) {
				continue;
			}
			IBond peptideBond = subfragment.getMolecule().getBond(
					peptideBondIndex);
			IAtom betaCarbon = null;
			IAtom backboneNitrogen = null;
			if (peptideBond.getAtom(0).getAtomTypeName().startsWith("N.")) {
				backboneNitrogen = peptideBond.getAtom(0);
				betaCarbon = peptideBond.getAtom(1);
			} else {
				backboneNitrogen = peptideBond.getAtom(1);
				betaCarbon = peptideBond.getAtom(0);
			}
			// For the double bond isomer, set the oxygen back to sp2 and double
			// bonded.
			if (peptideBond.getOrder() == IBond.Order.DOUBLE) {
				List<IAtom> connectedBetaAtoms = subfragment.getMolecule()
						.getConnectedAtomsList(betaCarbon);
				for (IAtom a : connectedBetaAtoms) {
					if (a.getAtomTypeName().startsWith("O.")) {
						// This is the carbonyl oxygen
						a.setAtomTypeName("O.sp2");
						subfragment.getMolecule().getBond(betaCarbon, a)
								.setOrder(IBond.Order.DOUBLE);
					}
				}
			}

			subfragment.getMolecule().removeBond(peptideBondIndex);

			// Add -O for the OH group on the carbon

			Atom hydroxideO = new Atom("O");
			hydroxideO.setAtomTypeName("O.sp3");
			subfragment.getMolecule().addAtom(hydroxideO);
			subfragment.getMolecule().addBond(new Bond(betaCarbon, hydroxideO));

			ArrayList<Fragment> fragments = subfragment
					.partitionIntoMonomerFragments();
			if (fragments.size() == 1) {
				boolean bondParticipatesInSmallRing = false;
				for(Set<IAtom> atomsInRing : smallRings) {
					if(atomsInRing.contains(peptideBond.getAtom(0)) &&
							atomsInRing.contains(peptideBond.getAtom(1))){
						bondParticipatesInSmallRing = true;
					}
				}
				if(!bondParticipatesInSmallRing) {
					subfragment.addAminoC(betaCarbon);
					subfragment.addAminoN(backboneNitrogen);
					subfragment.setAtomAfterCTerminus(backboneNitrogen);
					subfragment.setAtomAfterNTerminus(betaCarbon);
				}
				q.add(subfragment);
			} else {
				// there must be two fragments
				Fragment CFragment = null, NFragment = null;
				if (fragments.get(0).getMolecule().contains(betaCarbon)) {
					CFragment = fragments.get(0);
					NFragment = fragments.get(1);
				} else {
					CFragment = fragments.get(1);
					NFragment = fragments.get(0);
				}
				

				CFragment.addAminoC(betaCarbon);
				NFragment.addAminoN(backboneNitrogen);
				
				// If this peptide connectivity does not conflict with previously assigned connections
				if(CFragment.getAtomAfterCTerminus() == null &&
						NFragment.getAtomAfterNTerminus() == null) {	
					CFragment.setAtomAfterCTerminus(backboneNitrogen);
					NFragment.setAtomAfterNTerminus(betaCarbon);
				}
				
				for(IAtom aminoNtoAdd : subfragment.getAminoNs()){
					if(NFragment.getMolecule().contains(aminoNtoAdd)){
						if(!NFragment.getAminoNs().contains(aminoNtoAdd)) NFragment.addAminoN(aminoNtoAdd);
					}
				}

				int indexToReplace = monomerFragments.indexOf(subfragment);

				monomerFragments.remove(indexToReplace);
				monomerFragments.add(indexToReplace, CFragment);
				monomerFragments.add(indexToReplace + 1, NFragment);
				q.add(NFragment);
				q.add(CFragment);
			}

		}
		if (printSteps) {
			if (numStandardPeptideBonds > 0) {
				System.out.println("Processed " + numStandardPeptideBonds
						+ " standard peptide bonds");
			}
			if (numNonStandardPeptideBonds > 0) {
				System.out.println("Processed " + numNonStandardPeptideBonds
						+ " nonstandard peptide bonds");
			}
		}
	}

	/**
	 * Open lactone rings, labelling Fragment adjacent lactone fragments
	 * when possible. This method modifies monomerFragments and also updates the
	 * field lactoneCarboxylCList.
	 */
	private void openLactoneRings() {
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		int numLactoneRings = 0;
		// go through each nrp fragment. In current implementation there should
		// be only one.
		while (!q.isEmpty()) {

			Fragment fragment = q.poll();
			ArrayList<Fragment> subfragments = new ArrayList<Fragment>();
			subfragments.add(fragment);
			// Find all cyclic esters. First find the oxygen in the ring, check
			// that it has two carbon neighbours, check that one of the carbon
			// neighbours
			// is double bonded to another O, then check that cleaving the C-O
			// bond gives a connected product

			ArrayList<IAtom> hydroxylOList = new ArrayList<IAtom>();
			ArrayList<IAtom> hydroxylCList = new ArrayList<IAtom>();
			ArrayList<IAtom> carboxylCList = new ArrayList<IAtom>();

			for (int j = 0; j < fragment.getMolecule().getAtomCount(); j++) {

				if (!(fragment.getMolecule().getAtom(j).getAtomTypeName()
						.equals("O.sp3"))) {
					continue;
				}

				IAtom hydroxyl_O = fragment.getMolecule().getAtom(j);
				// check if it is linked to two carbons, one of which has a
				// double bond O.

				IAtom carboxylC = null;
				IAtom hydroxylC = null;

				for (IAtom c : fragment.getMolecule().getConnectedAtomsList(
						hydroxyl_O)) {
					// the carboxyl carbon must be sp2 hybridized and connected
					// to another oxygen via double bond
					if (c.getAtomTypeName().equals("C.sp2")) {
						for (IAtom candidateOtherO : fragment.getMolecule()
								.getConnectedAtomsList(c)) {
							if (candidateOtherO != hydroxyl_O
									&& fragment.getMolecule()
											.getBond(candidateOtherO, c)
											.getOrder() == IBond.Order.DOUBLE
									&& candidateOtherO.getAtomicNumber() == 8) {
								carboxylC = c;
							}
						}
					} else if (c.getAtomTypeName().equals("C.sp3")) {
						hydroxylC = c;
					}
				}

				if (carboxylC == null || hydroxylC == null) {
					continue;
				}

				// Then check if breaking the hydroxyl_O and carboxylC leaves
				// one piece.
				IBond removedBond = fragment.getMolecule().removeBond(
						hydroxyl_O, carboxylC);

				if (ConnectivityChecker.isConnected(fragment.getMolecule())) {
					hydroxylOList.add(hydroxyl_O);
					hydroxylCList.add(hydroxylC);
					carboxylCList.add(carboxylC);
					// Add this carboxyl C to the global list
					lactoneCarboxylCList.add(carboxylC);
				}
				fragment.getMolecule().addBond(removedBond);
			}

			if (showBonds) {
				// For testing: draw the molecule with detected bonds
				// highlighted
				ArrayList<IBond> bonds = new ArrayList<IBond>();
				for (int j = 0; j < hydroxylOList.size(); j++) {
					bonds.add(fragment.getMolecule().getBond(
							hydroxylOList.get(j), carboxylCList.get(j)));
				}
				SmilesIO.drawMoleculeHighlightingBonds(
						fragment.getMolecule(),
						SmilesIO.getCleanFileName(GrapeMain.currentName
								.replace(' ', '_') + "Lactone"), bonds);
			}

			// Break the cyclic esters

			for (int j = 0; j < hydroxylOList.size(); j++) {
				numLactoneRings++;
				int index = -1;
				Fragment currentFragment = null;
				for (int k = 0; k < subfragments.size(); k++) {
					if (subfragments.get(k).getMolecule()
							.contains(hydroxylOList.get(j))) {
						if (!subfragments.get(k).getMolecule()
								.contains(carboxylCList.get(j))) {
							// TODO: proper exception handling
							System.out
									.println("error: fragment does not contain both lactone O and C");
						}
						index = k;
						break;
					}
				}

				currentFragment = subfragments.get(index);
				IAtom currentHydroxylO = hydroxylOList.get(j);
				IAtom currentHydroxylC = hydroxylCList.get(j);
				IAtom currentCarboxylC = carboxylCList.get(j);
				currentFragment.getMolecule().removeBond(currentHydroxylO,
						currentCarboxylC);
				Atom atom = new Atom("O");
				atom.setAtomTypeName("O.sp3");
				currentFragment.getMolecule().addAtom(atom);
				currentFragment.getMolecule().addBond(
						new Bond(carboxylCList.get(j), atom));

				ArrayList<Fragment> partitions = currentFragment
						.partitionIntoMonomerFragments();
				if (partitions.size() == 1) {
					currentFragment.setLactoneHydroxylO(currentHydroxylO);
					currentFragment.setLactoneHydroxylC(currentHydroxylC);
					currentFragment.setLactoneCarboxylC(currentCarboxylC);
					currentFragment
							.setAtomAfterLactoneHydroxyl(currentCarboxylC);
					currentFragment
							.setAtomAfterLactoneCarboxyl(currentHydroxylO);
				} else {
					// there must be two partitions
					Fragment carboxylPartiton = null;
					Fragment hydroxylPartition = null;

					if (partitions.get(0).getMolecule()
							.contains(carboxylCList.get(j))) {
						carboxylPartiton = partitions.get(0);
						hydroxylPartition = partitions.get(1);
					} else {
						carboxylPartiton = partitions.get(1);
						hydroxylPartition = partitions.get(0);
					}

					carboxylPartiton.setLactoneCarboxylC(currentCarboxylC);
					hydroxylPartition.setLactoneHydroxylO(currentHydroxylO);
					hydroxylPartition.setLactoneHydroxylC(currentHydroxylC);
					carboxylPartiton
							.setAtomAfterLactoneCarboxyl(currentHydroxylO);
					hydroxylPartition
							.setAtomAfterLactoneHydroxyl(currentCarboxylC);

					subfragments.remove(index);
					subfragments.add(index, carboxylPartiton);
					subfragments.add(index + 1, hydroxylPartition);
				}
			}
			if (subfragments.size() > 1) {
				int indexToReplace = monomerFragments.indexOf(fragment);
				monomerFragments.remove(indexToReplace);
				for (int i = 0; i < subfragments.size(); i++) {
					monomerFragments.add(indexToReplace + i,
							subfragments.get(i));
				}
			}
		}
		if (printSteps) {
			if (numLactoneRings > 0) {
				System.out.println("Processed " + numLactoneRings
						+ " lactone rings");
			}
		}
	}

	/**
	 * Process thiazoles. Searches for thiazole rings and 'undoes' the reaction.
	 * This method modifies monomerFragments and updates the field
	 * thiazoleSList.
	 */
	private void processThiazoles() {
		// First, handle the situation corresponding to the five-membered S and
		// N containing ring
		// as in: CNC(=O)C1=CSC(CN)=N1
		boolean hasThiazole = false;
		for (Fragment m : monomerFragments) {
			IMolecule molecule = m.getMolecule();
			ArrayList<IAtom> carbonylCAtoms = new ArrayList<IAtom>();
			ArrayList<IAtom> thiazoneNAtoms = new ArrayList<IAtom>();
			ArrayList<IAtom> chainC1Atoms = new ArrayList<IAtom>();
			ArrayList<IAtom> chainC2Atoms = new ArrayList<IAtom>();
			ArrayList<IAtom> thiazoleSAtoms = new ArrayList<IAtom>();
			for (int i = 0; i < m.getMolecule().getAtomCount(); i++) {
				if (molecule.getAtom(i).getAtomicNumber() != 16) {
					continue;
				}
				// check if this sulfur is connected to a carbon (C) double
				// bonded with a nitrogen
				IAtom thiazoleS = null;
				IAtom carbonylC = null;
				IAtom thiazoleN = null;
				IAtom chainC_1 = null;
				IAtom chainC_2 = null;

				for (IAtom c : molecule.getConnectedAtomsList(molecule
						.getAtom(i))) {
					if (molecule.getBond(molecule.getAtom(i), c).getOrder() != IBond.Order.SINGLE) {
						continue;
					}
					if (!c.getAtomTypeName().equals("C.sp2")) {
						continue;
					}
					for (IAtom n : molecule.getConnectedAtomsList(c)) {
						if (n.getAtomTypeName().startsWith("C.")
								&& molecule.getBond(c, n).getOrder() == IBond.Order.DOUBLE) {
							chainC_1 = c;
							chainC_2 = n;
							continue;
						}
						if (n.getAtomTypeName().startsWith("N.")
								&& molecule.getBond(c, n).getOrder() == IBond.Order.DOUBLE) {
							carbonylC = c;
							thiazoleN = n;
						}
					}
				}
				if (carbonylC == null) {
					continue;
				}

				if (chainC_1 != null && chainC_2 != null) {
					molecule.getBond(chainC_1, chainC_2).setOrder(
							IBond.Order.SINGLE);
					chainC_1.setAtomTypeName("C.sp3");
					chainC_2.setAtomTypeName("C.sp3");
				}

				thiazoleS = molecule.getAtom(i);

				carbonylCAtoms.add(carbonylC);
				thiazoneNAtoms.add(thiazoleN);
				chainC1Atoms.add(chainC_1);
				chainC2Atoms.add(chainC_2);
				thiazoleSAtoms.add(thiazoleS);
			}
			for (int i = 0; i < thiazoleSAtoms.size(); i++) {
				IAtom thiazoleS = thiazoleSAtoms.get(i);
				IAtom carbonylC = carbonylCAtoms.get(i);
				IAtom thiazoleN = thiazoneNAtoms.get(i);
				IAtom chainC_1 = chainC1Atoms.get(i);
				IAtom chainC_2 = chainC2Atoms.get(i);

				// Cleave the C-S bond. If this alters connectivity, then it
				// wasn't part of a cycle so re-add it.
				IBond removedBond = molecule.removeBond(carbonylC, thiazoleS);
				if (!ConnectivityChecker.isConnected(molecule)) {
					molecule.addBond(removedBond);
					continue;
				}

				// Turn the C=N bond to C-N
				molecule.getBond(carbonylC, thiazoleN).setOrder(
						IBond.Order.SINGLE);

				// Add a =0 to that C
				Atom atom = new Atom("O");
				atom.setAtomTypeName("O.sp2");
				molecule.addAtom(atom);
				molecule.addBond(new Bond(carbonylC, atom, IBond.Order.DOUBLE));

				// Set carbonylC to sp2

				carbonylC.setAtomTypeName("C.sp2");
				thiazoleN.setAtomTypeName("N.amide");

				thiazoleSList.add(thiazoleS);
				hasThiazole = true;
			}
		}
		if (printSteps) {
			if (hasThiazole) {
				System.out.println("Processed thiazole");
			}
		}

		// Next, generally look for any cysteine whose sulfur is connected to a
		// carbon, with that S-C bond a part of a cycle.
		boolean hasCyclizedNonThiazoleCysteine = false;
		IMolecule template = null;
		try {
			template = SmilesIO.readSmiles("CSCC(N)C=O");
		} catch (IOException | CDKException e) {
			e.printStackTrace();
		}
		// Find the template bond
		IBond templateBond = null;
		for (int i = 0; i < template.getBondCount(); i++) {
			IBond bond = template.getBond(i);
			if (bond.getAtom(0).getAtomicNumber() == 16
					|| bond.getAtom(1).getAtomicNumber() == 16
					&& template.getConnectedAtomsCount(bond.getAtom(0)) == 1
					|| template.getConnectedAtomsCount(bond.getAtom(1)) == 1) {
				templateBond = bond;
			}
		}

		for (Fragment frag : monomerFragments) {
			List<IBond> matchingBonds = ChemicalUtilities.findMatchingBondsFromTemplate(template,
					templateBond, frag.getMolecule());
			for (IBond matchingBond : matchingBonds) {
				IAtom thiazoleS = null;
				IAtom carbonylC = null;
				if (matchingBond.getAtom(0).getAtomicNumber() == 16) {
					thiazoleS = matchingBond.getAtom(0);
					carbonylC = matchingBond.getAtom(1);
				} else {
					carbonylC = matchingBond.getAtom(0);
					thiazoleS = matchingBond.getAtom(1);
				}
				frag.getMolecule().removeBond(matchingBond);
				if (ConnectivityChecker.isConnected(frag.getMolecule())) {
					// Set all bonds connecting to carbonylC as single
					for (int i = 0; i < frag.getMolecule().getBondCount(); i++) {
						IBond b = frag.getMolecule().getBond(i);
						if (b.contains(carbonylC)) {
							b.setOrder(IBond.Order.SINGLE);
						}
					}
					Atom atom = new Atom("O");
					atom.setAtomTypeName("O.sp2");
					frag.getMolecule().addAtom(atom);
					frag.getMolecule().addBond(
							new Bond(carbonylC, atom, IBond.Order.DOUBLE));
					thiazoleSList.add(thiazoleS);
					hasCyclizedNonThiazoleCysteine = true;
				} else {
					frag.getMolecule().addBond(matchingBond);
				}
			}
		}
		if (printSteps) {
			if (hasCyclizedNonThiazoleCysteine) {
				System.out.println("Processed non-thiazole cyclized cysteine");
			}
		}
	}

	/**
	 * Process oxazoles. Searches for oxazole rings and 'undoes' the reaction.
	 * This method modifies monomerFragments and updates the field oxazoleOList.
	 */
	private void processOxazoles() {
		boolean hasOxazole = false;
		for (Fragment m : monomerFragments) {
			IMolecule molecule = m.getMolecule();
			IMolecule[] oxazoleTemplates = new IMolecule[2];
			
			try {
				oxazoleTemplates[0] = SmilesIO.readSmiles("O1C=CN=C1");
				oxazoleTemplates[1] = SmilesIO.readSmiles("C1CN=CO1");
			} catch(Exception e) {
				e.printStackTrace();
			}
			
			// ArrayList of template bonds corresponding to each oxazole template.
			//The first bond is the double bond with nitrogen, the second bond is the carbon-carbon bond, and the third bond is the oxygen-carbon bond which is broken.
			List<List<IBond>> templateBonds = new ArrayList<List<IBond>>();
			for(int i = 0; i < oxazoleTemplates.length; i++) {
				// Initialize the size of the template bonds
				templateBonds.add(new ArrayList<IBond>());
				for(int j = 0; j < 3; j++) {
					templateBonds.get(i).add(null);
				}
				for(int j = 0; j < oxazoleTemplates[i].getBondCount(); j++) {
					IBond currentBond = oxazoleTemplates[i].getBond(j);
					if(currentBond.getOrder() == IBond.Order.DOUBLE &&
						(currentBond.getAtom(0).getAtomicNumber() == 7 || currentBond.getAtom(1).getAtomicNumber() == 7)) {
						templateBonds.get(i).set(0, currentBond);
						continue;
					}
					if(currentBond.getAtom(0).getAtomicNumber() == 6 && currentBond.getAtom(1).getAtomicNumber() == 6) {
						templateBonds.get(i).set(1, currentBond);
						continue;
					}
					if(currentBond.getAtom(0).getAtomicNumber() == 8 || currentBond.getAtom(1).getAtomicNumber() == 8) {
						IAtom currentO = null;
						IAtom currentC = null;
						if(currentBond.getAtom(0).getAtomicNumber() == 8) {
							currentO = currentBond.getAtom(0);
							currentC = currentBond.getAtom(1);
						}
						else {
							currentO = currentBond.getAtom(1);
							currentC = currentBond.getAtom(0);
						}
						// Check that the C is connected to a nitrogen
						boolean connectedToNitrogen = false;
						for(IAtom connectedAtom : oxazoleTemplates[i].getConnectedAtomsList(currentC)) {
							if(connectedAtom.getAtomicNumber() == 7) {
								connectedToNitrogen = true;
							}
						}
						if(connectedToNitrogen) {
							templateBonds.get(i).set(2, currentBond);
						}
					}
				}
			}
			
			for(int i = 0; i < oxazoleTemplates.length; i++) {
				List<List<IBond>> matchingBondsList = ChemicalUtilities.findMatchingBondsFromTemplate(oxazoleTemplates[i], templateBonds.get(i), m.getMolecule());
				for(List<IBond> matchingBonds : matchingBondsList) {
					matchingBonds.get(0).setOrder(IBond.Order.SINGLE);
					if(matchingBonds.get(0).getAtom(0).getAtomicNumber() == 6) {
						matchingBonds.get(0).getAtom(0).setAtomTypeName("C.sp3");
						matchingBonds.get(0).getAtom(1).setAtomTypeName("N.sp3");
					}
					else {
						matchingBonds.get(0).getAtom(0).setAtomTypeName("N.sp3");
						matchingBonds.get(0).getAtom(1).setAtomTypeName("C.sp3");
					}
					matchingBonds.get(1).setOrder(IBond.Order.SINGLE);
					matchingBonds.get(1).getAtom(0).setAtomTypeName("C.sp3");
					matchingBonds.get(1).getAtom(1).setAtomTypeName("C.sp3");
					m.getMolecule().removeBond(matchingBonds.get(2));
					IAtom oxazoleO = null;
					IAtom newKetoneC = null;
					if(matchingBonds.get(2).getAtom(0).getAtomicNumber() == 8) {
						oxazoleO = matchingBonds.get(2).getAtom(0);
						newKetoneC = matchingBonds.get(2).getAtom(1);
					}
					else {
						oxazoleO = matchingBonds.get(2).getAtom(1);
						newKetoneC = matchingBonds.get(2).getAtom(0);
					}
					
					Atom ketoneO = new Atom("O");
					ketoneO.setAtomTypeName("O.sp2");
					newKetoneC.setAtomTypeName("C.sp2");
					molecule.addAtom(ketoneO);
					molecule.addBond(new Bond(newKetoneC, ketoneO, IBond.Order.DOUBLE));
					
					hasOxazole = true;
					oxazoleOList.add(oxazoleO);
				}
			}
			
			}
		if (printSteps) {
			//if (hazOxazole) {
				//System.out.println("Processed oxazole");
			//}
		}
	}

	/**
	 * Process ureido linkages
	 */
	private void processUreido() {
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}

		boolean hasUreido = false;
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();

			for (int i = 0; i < fragment.getMolecule().getAtomCount(); i++) {

				IAtom mainC = null;
				IAtom main_O = null;
				IAtom n_1 = null;
				IAtom c_1 = null;
				IAtom n_2 = null;
				IAtom c_2 = null;

				if (!fragment.getMolecule().getAtom(i).getAtomTypeName()
						.equals("C.sp2")) {
					continue;
				}
				mainC = fragment.getMolecule().getAtom(i);
				for (int j = 0; j < fragment.getMolecule()
						.getConnectedAtomsCount(mainC); j++) {
					IAtom connectedAtom = fragment.getMolecule()
							.getConnectedAtomsList(mainC).get(j);

					if (connectedAtom.getAtomTypeName().startsWith("O.")
							&& fragment.getMolecule()
									.getBond(mainC, connectedAtom).getOrder() == IBond.Order.DOUBLE) {
						main_O = connectedAtom;
					} else if (connectedAtom.getAtomTypeName().startsWith("N.")
							&& n_1 == null) {
						n_1 = connectedAtom;
					} else if (connectedAtom.getAtomTypeName().startsWith("N.")) {
						n_2 = connectedAtom;
					}

				}
				if (mainC == null || main_O == null || n_1 == null
						|| n_2 == null) {
					continue;
				}

				for (int j = 0; j < fragment.getMolecule()
						.getConnectedAtomsCount(n_1); j++) {
					IAtom connectedAtom = fragment.getMolecule()
							.getConnectedAtomsList(n_1).get(j);
					if (connectedAtom == mainC) {
						continue;
					} else if (connectedAtom.getAtomTypeName().startsWith("C.")) {
						c_1 = connectedAtom;
						break;
					}
				}
				for (int j = 0; j < fragment.getMolecule()
						.getConnectedAtomsCount(n_2); j++) {
					IAtom connectedAtom = fragment.getMolecule()
							.getConnectedAtomsList(n_2).get(j);
					if (connectedAtom == mainC) {
						continue;
					} else if (connectedAtom.getAtomTypeName().startsWith("C.")) {
						c_2 = connectedAtom;
						break;
					}
				}
				if (c_1 == null || c_2 == null) {
					continue;
				}

				// At this point we know that we have a -C-N-(C=O)-N-C-
				// structure). Now we break bonds.

				int indexToReplace = monomerFragments.indexOf(fragment);

				ArrayList<IBond> removedBonds = new ArrayList<IBond>();

				removedBonds.add(fragment.getMolecule().getBond(mainC, main_O));
				removedBonds.add(fragment.getMolecule().getBond(mainC, n_1));
				removedBonds.add(fragment.getMolecule().getBond(mainC, n_2));

				fragment.getMolecule().removeBond(
						fragment.getMolecule().getBond(mainC, main_O));
				fragment.getMolecule().removeAtom(main_O);
				fragment.getMolecule().removeBond(
						fragment.getMolecule().getBond(mainC, n_1));
				fragment.getMolecule().removeBond(
						fragment.getMolecule().getBond(mainC, n_2));
				fragment.getMolecule().removeAtom(mainC);

				IMoleculeSet fragments = ConnectivityChecker
						.partitionIntoMolecules(fragment.getMolecule());
				if (fragments.getAtomContainerCount() != 2) {
					// We must have had a ring, so re-add what we removed.
					fragment.getMolecule().addAtom(main_O);
					fragment.getMolecule().addAtom(mainC);
					for (IBond bond : removedBonds) {
						fragment.getMolecule().addBond(bond);
					}
					continue;
				}

				Fragment subfragment1 = new Fragment(
						fragments.getMolecule(0));
				Fragment subfragment2 = new Fragment(
						fragments.getMolecule(1));

				// In the future, additional data may be added to these
				// MonomerFragments to indicate the presence of this bond.
				// subfragment1.addOtherConnectedFragment(subfragment2);
				// subfragment2.addOtherConnectedFragment(subfragment1);

				monomerFragments.remove(indexToReplace);
				monomerFragments.add(indexToReplace, subfragment1);
				monomerFragments.add(indexToReplace + 1, subfragment2);
				q.add(subfragment1);
				q.add(subfragment2);
				hasUreido = true;
				break;
			}
		}
		if (printSteps) {
			if (hasUreido) {
				System.out.println("Processed ureido linkage");
			}
		}
	}

	/**
	 * Open up the rings of sulfur-containing beta-lactams: penams, penems, and
	 * cephems
	 */
	private void processBetaLactamLike() {
		// If this is a cephem, convert to a penam

		IMolecule cephemTemplate = null;
		try {
			cephemTemplate = SmilesIO.readSmiles("C1CN2C=CCSC12");
		} catch (Exception e) {
			e.printStackTrace();
		}
		// Store the double bond in the cephem template
		IBond cephemTemplateDoubleBond = null;
		// Store the cephem-carbon single bond that is to be broken
		IBond cephemTemplateSulfurBondToBreak = null;
		for (int i = 0; i < cephemTemplate.getBondCount(); i++) {
			if (cephemTemplate.getBond(i).getOrder() == IBond.Order.DOUBLE) {
				cephemTemplateDoubleBond = cephemTemplate.getBond(i);
				continue;
			}
			IBond candidateBond = cephemTemplate.getBond(i);
			if (candidateBond.getAtom(0).getAtomicNumber() != 16
					&& candidateBond.getAtom(1).getAtomicNumber() != 16) {
				continue;
			}
			IAtom currentSulfur = null;
			IAtom currentCarbon = null;
			if (candidateBond.getAtom(0).getAtomicNumber() == 16) {
				currentSulfur = candidateBond.getAtom(0);
				currentCarbon = candidateBond.getAtom(1);
			} else {
				currentSulfur = candidateBond.getAtom(1);
				currentCarbon = candidateBond.getAtom(0);
			}
			if (cephemTemplate.getConnectedAtomsCount(currentCarbon) != 2) {
				continue;
			}
			cephemTemplateSulfurBondToBreak = candidateBond;
		}
		ArrayList<IBond> cephemTemplateBonds = new ArrayList<IBond>();
		cephemTemplateBonds.add(cephemTemplateDoubleBond);
		cephemTemplateBonds.add(cephemTemplateSulfurBondToBreak);

		for (Fragment frag : monomerFragments) {
			List<List<IBond>> matchingBonds = ChemicalUtilities.findMatchingBondsFromTemplate(
					cephemTemplate, cephemTemplateBonds, frag.getMolecule());
			for (int i = 0; i < matchingBonds.size(); i++) {
				matchingBonds.get(i).get(0).setOrder(IBond.Order.SINGLE);
				matchingBonds.get(i).get(0).getAtom(0).setAtomTypeName("C.sp3");
				matchingBonds.get(i).get(0).getAtom(1).setAtomTypeName("C.sp3");
				frag.getMolecule().removeBond(matchingBonds.get(i).get(1));
				// Create new bond between the sulfur and the carbon in the
				// double bond to form a five-membered ring
				IAtom carbonPreviouslyConnectedToSulfur = null;
				IAtom carbonNewlyConnectingToSulfur = null;
				IAtom cephemSulfur = null;

				if (matchingBonds.get(i).get(1).getAtom(0).getAtomicNumber() == 13) {
					carbonPreviouslyConnectedToSulfur = matchingBonds.get(i)
							.get(1).getAtom(0);
					cephemSulfur = matchingBonds.get(i).get(1).getAtom(1);
				} else {
					carbonPreviouslyConnectedToSulfur = matchingBonds.get(i)
							.get(1).getAtom(1);
					cephemSulfur = matchingBonds.get(i).get(1).getAtom(0);
				}
				for (IAtom candidateNewlyConnectedCarbon : frag.getMolecule()
						.getConnectedAtomsList(
								carbonPreviouslyConnectedToSulfur)) {
					// If this candidate carbon is connected to one of the
					// previously double bonded carbons
					if (frag.getMolecule()
							.getConnectedAtomsList(
									candidateNewlyConnectedCarbon)
							.contains(matchingBonds.get(i).get(0).getAtom(0))
							|| frag.getMolecule()
									.getConnectedAtomsList(
											candidateNewlyConnectedCarbon)
									.contains(
											matchingBonds.get(i).get(0)
													.getAtom(1))) {
						carbonNewlyConnectingToSulfur = candidateNewlyConnectedCarbon;
					}
				}

				// Create a new bond
				IBond newCephemBond = new Bond(cephemSulfur,
						carbonNewlyConnectingToSulfur);

				frag.getMolecule().addBond(newCephemBond);
				frag.setBetaLactamS(cephemSulfur);
			}
		}

		IMolecule sulfurBetaLactamTemplate = null;
		try {
			sulfurBetaLactamTemplate = SmilesIO.readSmiles("CSC1CC(=O)N1");
		} catch (IOException | CDKException e) {
			e.printStackTrace();
		}
		// Store the bond in the four-membered ring that needs breaking
		IBond sulfurBetaLactamTemplateBond1 = null;
		// Store the bond in the five-membered ring that needs breaking
		IBond sulfurBetaLactamTemplateBond2 = null;
		for (int i = 0; i < sulfurBetaLactamTemplate.getBondCount(); i++) {
			IBond bond = sulfurBetaLactamTemplate.getBond(i);
			// Check if this is templateBond1
			if (bond.getAtom(0).getAtomicNumber() == 6
					&& bond.getAtom(1).getAtomicNumber() == 7) {
				if (ChemicalUtilities.hasNeighbourOfAtomicNumber(sulfurBetaLactamTemplate, bond.getAtom(0), 16)) {
					sulfurBetaLactamTemplateBond1 = bond;
				}
			} else if (bond.getAtom(1).getAtomicNumber() == 6
					&& bond.getAtom(0).getAtomicNumber() == 7) {
				if (ChemicalUtilities.hasNeighbourOfAtomicNumber(sulfurBetaLactamTemplate, bond.getAtom(1), 16)) {
					sulfurBetaLactamTemplateBond1 = bond;
				}
			} // Chcek if this is templateBond2
			else if (bond.getAtom(0).getAtomicNumber() == 6
					&& bond.getAtom(1).getAtomicNumber() == 16) { // Check if
																	// this is
																	// templateBond2
				if (!ChemicalUtilities.hasNeighbourOfAtomicNumber(sulfurBetaLactamTemplate, bond.getAtom(0), 7)) {
					sulfurBetaLactamTemplateBond2 = bond;
				}
			} else if (bond.getAtom(1).getAtomicNumber() == 6
					&& bond.getAtom(0).getAtomicNumber() == 16) { // Check if
																	// this is
																	// templateBond2
				if (!ChemicalUtilities.hasNeighbourOfAtomicNumber(sulfurBetaLactamTemplate, bond.getAtom(1), 7)) {
					sulfurBetaLactamTemplateBond2 = bond;
				}
			}
		}
		ArrayList<IBond> sulfurBetaLactamTemplateBonds = new ArrayList<IBond>();
		sulfurBetaLactamTemplateBonds.add(sulfurBetaLactamTemplateBond1);
		sulfurBetaLactamTemplateBonds.add(sulfurBetaLactamTemplateBond2);
		
		//check for clauvalnicAcidLike, these are the same as sulphur but oxygen replaces the sulfur
		IMolecule clavualnicAcidLikeTemplate = null;
		try {
			clavualnicAcidLikeTemplate = SmilesIO.readSmiles("COC1CC(=O)N1");
		} catch (IOException | CDKException e) {
			e.printStackTrace();
		}
		// Store the bond in the four-membered ring that needs breaking
		IBond clavualnicAcidLikeTemplateBond1 = null;
		// Store the bond in the five-membered ring that needs breaking
		IBond clavualnicAcidLikeTemplateBond2 = null;
		for (int i = 0; i < clavualnicAcidLikeTemplate.getBondCount(); i++) {
			IBond bond = clavualnicAcidLikeTemplate.getBond(i);
			// Check if this is templateBond1
			if (bond.getAtom(0).getAtomicNumber() == 6
					&& bond.getAtom(1).getAtomicNumber() == 7) {
				if (ChemicalUtilities.hasNeighbourOfAtomicNumber(clavualnicAcidLikeTemplate, bond.getAtom(0), 8)) {
					clavualnicAcidLikeTemplateBond1 = bond;
				}
			} else if (bond.getAtom(1).getAtomicNumber() == 6
					&& bond.getAtom(0).getAtomicNumber() == 7) {
				if (ChemicalUtilities.hasNeighbourOfAtomicNumber(clavualnicAcidLikeTemplate, bond.getAtom(1), 8)) {
					clavualnicAcidLikeTemplateBond1 = bond;
				}
			} // Chcek if this is templateBond2
			else if (bond.getAtom(0).getAtomicNumber() == 6
					&& bond.getAtom(1).getAtomicNumber() == 8) { // Check if
																	// this is
																	// templateBond2
				if (!ChemicalUtilities.hasNeighbourOfAtomicNumber(clavualnicAcidLikeTemplate, bond.getAtom(0), 7)) {
					clavualnicAcidLikeTemplateBond2 = bond;
				}
			} else if (bond.getAtom(1).getAtomicNumber() == 6
					&& bond.getAtom(0).getAtomicNumber() == 8) { // Check if
																	// this is
																	// templateBond2
				if (!ChemicalUtilities.hasNeighbourOfAtomicNumber(clavualnicAcidLikeTemplate, bond.getAtom(1), 7)) {
					clavualnicAcidLikeTemplateBond2 = bond;
				}
			}
		}
		ArrayList<IBond> clavualnicAcidLikeTemplateBonds = new ArrayList<IBond>();
		clavualnicAcidLikeTemplateBonds.add(clavualnicAcidLikeTemplateBond1);
		clavualnicAcidLikeTemplateBonds.add(clavualnicAcidLikeTemplateBond2);
		
		
		IMolecule gammaLactamBetaLactoneTemplate = null;
		try {
			gammaLactamBetaLactoneTemplate = SmilesIO.readSmiles("O=C1OC2CC(=O)NC12");
		} catch (IOException | CDKException e) {
			e.printStackTrace();
		}
		// Store the bond in the four-membered ring that needs breaking
		IBond gammaLactamBetaLactoneTemplateBond1 = gammaLactamBetaLactoneTemplate.getBond(gammaLactamBetaLactoneTemplate.getAtom(3),gammaLactamBetaLactoneTemplate.getAtom(8));
		// Store the bond in the five-membered ring that needs breaking
		IBond gammaLactamBetaLactoneTemplateBond2 = gammaLactamBetaLactoneTemplate.getBond(gammaLactamBetaLactoneTemplate.getAtom(3),gammaLactamBetaLactoneTemplate.getAtom(4));
		
		ArrayList<IBond> gammaLactamBetaLactoneTemplateBonds = new ArrayList<IBond>();
		gammaLactamBetaLactoneTemplateBonds.add(gammaLactamBetaLactoneTemplateBond1);
		gammaLactamBetaLactoneTemplateBonds.add(gammaLactamBetaLactoneTemplateBond2);		
		
		boolean foundMatch = false;
		for (Fragment frag : monomerFragments) {
			List<List<IBond>> matchingBonds = ChemicalUtilities.findMatchingBondsFromTemplate(
					sulfurBetaLactamTemplate, sulfurBetaLactamTemplateBonds, frag.getMolecule());
			if(matchingBonds.size() < 1){
				matchingBonds = ChemicalUtilities.findMatchingBondsFromTemplate(
						gammaLactamBetaLactoneTemplate, gammaLactamBetaLactoneTemplateBonds, frag.getMolecule());
			}
			if(matchingBonds.size() < 1){
				matchingBonds = ChemicalUtilities.findMatchingBondsFromTemplate(
						clavualnicAcidLikeTemplate, clavualnicAcidLikeTemplateBonds, frag.getMolecule());
			}	
			for (int i = 0; i < matchingBonds.size(); i++) {
				frag.getMolecule().removeBond(matchingBonds.get(i).get(0));
				frag.getMolecule().removeBond(matchingBonds.get(i).get(1));
				// If the molecule is no longer connected, then re-add these
				// bonds.
				if (!ConnectivityChecker.isConnected(frag.getMolecule())) {
					frag.getMolecule().addBond(matchingBonds.get(i).get(0));
					frag.getMolecule().addBond(matchingBonds.get(i).get(1));
				} else {
					foundMatch = true;
					//SmilesIO.drawMolecule(frag.getMolecule(), "afterbetalactam");
				}
			}
		}
		if (printSteps) {
			if (foundMatch) {
				System.out.println("Processed beta lactam like");
			}
		}
	}

	private void breakThioesters() {
		IMolecule template = null;
		try {
			template = SmilesIO.readSmiles("CSC(C)=O");
		} catch (IOException | CDKException e) {
			e.printStackTrace();
		}
		IBond templateBond = null;
		for (int i = 0; i < template.getBondCount(); i++) {
			// The bond of interest involves a sulfur and a carbon with three
			// connections
			if (template.getConnectedAtomsCount(template.getBond(i).getAtom(0)) == 3
					&& template.getBond(i).getAtom(1).getAtomicNumber() == 16
					|| template.getConnectedAtomsCount(template.getBond(i)
							.getAtom(1)) == 3
					&& template.getBond(i).getAtom(0).getAtomicNumber() == 16) {
				templateBond = template.getBond(i);
				break;
			}
		}

		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}

		while (!q.isEmpty()) {
			Fragment currentFragment = q.poll();
			List<IBond> matchedBonds = ChemicalUtilities.findMatchingBondsFromTemplate(template,
					templateBond, currentFragment.getMolecule());

			for (IBond currentBond : matchedBonds) {
				IAtom thioesterCarbon = null;
				if (currentBond.getAtom(0).getAtomicNumber() == 6) {
					thioesterCarbon = currentBond.getAtom(0);
				} else {
					thioesterCarbon = currentBond.getAtom(1);
				}
				currentFragment.getMolecule().removeBond(currentBond);
				Atom carboxylO = new Atom("O");
				carboxylO.setAtomTypeName("O.sp3");
				currentFragment.getMolecule().addAtom(carboxylO);
				Bond newBond = new Bond(thioesterCarbon, carboxylO);
				currentFragment.getMolecule().addBond(newBond);

				ArrayList<Fragment> fragments = currentFragment
						.partitionIntoMonomerFragments();

				if (fragments.size() == 2) {
					Fragment sFrag = null;
					Fragment cFrag = null;
					if (fragments.get(0).getMolecule()
							.contains(thioesterCarbon)) {
						cFrag = fragments.get(0);
						sFrag = fragments.get(1);
					} else {
						cFrag = fragments.get(1);
						sFrag = fragments.get(0);
					}
					int indexToReplace = monomerFragments
							.indexOf(currentFragment);
					monomerFragments.remove(indexToReplace);
					monomerFragments.add(indexToReplace, sFrag);
					monomerFragments.add(indexToReplace + 1, cFrag);
					q.add(sFrag);
					q.add(cFrag);
				}
			}
		}
	}

	/**
	 * Break bonds which connect two phenyl rings
	 */
	private void breakAdjoinedAromaticRings() {
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		boolean hasAdjoinedAromaticRings = false;
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IMolecule currentMolecule = fragment.getMolecule();
			IMolecule template = null;
			try {
				template = SmilesIO.readSmiles("C1=CC=C(C=C1)C1=CC=CC=C1");
			} catch (IOException | CDKException e) {
				e.printStackTrace();
			}
			IBond templateBond = null;
			for (int i = 0; i < template.getBondCount(); i++) {
				// The bond of interest is the bond connective two Carbons with
				// three connections
				if (template.getConnectedAtomsCount(template.getBond(i)
						.getAtom(0)) == 3
						&& template.getConnectedAtomsCount(template.getBond(i)
								.getAtom(1)) == 3) {
					templateBond = template.getBond(i);
					break;
				}
			}
			List<IBond> matchedBonds = ChemicalUtilities.findMatchingBondsFromTemplate(template,
					templateBond, currentMolecule);

			// At this point, conditions have been met. Break the bond if it
			// leads to one piece.
			for (IBond currentBond : matchedBonds) {
				currentMolecule.removeBond(currentBond);
				if (ConnectivityChecker.partitionIntoMolecules(currentMolecule)
						.getMoleculeCount() > 1) {
					currentMolecule.addBond(currentBond);
				} else {
					hasAdjoinedAromaticRings = true;
				}
			}
		}
		if (printSteps) {
			if (hasAdjoinedAromaticRings) {
				System.out.println("Processed adjoined aromatic rings");
			}
		}
	}

	/**
	 * Remove sugar groups
	 */
	private void breakSugarGroups() {
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		int numSugarBonds = 0;
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			ArrayList<Fragment> subfragments = new ArrayList<Fragment>();
			subfragments.add(fragment);
			IMolecule currentMolecule = fragment.getMolecule();
			// Prepare peptide template 1 (amino acids with C terminus) and
			// peptide template 2 (amino acids with a nitrogen attached to C
			// terminus)
			IMolecule[] sugarTemplates = new IMolecule[2];
			IBond[] sugarTemplateBonds = new IBond[2];
			try {
				sugarTemplates[0] = SmilesIO.readSmiles("COC1CCCCO1"); 
				sugarTemplates[1] = SmilesIO.readSmiles("COC1CCCO1");

			} catch (Exception e) {
				e.printStackTrace();
			}

			// Find the template bond
			for(int i = 0; i < sugarTemplateBonds.length; i++){
				IMolecule sugarTemplate = sugarTemplates[i];
				for (int j = 0; j < sugarTemplate.getBondCount(); j++) {
					if (sugarTemplate.getBond(j).getAtom(0).getAtomicNumber() == 6
							&& sugarTemplate.getBond(j).getAtom(1)
									.getAtomicNumber() == 8
							|| sugarTemplate.getBond(j).getAtom(0)
									.getAtomicNumber() == 8
							&& sugarTemplate.getBond(j).getAtom(1)
									.getAtomicNumber() == 6) {
						IAtom currentCarbon = null;
						IAtom currentOxygen = null;
						if (sugarTemplate.getBond(j).getAtom(0).getAtomicNumber() == 6) {
							currentCarbon = sugarTemplate.getBond(j).getAtom(0);
							currentOxygen = sugarTemplate.getBond(j).getAtom(1);
						} else {
							currentCarbon = sugarTemplate.getBond(j).getAtom(1);
							currentOxygen = sugarTemplate.getBond(j).getAtom(0);
						}

						if (sugarTemplate.getConnectedAtomsCount(currentCarbon) != 1) {
							continue;
						}
						sugarTemplateBonds[i] = sugarTemplate.getBond(j);
					}
				}
			}

			// Find all sugar bonds
			
			for(int i = 0; i < sugarTemplateBonds.length; i++){
				IMolecule sugarTemplate = sugarTemplates[i];
				IBond sugarTemplateBond = sugarTemplateBonds[i];
				List<IBond> sugarBonds = ChemicalUtilities.findMatchingBondsFromTemplate(
						sugarTemplate, sugarTemplateBond, currentMolecule);

				List<IAtom> carbonsNotInSugar = new ArrayList<IAtom>();
				for (IBond sugarBond : sugarBonds) {
					Fragment subfragment = null;
					for (Fragment possibleSubfragment : subfragments) {
						if (possibleSubfragment.getMolecule().contains(sugarBond)) {
							subfragment = possibleSubfragment;
						}
					}
					IAtom bondCarbon = null;
					IAtom bondOxygen = null;
					if (sugarBond.getAtom(0).getAtomTypeName().startsWith("C")) {
						bondCarbon = sugarBond.getAtom(0);
						bondOxygen = sugarBond.getAtom(1);
					} else {
						bondOxygen = sugarBond.getAtom(0);
						bondCarbon = sugarBond.getAtom(1);
					}

					subfragment.getMolecule().removeBond(sugarBond);
					// Check if this led to two pieces
					if (ConnectivityChecker.isConnected(subfragment.getMolecule())) {
						subfragment.getMolecule().addBond(sugarBond);
						continue;
					}

					sugarOxygens.add(bondOxygen);
					IAtom hydroxylO = new Atom("O");
					hydroxylO.setAtomTypeName("O.sp3");
					subfragment.getMolecule().addAtom(hydroxylO);
					subfragment.getMolecule().addBond(
							new Bond(hydroxylO, bondCarbon));

					ArrayList<Fragment> partitions = subfragment
							.partitionIntoMonomerFragments();

					// There must be two partitions
					Fragment sugarPart;
					Fragment nonSugarPart;
					if (partitions.get(0).getMolecule().contains(bondCarbon)) {
						nonSugarPart = partitions.get(0);
						sugarPart = partitions.get(1);
					} else {
						nonSugarPart = partitions.get(1);
						sugarPart = partitions.get(0);
					}
					nonSugarPart.setAtomInAttachedSugar(bondOxygen);
					sugarPart.setAtomOppositeSugarBond(bondCarbon);
					numSugarBonds++;
					int indexToReplace = subfragments.indexOf(subfragment);
					subfragments.remove(indexToReplace);
					subfragments.add(indexToReplace, nonSugarPart);
					subfragments.add(indexToReplace + 1, sugarPart);
				}
			}

			if(subfragments.size() > 1) {
				int indexToReplace = monomerFragments.indexOf(fragment);
				monomerFragments.remove(indexToReplace);
				for(int i = 0; i < subfragments.size(); i++) {
					monomerFragments.add(indexToReplace + i,
							subfragments.get(i));
				}
			}
		}
		if (printSteps) {
			if (numSugarBonds > 0) {
				System.out.println("Processed " + numSugarBonds
						+ " sugar bonds");
			}
		}
	}

	private void breakEsterLinkages() {
		IMolecule template = null;
		IBond templateBond = null;
		try {
			template = SmilesIO.readSmiles("COC=O");
		} catch (Exception e) {
			e.printStackTrace();
		}
		templateBond = template.getBond(template.getAtom(1),template.getAtom(2)); //ester bond
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		boolean foundEster = false;
		while (!q.isEmpty()) {
			Fragment frag = q.poll();
			List<IBond> esterBonds = ChemicalUtilities.findMatchingBondsFromTemplate(template,
					templateBond, frag.getMolecule());
			for (IBond esterBond : esterBonds) {
				foundEster = true;
				IAtom esterC = null;
				if (esterBond.getAtom(0).getAtomicNumber() == 6) {
					esterC = esterBond.getAtom(0);
				} else {
					esterC = esterBond.getAtom(1);
				}
				// Remove bond
				frag.getMolecule().removeBond(esterBond);
				IAtom carboxylO = new Atom("O");
				carboxylO.setAtomTypeName("O.sp3");
				frag.getMolecule().addAtom(carboxylO);
				frag.getMolecule().addBond(new Bond(esterC, carboxylO));
			}
			if (!ConnectivityChecker.isConnected(frag.getMolecule())) {
				ArrayList<Fragment> subfragments = frag
						.partitionIntoMonomerFragments();
				monomerFragments.remove(frag);
				for (Fragment subfrag : subfragments) {
					monomerFragments.add(subfrag);
				}
			}
		}
		if (printSteps) {
			if (foundEster) {
				System.out.println("Processed non-cyclic ester(s)");
			}
		}
	}

	/**
	 * Break off sulfate groups
	 */
	private void breakSulfateGroups() {
		IMolecule template = null;
		IBond templateBond = null;
		try {
			template = SmilesIO.readSmiles("COS");
		} catch (Exception e) {
			e.printStackTrace();
		}
		for (int i = 0; i < template.getBondCount(); i++) {
			IBond bond = template.getBond(i);
			if (bond.getAtom(0).getAtomicNumber() == 6
					|| bond.getAtom(1).getAtomicNumber() == 6) {
				templateBond = bond;
				break;
			}
		}
		boolean hasSulfate = false;
		for (Fragment frag : monomerFragments) {
			List<IBond> sulfateBonds = ChemicalUtilities.findMatchingBondsFromTemplate(template,
					templateBond, frag.getMolecule());
			for (IBond sulfateBond : sulfateBonds) {
				IAtom sulfateO = null;
				IAtom connectingCarbon = null;
				if (sulfateBond.getAtom(0).getAtomicNumber() == 8) {
					sulfateO = sulfateBond.getAtom(0);
					connectingCarbon = sulfateBond.getAtom(1);
				} else {
					connectingCarbon = sulfateBond.getAtom(0);
					sulfateO = sulfateBond.getAtom(1);
				}
				// Remove bond
				frag.getMolecule().removeBond(sulfateBond);
				// Check that the bond removal has led to two pieces, of of
				// which has exactly one sulfur and four oxygens
				IMoleculeSet partitions = ConnectivityChecker
						.partitionIntoMolecules(frag.getMolecule());
				if (partitions.getMoleculeCount() != 2) {
					frag.getMolecule().addBond(sulfateBond);
					continue;
				}
				IMolecule sulfatePiece = null;
				IMolecule nonSulfatePiece = null;
				if (partitions.getMolecule(0).contains(sulfateO)) {
					sulfatePiece = partitions.getMolecule(0);
					nonSulfatePiece = partitions.getMolecule(1);
				} else {
					sulfatePiece = partitions.getMolecule(1);
					nonSulfatePiece = partitions.getMolecule(0);
				}
				// Make sure there are no carbons and at least four oxygens
				int numOxygens = 0;
				boolean hasCarbon = false;
				for (int i = 0; i < sulfatePiece.getAtomCount(); i++) {
					if (sulfatePiece.getAtom(i).getAtomicNumber() == 8) {
						numOxygens++;
					}
					if (sulfatePiece.getAtom(i).getAtomicNumber() == 6) {
						hasCarbon = true;
					}
				}
				if (hasCarbon || numOxygens < 4) {
					frag.getMolecule().addBond(sulfateBond);
					continue;
				}
				hasSulfate = true;
				// At this point, we conclude that this is indeed a sulfate
				// piece. Set this fragment molecule the non-sulfur piece, set
				// enum.
				IAtom hydroxylO = new Atom("O");
				hydroxylO.setAtomTypeName("O.sp3");
				nonSulfatePiece.addAtom(hydroxylO);
				nonSulfatePiece.addBond(new Bond(connectingCarbon, hydroxylO));
				frag.setMolecule(nonSulfatePiece);
				frag.addTailoringDomain(TailoringDomainEnums.SULFOTRANSFERASE);
			}
			
			for (Fragment m : monomerFragments) { // check for *S(O)(O)O
				List<IAtom> sulfurAtoms = new ArrayList<IAtom>();
				for (int i = 0; i < m.getMolecule().getAtomCount(); i++) {
					if (m.getMolecule().getAtom(i).getAtomicNumber() == 16) {
						sulfurAtoms.add(m.getMolecule().getAtom(i));
					}
				}
				if (sulfurAtoms.size() == 0) {
					continue;
				}
				List<IAtom> oxygens = new ArrayList<IAtom>(); 
				for(IAtom sulfur : sulfurAtoms){
					for(IAtom connectedAtom : m.getMolecule().getConnectedAtomsList(sulfur)){
						if(connectedAtom.getAtomicNumber() == 8 && m.getMolecule().getConnectedAtomsList(connectedAtom).size() == 1){
							oxygens.add(connectedAtom);
						}
					}
					if(oxygens.size() == 3){
						m.getMolecule().removeAtomAndConnectedElectronContainers(sulfur);
						for(IAtom oxygen : oxygens){
							m.getMolecule().removeAtomAndConnectedElectronContainers(oxygen);
						}
						hasSulfate = true;
						m.getTailoringDomains().add(TailoringDomainEnums.SULFOTRANSFERASE);
					}
				}
			}
			
		}
		if (printSteps) {
			if (hasSulfate) {
				System.out.println("Processed sulfate");
			}
		}
	}

	/**
	 * Break aromatic carbon rings connected ba an ether
	 */
	private void breakAromaticEthers() {
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}

		int numBonds = 0;
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IMolecule currentMolecule = fragment.getMolecule();

			IMolecule[] templates = new IMolecule[2];
			IBond[] templateBonds = new IBond[2];
			try {
				// Benzene rings connected to each other
				templates[0] = SmilesIO.readSmiles("OC1=CC=CC=C1OC1=CC=CC=C1");
				templates[1] = SmilesIO
						.readSmiles("OC1=C(OC2=CC=CC=C2)C=CC=C1");
				// Benzene rings connected to another carbon chain
				// templates[2] = SmilesIO.readSmiles("CCOC1=CC=C(C)C=C1");
				// templates[3] = SmilesIO.readSmiles("CCOC1=CC(O)=CC(C)=C1");
				// templates[4] = SmilesIO.readSmiles("CCOC1=CC(C)=CC(O)=C1");

			} catch (Exception e) {
				e.printStackTrace();
			}
			// Find the template bonds for the first two templates
			for (int i = 0; i <= 1; i++) {
				// Walk from the hydroxyl oxygen to the ester bond
				IAtom hydroxylO = null;
				for (int j = 0; j < templates[i].getAtomCount(); j++) {
					if (templates[i].getAtom(j).getAtomicNumber() == 8
							|| templates[i].getConnectedAtomsCount(templates[i]
									.getAtom(j)) == 1) {
						hydroxylO = templates[i].getAtom(j);
						break;
					}
				}
				IAtom carbon1 = templates[i].getConnectedAtomsList(hydroxylO)
						.get(0);
				IAtom carbon2 = null;
				for (IAtom a : templates[i].getConnectedAtomsList(carbon1)) {
					if (templates[i].getConnectedAtomsCount(a) == 3) {
						carbon2 = a;
						break;
					}
				}
				for (IAtom a : templates[i].getConnectedAtomsList(carbon2)) {
					if (a.getAtomicNumber() == 8) {
						templateBonds[i] = templates[i].getBond(carbon2, a);
					}
				}
			}
			/*
			 * // Find the template bonds for templates 3-5 for(int i = 2; i <=
			 * 4; i++) { for(int j = 0; j < templates[i].getBondCount(); j++) {
			 * IBond bond = templates[i].getBond(j);
			 * if(bond.getAtom(0).getAtomicNumber() == 8 &&
			 * bond.getAtom(1).getAtomTypeName().equals("C.sp3")) {
			 * templateBonds[i] = bond; continue; }
			 * if(bond.getAtom(1).getAtomicNumber() == 8 &&
			 * bond.getAtom(0).getAtomTypeName().equals("C.sp3")) {
			 * templateBonds[i] = bond; continue; } } }
			 */
			for (int x = 0; x < templates.length; x++) {
				List<IBond> bondsToBreak = ChemicalUtilities.findMatchingBondsFromTemplate(
						templates[x], templateBonds[x], currentMolecule);
				for (IBond b : bondsToBreak) {
					if (currentMolecule.contains(b)) {
						currentMolecule.removeBond(b);
						numBonds++;
					}
				}
			}

			IMoleculeSet fragments = ConnectivityChecker
					.partitionIntoMolecules(currentMolecule);
			if (fragments.getMoleculeCount() > 1) {
				monomerFragments.remove(fragment);
				for (int i = 0; i < fragments.getMoleculeCount(); i++) {
					monomerFragments.add(new Fragment(fragments
							.getMolecule(i)));
				}
			}
		}
		if (printSteps) {
			if (numBonds > 0) {
				System.out.println("Processed " + numBonds
						+ " aromatic ether(s)");
			}
		}
	}

	/**
	 * Break sulfur-sulfur bonds
	 */
	private void breakDisulfideBridges() {
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		int numBonds = 0;
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IMolecule currentMolecule = fragment.getMolecule();

			IMolecule template = null;
			try {
				template = SmilesIO.readSmiles("CSSC");

			} catch (Exception e) {
				e.printStackTrace();
			}

			IBond templateBond = null;
			for (int i = 0; i < template.getBondCount(); i++) {
				if (template.getBond(i).getAtom(0).getAtomicNumber() == 16
						&& template.getBond(i).getAtom(1).getAtomicNumber() == 16) {
					templateBond = template.getBond(i);
					break;
				}
			}
			
			List<IBond> bondsToBreak = ChemicalUtilities.findMatchingBondsFromTemplate(template,
					templateBond, currentMolecule);
			
			if(bondsToBreak.size() > 0){
				if(fungal){
					IMolecule fungalTemplate = null;
					try {
						fungalTemplate = SmilesIO.readSmiles("O=C1NC2SSC11CC3=CC=CCC3N1C2=O");
					} catch(Exception e) {
						e.printStackTrace();
					}
					IBond fungalTemplateBond = fungalTemplate.getBond(fungalTemplate.getAtom(13),fungalTemplate.getAtom(14));

					List<IBond> fungalBondsToBreak = ChemicalUtilities.findMatchingBondsFromTemplate(fungalTemplate, fungalTemplateBond, currentMolecule);
					for(IBond bond : fungalBondsToBreak){
						currentMolecule.removeBond(bond);
					}
					for(IBond b : bondsToBreak){
						for(IAtom a : b.atoms()){
							currentMolecule.removeAtomAndConnectedElectronContainers(a);
						}
					}
				}else{
					for (IBond b : bondsToBreak) {
						if (currentMolecule.contains(b)) {
							currentMolecule.removeBond(b);
							numBonds++;
						}
					}
				}
				IMoleculeSet fragments = ConnectivityChecker
						.partitionIntoMolecules(currentMolecule);
				if (fragments.getMoleculeCount() > 1) {
					monomerFragments.remove(fragment);
					for (int i = 0; i < fragments.getMoleculeCount(); i++) {
						monomerFragments.add(new Fragment(fragments
								.getMolecule(i)));
					}
				}
			}
		}
		if (printSteps) {
			if (numBonds > 0) {
				System.out.println("Processed " + numBonds
						+ " disulfide bridge(s)");
			}
		}
	}

	/**
	 * Find, annotate, and remove N, C, and O methylations to amino acids
	 * TODO fix so will loop through multiple amino atoms
	 */
	private void processAAMethylations() {
		int numNMethylations = 0;
		int numOMethylations = 0;
		int numCMethylations = 0;
		for (Fragment m : monomerFragments) {
			for (int i = 0; i < m.getMolecule().getBondCount(); i++) {
				IBond currentBond = m.getMolecule().getBond(i);

				// N methylations
				if (m.getAminoNs().size() > 0) {
					IAtom candidateMethylC = null;
					if (currentBond.getAtom(0) == m.getAminoNs().get(0)
							&& currentBond.getAtom(1).getAtomicNumber() == 6) {
						candidateMethylC = currentBond.getAtom(1);
					} else if (currentBond.getAtom(1) == m.getAminoNs().get(0)
							&& currentBond.getAtom(0).getAtomicNumber() == 6) {
						candidateMethylC = currentBond.getAtom(0);
					}
					if (candidateMethylC != null
							&& m.getMolecule().getConnectedAtomsCount(
									candidateMethylC) == 1) {
						m.getMolecule().removeBond(currentBond);
						m.getMolecule().removeAtom(candidateMethylC);
						m.getTailoringDomains().add(
								TailoringDomainEnums.N_METHYLTRANSFERASE);
						numNMethylations++;
						continue;
					}
				}
				// C methylations
				if (m.getAminoNs().size() > 0) {
					IAtom candidateMethylC = null;
					IAtom candidateAlphaCarbon = null;
					if (currentBond.getAtom(0).getAtomicNumber() == 6
							&& currentBond.getAtom(1).getAtomicNumber() == 6) {
						if (m.getMolecule()
										.getConnectedAtomsList(
												currentBond.getAtom(0))
										.contains(m.getAminoNs().get(0))) {
							candidateAlphaCarbon = currentBond.getAtom(0);
							candidateMethylC = currentBond.getAtom(1);
						}
						if (m.getMolecule()
										.getConnectedAtomsList(
												currentBond.getAtom(1))
										.contains(m.getAminoNs().get(0))) {
							candidateAlphaCarbon = currentBond.getAtom(1);
							candidateMethylC = currentBond.getAtom(0);
						}
						if (m.getMolecule().getConnectedAtomsCount(
								candidateAlphaCarbon) == 4
								&& m.getMolecule().getConnectedAtomsCount(
										candidateMethylC) == 1) {
							m.getMolecule().removeBond(currentBond);
							m.getMolecule().removeAtom(candidateMethylC);
							m.getTailoringDomains().add(
									TailoringDomainEnums.C_METHYLTRANSFERASE);
							numCMethylations++;
							continue;
						}
					}
				}
				if (m.getAminoCs().size() > 0) {
					IAtom candidateMethylC = null;
					IAtom candidateAlphaCarbon = null;
					if (currentBond.getAtom(0).getAtomicNumber() == 6
							&& currentBond.getAtom(1).getAtomicNumber() == 6) {
						if (m.getMolecule()
								.getConnectedAtomsList(currentBond.getAtom(0))
								.contains(m.getAminoCs().get(0))
								) {
							candidateAlphaCarbon = currentBond.getAtom(0);
							candidateMethylC = currentBond.getAtom(1);
						}
						if (m.getMolecule()
								.getConnectedAtomsList(currentBond.getAtom(1))
								.contains(m.getAminoCs().get(0))
								) {
							candidateAlphaCarbon = currentBond.getAtom(1);
							candidateMethylC = currentBond.getAtom(0);
						}
						if (m.getMolecule().getConnectedAtomsCount(
								candidateAlphaCarbon) == 4
								&& m.getMolecule().getConnectedAtomsCount(
										candidateMethylC) == 1) {
							m.getMolecule().removeBond(currentBond);
							m.getMolecule().removeAtom(candidateMethylC);
							m.getTailoringDomains().add(
									TailoringDomainEnums.C_METHYLTRANSFERASE);
							numCMethylations++;
							continue;
						}
					}
				}
				// O methylations
				if (m.getAminoCs().size() > 0) {
					IAtom candidateMethylC = null;
					IAtom candidateO = null;
					if (currentBond.getAtom(0).getAtomicNumber() == 8
							&& currentBond.getAtom(1).getAtomicNumber() == 6) {
						candidateMethylC = currentBond.getAtom(1);
						candidateO = currentBond.getAtom(0);
					} else if (currentBond.getAtom(1).getAtomicNumber() == 8
							&& currentBond.getAtom(0).getAtomicNumber() == 6) {
						candidateMethylC = currentBond.getAtom(0);
						candidateO = currentBond.getAtom(1);
					}
					if (m.getMolecule().getConnectedAtomsList(candidateO)
							.contains(m.getAminoCs().get(0))
							&& m.getMolecule().getConnectedAtomsCount(
									candidateMethylC) == 1) {
						candidateO.setAtomTypeName("O.sp3");
						m.getMolecule().getBond(candidateO, m.getAminoCs().get(0))
								.setOrder(IBond.Order.DOUBLE);
						m.getMolecule().removeBond(currentBond);
						m.getMolecule().removeAtom(candidateMethylC);
						m.getTailoringDomains().add(
								TailoringDomainEnums.O_METHYLTRANSFERASE);
						numOMethylations++;
						continue;
					}
				}
			}
		}
		if (printSteps) {
			if (numNMethylations > 0) {
				System.out.println("Processed " + numNMethylations
						+ " N-methylated amino acids");
			}
			if (numOMethylations > 0) {
				System.out.println("Processed " + numOMethylations
						+ " O-methylated amino acids");
			}
			if (numCMethylations > 0) {
				System.out.println("Processed " + numCMethylations
						+ " C-methylated amino acids");
			}
		}
	}
	
	private void processAACHydroxyliations() {
		for (Fragment m : monomerFragments) {
			for (int i = 0; i < m.getMolecule().getBondCount(); i++) {
				IBond currentBond = m.getMolecule().getBond(i);
				if (m.getAminoNs().size() > 0) {
					IAtom candidateHydroxylO = null;
					IAtom candidateAlphaCarbon = null;
					if (currentBond.getAtom(0).getAtomicNumber() == 6
							&& currentBond.getAtom(1).getAtomicNumber() == 8) {
						candidateAlphaCarbon = currentBond.getAtom(0);
						candidateHydroxylO = currentBond.getAtom(1);
					}else if(currentBond.getAtom(0).getAtomicNumber() == 8
							&& currentBond.getAtom(1).getAtomicNumber() == 6){
						candidateAlphaCarbon = currentBond.getAtom(1);
						candidateHydroxylO = currentBond.getAtom(0);
					}
					if (m.getMolecule().getConnectedAtomsCount(
							candidateAlphaCarbon) == 4
							&& m.getMolecule().getConnectedAtomsCount(
									candidateHydroxylO) == 1) {
						m.getMolecule().removeBond(currentBond);
						m.getMolecule().removeAtom(candidateHydroxylO);
						m.getTailoringDomains().add(
								TailoringDomainEnums.C_HYDROXYLATION);
						continue;
					}
				}
				if (m.getAminoCs().size() > 0) {
					IAtom candidateHydroxylO = null;
					IAtom candidateAlphaCarbon = null;
					if (currentBond.getAtom(0).getAtomicNumber() == 6
							&& currentBond.getAtom(1).getAtomicNumber() == 8) {
						candidateAlphaCarbon = currentBond.getAtom(0);
						candidateHydroxylO = currentBond.getAtom(1);
					}else if(currentBond.getAtom(0).getAtomicNumber() == 8
							&& currentBond.getAtom(1).getAtomicNumber() == 6){
						candidateAlphaCarbon = currentBond.getAtom(1);
						candidateHydroxylO = currentBond.getAtom(0);
					}
					if (m.getMolecule().getConnectedAtomsCount(
							candidateAlphaCarbon) == 4
							&& m.getMolecule().getConnectedAtomsCount(
									candidateHydroxylO) == 1) {
						m.getMolecule().removeBond(currentBond);
						m.getMolecule().removeAtom(candidateHydroxylO);
						m.getTailoringDomains().add(
								TailoringDomainEnums.C_HYDROXYLATION);
						continue;
					}
				}
			}
		}
	}

	/**
	 * Find, annotate, and remove chlorine groups
	 */
	private void processChlorinations() {
		boolean hasChlorine = false;
		for (Fragment m : monomerFragments) {
			List<IAtom> chlorineAtoms = new ArrayList<IAtom>();
			for (int i = 0; i < m.getMolecule().getAtomCount(); i++) {
				if (m.getMolecule().getAtom(i).getAtomicNumber() == 17) {
					chlorineAtoms.add(m.getMolecule().getAtom(i));
				}
			}
			if (chlorineAtoms.size() == 0) {
				continue;
			}
			m.getTailoringDomains().add(TailoringDomainEnums.CHLORINATION);
			
			if(chlorineAtoms.size() >= 3) {
				for(IAtom chlorine : chlorineAtoms){
					IAtom connectedToChlorine = null;
					int numConnectedToSameAtom = 0;
					if(m.getMolecule().getConnectedAtomsCount(chlorine) != 1) continue;
					connectedToChlorine = m.getMolecule().getConnectedAtomsList(chlorine).get(0);
					for(IAtom atom : m.getMolecule().getConnectedAtomsList(connectedToChlorine)){
						if(atom.getAtomicNumber() == 17){
							numConnectedToSameAtom ++;
						}
					}
					if(numConnectedToSameAtom == 3){
						Atom firstO = new Atom("O");
						firstO.setAtomTypeName("O.sp3");
						Atom secondO = new Atom("O");
						firstO.setAtomTypeName("O.sp2");
						Atom carbonylCarbon = new Atom("C");
						carbonylCarbon.setAtomTypeName("C.sp2");
						m.getMolecule().addAtom(firstO);
						m.getMolecule().addAtom(secondO);
						m.getMolecule().addAtom(carbonylCarbon);
						m.getMolecule().addBond(new Bond(carbonylCarbon, firstO, IBond.Order.SINGLE));
						m.getMolecule().addBond(new Bond(carbonylCarbon, secondO, IBond.Order.DOUBLE));
						m.getMolecule().addBond(new Bond(carbonylCarbon, connectedToChlorine, IBond.Order.SINGLE));
						break;
					}
				}
			}
			for (IAtom chlorine : chlorineAtoms) {
				// Check that there is exactly one connection
				if (m.getMolecule().getConnectedAtomsCount(chlorine) != 1) {
					continue;
				}
				hasChlorine = true;
				m.getMolecule().removeBond(
						m.getMolecule().getBond(
								chlorine,
								m.getMolecule().getConnectedAtomsList(chlorine)
										.get(0)));
				m.getMolecule().removeAtom(chlorine);
			}
		}
		if (printSteps) {
			if (hasChlorine) {
				System.out.println("Processed chlorination");
			}
		}
	}

	/**
	 * Find and annotate epoxide groups
	 */
	private void findEpoxides() {
		boolean hasEpoxide = false;
		for (Fragment m : monomerFragments) {
			for (int i = 0; i < m.getMolecule().getAtomCount(); i++) {
				IAtom o = m.getMolecule().getAtom(i);
				if (o.getAtomicNumber() != 8)
					continue;
				if (m.getMolecule().getConnectedAtomsCount(o) != 2)
					continue;
				IAtom c1 = m.getMolecule().getConnectedAtomsList(o).get(0);
				IAtom c2 = m.getMolecule().getConnectedAtomsList(o).get(1);
				if (c1.getAtomicNumber() != 6 || c2.getAtomicNumber() != 6)
					continue;
				if (!m.getMolecule().getConnectedAtomsList(c1).contains(c2))
					continue;
				hasEpoxide = true;
				m.getTailoringDomains().add(TailoringDomainEnums.P450);
			}
		}
		if (printSteps) {
			if (hasEpoxide) {
				System.out.println("Found epoxide");
			}
		}
	}
}
