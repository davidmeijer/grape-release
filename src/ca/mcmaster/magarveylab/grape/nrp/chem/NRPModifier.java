package ca.mcmaster.magarveylab.grape.nrp.chem;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;

import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Element;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.isomorphism.mcss.RMap;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;
import org.openscience.cdk.smiles.smarts.parser.SMARTSParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

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

	IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
	UniversalIsomorphismTester uit = new UniversalIsomorphismTester();
	private final static boolean imageDump = false;
	private final static boolean showBonds = false;
	private final static boolean printSteps = false;

	private List<Fragment> monomerFragments;

	// This global list is necessary for determining the presence and
	// directionality of lactone bonds in the case that there is one lactone
	// whose adjacent monomers are separated in the peptide cleavage stage.
	private List<IAtom> lactoneCarboxylCList;
	private List<IAtom> thiazoleSList;
	private List<IAtom> thiazolineSList;
	private List<IAtom> methylatedCarbonsList;
	private List<IAtom> oxazoleOList;
	private List<IAtom> oxazolineOList;
	private List<IAtom> sugarCarbons;
	private ChemicalAbstraction chemicalAbstraction;

	/**
	 * Constructor for this NRPModifier. Typical usage: construct an NRPModifier
	 * from an original NRP IAtomContainer, and get the output from method
	 * getMonomerFragments.
	 * 
	 * @param originalAtomContainer
	 */
	public NRPModifier(IAtomContainer originalAtomContainer, ChemicalAbstraction chemicalAbstraction) {
		lactoneCarboxylCList = new ArrayList<IAtom>();
		methylatedCarbonsList = new ArrayList<IAtom>();
		oxazoleOList = new ArrayList<IAtom>();
		thiazoleSList = new ArrayList<IAtom>();
		thiazolineSList = new ArrayList<IAtom>();
		oxazoleOList = new ArrayList<IAtom>();
		oxazolineOList = new ArrayList<IAtom>();
		sugarCarbons = new ArrayList<IAtom>();
		monomerFragments = new ArrayList<Fragment>();
		monomerFragments.add(new Fragment(originalAtomContainer));
		this.chemicalAbstraction = chemicalAbstraction;
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
	public void performAllNrpModifications() { //TODO start splitting these up into different classes
		// if(false)
		if (!imageDump) {
			breakDisulfideBridges(); 
			breakAromaticEthers();
			breakAdjoinedAromaticRings();
			processBetaLactamLike();
			processMonobactams();
			breakLinearEthers(); 
			processUniqueSubstructures();
			processThiazs();
			processOxazs();
			modifyImine(); 
			breakThioesters();
			openLactoneRings();
			processUreido();
			breakEsterLinkages();
			breakPeptideBonds(); // do a check here and see why in prenly the amino acids are not being cleaved
			breakSugarGroups();
			breakSulfateGroups();
			processAAMethylations();
			processAACHydroxyliations();
			processHalogenations();
			processSecondaryAmineExtensions();
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
			processThiazs();
			processOxazs();
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
			processHalogenations();
			findEpoxides();
			Fragment.drawMonomerFragments(monomerFragments, GrapeMain.currentName
					+ "_04-Tailoring");
		}

		// Update tailoring fields fields
		for (Fragment m : monomerFragments) {
			for (IAtom t : methylatedCarbonsList) {
				if (m.getAtomContainer().contains(t)) {
					m.addTailoringDomain(TailoringDomainEnums.C_METHYLTRANSFERASE);
				}
			}
			for (IAtom t : thiazoleSList) {
				if (m.getAtomContainer().contains(t)) {
					m.setThiazoleS(t);
					m.addTailoringDomain(TailoringDomainEnums.THIAZOLE);
				}
			}
			for (IAtom t : thiazolineSList) {
				if (m.getAtomContainer().contains(t)) {
					m.setThiazolineS(t);
					m.addTailoringDomain(TailoringDomainEnums.THIAZOLINE);
				}
			}
			for (IAtom o : oxazoleOList) {
				if (m.getAtomContainer().contains(o)) {
					m.setOxazoleO(o);
					m.addTailoringDomain(TailoringDomainEnums.OXAZOLE);
				}
			}
			for (IAtom o : oxazolineOList) {
				if (m.getAtomContainer().contains(o)) {
					m.setOxazolineO(o);
					m.addTailoringDomain(TailoringDomainEnums.OXAZOLINE);
				}
			}
			for (IAtom o : sugarCarbons) {
				if (m.getAtomContainer().contains(o)) {
					m.setFragmentType(FragmentType.SUGAR);
				}
			}
			if(m.getBetaLactamS() != null) {
				m.addTailoringDomain(TailoringDomainEnums.SULFUR_BETA_LACTAM);
			}
			for(IAtom a : m.getLactamAtoms()){
				if (m.getAtomContainer().contains(a)) {
					m.addTailoringDomain(TailoringDomainEnums.LACTAM);
				}
			}
		}
	}
	
	private void processUniqueSubstructures() {
		processKendoLikeSubstructures();
		processPieriLikeSubstructures();
		processAnthramycinLikeSubstructures();
		breakHybridDiCystines();
		processSpectinomycinLikedoubleGlycosidicBond();
		processHygromycinLikeGlycosidicBond();
		
	}

	private void processHygromycinLikeGlycosidicBond() { //CN[C@H]1C[C@@H](N)[C@H](O)[C@@H](O[C@@H]2O[C@H](CO)[C@H](O)[C@@H]3OC4(O[C@H]23)O[C@H]([C@@H](N)CO)[C@H](O)[C@@H](O)[C@H]4O)[C@@H]1O
IAtomContainer tempalte = null;
		
		try {
			tempalte = SmilesIO.readSmilesTemplates("COC1(C)OCCO1");
		} catch(Exception e) {
			e.printStackTrace();
		}
		IBond templateBondToBreak = tempalte.getBond(tempalte.getAtom(2),tempalte.getAtom(4));
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IAtomContainer mol = fragment.getAtomContainer();
			List<IBond> templateBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(tempalte, templateBondToBreak, mol);
			for(IBond bond : templateBondMatches){
				mol.removeBond(bond);
				IAtomContainerSet molFrags = ConnectivityChecker.partitionIntoMolecules(mol);
				if (molFrags.getAtomContainerCount() > 1){
					mol.addBond(bond);
					break;
				}
			}
		}
	}

	private void processSpectinomycinLikedoubleGlycosidicBond() { //CN[C@@H]1[C@H](O)[C@H](NC)[C@H]2O[C@]3(O)[C@@H](O[C@H](C)CC3=O)O[C@@H]2[C@H]1O
		IAtomContainer tempalte = null;
		
		try {
			tempalte = SmilesIO.readSmilesTemplates("OC1COCCO1");
		} catch(Exception e) {
			e.printStackTrace();
		}
		IBond templateBondToBreak = tempalte.getBond(tempalte.getAtom(1),tempalte.getAtom(6));
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IAtomContainer mol = fragment.getAtomContainer();
			List<IBond> templateBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(tempalte, templateBondToBreak, mol);
			for(IBond bond : templateBondMatches){
				mol.removeBond(bond);
				IAtomContainerSet molFrags = ConnectivityChecker.partitionIntoMolecules(mol);
				if (molFrags.getAtomContainerCount() > 1){
					mol.addBond(bond);
					break;
				}
			}
		}
	}

	private void processPeptideEpoxiKetones() {
		IAtomContainer peptideEpoxiKetoneTemplate = null;
				
		try {
			peptideEpoxiKetoneTemplate = SmilesIO.readSmilesTemplates("CC1(CO1)C(=O)CN");
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
			IAtomContainer mol = fragment.getAtomContainer();
			List<IBond> templateBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(peptideEpoxiKetoneTemplate, templateBondToBreak, mol);
			for(IBond bond : templateBondMatches){
				mol.removeBond(bond);
				for(IAtom atom : bond.atoms()){
					if(mol.getConnectedAtomsCount(atom) == 2){
						
						Atom oxygen = new Atom("O");
						oxygen.setAtomTypeName("O.sp3");
						mol.addAtom(oxygen);
						mol.addBond(
								fragment.getAtomContainer().getAtomNumber(atom),
								fragment.getAtomContainer().getAtomNumber(oxygen),
								IBond.Order.SINGLE);
					}
				}
			}

			IAtomContainerSet molFrags = ConnectivityChecker.partitionIntoMolecules(mol);
			if (molFrags.getAtomContainerCount() > 1){
				int indexToReplace = monomerFragments.indexOf(fragment);
				monomerFragments.remove(indexToReplace);
				IAtomContainer[] molPiece = new IAtomContainer[2];
				try {
					molPiece[0] = SmilesIO.readSmilesTemplates("CC1CO1");
					molPiece[1] = SmilesIO.readSmilesTemplates("OCC1CO1");
				} catch(Exception e) {
					e.printStackTrace();
				}
				 int addCount = 0;
				for(Fragment fragmentToAdd : fragment.partitionIntoMonomerFragments()){
					IAtomContainer molFrag = fragmentToAdd.getAtomContainer();
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
					}
					monomerFragments.add(indexToReplace + addCount, fragmentToAdd);
					addCount += 1;
				}
			}
		}
	}

	private void processMonobactams() {
		
		IAtomContainer monobactamTemplate = null;
		
		try {
			monobactamTemplate = SmilesIO.readSmilesTemplates("O=C1CCN1");
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
			IAtomContainer mol = fragment.getAtomContainer();
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
									fragment.getAtomContainer().getAtomNumber(atom),
									fragment.getAtomContainer().getAtomNumber(oxygen),
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
								fragment.getAtomContainer().getAtomNumber(atom),
								fragment.getAtomContainer().getAtomNumber(oxygen),
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
			IAtomContainer mol = fragment.getAtomContainer();
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
					if(connectedC.getAtomicNumber() != 6 //Must be connected to two carbons
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
							fragment.getAtomContainer().getAtomNumber(connectedC),
							fragment.getAtomContainer().getAtomNumber(oxygen),
							IBond.Order.SINGLE);
				}
			}
			//set new fragments if not connected
			IAtomContainerSet fragments = ConnectivityChecker
					.partitionIntoMolecules(mol);
			if (fragments.getAtomContainerCount() > 1) { //add bridge tailor
				monomerFragments.remove(fragment);
				for (Fragment newFragment : fragment.partitionIntoMonomerFragments()) {
					monomerFragments.add(newFragment);
				}
			}
		}
	}

	private void breakHybridDiCystines() {
		IAtomContainer template = null;
		try {template = SmilesIO.readSmilesTemplates("CC1=NC(=CS1)C1=NC(C)=CS1");} catch (Exception e) {}
		IBond templateToBreakBond = template.getBond(template.getAtom(8),template.getAtom(9));
		IAtomContainer template2 = null;
		try {template2 = SmilesIO.readSmilesTemplates("CC1CSC(N1)C1CSC(C)=N1");} catch (Exception e) {}
		IBond templateToBreakBond2 = template2.getBond(template2.getAtom(1),template2.getAtom(0));
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {

			Fragment fragment = q.poll();
			List<IBond> templateDoubleBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateToBreakBond, fragment.getAtomContainer());
			templateDoubleBondMatches.addAll(ChemicalUtilities.findMatchingBondsFromTemplate(template2, templateToBreakBond2, fragment.getAtomContainer()));
			if(templateDoubleBondMatches.size() < 1) continue;
			IAtomContainer mol = fragment.getAtomContainer();
			for(IBond bond : templateDoubleBondMatches){
				mol.removeBond(bond);
				for(IAtom atom : bond.atoms()){
					IAtomContainer carboxylicAcid = null;
					try {carboxylicAcid = SmilesIO.readSmilesTemplates("C(O)=O");} catch (Exception e) {}
					IAtom carboxilicC = carboxylicAcid.getAtom(0);
					mol.add(carboxylicAcid);
					mol.addBond(
							mol.getAtomNumber(atom),
							mol.getAtomNumber(carboxilicC),
							IBond.Order.SINGLE);
					try {
						AtomContainerManipulator.percieveAtomTypesAndConfigureUnsetProperties(mol);
					} catch (CDKException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			//set new fragments if not connected
			IAtomContainerSet fragments = ConnectivityChecker
					.partitionIntoMolecules(mol);
			if (fragments.getAtomContainerCount() > 1) { //add bridge tailor
				monomerFragments.remove(fragment);
				for (Fragment newFragment : fragment.partitionIntoMonomerFragments()) {
					monomerFragments.add(newFragment);
				}
			}
		}
	}
	
	private void processAnthramycinLikeSubstructures() {
		IAtomContainer[] templates = new IAtomContainer[2];
		try {
			templates[0] = SmilesIO.readSmilesTemplates("O=C1NCC=NC=C1");
			templates[1] = SmilesIO.readSmilesTemplates("O=C1CCN=CCN1");
		} catch(Exception e) {
			e.printStackTrace();
		}
		IBond[] templatesBond = new IBond[]{
				templates[0].getBond(templates[0].getAtom(4),templates[0].getAtom(5)),
				templates[1].getBond(templates[1].getAtom(4),templates[1].getAtom(5))
		};

		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IAtomContainer mol = fragment.getAtomContainer();
			List<IBond> templateBondMatches = new ArrayList<IBond>();
			int templateIndex = 0; //template index to check
			boolean match = false; //change to true if a match comes up in the next statement
			while(!match && templates.length > templateIndex){
				templateBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(templates[templateIndex], templatesBond[templateIndex], mol);
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
		IAtomContainer template = null;
		try {template = SmilesIO.readSmilesTemplates("CC1=C(C)C(=O)C(O)=C(O)N1");} catch (Exception e) {}
		IBond templateToBreakBond = template.getBond(template.getAtom(0),template.getAtom(1));
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		// go through each fragment.
		while (!q.isEmpty()) {

			Fragment fragment = q.poll();
			IAtomContainer molClone = null;
			try {molClone = fragment.getAtomContainer().clone();}catch (CloneNotSupportedException e) {}
			List<IBond> templateDoubleBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateToBreakBond, molClone);
			
			if(templateDoubleBondMatches.size() < 1) continue;
			
			for(IBond bond : templateDoubleBondMatches){
				molClone.removeBond(bond);
				
				if(!ConnectivityChecker.isConnected(molClone)){
					IAtomContainerSet splitFragments = ConnectivityChecker.partitionIntoMolecules(molClone);
					if(splitFragments.getAtomContainerCount() == 2){
						if(splitFragments.getAtomContainer(0).getAtomCount() > splitFragments.getAtomContainer(1).getAtomCount()){
							molClone.addBond(bond);
							continue;
						}
						if(splitFragments.getAtomContainer(0).getAtomCount() > splitFragments.getAtomContainer(1).getAtomCount()){
							molClone = splitFragments.getAtomContainer(0);
						}else{
							molClone = splitFragments.getAtomContainer(1);
						}
						IAtom terminalCarbon = null;
						if(molClone.contains(bond.getAtom(0))){
							terminalCarbon = bond.getAtom(0);
						}else if(molClone.contains(bond.getAtom(1))){
							terminalCarbon = bond.getAtom(1);
						}else{
							System.err.println("Something went wrong with method processPieriLikeSubstructures -- This should never happen, this method was ignored for the final breakdown of this compound");
							break;
						}
						IAtomContainer pieceToAdd = null;
						try {pieceToAdd = SmilesIO.readSmilesTemplates("CC(C=O)C(=O)CC(O)=O");} catch (Exception e) {}
						IAtom atomToConnect = pieceToAdd.getAtom(2);
						molClone.add(pieceToAdd);
						molClone.addBond(molClone.getAtomNumber(terminalCarbon), molClone.getAtomNumber(atomToConnect), IBond.Order.SINGLE);
						fragment.setAtomContainer(molClone);
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

		IAtomContainer template = null;
		try {template = SmilesIO.readSmilesTemplates("CC1=C2OC(C)(O)C=C2C(C)=C(O)C1=O");} catch (Exception e) {}
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
			List<IBond> templateDoubleBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateDoubleBond, fragment.getAtomContainer());
			
			if(templateDoubleBondMatches.size() < 1) continue;
			
			List<IBond> templateEtherBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateEtherBond, fragment.getAtomContainer());
			List<IBond> templateHydroxylBondMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateHydroxylBond, fragment.getAtomContainer());
			
			for(IBond bond : templateDoubleBondMatches){
				fragment.getAtomContainer().removeBond(bond);
				for(IAtom atom : bond.atoms()){
					if(fragment.getAtomContainer().getConnectedAtomsCount(atom) == 2){ //Add ketone to carbon in the starter unit
						Atom ketone = new Atom("O");
						ketone.setAtomTypeName("O.sp2");
						fragment.getAtomContainer().addAtom(ketone);
						fragment.getAtomContainer().addBond(
								fragment.getAtomContainer().getAtomNumber(ketone),
								fragment.getAtomContainer().getAtomNumber(atom),
								IBond.Order.DOUBLE);
					}else if(fragment.getAtomContainer().getConnectedAtomsCount(atom) == 1){ //Add carboxcylic acid to carbon not in starter unit
						IAtomContainer carboxylicAcid = null;
						try {carboxylicAcid = SmilesIO.readSmilesTemplates("C(O)=O");} catch (Exception e) {}
						IAtom carboxilicC = carboxylicAcid.getAtom(0);
						fragment.getAtomContainer().add(carboxylicAcid);
						fragment.getAtomContainer().addBond(
								fragment.getAtomContainer().getAtomNumber(atom),
								fragment.getAtomContainer().getAtomNumber(carboxilicC),
								IBond.Order.SINGLE);
						
					}else{
						System.err.println("Something went wrong with method processKendoLikeSubstructures (Double bond match doesn't have appropriate carbons): reattached broken bond"); //TODO make exception
						fragment.getAtomContainer().addBond(bond);
					}				
				}		
			}
			for(IBond bond : templateEtherBondMatches){
				fragment.getAtomContainer().removeBond(bond);
				if(!ConnectivityChecker.isConnected(fragment.getAtomContainer())){
					fragment.getAtomContainer().addBond(bond);
					System.err.println("Something went wrong with method processKendoLikeSubstructures (Breaking ether makes AtomContainer not connected): reattached broken bond");
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

			for (int j = 0; j < fragment.getAtomContainer().getAtomCount(); j++) {
				if (!(fragment.getAtomContainer().getAtom(j).getAtomTypeName()
						.equals("N.sp2"))) {
					continue;
				}
				IAtom imineN = fragment.getAtomContainer().getAtom(j);
				List<Set<IAtom>> rings = ChemicalUtilities.getSmallestRings(fragment.getAtomContainer());
				boolean inSmallRing = false;
				for(Set<IAtom> ringAtoms : rings){
					if(ringAtoms.contains(imineN)){
						inSmallRing = true;
						break;
					}
				}
				if(inSmallRing) continue;
				
				IAtom imineC = null;
				for (IAtom c : fragment.getAtomContainer().getConnectedAtomsList(imineN)) {
					if (c.getAtomTypeName().equals("C.sp2")) {

						if (fragment.getAtomContainer()
										.getBond(imineN, c)
										.getOrder() == IBond.Order.DOUBLE
										&& fragment.getAtomContainer().getConnectedAtomsList(c).size() == 2){
							imineC = c;
							
						}
					}
				}
				if(imineC == null) continue;
				// Then check if breaking the imineN and imineC leaves
				// one piece.
				IBond removedBond = fragment.getAtomContainer().removeBond(
						imineN, imineC);

				if (ConnectivityChecker.isConnected(fragment.getAtomContainer())) {
					//Do the removal of double bond and create single
					fragment.getAtomContainer().addBond(
							fragment.getAtomContainer().getAtomNumber(imineN),
							fragment.getAtomContainer().getAtomNumber(imineC),
							IBond.Order.SINGLE);
					Atom ketone = new Atom("O");
					ketone.setAtomTypeName("O.sp2");
					fragment.getAtomContainer().addAtom(ketone);
					fragment.getAtomContainer().addBond(
							fragment.getAtomContainer().getAtomNumber(ketone),
							fragment.getAtomContainer().getAtomNumber(imineC),
							IBond.Order.DOUBLE);
					
				}else{
					//re-add the original bond
					fragment.getAtomContainer().addBond(removedBond);
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
		List<IBond> peptideBonds = new ArrayList<IBond>();

		List<Set<IAtom>> smallRings = new ArrayList<Set<IAtom>>();
		List<IAtom> atomsInRings = new ArrayList<IAtom>();
		
		// Prepare peptide template 1 (amino acids with C terminus) and
					// peptide template 2 (amino acids with a nitrogen attached to C
					// terminus)
					IAtomContainer[] peptideTemplates = new IAtomContainer[5];
					IBond[] peptideTemplateBonds = new IBond[5];
					try {
						peptideTemplates[0] = SmilesIO.readSmilesTemplates("CCNC(C)=O");
						peptideTemplates[1] = SmilesIO.readSmilesTemplates("CC(=O)NC=C");
						peptideTemplates[2] = SmilesIO.readSmilesTemplates("CC\\N=C(/C)O");
						peptideTemplates[3] = SmilesIO.readSmilesTemplates("CCNC(C)OC");
						peptideTemplates[4] = SmilesIO.readSmilesTemplates("C\\C(O)=N/C=C");
						// peptideTemplates[5] = SmilesIO.readSmiles("CC(=O)NC=C");
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
						IAtomContainer peptideTemplate = peptideTemplates[i];
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
		
		for (Fragment frag : monomerFragments) {
			IAtomContainer m = frag.getAtomContainer();
			// Keep track of rings with six or fewer atoms. If the C-N or C=N bond participates in this ring, break it but do not
			// change the connectivity annotation
			smallRings.addAll(ChemicalUtilities.getSmallestRings(frag.getAtomContainer()));
			
			// Keep track of atoms in rings of any size
			atomsInRings.addAll(ChemicalUtilities.getAtomsInRings(m));
			
			for (int x = 0; x < peptideTemplates.length; x++) {
				IAtomContainer peptideTemplate = peptideTemplates[x];
				IBond peptideTemplateBond = peptideTemplateBonds[x];
				List<List<RMap>> templateMatchMap = null;
				try {
					templateMatchMap = uit
							.getSubgraphMaps(m, peptideTemplate);
				} catch (CDKException e) {
					e.printStackTrace();
				}
				for (int i = 0; i < templateMatchMap.size(); i++) {
					// None of the carbons attached to another carbon can be
					// terminal.
//					boolean foundWrongTerminalCarbon = false;
//					for (int j = 0; j < templateMatchMap.get(i).size(); j++) {
//						IBond currentTemplateMatchBond = peptideTemplate
//								.getBond(templateMatchMap.get(i).get(j)
//										.getId2());
//						for (int k = 0; k < 2; k++) {
//							// For each bond in the AtomContainer structure
//							IAtom a = m.getBond(
//									templateMatchMap.get(i).get(j).getId1())
//									.getAtom(k);
//							if (m.getConnectedAtomsCount(a) == 1
//									&& a.getAtomicNumber() == 6) {
//								for (IAtom connectedAtom : m
//										.getConnectedAtomsList(a)) {
//									if (connectedAtom.getAtomicNumber() == 6) {
//										foundWrongTerminalCarbon = true;
//									}
//								}
//							}
//						}
//					}
//					if (foundWrongTerminalCarbon) {
//						continue;
//					}
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
							if (nitrogenAlreadyPresent && frag.getAtomContainer().getConnectedAtomsCount(currentNitrogen) < 3) {
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
		
		// Swap the order of peptide bonds so that bonds with the carbon in a cycle are processed first.
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
		IBond lactamBond = null;
		if(lastUnvisitedIndex > 0){
			lactamBond = peptideBonds.get(lastUnvisitedIndex);
		}
		
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		while (!q.isEmpty()) {
			Fragment subfragment = q.poll();
			// get peptide bond index or move on if there are no peptide bonds
			int peptideBondIndex = -1;
			int bondCount = subfragment.getAtomContainer().getBondCount();
			for (int i = 0; i < bondCount; i++) {
				if (peptideBonds.contains(subfragment.getAtomContainer().getBond(i))) {
					peptideBondIndex = i;
					break;
				}
			}
			if (peptideBondIndex == -1) {
				continue;
			}
			IBond peptideBond = subfragment.getAtomContainer().getBond(
					peptideBondIndex);
			if(peptideBond.equals(lactamBond)){
				for(IAtom atom : peptideBond.atoms()){
					subfragment.addLactamAtom(atom);
				}
			}
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
				List<IAtom> connectedBetaAtoms = subfragment.getAtomContainer()
						.getConnectedAtomsList(betaCarbon);
				for (IAtom a : connectedBetaAtoms) {
					if (a.getAtomTypeName().startsWith("O.")) {
						// This is the carbonyl oxygen
						a.setAtomTypeName("O.sp2");
						subfragment.getAtomContainer().getBond(betaCarbon, a)
								.setOrder(IBond.Order.DOUBLE);
					}
				}
			}

			subfragment.getAtomContainer().removeBond(peptideBondIndex);

			// Add -O for the OH group on the carbon

			Atom hydroxideO = new Atom("O");
			hydroxideO.setAtomTypeName("O.sp3");
			subfragment.getAtomContainer().addAtom(hydroxideO);
			subfragment.getAtomContainer().addBond(new Bond(betaCarbon, hydroxideO));

			List<Fragment> fragments = subfragment
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
				if (fragments.get(0).getAtomContainer().contains(betaCarbon)) {
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
					if(NFragment.getAtomContainer().contains(aminoNtoAdd)){
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

			for (int j = 0; j < fragment.getAtomContainer().getAtomCount(); j++) {

				if (!(fragment.getAtomContainer().getAtom(j).getAtomTypeName()
						.equals("O.sp3"))) {
					continue;
				}

				IAtom hydroxyl_O = fragment.getAtomContainer().getAtom(j);
				// check if it is linked to two carbons, one of which has a
				// double bond O.

				IAtom carboxylC = null;
				IAtom hydroxylC = null;

				for (IAtom c : fragment.getAtomContainer().getConnectedAtomsList(
						hydroxyl_O)) {
					// the carboxyl carbon must be sp2 hybridized and connected
					// to another oxygen via double bond
					if (c.getAtomTypeName().equals("C.sp2")) {
						for (IAtom candidateOtherO : fragment.getAtomContainer()
								.getConnectedAtomsList(c)) {
							if (candidateOtherO != hydroxyl_O
									&& fragment.getAtomContainer()
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
				IBond removedBond = fragment.getAtomContainer().removeBond(
						hydroxyl_O, carboxylC);

				if (ConnectivityChecker.isConnected(fragment.getAtomContainer())) {
					hydroxylOList.add(hydroxyl_O);
					hydroxylCList.add(hydroxylC);
					carboxylCList.add(carboxylC);
					// Add this carboxyl C to the global list
					lactoneCarboxylCList.add(carboxylC);
				}
				fragment.getAtomContainer().addBond(removedBond);
			}

			if (showBonds) {
				// For testing: draw the AtomContainer with detected bonds
				// highlighted
				ArrayList<IBond> bonds = new ArrayList<IBond>();
				for (int j = 0; j < hydroxylOList.size(); j++) {
					bonds.add(fragment.getAtomContainer().getBond(
							hydroxylOList.get(j), carboxylCList.get(j)));
				}
				SmilesIO.drawMoleculeHighlightingBonds(
						fragment.getAtomContainer(),
						SmilesIO.getCleanFileName(GrapeMain.currentName
								.replace(' ', '_') + "Lactone"), bonds);
			}

			// Break the cyclic esters

			for (int j = 0; j < hydroxylOList.size(); j++) {
				numLactoneRings++;
				int index = -1;
				Fragment currentFragment = null;
				for (int k = 0; k < subfragments.size(); k++) {
					if (subfragments.get(k).getAtomContainer()
							.contains(hydroxylOList.get(j))) {
						if (!subfragments.get(k).getAtomContainer()
								.contains(carboxylCList.get(j))) {
							// TODO: proper exception handling
							System.err
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
				currentFragment.getAtomContainer().removeBond(currentHydroxylO,
						currentCarboxylC);
				Atom atom = new Atom("O");
				atom.setAtomTypeName("O.sp3");
				currentFragment.getAtomContainer().addAtom(atom);
				currentFragment.getAtomContainer().addBond(
						new Bond(carboxylCList.get(j), atom));

				List<Fragment> partitions = currentFragment
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

					if (partitions.get(0).getAtomContainer()
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
	 * Process thiazoles and thiazolines. Searches for thiazole/thizoline rings and 'undoes' the reaction.
	 * This method modifies monomerFragments and updates the field
	 * thiazoleSList.
	 */
	private void processThiazs() {
		// First, handle the situation corresponding to the five-membered S and
		// N containing ring
		// as in: CNC(=O)C1=CSC(CN)=N1
		List<Set<IAtom>> smallestRings = null;
		boolean hasThiaz = false;
		for (Fragment m : monomerFragments) {
			IAtomContainer mol = m.getAtomContainer();
			List<IAtom> carbonylCAtoms = new ArrayList<IAtom>();
			List<IAtom> thiazoneNAtoms = new ArrayList<IAtom>();
			List<IAtom> chainC1Atoms = new ArrayList<IAtom>();
			List<IAtom> chainC2Atoms = new ArrayList<IAtom>();
			List<IAtom> thiazSAtoms = new ArrayList<IAtom>();
			List<Boolean> thiazoles = new ArrayList<Boolean>();
			for (int i = 0; i < m.getAtomContainer().getAtomCount(); i++) {
				if (mol.getAtom(i).getAtomicNumber() != 16) {
					continue;
				}
				// check if this sulfur is connected to a carbon (C) double
				// bonded with a nitrogen
				IAtom thiazS = null;
				IAtom carbonylC = null;
				IAtom thiazN = null;
				IAtom chainC_1 = null;
				IAtom chainC_2 = null;
				boolean thiazole = false;

				for (IAtom c : mol.getConnectedAtomsList(mol
						.getAtom(i))) {
					if (mol.getBond(mol.getAtom(i), c).getOrder() != IBond.Order.SINGLE) {
						continue;
					}
					if (!c.getAtomTypeName().equals("C.sp2")) {
						continue;
					}
					for (IAtom n : mol.getConnectedAtomsList(c)) {
						if (n.getAtomTypeName().startsWith("C.")
								&& mol.getBond(c, n).getOrder() == IBond.Order.DOUBLE) {
							chainC_1 = c;
							chainC_2 = n;
							thiazole = true;
							continue;
						}
						if (n.getAtomTypeName().startsWith("N.")
								&& mol.getBond(c, n).getOrder() == IBond.Order.DOUBLE) {
							carbonylC = c;
							thiazN = n;
						}
					}
				}
				
				if(smallestRings == null){
					smallestRings = ChemicalUtilities.getSmallestRings(mol); //make sure if in ring it's 6+ size
				}
				
				boolean inRingSizeFive = false; // must be in a ring size 5
				
				for(Set<IAtom> ring : smallestRings){
					if(ring.contains(thiazN) && ring.size() == 5){
						inRingSizeFive = true;
						break;
					}
				}				
				if (carbonylC == null || !inRingSizeFive) {
					continue;
				}
				
				if (chainC_1 != null && chainC_2 != null) {
					IBond bond = mol.getBond(chainC_1, chainC_2);
					bond.setOrder(IBond.Order.SINGLE);
					chainC_1.setAtomTypeName("C.sp3");
					chainC_2.setAtomTypeName("C.sp3");
				}

				thiazS = mol.getAtom(i);

				carbonylCAtoms.add(carbonylC);
				thiazoneNAtoms.add(thiazN);
				chainC1Atoms.add(chainC_1);
				chainC2Atoms.add(chainC_2);
				thiazSAtoms.add(thiazS);
				thiazoles.add(thiazole);
			}
			for (int i = 0; i < thiazSAtoms.size(); i++) {
				IAtom thiazS = thiazSAtoms.get(i);
				IAtom carbonylC = carbonylCAtoms.get(i);
				IAtom thiazN = thiazoneNAtoms.get(i);
				IAtom chainC_1 = chainC1Atoms.get(i);
				IAtom chainC_2 = chainC2Atoms.get(i);
				boolean thiazole = thiazoles.get(i);

				// Cleave the C-S bond. If this alters connectivity, then it
				// wasn't part of a cycle so re-add it.
				IBond removedBond = mol.removeBond(carbonylC, thiazS);
				if (!ConnectivityChecker.isConnected(mol)) {
					mol.addBond(removedBond);
					continue;
				}

				// Turn the C=N bond to C-N
				mol.getBond(carbonylC, thiazN).setOrder(
						IBond.Order.SINGLE);

				// Add a =0 to that C
				Atom atom = new Atom("O");
				atom.setAtomTypeName("O.sp2");
				mol.addAtom(atom);
				mol.addBond(new Bond(carbonylC, atom, IBond.Order.DOUBLE));

				// Set carbonylC to sp2

				carbonylC.setAtomTypeName("C.sp2");
				thiazN.setAtomTypeName("N.amide");
				
				if(thiazole){
					thiazoleSList.add(thiazS);
				}else{
					thiazolineSList.add(thiazS);
				}

				hasThiaz = true;
			}
		}
		if (printSteps) {
			if (hasThiaz) {
				System.out.println("Processed thiaz");
			}
		}

		// Next, generally look for any cysteine whose sulfur is connected to a
		// carbon, with that S-C bond a part of a cycle.
		boolean hasCyclizedNonThiazoleCysteine = false;
		IAtomContainer template = null;
		try {
			template = SmilesIO.readSmilesTemplates("CSCC(N)C=O");
		} catch (IOException | CDKException e) {
			e.printStackTrace();
		}
		
		IBond templateBond = template.getBond(template.getAtom(0), template.getAtom(1));

		for (Fragment frag : monomerFragments) {
			List<IBond> matchingBonds = ChemicalUtilities.findMatchingBondsFromTemplate(template,
					templateBond, frag.getAtomContainer());
			for (IBond matchingBond : matchingBonds) {
				IAtom thiazS = null;
				IAtom carbonylC = null;
				if (matchingBond.getAtom(0).getAtomicNumber() == 16) {
					thiazS = matchingBond.getAtom(0);
					carbonylC = matchingBond.getAtom(1);
				} else {
					carbonylC = matchingBond.getAtom(0);
					thiazS = matchingBond.getAtom(1);
				}
				
				boolean inRingSizeFive = false; // must be in a ring size 5
				
				for(Set<IAtom> ring : smallestRings){
					if(ring.contains(thiazS) && ring.size() == 5){
						inRingSizeFive = true;
						break;
					}
				}
				if(!inRingSizeFive){
					continue;
				}
				
				frag.getAtomContainer().removeBond(matchingBond);
				if (ConnectivityChecker.isConnected(frag.getAtomContainer())) {
					// Set all bonds connecting to carbonylC as single
					for (int i = 0; i < frag.getAtomContainer().getBondCount(); i++) {
						IBond b = frag.getAtomContainer().getBond(i);
						if (b.contains(carbonylC)) {
							b.setOrder(IBond.Order.SINGLE);
						}
					}
					Atom atom = new Atom("O");
					atom.setAtomTypeName("O.sp2");
					frag.getAtomContainer().addAtom(atom);
					frag.getAtomContainer().addBond(
							new Bond(carbonylC, atom, IBond.Order.DOUBLE));
					thiazoleSList.add(thiazS);
					hasCyclizedNonThiazoleCysteine = true;
				} else {
					frag.getAtomContainer().addBond(matchingBond);
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
	 * Process oxazoles and oxazolines. Searches for oxazole/oxazoline rings and 'undoes' the reaction.
	 * Looks for CMT on the hydroxyl carbon
	 * This method modifies monomerFragments and updates the field oxazoleOList.
	 */
	private void processOxazs() {
		boolean hasOxazole = false;
		for (Fragment m : monomerFragments) {
			IAtomContainer mol = m.getAtomContainer();
			IAtomContainer[] oxazoleTemplates = new IAtomContainer[2];
			
			try {
				oxazoleTemplates[0] = SmilesIO.readSmilesTemplates("O1C=CN=C1");
				oxazoleTemplates[1] = SmilesIO.readSmilesTemplates("C1CN=CO1");
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
				List<List<IBond>> matchingBondsList = ChemicalUtilities.findMatchingBondsFromTemplate(oxazoleTemplates[i], templateBonds.get(i), m.getAtomContainer());
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
					m.getAtomContainer().removeBond(matchingBonds.get(2));
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
					
					//check if the oxazole carbon is methylated
					for(IAtom atom : mol.getConnectedAtomsList(oxazoleO)){
						if(mol.getConnectedBondsCount(atom) == 3){
							for(IAtom atom2 : mol.getConnectedAtomsList(atom)){
								if(atom2.getAtomicNumber() == 6 
										&& ChemicalUtilities.getConnectedAtomsCountNonHydrogen(mol, atom2) == 1){
									mol.removeAtomAndConnectedElectronContainers(atom2);
									methylatedCarbonsList.add(atom);
								}
							}
						}
					}
					
					Atom ketoneO = new Atom("O");
					ketoneO.setAtomTypeName("O.sp2");
					newKetoneC.setAtomTypeName("C.sp2");
					mol.addAtom(ketoneO);
					mol.addBond(new Bond(newKetoneC, ketoneO, IBond.Order.DOUBLE));
					
					hasOxazole = true;
					if(i == 0) { //oxazole
						oxazoleOList.add(oxazoleO);
					}else { //oxazoline
						oxazolineOList.add(oxazoleO);
					}
				}
			}
			
			}
		if (printSteps) {
			if (hasOxazole) {
				System.out.println("Processed oxazole");
			}
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

			for (int i = 0; i < fragment.getAtomContainer().getAtomCount(); i++) {

				IAtom mainC = null;
				IAtom main_O = null;
				IAtom n_1 = null;
				IAtom c_1 = null;
				IAtom n_2 = null;
				IAtom c_2 = null;

				if (!fragment.getAtomContainer().getAtom(i).getAtomTypeName()
						.equals("C.sp2")) {
					continue;
				}
				mainC = fragment.getAtomContainer().getAtom(i);
				for (int j = 0; j < fragment.getAtomContainer()
						.getConnectedAtomsCount(mainC); j++) {
					IAtom connectedAtom = fragment.getAtomContainer()
							.getConnectedAtomsList(mainC).get(j);

					if (connectedAtom.getAtomTypeName().startsWith("O.")
							&& fragment.getAtomContainer()
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

				for (int j = 0; j < fragment.getAtomContainer()
						.getConnectedAtomsCount(n_1); j++) {
					IAtom connectedAtom = fragment.getAtomContainer()
							.getConnectedAtomsList(n_1).get(j);
					if (connectedAtom == mainC) {
						continue;
					} else if (connectedAtom.getAtomTypeName().startsWith("C.")) {
						c_1 = connectedAtom;
						break;
					}
				}
				for (int j = 0; j < fragment.getAtomContainer()
						.getConnectedAtomsCount(n_2); j++) {
					IAtom connectedAtom = fragment.getAtomContainer()
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

				removedBonds.add(fragment.getAtomContainer().getBond(mainC, main_O));
				removedBonds.add(fragment.getAtomContainer().getBond(mainC, n_1));
				removedBonds.add(fragment.getAtomContainer().getBond(mainC, n_2));

				fragment.getAtomContainer().removeBond(
						fragment.getAtomContainer().getBond(mainC, main_O));
				fragment.getAtomContainer().removeAtom(main_O);
				fragment.getAtomContainer().removeBond(
						fragment.getAtomContainer().getBond(mainC, n_1));
				fragment.getAtomContainer().removeBond(
						fragment.getAtomContainer().getBond(mainC, n_2));
				fragment.getAtomContainer().removeAtom(mainC);

				IAtomContainerSet fragments = ConnectivityChecker
						.partitionIntoMolecules(fragment.getAtomContainer());
				if (fragments.getAtomContainerCount() != 2) {
					// We must have had a ring, so re-add what we removed.
					fragment.getAtomContainer().addAtom(main_O);
					fragment.getAtomContainer().addAtom(mainC);
					for (IBond bond : removedBonds) {
						fragment.getAtomContainer().addBond(bond);
					}
					continue;
				}

				if (fragments.getAtomContainerCount() > 1) { //add bridge tailor
					monomerFragments.remove(fragment);
					for (Fragment newFragment : fragment.partitionIntoMonomerFragments()) {
						monomerFragments.add(newFragment);
					}
				}
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

		IAtomContainer cephemTemplate = null;
		try {
			cephemTemplate = SmilesIO.readSmilesTemplates("C1CN2C=CCSC12");
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
					cephemTemplate, cephemTemplateBonds, frag.getAtomContainer());
			for (int i = 0; i < matchingBonds.size(); i++) {
				matchingBonds.get(i).get(0).setOrder(IBond.Order.SINGLE);
				matchingBonds.get(i).get(0).getAtom(0).setAtomTypeName("C.sp3");
				matchingBonds.get(i).get(0).getAtom(1).setAtomTypeName("C.sp3");
				frag.getAtomContainer().removeBond(matchingBonds.get(i).get(1));
				// Create new bond between the sulfur and the carbon in the
				// double bond to form a five-membered ring
				IAtom carbonPreviouslyConnectedToSulfur = null;
				IAtom carbonNewlyConnectingToSulfur = null;
				IAtom cephemSulfur = null;

				if (matchingBonds.get(i).get(1).getAtom(0).getAtomicNumber() == 6) {
					carbonPreviouslyConnectedToSulfur = matchingBonds.get(i)
							.get(1).getAtom(0);
					cephemSulfur = matchingBonds.get(i).get(1).getAtom(1);
				} else {
					carbonPreviouslyConnectedToSulfur = matchingBonds.get(i)
							.get(1).getAtom(1);
					cephemSulfur = matchingBonds.get(i).get(1).getAtom(0);
				}
				for (IAtom candidateNewlyConnectedCarbon : frag.getAtomContainer()
						.getConnectedAtomsList(
								carbonPreviouslyConnectedToSulfur)) {
					// If this candidate carbon is connected to one of the
					// previously double bonded carbons
					if (frag.getAtomContainer()
							.getConnectedAtomsList(
									candidateNewlyConnectedCarbon)
							.contains(matchingBonds.get(i).get(0).getAtom(0))
							|| frag.getAtomContainer()
									.getConnectedAtomsList(
											candidateNewlyConnectedCarbon)
									.contains(
											matchingBonds.get(i).get(0)
													.getAtom(1))) {
						carbonNewlyConnectingToSulfur = candidateNewlyConnectedCarbon;
					}
				}
				IBond newCephemBond = new Bond(cephemSulfur,
						carbonNewlyConnectingToSulfur);
				frag.getAtomContainer().addBond(newCephemBond);
				frag.setBetaLactamS(cephemSulfur);
			}
		}

		IAtomContainer sulfurBetaLactamTemplate = null;
		try {
			sulfurBetaLactamTemplate = SmilesIO.readSmilesTemplates("CSC1CC(=O)N1");
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
		IAtomContainer clavualnicAcidLikeTemplate = null;
		try {
			clavualnicAcidLikeTemplate = SmilesIO.readSmilesTemplates("COC1CC(=O)N1");
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
		
		
		IAtomContainer gammaLactamBetaLactoneTemplate = null;
		try {
			gammaLactamBetaLactoneTemplate = SmilesIO.readSmilesTemplates("O=C1OC2CC(=O)NC12");
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
					sulfurBetaLactamTemplate, sulfurBetaLactamTemplateBonds, frag.getAtomContainer());
			if(matchingBonds.size() < 1){
				matchingBonds = ChemicalUtilities.findMatchingBondsFromTemplate(
						gammaLactamBetaLactoneTemplate, gammaLactamBetaLactoneTemplateBonds, frag.getAtomContainer());
			}
			if(matchingBonds.size() < 1){
				matchingBonds = ChemicalUtilities.findMatchingBondsFromTemplate(
						clavualnicAcidLikeTemplate, clavualnicAcidLikeTemplateBonds, frag.getAtomContainer());
			}
			for (int i = 0; i < matchingBonds.size(); i++) {
				frag.getAtomContainer().removeBond(matchingBonds.get(i).get(0));
				frag.getAtomContainer().removeBond(matchingBonds.get(i).get(1));
				// If the AtomContainer is no longer connected, then re-add these
				// bonds.
				if (!ConnectivityChecker.isConnected(frag.getAtomContainer())) {
					frag.getAtomContainer().addBond(matchingBonds.get(i).get(0));
					frag.getAtomContainer().addBond(matchingBonds.get(i).get(1));
				} else {
					foundMatch = true;
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
		IAtomContainer template = null;
		try {
			template = SmilesIO.readSmilesTemplates("CSC(C)=O");
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
					templateBond, currentFragment.getAtomContainer());

			for (IBond currentBond : matchedBonds) {
				IAtom thioesterCarbon = null;
				if (currentBond.getAtom(0).getAtomicNumber() == 6) {
					thioesterCarbon = currentBond.getAtom(0);
				} else {
					thioesterCarbon = currentBond.getAtom(1);
				}
				currentFragment.addThioEsterAtom(currentBond.getAtom(0));
				currentFragment.addThioEsterAtom(currentBond.getAtom(1));
				currentFragment.getAtomContainer().removeBond(currentBond);
				Atom carboxylO = new Atom("O");
				carboxylO.setAtomTypeName("O.sp3");
				currentFragment.getAtomContainer().addAtom(carboxylO);
				Bond newBond = new Bond(thioesterCarbon, carboxylO);
				currentFragment.getAtomContainer().addBond(newBond);

				List<Fragment> fragments = currentFragment
						.partitionIntoMonomerFragments();

				if (fragments.size() == 2) {
					Fragment sFrag = null;
					Fragment cFrag = null;
					if (fragments.get(0).getAtomContainer()
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
			IAtomContainer currentAtomContainer = fragment.getAtomContainer();
			IAtomContainer template = null;
			try {
				template = SmilesIO.readSmilesTemplates("C1=CC=C(C=C1)C1=CC=CC=C1");
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
					templateBond, currentAtomContainer);

			// At this point, conditions have been met. Break the bond if it
			// leads to one piece.
			for (IBond currentBond : matchedBonds) {
				currentAtomContainer.removeBond(currentBond);
				if (ConnectivityChecker.partitionIntoMolecules(currentAtomContainer)
						.getAtomContainerCount() > 1) {
					currentAtomContainer.addBond(currentBond);
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
			IAtomContainer currentAtomContainer = fragment.getAtomContainer();
			//prepare sugar templates, and find sugar bonds
			IAtomContainer[] sugarTemplates = new IAtomContainer[4];
			IBond[] sugarTemplateBonds = new IBond[4];
			try {
				sugarTemplates[0] = SmilesIO.readSmilesTemplates("COC1CCCCO1"); 
				sugarTemplates[1] = SmilesIO.readSmilesTemplates("COC1CCCO1");
				sugarTemplates[2] = SmilesIO.readSmilesTemplates("CSC1CCCO1");
				sugarTemplates[3] = SmilesIO.readSmilesTemplates("CSC1CCCCO1");

			} catch (Exception e) {
				e.printStackTrace();
			}

			// Set template bonds
			sugarTemplateBonds[0] = sugarTemplates[0].getBond(sugarTemplates[0].getAtom(2), sugarTemplates[0].getAtom(1));
			sugarTemplateBonds[1] = sugarTemplates[1].getBond(sugarTemplates[1].getAtom(2), sugarTemplates[1].getAtom(1));
			sugarTemplateBonds[2] = sugarTemplates[2].getBond(sugarTemplates[2].getAtom(2), sugarTemplates[2].getAtom(1));
			sugarTemplateBonds[3] = sugarTemplates[3].getBond(sugarTemplates[3].getAtom(2), sugarTemplates[3].getAtom(1));
	

			// Find all sugar bonds
			
			for(int i = 0; i < sugarTemplateBonds.length; i++){
				IAtomContainer sugarTemplate = sugarTemplates[i];
				IBond sugarTemplateBond = sugarTemplateBonds[i];
				List<IBond> sugarBonds = ChemicalUtilities.findMatchingBondsFromTemplate(
						sugarTemplate, sugarTemplateBond, currentAtomContainer);

				for (IBond sugarBond : sugarBonds) {
					Fragment subfragment = null;
					for (Fragment possibleSubfragment : subfragments) {
						if (possibleSubfragment.getAtomContainer().contains(sugarBond)) {
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
					
					//add check where the oxygen cannot be part of a small cycle

					subfragment.getAtomContainer().removeBond(sugarBond);
					// Check if this led to two pieces
					if (ConnectivityChecker.isConnected(subfragment.getAtomContainer())) {
						subfragment.getAtomContainer().addBond(sugarBond);
						continue;
					}

					sugarCarbons.add(bondCarbon);
					IAtom hydroxylO = new Atom("O");
					hydroxylO.setAtomTypeName("O.sp3");
					subfragment.getAtomContainer().addAtom(hydroxylO);
					subfragment.getAtomContainer().addBond(
							new Bond(hydroxylO, bondCarbon));

					List<Fragment> partitions = subfragment
							.partitionIntoMonomerFragments();

					// There must be two partitions
					Fragment sugarPart;
					Fragment nonSugarPart;
					if (partitions.get(0).getAtomContainer().contains(bondCarbon)) {
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
		IAtomContainer template = null;
		IBond templateBond = null;
		try {
			template = SmilesIO.readSmilesTemplates("COC=O");
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
					templateBond, frag.getAtomContainer());
			for (IBond esterBond : esterBonds) {
				foundEster = true;
				IAtom esterC = null;
				if (esterBond.getAtom(0).getAtomicNumber() == 6) {
					esterC = esterBond.getAtom(0);
				} else {
					esterC = esterBond.getAtom(1);
				}
				// Remove bond
				frag.getAtomContainer().removeBond(esterBond);
				IAtom carboxylO = new Atom("O");
				carboxylO.setAtomTypeName("O.sp3");
				frag.getAtomContainer().addAtom(carboxylO);
				frag.getAtomContainer().addBond(new Bond(esterC, carboxylO));
			}
			if (!ConnectivityChecker.isConnected(frag.getAtomContainer())) {
				List<Fragment> subfragments = frag
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
		IAtomContainer template = null;
		IBond templateBond = null;
		try {
			template = SmilesIO.readSmilesTemplates("COS");
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
					templateBond, frag.getAtomContainer());
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
				frag.getAtomContainer().removeBond(sulfateBond);
				// Check that the bond removal has led to two pieces, of of
				// which has exactly one sulfur and four oxygens
				IAtomContainerSet partitions = ConnectivityChecker
						.partitionIntoMolecules(frag.getAtomContainer());
				if (partitions.getAtomContainerCount() != 2) {
					frag.getAtomContainer().addBond(sulfateBond);
					continue;
				}
				IAtomContainer sulfatePiece = null;
				IAtomContainer nonSulfatePiece = null;
				if (partitions.getAtomContainer(0).contains(sulfateO)) {
					sulfatePiece = partitions.getAtomContainer(0);
					nonSulfatePiece = partitions.getAtomContainer(1);
				} else {
					sulfatePiece = partitions.getAtomContainer(1);
					nonSulfatePiece = partitions.getAtomContainer(0);
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
					frag.getAtomContainer().addBond(sulfateBond);
					continue;
				}
				hasSulfate = true;
				// At this point, we conclude that this is indeed a sulfate
				// piece. Set this fragment AtomContainer the non-sulfur piece, set
				// enum.
				IAtom hydroxylO = new Atom("O");
				hydroxylO.setAtomTypeName("O.sp3");
				nonSulfatePiece.addAtom(hydroxylO);
				nonSulfatePiece.addBond(new Bond(connectingCarbon, hydroxylO));
				frag.setAtomContainer(nonSulfatePiece);
				frag.addTailoringDomain(TailoringDomainEnums.SULFOTRANSFERASE);
			}
			
			for (Fragment m : monomerFragments) { // check for *S(O)(O)O
				List<IAtom> sulfurAtoms = new ArrayList<IAtom>();
				for (int i = 0; i < m.getAtomContainer().getAtomCount(); i++) {
					if (m.getAtomContainer().getAtom(i).getAtomicNumber() == 16) {
						sulfurAtoms.add(m.getAtomContainer().getAtom(i));
					}
				}
				if (sulfurAtoms.size() == 0) {
					continue;
				}
				List<IAtom> oxygens = new ArrayList<IAtom>(); 
				for(IAtom sulfur : sulfurAtoms){
					for(IAtom connectedAtom : m.getAtomContainer().getConnectedAtomsList(sulfur)){
						if(connectedAtom.getAtomicNumber() == 8 && m.getAtomContainer().getConnectedAtomsList(connectedAtom).size() == 1){
							oxygens.add(connectedAtom);
						}
					}
					if(oxygens.size() == 3){
						m.getAtomContainer().removeAtomAndConnectedElectronContainers(sulfur);
						for(IAtom oxygen : oxygens){
							m.getAtomContainer().removeAtomAndConnectedElectronContainers(oxygen);
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
			IAtomContainer currentAtomContainer = fragment.getAtomContainer();

			IAtomContainer[] templates = new IAtomContainer[2];
			IBond[] templateBonds = new IBond[2];
			try {
				// Benzene rings connected to each other
				templates[0] = SmilesIO.readSmilesTemplates("OC1=CC=CC=C1OC1=CC=CC=C1");
				templates[1] = SmilesIO
						.readSmilesTemplates("OC1=C(OC2=CC=CC=C2)C=CC=C1");
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
						templates[x], templateBonds[x], currentAtomContainer);
				for (IBond b : bondsToBreak) {
					if (currentAtomContainer.contains(b)) {
						currentAtomContainer.removeBond(b);
						numBonds++;
					}
				}
			}

			IAtomContainerSet fragments = ConnectivityChecker
					.partitionIntoMolecules(currentAtomContainer);
			if (fragments.getAtomContainerCount() > 1) {
				monomerFragments.remove(fragment);
				for (int i = 0; i < fragments.getAtomContainerCount(); i++) {
					monomerFragments.add(new Fragment(fragments
							.getAtomContainer(i)));
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
			IAtomContainer currentAtomContainer = fragment.getAtomContainer();

			IAtomContainer template = null;
			try {
				template = SmilesIO.readSmilesTemplates("CSSC");

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
					templateBond, currentAtomContainer);
			
			if(bondsToBreak.size() > 0){
				for (IBond b : bondsToBreak) {
					if (currentAtomContainer.contains(b)) {
						currentAtomContainer.removeBond(b);
						numBonds++;
					}
				}
				IAtomContainerSet fragments = ConnectivityChecker
						.partitionIntoMolecules(currentAtomContainer);
				if (fragments.getAtomContainerCount() > 1) { //add bridge tailor
					monomerFragments.remove(fragment);
					for (int i = 0; i < fragments.getAtomContainerCount(); i++) {
						monomerFragments.add(new Fragment(fragments
								.getAtomContainer(i)));
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
			
			//guess amino N if it wasn't cleaved (ie N terminus)
			IAtomContainer mol = m.getAtomContainer();
			if(m.getAminoNs().size() == 0 && m.getAminoCs().size() == 1){
				IAtom aminoC =  m.getAminoCs().get(0);
				for(IAtom atom : mol.getConnectedAtomsList(aminoC)){
					if(atom.getAtomicNumber() == 6){
						for(IAtom atom2 : mol.getConnectedAtomsList(atom)){
							if(atom2.getAtomicNumber() == 7){
								m.addAminoN(atom2);
							}
						}
					}
				}
			}
			
			for (int i = 0; i < m.getAtomContainer().getBondCount(); i++) {
				IBond currentBond = m.getAtomContainer().getBond(i);
				
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
							&& m.getAtomContainer().getConnectedAtomsCount(
									candidateMethylC) == 1) {
						m.getAtomContainer().removeBond(currentBond);
						m.getAtomContainer().removeAtom(candidateMethylC);
						m.getTailoringDomains().add(
								TailoringDomainEnums.N_METHYLTRANSFERASE);
						numNMethylations++;
						i --;
						continue;
					}
				}
				// C methylations
				if (m.getAminoNs().size() > 0) {
					IAtom candidateMethylC = null;
					IAtom candidateAlphaCarbon = null;
					if (currentBond.getAtom(0).getAtomicNumber() == 6
							&& currentBond.getAtom(1).getAtomicNumber() == 6) {
						if (m.getAtomContainer()
										.getConnectedAtomsList(
												currentBond.getAtom(0))
										.contains(m.getAminoNs().get(0))) {
							candidateAlphaCarbon = currentBond.getAtom(0);
							candidateMethylC = currentBond.getAtom(1);
						}
						if (m.getAtomContainer()
										.getConnectedAtomsList(
												currentBond.getAtom(1))
										.contains(m.getAminoNs().get(0))) {
							candidateAlphaCarbon = currentBond.getAtom(1);
							candidateMethylC = currentBond.getAtom(0);
						}
						if (m.getAtomContainer().getConnectedAtomsCount(
								candidateAlphaCarbon) == 4
								&& m.getAtomContainer().getConnectedAtomsCount(
										candidateMethylC) == 1) {
							m.getAtomContainer().removeBond(currentBond);
							m.getAtomContainer().removeAtom(candidateMethylC);
							m.getTailoringDomains().add(
									TailoringDomainEnums.C_METHYLTRANSFERASE);
							numCMethylations++;
							i --;
							continue;
						}
					}
				}
				if (m.getAminoCs().size() > 0) {
					IAtom candidateMethylC = null;
					IAtom candidateAlphaCarbon = null;
					if (currentBond.getAtom(0).getAtomicNumber() == 6
							&& currentBond.getAtom(1).getAtomicNumber() == 6) {
						if (m.getAtomContainer()
								.getConnectedAtomsList(currentBond.getAtom(0))
								.contains(m.getAminoCs().get(0))
								) {
							candidateAlphaCarbon = currentBond.getAtom(0);
							candidateMethylC = currentBond.getAtom(1);
						}
						if (m.getAtomContainer()
								.getConnectedAtomsList(currentBond.getAtom(1))
								.contains(m.getAminoCs().get(0))
								) {
							candidateAlphaCarbon = currentBond.getAtom(1);
							candidateMethylC = currentBond.getAtom(0);
						}
						if (m.getAtomContainer().getConnectedAtomsCount(
								candidateAlphaCarbon) == 4
								&& m.getAtomContainer().getConnectedAtomsCount(
										candidateMethylC) == 1) {
							m.getAtomContainer().removeBond(currentBond);
							m.getAtomContainer().removeAtom(candidateMethylC);
							m.getTailoringDomains().add(
									TailoringDomainEnums.C_METHYLTRANSFERASE);
							numCMethylations++;
							i --;
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
					if (m.getAtomContainer().getConnectedAtomsList(candidateO)
							.contains(m.getAminoCs().get(0))
							&& m.getAtomContainer().getConnectedAtomsCount(
									candidateMethylC) == 1) {
						candidateO.setAtomTypeName("O.sp3");
						m.getAtomContainer().getBond(candidateO, m.getAminoCs().get(0))
								.setOrder(IBond.Order.DOUBLE);
						m.getAtomContainer().removeBond(currentBond);
						m.getAtomContainer().removeAtom(candidateMethylC);
						m.getTailoringDomains().add(
								TailoringDomainEnums.O_METHYLTRANSFERASE);
						numOMethylations++;
						i --;
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
			for (int i = 0; i < m.getAtomContainer().getBondCount(); i++) {
				IBond currentBond = m.getAtomContainer().getBond(i);
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
					if (m.getAtomContainer().getConnectedAtomsCount(
							candidateAlphaCarbon) == 4
							&& m.getAtomContainer().getConnectedAtomsCount(
									candidateHydroxylO) == 1) {
						m.getAtomContainer().removeBond(currentBond);
						m.getAtomContainer().removeAtom(candidateHydroxylO);
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
					if (m.getAtomContainer().getConnectedAtomsCount(
							candidateAlphaCarbon) == 4
							&& m.getAtomContainer().getConnectedAtomsCount(
									candidateHydroxylO) == 1) {
						m.getAtomContainer().removeBond(currentBond);
						m.getAtomContainer().removeAtom(candidateHydroxylO);
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
	private void processHalogenations() {
		boolean hasHalogen = false;
		for (Fragment m : monomerFragments) {
			List<IAtom> halogenAtoms = new ArrayList<IAtom>();
			for (int i = 0; i < m.getAtomContainer().getAtomCount(); i++) {
				if (m.getAtomContainer().getAtom(i).getAtomicNumber() == 17
						|| m.getAtomContainer().getAtom(i).getAtomicNumber() == 35) {
					halogenAtoms.add(m.getAtomContainer().getAtom(i));
				}
			}
			if (halogenAtoms.size() == 0) {
				continue;
			}
			m.getTailoringDomains().add(TailoringDomainEnums.HALOGENATION);
			
			if(halogenAtoms.size() >= 3) {
				for(IAtom chlorine : halogenAtoms){
					IAtom connectedToChlorine = null;
					int numConnectedToSameAtom = 0;
					if(m.getAtomContainer().getConnectedAtomsCount(chlorine) != 1) continue;
					connectedToChlorine = m.getAtomContainer().getConnectedAtomsList(chlorine).get(0);
					for(IAtom atom : m.getAtomContainer().getConnectedAtomsList(connectedToChlorine)){
						if(atom.getAtomicNumber() == 17
								|| atom.getAtomicNumber() == 35){
							numConnectedToSameAtom ++;
						}
					}
					if(numConnectedToSameAtom == 3){
						Atom firstO = new Atom("O");
						firstO.setAtomTypeName("O.sp3");
						Atom secondO = new Atom("O");
						secondO.setAtomTypeName("O.sp2");
						Atom carbonylCarbon = new Atom("C");
						carbonylCarbon.setAtomTypeName("C.sp2");
						m.getAtomContainer().addAtom(firstO);
						m.getAtomContainer().addAtom(secondO);
						m.getAtomContainer().addAtom(carbonylCarbon);
						m.getAtomContainer().addBond(new Bond(carbonylCarbon, firstO, IBond.Order.SINGLE));
						m.getAtomContainer().addBond(new Bond(carbonylCarbon, secondO, IBond.Order.DOUBLE));
						m.getAtomContainer().addBond(new Bond(carbonylCarbon, connectedToChlorine, IBond.Order.SINGLE));
						break;
					}
				}
			}
			for (IAtom halogen : halogenAtoms) {
				// Check that there is exactly one connection
				if (m.getAtomContainer().getConnectedAtomsCount(halogen) != 1) {
					continue;
				}
				hasHalogen = true;
				m.getAtomContainer().removeBond(
						m.getAtomContainer().getBond(
								halogen,
								m.getAtomContainer().getConnectedAtomsList(halogen)
										.get(0)));
				m.getAtomContainer().removeAtom(halogen);
			}
		}
		if (printSteps) {
			if (hasHalogen) {
				System.out.println("Processed halogenation");
			}
		}
	}
	
	private void processSecondaryAmineExtensions() {
	
		Queue<Fragment> q = new LinkedList<Fragment>();
		for (Fragment m : monomerFragments) {
			q.add(m);
		}
		
		while (!q.isEmpty()) {
			Fragment fragment = q.poll();
			IAtomContainer mol = fragment.getAtomContainer();
			Set<IBond> bondsToBreak = new HashSet<IBond>();
			List<Set<IAtom>> smallestRings = null;
			//find extended secondary amines
			for(IAtom atom : mol.atoms()){
				
				if(smallestRings == null){
					smallestRings = ChemicalUtilities.getSmallestRings(mol); //make sure if in ring it's 6+ size
				}

				boolean inSmallRing = false; // cannot be in a ring < 7 (then is potentially a thiazol)
				for(Set<IAtom> ring : smallestRings){
					if(ring.contains(atom) && ring.size() < 7){
						inSmallRing = true;
						break;
					}
				}
				if(inSmallRing == true){
					continue;
				}
				
				//must be a nitrogen
				boolean possible = true;
				if(atom.getAtomicNumber() != 7){
					continue;
				}
				List<IAtom> connectedCarbons = new ArrayList<IAtom>();
				for(IAtom connectedAtom : mol.getConnectedAtomsList(atom)){
					//can only be connected to carbons
					if(connectedAtom.getAtomicNumber() == 6){
						//carbons cannot be a carbonyl
						for(IAtom connectedToCarbon : mol.getConnectedAtomsList(connectedAtom)){
							if(connectedToCarbon.getAtomicNumber() == 8 && mol.getBond(connectedAtom, connectedToCarbon).getOrder().equals(Order.DOUBLE)){
								possible = false;
								break;
							}
						}
						if(!possible){
							break;
						}else{
							connectedCarbons.add(connectedAtom);
						}
					}
				}
				//make sure it's still possible and is connected to 2 carbons
				if(!possible || connectedCarbons.size() != 2){
					continue;
				}
				IAtom aminoCarbon = null;
				for(IAtom carbon : connectedCarbons){
					//determine which is connected to a COOH
					int numConnections = 0;
					for(IAtom connectedCarbon : mol.getConnectedAtomsList(carbon)){
						//ensure it's a carbon
						if(connectedCarbon.getAtomicNumber() != 1){
							numConnections ++;
						}
						
						if(connectedCarbon.getAtomicNumber() != 6){
							continue;
						}
						
						boolean connectedToOH = false; // OH
						boolean connectedToDO = false; // =O
						for(IAtom connectedAtom : mol.getConnectedAtomsList(connectedCarbon)){
							if(connectedAtom.getAtomicNumber() == 8){
								if(mol.getBond(connectedCarbon, connectedAtom).getOrder().equals(Order.DOUBLE)){
									connectedToDO = true;
								}else if(mol.getBond(connectedCarbon, connectedAtom).getOrder().equals(Order.SINGLE)){
									connectedToOH = true;
								}
							}
						}
						
						if(connectedToDO && connectedToOH){
							aminoCarbon = carbon;
							break;
						}						
					}
					if(numConnections < 2){
						possible = false;
						break;
					}
					if(aminoCarbon != null){
						break;
					}
					
				}
				if(aminoCarbon != null && possible){
					for(IAtom carbon : connectedCarbons){
						if(!carbon.equals(aminoCarbon)){
							bondsToBreak.add(mol.getBond(carbon, atom));
						}
					}
					
				}
			}
			//break the secondary amines
			for(IBond bondToBreak : bondsToBreak){
				mol.removeBond(bondToBreak);
				
				//get the carbon
				IAtom carbon = null;
				for(IAtom atom : bondToBreak.atoms()){
					
					if(atom.getAtomicNumber() == 6){
						carbon = atom;
						break;
					}					
				}
				
				for(IAtom atom : mol.getConnectedAtomsList(carbon)){
					mol.removeBond(atom, carbon);
					IBond bond = new Bond(atom, carbon, Order.SINGLE);
					mol.addBond(bond);
				}
				
				//add carboxylic acid to the carbon
				Atom firstO = new Atom("O");
				firstO.setAtomTypeName("O.sp3");
				Atom secondO = new Atom("O");
				secondO.setAtomTypeName("O.sp2");
				mol.addAtom(firstO);
				mol.addAtom(secondO);
				mol.addAtom(carbon);
				mol.addBond(new Bond(carbon, firstO, IBond.Order.SINGLE));
				mol.addBond(new Bond(carbon, secondO, IBond.Order.DOUBLE));
				
			}
			
			//set new fragments if not connected
			IAtomContainerSet fragments = ConnectivityChecker
					.partitionIntoMolecules(mol);
			if (fragments.getAtomContainerCount() > 1) { //add bridge tailor
				monomerFragments.remove(fragment);
				for (Fragment newFragment : fragment.partitionIntoMonomerFragments()) {
					monomerFragments.add(newFragment);
				}
			}		
		}
	}

	/**
	 * Find and annotate epoxide groups
	 */
	private void findEpoxides() {
		boolean hasEpoxide = false;
		for (Fragment m : monomerFragments) {
			for (int i = 0; i < m.getAtomContainer().getAtomCount(); i++) {
				IAtom o = m.getAtomContainer().getAtom(i);
				if (o.getAtomicNumber() != 8)
					continue;
				if (m.getAtomContainer().getConnectedAtomsCount(o) != 2)
					continue;
				IAtom c1 = m.getAtomContainer().getConnectedAtomsList(o).get(0);
				IAtom c2 = m.getAtomContainer().getConnectedAtomsList(o).get(1);
				if (c1.getAtomicNumber() != 6 || c2.getAtomicNumber() != 6)
					continue;
				if (!m.getAtomContainer().getConnectedAtomsList(c1).contains(c2))
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
