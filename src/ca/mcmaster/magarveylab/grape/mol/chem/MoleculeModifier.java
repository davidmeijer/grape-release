package ca.mcmaster.magarveylab.grape.mol.chem;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import org.openscience.cdk.Atom;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.PolyKetideDomainEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.TailoringDomainEnums;
import ca.mcmaster.magarveylab.grape.mol.chem.Fragment.FragmentType;
import ca.mcmaster.magarveylab.grape.mol.chem.modifications.BridgeBreakages;
import ca.mcmaster.magarveylab.grape.mol.chem.modifications.SubstrateBreakages;
import ca.mcmaster.magarveylab.grape.mol.chem.modifications.UniqueSubstructureBreakages;
import ca.mcmaster.magarveylab.grape.pk.modules.PKsubstrate;
import ca.mcmaster.magarveylab.grape.util.ChemicalUtilities;
import ca.mcmaster.magarveylab.grape.util.io.SmilesIO;

/**
 * Perform modifications to NRPs, breaking them into constitudent monomer units
 * 
 * @author gmchen cDejong
 */
public class MoleculeModifier {

	public final static boolean imageDump = false;
	public final static boolean showBonds = false;
	public final static boolean printSteps = false;

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
	private String name;
	private List<Fragment> fragments;

	//Builders and common objects used in modifications
	IChemObjectBuilder cob = DefaultChemObjectBuilder.getInstance();
	UniversalIsomorphismTester uit = new UniversalIsomorphismTester();

	
	/**
	 * Constructor for this NRPModifier. Typical usage: construct an NRPModifier
	 * from an original NRP IAtomContainer, and get the output from method
	 * getMonomerFragments.
	 * 
	 * @param originalAtomContainer
	 * @param chemicalAbstraction
	 * @param name
	 */
	public MoleculeModifier(IAtomContainer originalAtomContainer, ChemicalAbstraction chemicalAbstraction) {
		lactoneCarboxylCList = new ArrayList<IAtom>();
		methylatedCarbonsList = new ArrayList<IAtom>();
		oxazoleOList = new ArrayList<IAtom>();
		thiazoleSList = new ArrayList<IAtom>();
		thiazolineSList = new ArrayList<IAtom>();
		oxazoleOList = new ArrayList<IAtom>();
		oxazolineOList = new ArrayList<IAtom>();
		sugarCarbons = new ArrayList<IAtom>();
		fragments = new ArrayList<Fragment>();
		fragments.add(new Fragment(originalAtomContainer));
		this.name = chemicalAbstraction.getName();
	}
	
	public void addThiazoleS(IAtom atom){
		thiazoleSList.add(atom);
	}
	
	public void addThiazolineS(IAtom atom){
		thiazolineSList.add(atom);
	}
	
	public void addMethylatedCarbon(IAtom atom){
		methylatedCarbonsList.add(atom);
	}
	
	public void addOxazoleO(IAtom atom){
		oxazoleOList.add(atom);
	}
	
	public void addOxazolineO(IAtom atom){
		oxazolineOList.add(atom);
	}
	
	public void addSugarCarbon(IAtom atom){
		sugarCarbons.add(atom);
	}
	
	public void addLactoneCarboxylicC(IAtom atom){
		lactoneCarboxylCList.add(atom);
	}
	
	public String getName(){
		return name;
	}

	/**
	 * Get the monomer fragments from the modifications performed by the
	 * NRPModifier
	 * 
	 * @return
	 */
	public List<Fragment> getFragments(){
		return fragments;
	}
	
	public UniversalIsomorphismTester getUIT(){
		return uit;
	}
	
	public IChemObjectBuilder getCOB(){
		return cob;
	}


	/**
	 * Perform all NRP modifications in the intended order
	 */
	public void performAllNrpModifications() {
		BridgeBreakages.process(this);
		UniqueSubstructureBreakages.process(this);
		SubstrateBreakages.process(this);
		findEpoxides();
		processPeptideEpoxiKetones();

		// Update tailoring fields fields
		for (Fragment m : fragments) {
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
			for(IAtom a : m.getThioEsterAtoms()){
				if (m.getAtomContainer().contains(a)) {
					m.addTailoringDomain(TailoringDomainEnums.THIOL_ESTHER);
				}
			}
			for(IAtom a : m.getLactamAtoms()){
				if (m.getAtomContainer().contains(a)) {
					m.addTailoringDomain(TailoringDomainEnums.LACTAM);
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
		for (Fragment m : fragments) {
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
				int indexToReplace = fragments.indexOf(fragment);
				fragments.remove(indexToReplace);
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
					fragments.add(indexToReplace + addCount, fragmentToAdd);
					addCount += 1;
				}
			}
		}
	}

	/**
	 * Find and annotate epoxide groups
	 */
	private void findEpoxides() {
		boolean hasEpoxide = false;
		for (Fragment m : fragments) {
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
