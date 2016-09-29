package ca.mcmaster.magarveylab.grape.pk.chem;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;

import ca.mcmaster.magarveylab.grape.enums.DomainEnums.PolyKetideDomainEnums;
import ca.mcmaster.magarveylab.grape.pk.modules.PKsubstrate;
import ca.mcmaster.magarveylab.grape.util.ChemicalUtilities;
import ca.mcmaster.magarveylab.grape.util.io.SmilesIO;

/**
 * @author cDejong
 *
 */
public class PolyketideModulePredictor {
	private IAtomContainer pk; //potential polyketide to be analyzed
	private IAtom startCarbon; //the start of the polyketide (actually the biosynthetic end)
	private IAtom endCarbon; //the end of the polyketide (actually the biosynthetic start)
	private Integer pkType = null;
	private List<PKsubstrate> loadingUnits = new ArrayList<PKsubstrate>(); //the type of molecule [acetate] loaded
	private List<List<PolyKetideDomainEnums>> allDomains = new ArrayList<List<PolyKetideDomainEnums>>(); //the resulting module prediction
	private boolean isPK = false;
	
	/**
	 * @param pk polyketide fragment to be analyzed
	 * @param startCarbon the start of the polyketide, 
	 * @param endCarbon //the end of the polyketide generally a carboxylic acid
	 */
	public PolyketideModulePredictor(IAtomContainer pk, IAtom endCarbon, IAtom startCarbon, Integer pkType){
		IAtomContainer pkClone = null;
		int numStart = pk.getAtomNumber(startCarbon);
		int numEnd = pk.getAtomNumber(endCarbon);
		try {
			pkClone = pk.clone();
		} catch (CloneNotSupportedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		this.pk = pkClone;
		this.pkType = pkType;
		this.startCarbon = this.pk.getAtom(numStart);
		this.endCarbon = this.pk.getAtom(numEnd);
		analyzePK();
		isPKfrag();
	}
	
	/**
	 * @param pk polyketide fragment to be analyzed
	 * @param startCarbon //the start of the polyketide 
	 * @throws CDKException 
	 */
	public PolyketideModulePredictor(IAtomContainer pk, IAtom endCarbon, Integer pkType) throws CDKException{
		IAtomContainer pkClone = null;
		int numEnd = pk.getAtomNumber(endCarbon);
		try {
			pkClone = pk.clone();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		this.pk = pkClone;
		this.pkType = pkType;
		this.endCarbon = this.pk.getAtom(numEnd);
		predictStartCarbon();
		analyzePK();
		isPKfrag();
	}
	
	/**
	 * @param pk polyketide fragment to be analyzed
	 */
	public PolyketideModulePredictor(IAtomContainer pk, Integer pkType){
		IAtomContainer pkClone = null;
		try {
			pkClone = pk.clone();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		};
		this.pk = pkClone;
		this.pkType = pkType;
		endCarbon = predictEndCarbon(this.pk);
		if(endCarbon != null){
			predictStartCarbon();
			analyzePK();	
		}
		isPKfrag();
	}
	
	/**
	 * This method takes the polyketide and attempts to analyze it using BackboneAnalyser
	 * It will populate the loadingUnits and allDomains variables, while also setting the appropriate start domain.
	 * @param pk the polyketide being analyzed
	 * @throws CloneNotSupportedException 
	 */
	private void analyzePK() {
		if(startCarbon == null || endCarbon == null || startCarbon == endCarbon){
			return;
		}
		boolean isLinear = true;
		if(pkType > 0){
			isLinear = false;
		}
		BackboneAnalyser ba = new BackboneAnalyser(pk, startCarbon, endCarbon, isLinear);
		if(ConnectivityChecker.isConnected(pk)){
			loadingUnits = ba.getLoadingUnits();
			allDomains =  ba.getAllDomains();
		}
	}
	
	/**
	 * Predicts the start carbon. The start carbon is predicted to be the carboxylic carbon if there is only one carboxylic carbon. Else this returns no prediction
	 */
	public static IAtom predictEndCarbon(IAtomContainer molecule) {
		IAtom endCarbon = null;
		List<IAtom> endCarbonCandidates = new ArrayList<IAtom>(); 
		for(IAtom carbon : molecule.atoms()){
			if(isCarboxylicCarbon(molecule,carbon)) {
				endCarbonCandidates.add(carbon);
			}
		}
		if(endCarbonCandidates.size() == 1){
			endCarbon = endCarbonCandidates.get(0);
		}else if(endCarbonCandidates.size() > 1){
			int longestChain = -1;
			for(IAtom candidate : endCarbonCandidates){
				for(IAtom atom : molecule.atoms()){
					if(atom == candidate || atom.getAtomicNumber() != 6) continue;
					int chainLength = getCarbonOnlyPath(molecule, atom, candidate).size();
					if(chainLength > longestChain){
						longestChain = chainLength;
						endCarbon = candidate;
					}
				}
			}			
		}
		if(endCarbon == null){
			for(IAtom carbon : molecule.atoms()){
				if(carbon.getAtomicNumber() != 6) continue;
				if(isBesideTerminalDoubleBondCarbon(molecule, carbon)){
					for(IAtom terminalCarbonCandidate : molecule.getConnectedAtomsList(carbon)){
						if(ChemicalUtilities.getConnectedAtomsCountNonHydrogen(molecule,terminalCarbonCandidate) != 1) continue;
						IBond bondToChange = molecule.getBond(terminalCarbonCandidate,
								molecule.getConnectedAtomsList(terminalCarbonCandidate).get(0)); //change to single bond, it is only connected to a single atom
						bondToChange.setOrder(Order.SINGLE);
						
						Atom firstO = new Atom("O");
						firstO.setAtomTypeName("O.sp3");
						Atom secondO = new Atom("O");
						firstO.setAtomTypeName("O.sp2");
						Atom carboxylCarbon = new Atom("C");
						carboxylCarbon.setAtomTypeName("C.sp2");
						molecule.addAtom(firstO);
						molecule.addAtom(secondO);
						molecule.addAtom(carboxylCarbon);
						molecule.addBond(new Bond(carboxylCarbon, firstO, IBond.Order.SINGLE));
						molecule.addBond(new Bond(carboxylCarbon, secondO, IBond.Order.DOUBLE));
						molecule.addBond(new Bond(carboxylCarbon, terminalCarbonCandidate, IBond.Order.SINGLE));
						endCarbon = carboxylCarbon;
						break;
					}
				}
			}
		}
		return endCarbon;
	}

	private static boolean isBesideTerminalDoubleBondCarbon(IAtomContainer pk2, IAtom carbon) {
		boolean passes = false;
		boolean besideTerminal = false;
		int orderSum = 0;
		int carbonSum = 0;
		if(ChemicalUtilities.getConnectedAtomsCountNonHydrogen(pk2, carbon) == 2){
			for(IAtom atom : pk2.getConnectedAtomsList(carbon)){
				if(atom.getAtomicNumber() == 6){
					carbonSum ++;
					orderSum += pk2.getBond(carbon, atom).getOrder().numeric();
					if(ChemicalUtilities.getConnectedAtomsCountNonHydrogen(pk2, atom) == 1) besideTerminal = true;
				}
			}
		}
		if(orderSum == 3 && carbonSum == 2 && besideTerminal){
			passes = true;
		}
		return passes;
	}

	/**
	 * Predicts the end carbon, requires the start carbon to be known. It looks for the farthest connected carbon with only carbons (continuous backbon)
	 */
	private void predictStartCarbon() {
		List<IAtom> atomsCarbonOnlyPath = new ArrayList<IAtom>();
		List<IAtom> cyclicAtoms = ChemicalUtilities.getAtomsInRings(pk);
		List<IAtom> ringsWithoutOxygenAtoms = ChemicalUtilities.getAtomsInRingsWithoutOxygen(pk); //change to multi oxygen as well

		for(IAtom atom : pk.atoms()){
			if(endCarbon == atom || atom.getAtomicNumber() != 6 || numConnectedCarbons(atom) < 2) continue; // can't be the same atom at the carboxylic and must be a carbon must be connected to two other carbons to be part of a pk chain
			
			if(ringsWithoutOxygenAtoms.contains(atom)) continue; // can't end with a ring unless it has oxygen in it
			
			List<IAtom> currentAtomsCarbonOnlyPath = getCarbonOnlyPath(pk, atom, endCarbon);
			if(currentAtomsCarbonOnlyPath.size() % 2 == 0) continue; //Must be odd number
			if(currentAtomsCarbonOnlyPath.size() < 3) continue; //must be larger than 3 carbons
			if(currentAtomsCarbonOnlyPath.size() < 8 && cyclicAtoms.contains(atom)) continue; //if chain is < 8 carbons but is in ring then is almost certainly not part of pk
			
			boolean alphaInCarbonOnlyRing = false;
			
			for(IAtom connectedAtom : pk.getConnectedAtomsList(atom)){ //Check if alpha carbon is in a carbon only ring, if so, this is not correct
				if(!currentAtomsCarbonOnlyPath.contains(connectedAtom) && connectedAtom.getAtomTypeName().contains("C.")){
					if(ringsWithoutOxygenAtoms.contains(connectedAtom)){
						alphaInCarbonOnlyRing = true;
						continue;
					}
				}
			}
			if(alphaInCarbonOnlyRing) continue;
			
			if(currentAtomsCarbonOnlyPath.size() > atomsCarbonOnlyPath.size()){
				atomsCarbonOnlyPath = currentAtomsCarbonOnlyPath;
			}
		}
		if(atomsCarbonOnlyPath.size() > 0){
			startCarbon = atomsCarbonOnlyPath.get(0);
		}
	}
	private int numConnectedCarbons(IAtom atom) {
		int connectedCarbons = 0;
		for(IAtom connected : pk.getConnectedAtomsList(atom)){
			if(connected.getAtomicNumber() == 6){
				connectedCarbons ++;
			}
		}
		return connectedCarbons;
	}
	/**
	 * Finds the carbon only path between two carbon atoms by removing all other bonds. If there is no connection nothing will be returned. 
	 * @param molecule for which the path is being looked for between the two carbons (carbonOne and carbonTwo)
	 * @param carbonOne start point
	 * @param carbonTwo end point
	 * @return the list of carbons in the carbon only path between two carbons (from carbonOne to carbonTwo).  If there is no path this will return null
	 */
	public static List<IAtom> getCarbonOnlyPath(IAtomContainer molecule, IAtom carbonOne, IAtom carbonTwo) {
		
		int numAtom = molecule.getAtomNumber(carbonOne);
		int numStartCarbon = molecule.getAtomNumber(carbonTwo);
		
		IAtomContainer moleculeClone = null;
		try {
			moleculeClone = molecule.clone();
		} catch (CloneNotSupportedException e1) {
			e1.printStackTrace();
		}
		IAtom cloneAtom = moleculeClone.getAtom(numAtom);
		IAtom cloneStartCarbon = moleculeClone.getAtom(numStartCarbon);
		
		List<IAtom> cloneShortestCarbonPath = new ArrayList<IAtom>();
		IAtomContainerSet carbonOnlyMolecules = removeOxygens(moleculeClone); 
		
		IAtomContainer carbonOnlyMoleculeWithCarboxyclicCarbon = null;
		for(int i = 0; carbonOnlyMolecules.getAtomContainerCount() > i; i ++){
			IAtomContainer mol = carbonOnlyMolecules.getAtomContainer(i);
			if(mol.contains(cloneStartCarbon)){
				carbonOnlyMoleculeWithCarboxyclicCarbon = mol;
			}
		}
		
		if(carbonOnlyMoleculeWithCarboxyclicCarbon.contains(cloneAtom)){
			cloneShortestCarbonPath = PathTools.getShortestPath(carbonOnlyMoleculeWithCarboxyclicCarbon, cloneAtom, cloneStartCarbon);
		}
		List<IAtom> shortestCarbonPath = new ArrayList<IAtom>();
		if(cloneShortestCarbonPath.size() > 0){
			shortestCarbonPath = getPathFromClone(molecule, moleculeClone,cloneShortestCarbonPath);
		}
		
		return shortestCarbonPath;
	}

	/**
	 * Function that returns the proper path from the parent (non clone molecule) from the clone.
	 * @param molecule parent
	 * @param moleculeClone clone
	 * @param cloneShortestCarbonPath shortest path of clone
	 * @return shortest path of parent that is the parents atoms, not the clones atoms
	 */
	private static List<IAtom> getPathFromClone(IAtomContainer molecule, IAtomContainer moleculeClone, List<IAtom> cloneShortestCarbonPath) {
		List<Integer> atomNumbers = new ArrayList<Integer>();
		for(IAtom atom : cloneShortestCarbonPath){
			atomNumbers.add(moleculeClone.getAtomNumber(atom));
		}
		List<IAtom> shortestCarbonPath = new ArrayList<IAtom>();
		for(Integer num : atomNumbers){
			shortestCarbonPath.add(molecule.getAtom(num));
		}
		
		return shortestCarbonPath;
	}
	/**
	 * removes all non carbon atoms from a molecule
	 * @param molecule
	 * @return molecules without oxygens
	 */
	private static IAtomContainerSet removeOxygens(IAtomContainer molecule) {
		for(IAtom atom : molecule.atoms()){
			if(atom.getAtomicNumber() == 8){
				for(IAtom atom2 : molecule.getConnectedAtomsList(atom)){
					molecule.removeBond(molecule.getBond(atom,atom2));
				}
			}
		}		
		
		IAtomContainerSet partitions = ConnectivityChecker.partitionIntoMolecules(molecule);
		return partitions;
	}

	/**
	 * Gets the number of carboxylic acids from a molecule
	 * @param molecule to be checked for carboxylic acids
	 * @return number of carboxylic acids contained in molecule
	 */
	public static int getNumCarboxylicAcids(IAtomContainer molecule) {
		int numCarboxylicAcids = 0;
		for(IAtom carbon : molecule.atoms()){
			if(isCarboxylicCarbon(molecule,carbon)){
				numCarboxylicAcids ++;
			}
		}
		return numCarboxylicAcids;
	}
	
	/**
	 * Checks to see if the carbon in question is a carboxylic carbon
	 * @param molecule containing the carbon
	 * @param carbon to check
	 * @return true if carbon is the carboxcylic carbon in molecule
	 */
	public static boolean isCarboxylicCarbon(IAtomContainer molecule, IAtom carbon){
		if(carbon.getAtomicNumber() != 6){
			return false;
		}
		
		boolean isCarboyxlicCarbon = false;
		
		int singleBondO = 0;
		int doubleBondO = 0;
		int singleBondC = 0;
		for(IAtom atom : molecule.getConnectedAtomsList(carbon)){
			if(atom.getAtomicNumber() == 1 || atom.getAtomTypeName() == null){
				continue;
			}
			if(atom.getAtomTypeName().equals("O.sp3")){
				int numCarbon = 0;
				for(IAtom connectedToOxygen : molecule.getConnectedAtomsList(atom)){
					if(connectedToOxygen.getAtomicNumber() == 6){
						numCarbon++;
					}
				}
				if(numCarbon == 1){
					singleBondO ++;
				}
			}
			if(atom.getAtomTypeName().equals("O.sp2")){
				doubleBondO ++;
			}
			if(atom.getAtomicNumber() == 6){
				singleBondC ++;
			}
		}
		if(singleBondO == 1 && doubleBondO == 1 && singleBondC == 1){
			isCarboyxlicCarbon = true;
		}
		return isCarboyxlicCarbon;
	}
		
	/**
	 * Does a check on pk using allDomains to see whether there are domains, and whether it is likely a fatty acid
	 * fatty acids are predicted when there are >75% of the domains being single bonded without oxygen (oxidative state of ZERO)
	 * 	 
	 */
	private void isPKfrag(){
		int numModules = allDomains.size();
		//int numDHdomains = BackboneAnalyser.getDomainCount(allDomains,PolyKetideDomainEnums.DOUBLEBOND) + BackboneAnalyser.getDomainCount(allDomains,PolyKetideDomainEnums.SINGLEBOND); 
		if(numModules <= 2 && numModules > 0){
			isPK = true;
		}else if(numModules > 2){
			isPK = true;
		}
	}
	
	/**
	 * @return true if PolyketideModulePredictor assesses the IAtomContainer pk to be a polyketide or a fragment
	 */
	public boolean isPK(){
		return isPK;
	}
	
	/**
	 * @return the modules
	 */
	public List<List<PolyKetideDomainEnums>> getDomains(){
		return allDomains;
	}
	/**
	 * @return the list of loading units
	 */
	public List<PKsubstrate> getLoadingUnits(){
		return loadingUnits;
	}

	public IAtom getEndCarbon() {
		return endCarbon;
	}
}
