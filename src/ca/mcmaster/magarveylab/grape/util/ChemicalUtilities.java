package ca.mcmaster.magarveylab.grape.util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.Stack;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.isomorphism.mcss.RMap;
import org.openscience.cdk.smsd.Isomorphism;
import org.openscience.cdk.smsd.interfaces.Algorithm;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ca.mcmaster.magarveylab.grape.enums.DomainEnums.CStarterSubstrate;
import ca.mcmaster.magarveylab.grape.nrp.chem.ChemicalAbstraction;
import ca.mcmaster.magarveylab.grape.nrp.chem.Fragment;
import ca.mcmaster.magarveylab.grape.nrp.chem.Fragment.FragmentType;
import ca.mcmaster.magarveylab.grape.nrp.chem.MoleculePredictor;
import ca.mcmaster.magarveylab.grape.util.io.SmilesIO;

public class ChemicalUtilities {

	static UniversalIsomorphismTester uit = new UniversalIsomorphismTester();
	/**
	 * Get a list of all smallest rings with size less than or equal to 6 atoms
	 * Uses a breadth-first search starting from each atom.
	 * 
	 * @param molecule
	 * @return
	 */
	public static List<Set<IAtom>> getSmallestRings(IAtomContainer molecule) {
		List<Set<IAtom>> rings = new ArrayList<Set<IAtom>>();
	
		class AtomNode {
			IAtom atom;
			AtomNode parent;
			int depth;
	
			AtomNode(IAtom atom, AtomNode parent, int depth) {
				this.atom = atom;
				this.parent = parent;
				this.depth = depth;
			}
		}
	
		for (int i = 0; i < molecule.getAtomCount(); i++) {
			IAtom startAtom = molecule.getAtom(i);
			// For further speed optimization, a visited list can be tracked and
			// a traceback of two ends of depth 3 can be used
			// Set<IAtom> visited = new HashSet<IAtom>();
			Queue<AtomNode> q = new LinkedList<AtomNode>();
	
			AtomNode firstNode = new AtomNode(startAtom, null, 0);
			q.add(firstNode);
			// visited.add(startAtom);
			AtomNode lastNodeInSmallestRing = null;
	
			while (!q.isEmpty()) {
				AtomNode currentNode = q.poll();
				for (IAtom connectedAtom : molecule
						.getConnectedAtomsList(currentNode.atom)) {
					if (currentNode.parent != null
							&& connectedAtom == currentNode.parent.atom) {
						continue;
					}
					if (connectedAtom == firstNode.atom) {
						lastNodeInSmallestRing = new AtomNode(connectedAtom,
								currentNode, currentNode.depth + 1);
						break;
					}
					// Make sure that this atom is not repeated in the current
					// chain
	
					AtomNode tracebackNode = currentNode;
					boolean repeat = false;
					while (tracebackNode != null) {
						if (tracebackNode.atom == connectedAtom) {
							repeat = true;
							break;
						}
						tracebackNode = tracebackNode.parent;
					}
					if (repeat) {
						continue;
					}
					if (currentNode.depth + 1 <= 6) {
						q.add(new AtomNode(connectedAtom, currentNode,
								currentNode.depth + 1));
						// visited.add(connectedAtom);
					}
				}
				if (lastNodeInSmallestRing != null) {
					break;
				}
			}
			if (lastNodeInSmallestRing != null) {
				// Trace back to the first node to get all atoms in the ring
	
				Set<IAtom> ring = new HashSet<IAtom>();
				AtomNode currentNode = lastNodeInSmallestRing;
				while (currentNode != null) {
					ring.add(currentNode.atom);
					currentNode = currentNode.parent;
				}
				boolean repeat = false;
				for (Set<IAtom> otherRing : rings) {
					if (otherRing.equals(ring)) {
						repeat = true;
						break;
					}
				}
				if (repeat == false) {
					rings.add(ring);
				}
			}
		}
	
		return rings;
	}

	/**
	 * From a molecule, find all atoms which participate in a cycle. This method
	 * performs a non-recursive depth-first search and find atoms in cycles by
	 * tracing back-edges of the DFS search tree as it is made.
	 * 
	 * @param molecule
	 *            A molecule object which may be connected or disconnected.
	 * @return
	 */
	public static ArrayList<IAtom> getAtomsInRings(IAtomContainer molecule) {
		Set<IAtom> atomsInRings = new HashSet<IAtom>();
	
		// Track what atoms are visited
		Set<IAtom> visited = new HashSet<IAtom>();
	
		// Iterate through atoms to use as root. For a connected molecule, this
		// look occurs once since all atoms will be discovered.
		for (int i = 0; i < molecule.getAtomCount(); i++) {
	
			if (visited.contains(molecule.getAtom(i))) {
				continue;
			}
	
			Stack<IAtom> nodesToVisit = new Stack<IAtom>();
	
			HashMap<IAtom, IAtom> parentInDFSTree = new HashMap<IAtom, IAtom>();
	
			// Use the zero-index atom as the root
			nodesToVisit.push(molecule.getAtom(i));
	
			while (!nodesToVisit.empty()) {
				IAtom currentAtom = nodesToVisit.peek();
				
				boolean firstVisit = false;
				if (!visited.contains(currentAtom)) {
					firstVisit = true;
					visited.add(currentAtom);
				}
				boolean hasNext = false;
				for (IAtom nextAtom : molecule.getConnectedAtomsList(currentAtom)) {
					if (visited.contains(nextAtom)) {
						if (firstVisit == true
								&& nextAtom != parentInDFSTree.get(currentAtom)) {
							// This is a cycle. Trace all atoms in this ring and
							// add to the atoms-in-rings set.
							// Iterate through the stack, find nextAtom, and add
							// it and all subsequent atoms.
							IAtom traceAtom = currentAtom;
							atomsInRings.add(traceAtom);
							do {
								traceAtom = parentInDFSTree.get(traceAtom);
								atomsInRings.add(traceAtom);
							} while(traceAtom != nextAtom);
						}
						continue;
					}
					// Push this atom to the top of the stack
					nodesToVisit.push(nextAtom);
					parentInDFSTree.put(nextAtom, currentAtom);
					hasNext = true;
				}
				if (!hasNext) {
					// Remove the object at the top of the stack
					nodesToVisit.pop();
				}
			}
		}
	
		return new ArrayList<IAtom>(atomsInRings);
	}

	public static boolean hasCarbonPath(IAtomContainer mol, IAtom carbon1,
			IAtom carbon2, ArrayList<IAtom> visited) {
		if (carbon1.getAtomicNumber() != 6) {
			return false;
		}
		// Base case
		if (carbon1.equals(carbon2)) {
			return true;
		}
		visited.add(carbon1);
		for (IAtom a : mol.getConnectedAtomsList(carbon1)) {
			if (!visited.contains(a) && a.getAtomicNumber() == 6) {
				if (hasCarbonPath(mol, a, carbon2, visited)) {
					return true;
				}
			}
		}
		return false;
	}

	public static boolean hasCarbonPath(IAtomContainer mol, IAtom carbon1,
			IAtom carbon2) {
		ArrayList<IAtom> visited = new ArrayList<IAtom>();
		return hasCarbonPath(mol, carbon1, carbon2, visited);
	}

	/**
	 * From a molecule, find all atoms which participate in a cycle that does
	 * not contain an oxygen.
	 * 
	 * @param molecule
	 * @return
	 */
	public static ArrayList<IAtom> getAtomsInRingsWithoutOxygen(
			IAtomContainer molecule) {
		// List<IBond> bondsToRemove = new ArrayList<IBond>
	
		// Remove all bonds which contain an oxygen
	
		ArrayList<IBond> removedBonds = new ArrayList<IBond>();
		ArrayList<IAtom> removedAtoms = new ArrayList<IAtom>();
	
		for (int i = 0; i < molecule.getBondCount(); i++) {
			if (molecule.getBond(i).getAtom(0).getAtomicNumber() == 8
					|| molecule.getBond(i).getAtom(1).getAtomicNumber() == 8) {
				removedBonds.add(molecule.getBond(i));
			}
		}
		for(IBond bond : removedBonds) {
			molecule.removeBond(bond);
		}
	
		// Remove all oxygens
		for (int i = 0; i < molecule.getAtomCount(); i++) {
			if (molecule.getAtom(i).getAtomicNumber() == 8) {
				removedAtoms.add(molecule.getAtom(i));
			}
		}
		for(IAtom atom : removedAtoms) {
			molecule.removeAtom(atom);
		}
	
		ArrayList<IAtom> atomsInRing = getAtomsInRings(molecule);
		
		// Replace atoms and bonds
		for (IAtom a : removedAtoms)
			molecule.addAtom(a);
		for (IBond b : removedBonds)
			molecule.addBond(b);
		
		return atomsInRing;
	}

	/**
	 * From a molecule, find all atoms which participate in a carbon cycle.
	 * 
	 * @param molecule
	 * @return
	 */
	public static ArrayList<IAtom> getAtomsInCarbonRings(IAtomContainer molecule) {
		// List<IBond> bondsToRemove = new ArrayList<IBond>
	
		// Remove all bonds which do not contain a carbon
	
		ArrayList<IBond> removedBonds = new ArrayList<IBond>();
		ArrayList<IAtom> removedAtoms = new ArrayList<IAtom>();
	
		for (int i = 0; i < molecule.getBondCount(); i++) {
			if (molecule.getBond(i).getAtom(0).getAtomicNumber() != 6
					|| molecule.getBond(i).getAtom(1).getAtomicNumber() != 6) {
				removedBonds.add(molecule.getBond(i));
				molecule.removeBond(i);
			}
		}
	
		// Remove all non-carbons
		for (int i = 0; i < molecule.getAtomCount(); i++) {
			if (molecule.getAtom(i).getAtomicNumber() != 6) {
				removedAtoms.add(molecule.getAtom(i));
				molecule.removeAtom(i);
				i--;
			}
		}
	
		ArrayList<IAtom> carbonsInRing = getAtomsInRings(molecule);
	
		// Replace atoms and bonds
		for (IAtom a : removedAtoms)
			molecule.addAtom(a);
		for (IBond b : removedBonds)
			molecule.addBond(b);
		return carbonsInRing;
	}

	public static boolean hasConnectedAtomOfAtomicNumber(IAtomContainer mol,
			IAtom atom, int atomicNumber) {
		for (IAtom a : mol.getConnectedAtomsList(atom)) {
			if (a.getAtomicNumber() == atomicNumber) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Check if an atom in a molecule has a neighbour of a given atomic number
	 * 
	 * @param molecule
	 * @param atom
	 * @param neighbourAtomicNumber
	 * @return
	 */
	public static boolean hasNeighbourOfAtomicNumber(IAtomContainer molecule, IAtom atom,
			int neighbourAtomicNumber) {
		List<IAtom> connectedAtomsList = molecule.getConnectedAtomsList(atom);
		for (IAtom neighbour : connectedAtomsList) {
			if (neighbour.getAtomicNumber() == neighbourAtomicNumber) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Given a template structure with a specified bond, identify all instances
	 * of this structure as a substructure in another molecule The outer list
	 * corresponds to a substructure match, and the inner list contains the
	 * (ordered) matching bonds in the molecule
	 * 
	 * @param template
	 * @param templateBond
	 * @param molecule
	 * @return
	 */
	public static List<List<IBond>> findMatchingBondsFromTemplate(IAtomContainer template,
			List<IBond> templateBonds, IAtomContainer molecule) {
		List<List<RMap>> templateMatchMap = null;
		try {
			templateMatchMap = uit.getSubgraphMaps(
					molecule, template);
		} catch (Exception e) {
			System.err.println("Error in UniversalIsomorphismTester");
			return new ArrayList<List<IBond>>();
		}
	
		List<List<IBond>> allMatchingBonds = new ArrayList<List<IBond>>();
	
		for (int i = 0; i < templateMatchMap.size(); i++) {
			List<IBond> matchingBonds = new ArrayList<IBond>();
			for (IBond templateBond : templateBonds) {
				for (int j = 0; j < templateMatchMap.get(i).size(); j++) {
					IBond currentTemplateMatchBond = template
							.getBond(templateMatchMap.get(i).get(j).getId2());
					if (currentTemplateMatchBond == templateBond) {
						IBond matchingBond = molecule.getBond(templateMatchMap
								.get(i).get(j).getId1());
						matchingBonds.add(matchingBond);
						break;
					}
				}
			}
			allMatchingBonds.add(matchingBonds);
		}
		return allMatchingBonds;
	}

	/**
	 * Given a template structure with a specified bond, identify all instances
	 * of this structure as a substructure in another molecule Each bond will
	 * appear at most once in the final list, even if multiple substructures
	 * with this specified bond exist.
	 * 
	 * @param template
	 * @param templateBond
	 * @param molecule
	 * @return
	 */
	public static  List<IBond> findMatchingBondsFromTemplate(IAtomContainer template,
			IBond templateBond, IAtomContainer molecule) {
		ArrayList<IBond> templateBonds = new ArrayList<IBond>();
		templateBonds.add(templateBond);
		List<List<IBond>> matchingBonds = findMatchingBondsFromTemplate(
				template, templateBonds, molecule);
		HashSet<IBond> matchesWithoutDuplicates = new HashSet<IBond>();
		for (int i = 0; i < matchingBonds.size(); i++) {
			matchesWithoutDuplicates.add(matchingBonds.get(i).get(0));
		}
		return new ArrayList<IBond>(matchesWithoutDuplicates);
	}

	public static double getTanimotoScore(IAtomContainer molecule1, IAtomContainer molecule2) {
		Isomorphism iso = new Isomorphism(Algorithm.DEFAULT, true);
		try {
			iso.init(molecule1, molecule2, true, true);
		} catch (CDKException e) {
			e.printStackTrace();
		}
		iso.setChemFilters(true, true, true);
		double tanimoto = 0;
		try {
			tanimoto = iso.getTanimotoSimilarity();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return tanimoto;
	}

	/**
	 * Implement the smsd Isomorphism algorithm to get Tanimoto scores comparing a given structure with a list of possible matches
	 * @param molecule The molecule to compare
	 * @param possibleMatches The list of possible matches
	 * @return An ArrayList of same size of possibleMatches, corresponding with the Tanimoto score of molecule to each element of possibleMatches
	 */
	public static ArrayList<Double> getTanimotoScores(IAtomContainer molecule, List<IAtomContainer> possibleMatches) {
		
		// Get the max values
		ArrayList<Double> tanimotoScores = new ArrayList<Double>();
		for(int i = 0; i < possibleMatches.size(); i++) {
			IAtomContainer p = possibleMatches.get(i);
			double score = getTanimotoScore(molecule, p);
			tanimotoScores.add(score);
		}
		return(tanimotoScores);
	}

	/**
	 * Find highest Tanimoto score achieved by a molecule matched against a list of possible matches
	 * @param molecule
	 * @param possibleMatches
	 * @return
	 */
	public static double getHighestTanimotoScore(IAtomContainer molecule, List<IAtomContainer> possibleMatches) {
		ArrayList<Double> tanimotoScores = getTanimotoScores(molecule, possibleMatches);
		double maxval = 0;
		for(double score:tanimotoScores) {
			if(score > maxval) {
				maxval = score;
			}
		}
		return(maxval);
	}

	public static int getHighestTanimotoScoreIndex(ArrayList<Double> tanimotoScores) {
		double maxval = 0;
		int index = 0;
		for(int i = 0; i < tanimotoScores.size(); i++) {
			if(tanimotoScores.get(i) > maxval) {
				maxval = tanimotoScores.get(i);
				index = i;
			}
		}
		return(index);
	}

	public static CStarterSubstrate matchCStarter(IAtomContainer mol, Map<CStarterSubstrate, IAtomContainer> cStarterSubstrates) {
		for(java.util.Map.Entry<CStarterSubstrate, IAtomContainer> matching : cStarterSubstrates.entrySet()){
			try {
				UniversalIsomorphismTester uit = new UniversalIsomorphismTester();
				if(uit.isSubgraph(mol, matching.getValue()))
					return matching.getKey();
			} catch (CDKException e) {
				e.printStackTrace();
			}
		}
		return null;
	}

	/**
	 * Get monomer fragments from SMILES string
	 * @param smiles
	 * @param name 
	 * @throws ExceptionSmilesTooLarge 
	 * @throws ExceptionSmilesNotConnected 
	 * @throws InvalidSmilesException 
	 */
	public static ChemicalAbstraction getChemicalAbstractionFromSmiles(String smiles, String name, MoleculePredictor predictor){
		IAtomContainer currentMolecule = null;
		ChemicalAbstraction chemicalAbstraction = new ChemicalAbstraction(name);
		try {
			currentMolecule = SmilesIO.readSmiles(smiles);
		} catch (IllegalArgumentException | IOException | CDKException e) {
			e.printStackTrace();
			chemicalAbstraction.setErrorMessage("Bad smiles");
			return chemicalAbstraction;
		}
		if(!ConnectivityChecker.isConnected(currentMolecule)) {
			chemicalAbstraction.setErrorMessage("Molecule not connected");
		}else{
			try {
				if(AtomContainerManipulator.getNaturalExactMass(currentMolecule) > 2500){ //TODO: currently no overrride
					System.out.println("Molecule is large... this may take a while");
				}
				chemicalAbstraction = predictor.getChemicalAbstraction(currentMolecule, name);
			} catch(Exception e) {
				e.printStackTrace();
				chemicalAbstraction.setErrorMessage("Could not be broken down");			
			}
		}
		return chemicalAbstraction;
	}

	/**
	 * Given an unordered list of monomer fragments, get a list of sequences of monomerFragments ordered in a NC - NC - NC fashion
	 * @param unorderedMonomerFragmentList
	 * @return
	 */
	public static List<ArrayList<Fragment>> getOrderedMonomerFragmentSequences(List<Fragment> unorderedMonomerFragments) {
		List<ArrayList<Fragment>> orderedMonomerFragmentsList = new ArrayList<ArrayList<Fragment>>();
		
		// Track which indices of unorderedMonomerFragments have been added to the ordered list of lists
		HashMap<Fragment, Boolean> visited = new HashMap<Fragment, Boolean>();
		for(Fragment m : unorderedMonomerFragments) {
			visited.put(m, false);
		}
		// Find fragments which have no other fragment in the N-direction, then trace in an NC-NC-NC fashion
		for(Fragment m : unorderedMonomerFragments) {
			if(visited.get(m))
				continue;
		
			if(m.getAtomAfterNTerminus() == null) {
				ArrayList<Fragment> orderedMonomerFragments = new ArrayList<Fragment>();
				// Trace the fragments from the starting fragment in a NC - NC - NC direction
				Fragment currentFragment = m;
				while (currentFragment != null && visited.get(currentFragment) == false){
					orderedMonomerFragments.add(currentFragment);
					visited.put(currentFragment, true);
					currentFragment = currentFragment.getFragmentAfterCTerminus(unorderedMonomerFragments);
				}
				orderedMonomerFragmentsList.add(orderedMonomerFragments);
			}			
		}
		// Find unvisited fragments that are not amino acids, then trace in an NC-NC-NC fashion
		for(Fragment m : unorderedMonomerFragments) {
			if(visited.get(m))
				continue;
			if(m.getFragmentType() != FragmentType.AMINO_ACID) {
				ArrayList<Fragment> orderedMonomerFragments = new ArrayList<Fragment>();
				// Trace the fragments from the starting fragment in a NC - NC - NC direction
				Fragment currentFragment = m;
				while (currentFragment != null && visited.get(currentFragment) == false){
					orderedMonomerFragments.add(currentFragment);
					visited.put(currentFragment, true);
					currentFragment = currentFragment.getFragmentAfterCTerminus(unorderedMonomerFragments);
				}
				orderedMonomerFragmentsList.add(orderedMonomerFragments);
			}
		}
		// What remains must be amino acid cycles.
		for(Fragment m : unorderedMonomerFragments) {
			//System.out.println("here");
			if(visited.get(m))
				continue;
			//System.out.println("here2");
			ArrayList<Fragment> orderedMonomerFragments = new ArrayList<Fragment>();
			// Trace the fragments from the starting fragment in a NC - NC - NC direction
			Fragment currentFragment = m;
			while (currentFragment != null && visited.get(currentFragment) == false){
				if(currentFragment.getFragmentType() == FragmentType.AMINO_ACID){
					currentFragment.setKnownStart(false);
				}
				orderedMonomerFragments.add(currentFragment);
				visited.put(currentFragment, true);
				currentFragment = currentFragment.getFragmentAfterCTerminus(unorderedMonomerFragments);				
			}
			orderedMonomerFragmentsList.add(orderedMonomerFragments);
		}
		return(orderedMonomerFragmentsList);
	}

	public static int getConnectedAtomsCountNonHydrogen(IAtomContainer molecule, IAtom terminalCarbonCandidate) {
		int count = 0;
		for(IAtom connectedAtom : molecule.getConnectedAtomsList(terminalCarbonCandidate)){
			if(connectedAtom.getAtomicNumber() != 1){
				count ++;
			}
		}
		return count;
	}
	
	public static int getTotalBondOrderExcludingHydrogen(IAtomContainer molecule, IAtom atom) {
		int order = 0;
		for(IAtom connectedAtom : molecule.getConnectedAtomsList(atom)){
			if(connectedAtom.getAtomicNumber() != 1){
				IBond bond = molecule.getBond(atom, connectedAtom);
				order += bond.getOrder().ordinal() + 1;
			}
		}
		return order;
	}
	
	public static boolean connectedToAtomicNumber(IAtomContainer mol, IAtom atom, int i) {
		for(IAtom connectedAtom : mol.getConnectedAtomsList(atom)){
			if(connectedAtom.getAtomicNumber() == i){
				return true;
			}
		}
		return false;
	}

}
