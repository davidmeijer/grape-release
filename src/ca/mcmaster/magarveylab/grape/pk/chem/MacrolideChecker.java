package ca.mcmaster.magarveylab.grape.pk.chem;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.smsd.ring.HanserRingFinder;
import org.openscience.cdk.tools.manipulator.RingSetManipulator;

/**
 * Checks whether the input molecule is a macrolide/macrolactam or not. Returns 0 if not a macrolide/lactam 1 if macrolide 2 if macrolactam
 * @author cDejong
 * 
 */

public class MacrolideChecker {
	private static List<IAtomContainer> allRingsSplit = null;
	private static int maxNumberOfRings = 500; //Higher the number the more that will be parsed but the slower it will be. Really slows down over ~500
	/**
	 * @param molecule being checked if a macrolide/lactam
	 * @return 0 if not macrolide/lactam, 1 if macrolide, 2 if macrolactam
	 */
	public static int isMacrolide(IMolecule molecule){
		int typeOfMolecule = 0; //0 for not proper, 1 for lactone, 2 for lactam
		
		HanserRingFinder ringSearch = new HanserRingFinder();
		
		IRingSet Rings = null;
		try {
			Rings = ringSearch.getRingSet(molecule);
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		RingSetManipulator.sort(Rings);

		allRingsSplit = RingSetManipulator.getAllAtomContainers(Rings);
		
		if(allRingsSplit.size() > maxNumberOfRings){//tweakable -- higher values will slow down run, but process more complex molecules
			System.out.println("Too many rings to parse out in a decent amount of time (over " + maxNumberOfRings + ")" + allRingsSplit.size());
			return 0;
		}
		int lactoneCount = MacrocycleFinder.getLactoneCount(molecule); 
		int lactamCount = MacrocycleFinder.getLactamCount(molecule);
		if(lactoneCount > 0){ //for lactone macrolides, cannot contain a lactam, and may have multiple lactones
			if(isProperLactone(molecule)){
				//System.out.print("LACTONE: ");
				typeOfMolecule = 1;
			}
		}if (lactamCount > 0 && typeOfMolecule == 0){ //for macrolactams
			if(isProperLactam(molecule)){
				//System.out.print("LACTAM: ");
				typeOfMolecule = 2;
			}
		}
		return typeOfMolecule;
	}


	/**
	 * @param molecule being checked
	 * @return true if a proper lactone for a macrolide
	 */
	private static boolean isProperLactone(IMolecule molecule) {
		List<List<IAtom>> lactoneAtomsList = MacrocycleFinder.getLactoneCandO(molecule);
		List<IAtom> lactoneAtoms = new ArrayList<IAtom>();
		for(List<IAtom>singleLactoneAtoms:lactoneAtomsList){
			lactoneAtoms.addAll(singleLactoneAtoms);			
		}
		
		List<IAtomContainer> ringsContainingLactone = new ArrayList<IAtomContainer>();
		for(IAtomContainer ring: allRingsSplit){
			if(hasAllAtoms(ring, lactoneAtoms)){
				ringsContainingLactone.add(ring);
			}
		}
		
//		if(singleLactoneMacrocycle == false){
//			System.out.println("Lactones are in two different macrocycles, cannot be parsed");
//			return false;
//		}
		if(ringsContainingLactone.size() == 0){
			return false;
		}
		IAtomContainer smallestRing = ringsContainingLactone.get(0);
		int nitrogenCounter = 0;
		for(IAtom atom: smallestRing.atoms()){
			if(atom.getAtomTypeName().contains("N.")){
				nitrogenCounter ++;
			}
		}
		if(smallestRing.getAtomCount() > 8
		&& nitrogenCounter < 2
		//&& !containsAtomType(smallestRing, "S.")
		){
			return true;
		}
		return false;
	}


	/**
	 * @param molecule being checked
	 * @return true if a proper lactam for a macrolactam
	 */
	private static boolean isProperLactam(IMolecule molecule) {
		List<IAtom> lactamAtoms = MacrocycleFinder.getLactamCandN(molecule);
			
		
		List<IAtomContainer> ringsContainingLactam = new ArrayList<IAtomContainer>();
		for(IAtomContainer ring: allRingsSplit){
			if(hasAllAtoms(ring, lactamAtoms)){
				ringsContainingLactam.add(ring);
			}
		}
		
//		if(singleLactoneMacrocycle == false){
//			System.out.println("Lactones are in two different macrocycles, cannot be parsed");
//			return false;
//		}
		if(ringsContainingLactam.size() == 0){
			return false;
		}
		
		IAtomContainer smallestRing = ringsContainingLactam.get(0);
		int nitrogenCounter = 0;
		for(IAtom atom: smallestRing.atoms()){
			if(atom.getAtomTypeName().contains("N.")){
				nitrogenCounter ++;
			}
		}
		
		if(smallestRing.getAtomCount() > 8
		&& nitrogenCounter < 2
		//&& !containsAtomType(ringContainingLacam, "S.")
		){
			return true;
		}
		return false;
	}
	

	/**
	 * Utility method that may be moved to a more general utility class
	 * @param ring being checked for atoms
	 * @param atomsToMatch atoms being checked if in the ring
	 * @return true if atomsToMatch are all in ring
	 */
	private static boolean hasAllAtoms(IAtomContainer ring, List<IAtom> atomsToMatch) {
		int numberOfAtoms = atomsToMatch.size();
		int matchedNumberOfAtoms = 0;
		for (IAtom atom: atomsToMatch){
			if(ring.contains(atom)){
				matchedNumberOfAtoms ++;
			}
		}
				
		if (numberOfAtoms == matchedNumberOfAtoms){
			return true;
		}			
		return false;
	}
}