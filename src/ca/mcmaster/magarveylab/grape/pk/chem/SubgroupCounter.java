package ca.mcmaster.magarveylab.grape.pk.chem;

import java.util.List;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.smsd.ring.HanserRingFinder;
import org.openscience.cdk.tools.manipulator.RingSetManipulator;

import ca.mcmaster.magarveylab.grape.pk.modules.Modifications;

/**
 * This class determines the subgroups, or 'tailors' in the macrolide
 * @author cDejong
*/

public class SubgroupCounter {

	/**
	 * @param molecule being analyzed for subgroups
	 * @return which subgroups are present and how many of each
	 */
	
	public static Modifications getSubgroupCount(IAtomContainer molecule){
		
		Modifications modifications = new Modifications(); //May be better to just be an Array of strings? Don't currently utilize the hashMap functionality
		
		HanserRingFinder ringSearch = new HanserRingFinder();
		IRingSet ringAtoms = null;
		try {
			ringAtoms = ringSearch.getRingSet(molecule);
		} catch (CDKException e) {
			return modifications;
		}
		List<IAtomContainer> allRings = RingSetManipulator.getAllAtomContainers(ringAtoms);	
		
		modifications.addChlorineCount(atomType(molecule,"Cl"));
		
		for(IAtomContainer ring: allRings){
			
			if(isSugar(ring, allRings)){
				if(isDeoxySugar(ring,molecule)){
					modifications.addDeoxySugar();
				}else{
					modifications.addSugar();
				}
			}
			if(isEpoxide(ring, allRings)){
				modifications.addEpoxide();
			}
		}
		/*
		String starterUnit = null;//checkForStarterUnit(molecule); TO ADD BACK WHEN HAVE STARTERUNITS
		
		if(starterUnit != null){
			subgroupCounts.put("Starter unit:", starterUnit);
		}
		*/
		
		return modifications;
	}
	
	/**
	 * @param ring being checked 
	 * @param molecule being analyzed
	 * @return weather the ring is the cyclic part of a deoxysugar
	 */
	private static boolean isDeoxySugar(IAtomContainer ring, IAtomContainer molecule) {
		int oxygenCount = 0;
		for (IAtom atom: ring.atoms()){ //Count how many oxygens are connected to the carbons in the sugar ring
			if (atom.getAtomTypeName().contains("C.")){
				for (IAtom connectedAtom: molecule.getConnectedAtomsList(atom)){
					if (connectedAtom.getAtomTypeName().contains("O.")){
						oxygenCount ++;
					}
				}
			}
		}
		if(oxygenCount < 5){
			return true;
		}
		return false;
	}
	/**
	 * @param ring being checked
	 * @param allRings of entire molecule being analyzed
	 * @return weather the ring is an epoxide
	 */
	private static boolean isEpoxide(IAtomContainer ring,
			List<IAtomContainer> allRings) {
		if(ring.getAtomCount() == 3){
			for(IAtom atom: ring.atoms()){
				if(atom.getAtomTypeName().contains("O.sp3")){
					return true;
				}
			}
		}
		return false;
	}

	/**
	 * @param ring
	 * @param allRings
	 * @return
	 */
	/**
	 * @param ring being checked
	 * @param allRings of entire molecule being analyzed
	 * @return weather the ring is the cyclic part of a sugar (may be deoxy, is non specific)
	 */
	private static boolean isSugar(IAtomContainer ring, 
			List<IAtomContainer> allRings){
		
		for(IAtom atom: ring.atoms()){
			if(atom.getAtomTypeName().contains("O.")
			&& getNumberOfRingsContainingAtom(atom, allRings) == 1
			&& ring.getAtomCount() == 6
			){
				for(IAtom atom2: ring.atoms()){
					if(getNumberOfRingsContainingAtom(atom2, allRings) != 1){
						return false;
					}
				}
				return true;
			}			
		}
		return false;
	}
	/**
	 * Should potentially be in a utility class
	 * @param atom being checked
	 * @param allRings in molecule
	 * @return number of times the specific atom is in a different ring from allRings
	 */
	private static int getNumberOfRingsContainingAtom(IAtom atom, 
			List<IAtomContainer> allRings){
		int counter = 0;
		for(IAtomContainer ring: allRings){
			if(ring.contains(atom)){counter++;}
		}
		return counter;
	}
	
	/**
	 * Should potentially be in a utility class
	 * @param moleule having it's atoms counted for type
	 * @param type of atom
	 * @return number of type in molecule
	 */
	public static int atomType(IAtomContainer moleule, String type){
		int atomCount = 0;
		for(IAtom atom: moleule.atoms()){
			if(atom.getAtomTypeName().contains(type)){
				atomCount ++;
			}
		}
		
		return atomCount;
	}
}


