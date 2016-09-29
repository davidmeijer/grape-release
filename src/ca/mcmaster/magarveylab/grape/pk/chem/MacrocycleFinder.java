package ca.mcmaster.magarveylab.grape.pk.chem;

import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;

/**
 * Finds macrocycles of lactones/lactams, can find counts of each, and detemrine the key molecules involved (hydroxyl and carboxyl carbons and hydroxyl oxygen)
 * @author cDejong &
 * adapted from gchen lactone break in rprism (NRPSmodifications.java)
 */

public class MacrocycleFinder {
	
	/**
	 * @param molecule being analyzed
	 * @return the lactone carbons and hydroxyl oxygen
	 */
	
	public static List<List<IAtom>> getLactoneCandO(IMolecule molecule){
		
		List<List<IAtom>> lactoneAtoms = getLactCandO(molecule, "O.sp3");
		
		return lactoneAtoms;
	}
	
	/**
	 * @param molecule being analyzed
	 * @return the lactam carbons and hydroxyl oxygen
	 */
	public static List<IAtom> getLactamCandN(IMolecule molecule){
		
		List<List<IAtom>> lactamAtoms = getLactCandO(molecule, "N.amide");
		return lactamAtoms.get(0);
		
	}
	
	/**
	 * @param molecule being analyzed
	 * @param type of atom being looked for, either N.amide (lactam) or O.sp3 (lactone)
	 * @return the lact carbons and hydroxyl oxygen
	 */
	private static List<List<IAtom>> getLactCandO(IMolecule molecule, String type) {
		List<List<IAtom>> allCandOs = new ArrayList<List<IAtom>>();		
		
		for(int index = 0; index < molecule.getAtomCount(); index++) {
			if(!(molecule.getAtom(index).getAtomTypeName().equals(type))) {
				continue;
			}
			 List<IAtom> CandOs = new ArrayList<IAtom>();
			IAtom hydroxylO = molecule.getAtom(index);
			// check if it is linked to two carbons, one of which has a double bond O.
			
			IAtom carboxylC = null;
			IAtom hydroxylC = null;
			
			for(IAtom c:molecule.getConnectedAtomsList(hydroxylO)) {
				// the carboxyl carbon must be sp2 hybridized and connected to another oxygen via double bond
				if(c.getAtomTypeName().equals("C.sp2")) {
					for(IAtom candidateOtherO:molecule.getConnectedAtomsList(c)) {
						if(candidateOtherO != hydroxylO 
						&& candidateOtherO.getAtomTypeName().contains("O.")
						&& molecule.getBond(candidateOtherO, c).getOrder() == IBond.Order.DOUBLE
						){
							carboxylC = c;
						}
					}
				}if(c.getAtomTypeName().contains("C.") && c != carboxylC){
					hydroxylC = c; // want to save this carbon
				}
			}
			
			if(carboxylC == null || hydroxylC == null) {
				continue;
			}
			IBond removedBond = molecule.removeBond(hydroxylO, carboxylC);
			//Then check if breaking the hydroxylO and carboxylC leaves one piece.					
			if(ConnectivityChecker.isConnected(molecule)) {
				CandOs.add(hydroxylO);
				CandOs.add(hydroxylC);
				CandOs.add(carboxylC);
				allCandOs.add(CandOs);
			}
			molecule.addBond(removedBond);
		}
		return allCandOs;		
	}

	/**
	 * @param molecule being analyzed
	 * @return number of lactams
	 */
	public static int getLactamCount(IMolecule molecule) {
		int lactamCount = getLactCount(molecule, "N.amide");
		return lactamCount;
	}
	/**
	 * @param molecule being analyzed
	 * @return number of lactones
	 */
	public static int getLactoneCount(IAtomContainer molecule) {
		
		int lactoneCount = getLactCount(molecule, "O.sp3");
		return lactoneCount;
	}
	
	/**
	 * @param molecule being analyzed
	 * @param type of atom being looked for, either N.amide (lactam) or O.sp3 (lactone)
	 * @return number of lacts 
	 */
	private static int getLactCount(IAtomContainer molecule, String type) {
		int lactCount = 0;
		for(int index = 0; index < molecule.getAtomCount(); index++) {
			
			if(!(molecule.getAtom(index).getAtomTypeName().equals(type))) {
				continue;
			}
			
			IAtom hydroxylO = molecule.getAtom(index);
			// check if it is linked to two carbons, one of which has a double bond O.
			
			IAtom carboxylC = null;
			IAtom hydroxylC = null;
			for(IAtom c:molecule.getConnectedAtomsList(hydroxylO)) {
				// the carboxyl carbon must be sp2 hybridized and connected to another oxygen via double bond
				if(c.getAtomTypeName().contains("C.sp2")) {
					for(IAtom candidateOtherO:molecule.getConnectedAtomsList(c)) {
						if(candidateOtherO != hydroxylO 
						&& candidateOtherO.getAtomTypeName().contains("O.")
						&& molecule.getBond(candidateOtherO, c).getOrder() == IBond.Order.DOUBLE
						){
							carboxylC = c;
						}
					}
				}if(c.getAtomTypeName().contains("C.") && c != carboxylC){
					hydroxylC = c; // want to save this carbon
				}
			}
			
			if(carboxylC == null || hydroxylC == null) {
				continue;
			}
			IBond removedBond = molecule.removeBond(hydroxylO, carboxylC);
			//Then check if breaking the hydroxylO and carboxylC leaves one piece.					
			if(ConnectivityChecker.isConnected(molecule)) {
				molecule.addBond(removedBond);
				lactCount ++;
			}
			molecule.addBond(removedBond);
		}
		return lactCount;		
	}
}
