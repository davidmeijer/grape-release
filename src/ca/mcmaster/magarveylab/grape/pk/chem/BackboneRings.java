/**
 * 
 */
package ca.mcmaster.magarveylab.grape.pk.chem;

import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;

import ca.mcmaster.magarveylab.grape.enums.DomainEnums.*;
import ca.mcmaster.magarveylab.grape.pk.chem.BackboneAnalyser;


import org.openscience.cdk.Atom;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.PathTools;

import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.tools.manipulator.RingSetManipulator;

import ca.mcmaster.magarveylab.grape.util.ChemicalUtilities;
import ca.mcmaster.magarveylab.grape.util.io.SmilesIO;

import org.openscience.cdk.smsd.ring.HanserRingFinder;



/**
 * Open ring structures that are in the carbon backbone
 * As well as record the two possible states of the carbons' oxidation state prior to cyclization
 * @author cDejong
 * 
 */

public class BackboneRings{

private IAtomContainer pk; //The polyketide being analyzed
private IAtom startC; // the "start" carbon of the scaffold. is actually the end (of the macrocycle) during synthesis
private IAtom endC; // the "end" carbon of the scaffold. is actualy the start during synthesis
private HashMap<IAtom, List<PolyKetideDomainEnums>> variableOxiState = new HashMap<IAtom,List<PolyKetideDomainEnums>>();
private boolean isLinear;

	/**
	 * @param pk polyketide being analyzed
	 * @param startC first carbon of the scaffold
	 * @param endC last carbon of the scaffold
	 * @param isLinear 
	 */
	public BackboneRings(IAtomContainer pk, 
			IAtom startC, IAtom endC, boolean isLinear){
		
		this.pk = pk;
		this.startC = startC;
		this.endC = endC;
		this.isLinear = isLinear;
		open();	
	}
	
	/**
	 * checks what type of ring and call the appropriate method to open it
	 */
	private void open(){ 
		HanserRingFinder ringSearch = new HanserRingFinder();
		
		openAvermectinLikeRings(); // Change this to look for substeructuees outside of this forloop.
		
		List<IAtom> backboneCs = BackboneAnalyser.getBackbone(pk, startC, endC);
		
		for (int index = 0; index < backboneCs.size(); index ++){ //go through each backbone carbon and see if it is in a ring
			IRingSet RingAtoms; // holds the atoms that are in rings
			
			try {
				RingAtoms = ringSearch.getRingSet(pk); // get all atoms in rings
				
			} catch (CDKException e) {
				e.printStackTrace();
				break;
				
			}
			if(RingAtoms.isEmpty() == true){ // if there are no rings in the pk stop looking for carbons in (nonExistant) rings
				break;
			}
			while (index < backboneCs.size() - 1 && !RingAtoms.contains(backboneCs.get(index))){ //cycles through the carbons, this is done so there's no reRun of the ringfinder if nothing is changed
				index ++;
			}
			
			IAtom currentCarbon = backboneCs.get(index);
			
			if(RingAtoms.contains(currentCarbon)){ // if the current carbon is in a ring then go through the break
				IRingSet ringContainingCurrentCarbon = RingAtoms.getRings(currentCarbon);
				if(index >= backboneCs.size()) continue;
				//TODO:if it was one of the following, don't do rest
				epoxideCheckForAndBreak(pk, backboneCs, ringContainingCurrentCarbon, currentCarbon, index);
				oxygenInBackboneCheckForAndBreak(pk, backboneCs, ringContainingCurrentCarbon, currentCarbon, index);
				carbonCarbonRingCheckForAndBreak(pk, backboneCs, ringContainingCurrentCarbon, currentCarbon, index);
				
			}
			backboneCs = BackboneAnalyser.getBackbone(pk, startC, endC);	
		}
	}
	/**
	 * @param pk the polyketide molecule
	 * @param backboneCs list of backbone carbons, ordered
	 * @param ringContainingCurrentCarbon the ring with the backbone carbon to be opened
	 * @param currentCarbon 
	 * @param currentCarbonPosition the currentCarbons number in the backboneCs list
	 */
	private void carbonCarbonRingCheckForAndBreak(IAtomContainer pk, List<IAtom> backboneCs, IRingSet ringContainingCurrentCarbon, IAtom currentCarbon, int currentCarbonPosition) {
		
		IAtom nextAtom = backboneCs.get(currentCarbonPosition+1);
		IAtom previousAtom = null;
		try {previousAtom = backboneCs.get(currentCarbonPosition-1);}catch(IndexOutOfBoundsException e){return;}
		IAtom next2Atom = null;
		try {next2Atom = backboneCs.get(currentCarbonPosition+2);}catch(IndexOutOfBoundsException e){return;}
		
		RingSetManipulator.sort(ringContainingCurrentCarbon);
		List<IAtomContainer> allRings = RingSetManipulator.getAllAtomContainers(ringContainingCurrentCarbon);
		IAtomContainer allRingsTogether = RingSetManipulator.getAllInOneContainer(ringContainingCurrentCarbon);
		IAtomContainer smallestRingWithCurrentCarbon = getSmallestRingWithAtom(allRings, currentCarbon);
		
		try {
			if (smallestRingWithCurrentCarbon.getAtomCount() == 6 
			&& CDKHueckelAromaticityDetector.detectAromaticity(smallestRingWithCurrentCarbon)){
				//begin steps to open and restore the ketones in an aromatic ring
				pk.removeBond(currentCarbon,nextAtom);				
				Atom atom = new Atom("O");
				atom.setAtomTypeName("O.sp2");
				pk.addAtom(atom);
				
				if (isBetaCarbon(backboneCs, currentCarbon)){
					pk.addBond(pk.getAtomNumber(currentCarbon),pk.getAtomNumber(atom), IBond.Order.DOUBLE);
				}
				else{
					pk.addBond(pk.getAtomNumber(nextAtom),pk.getAtomNumber(atom), IBond.Order.DOUBLE);
				}
				
				List<IAtom> orderedRing = PathTools.getShortestPath(pk, currentCarbon, nextAtom);
				for(int index = 1; orderedRing.size() > index + 1; index++){
					IAtom ringCarbon = orderedRing.get(index);
					
						for(IAtom connectedAtom: pk.getConnectedAtomsList(ringCarbon)){
							if (connectedAtom.getAtomTypeName().contains("O.")){
								pk.removeBond(ringCarbon,connectedAtom);
								pk.addBond(pk.getAtomNumber(ringCarbon),pk.getAtomNumber(connectedAtom), IBond.Order.DOUBLE);
							}
						}
				}
				return;
			}
		} catch (CDKException e) {}
		
		if (nextAtom.getAtomTypeName().contains("C.")
		&& previousAtom.getAtomTypeName().contains("C.")
		&& allRingsTogether.contains(nextAtom)
		&& !allRingsTogether.contains(previousAtom)
		&& !allRingsTogether.contains(next2Atom)
		&& allRingsTogether.getAtomCount() > 3
		){
			pk.removeBond(currentCarbon,nextAtom);
			Atom atom = new Atom("O");
			atom.setAtomTypeName("O.sp2");
			pk.addAtom(atom);
			
			if(allRingsTogether.getAtomCount() % 2 == 0){//checks for evenness see if we know where the ketone is or not
				if (isBetaCarbon(backboneCs, currentCarbon)){
					pk.addBond(pk.getAtomNumber(currentCarbon),pk.getAtomNumber(atom), IBond.Order.DOUBLE);
				}else{
					pk.addBond(pk.getAtomNumber(nextAtom),pk.getAtomNumber(atom), IBond.Order.DOUBLE);
				}
			}else{
				pk.addBond(pk.getAtomNumber(currentCarbon),pk.getAtomNumber(atom), IBond.Order.DOUBLE);
				List<PolyKetideDomainEnums> domains = new ArrayList<PolyKetideDomainEnums>();
				domains.add(PolyKetideDomainEnums.SINGLEBOND);
				domains.add(PolyKetideDomainEnums.KETONE);
				variableOxiState.put(currentCarbon, domains);
				variableOxiState.put(nextAtom, domains);
			}
		}	
	}

	/**
	 * @param allRings of the molecule ordered by size, smallest first
	 * @param currentCarbon that contains a ring
	 * @return smallest ring
	 */
	private IAtomContainer getSmallestRingWithAtom(List<IAtomContainer> allRings, 
			IAtom currentCarbon) { // add exception if not contained?
		
		IAtomContainer smallestRing = null;
		
		for(IAtomContainer ring: allRings){
			if (ring.contains(currentCarbon)){
				smallestRing = ring;
				break;
			}
		}
		return smallestRing;
	}

	/**
	 * @param pk the polyketide molecule
	 */
	private void openAvermectinLikeRings() {
		openDoubleRingOxygen();
		openAvermectinLikeEnd();
	}

	private void openDoubleRingOxygen() {
		IAtomContainer template = null;
		try {template = SmilesIO.readSmilesTemplates("CC(C)COC1(C)CC(O)CC(C)O1");} catch (Exception e) {} //doubleRing template to match
		IBond templateCarbonOxygenBond1 = template.getBond(template.getAtom(5),template.getAtom(13));
		IBond templateCarbonOxygenBond2 = template.getBond(template.getAtom(5),template.getAtom(4));
		List<IBond> templateCarbonOxygenBond1Matches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateCarbonOxygenBond1, (IAtomContainer)pk);
		if(templateCarbonOxygenBond1Matches.size() < 1) return;
		List<IBond> templateCarbonOxygenBond2Matches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateCarbonOxygenBond2, (IAtomContainer)pk);
		
		for(IBond bond : templateCarbonOxygenBond1Matches){
			pk.removeBond(bond);
		}
		for(IBond bond : templateCarbonOxygenBond2Matches){
			pk.removeBond(bond);
		}
		for(IBond bond : templateCarbonOxygenBond1Matches){
			for(IAtom atom : bond.atoms()){
				if(atom.getAtomicNumber() != 6) continue;
				Atom o = new Atom("O");
				atom.setAtomTypeName("O.sp2");
				pk.addAtom(o);
				pk.addBond(pk.getAtomNumber(o), pk.getAtomNumber(atom), IBond.Order.DOUBLE);
			}
		}
	}

	private void openAvermectinLikeEnd() {
		IAtomContainer template = null;
		try {template = SmilesIO.readSmilesTemplates("OC1C=CC(C(O)=O)[C@]2(O)CCO[C@H]12");} catch (Exception e) {}
		IBond templateHydroxyl = template.getBond(template.getAtom(8),template.getAtom(9)); //change to double bond
		IBond templateCarbonRingBreak = template.getBond(template.getAtom(4),template.getAtom(8)); //break this bond
		IBond templateOxygenRemove = template.getBond(template.getAtom(12),template.getAtom(13)); //Remove this oxygen
		
		List<IBond> templateHydroxylMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateHydroxyl, (IAtomContainer)pk);
		if(templateHydroxylMatches.size() < 1) return;
		List<IBond> templateCarbonRingBreakMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateCarbonRingBreak, (IAtomContainer)pk);
		List<IBond> templateOxygenRemoveMatches = ChemicalUtilities.findMatchingBondsFromTemplate(template, templateOxygenRemove, (IAtomContainer)pk);
		
		for(IBond bond : templateHydroxylMatches){
			bond.setOrder(IBond.Order.DOUBLE);
		}
		for(IBond bond : templateOxygenRemoveMatches){
			for(IAtom atom : bond.atoms()){
				if(atom.getAtomicNumber() == 8){
					pk.removeAtomAndConnectedElectronContainers(atom);
				}
			}
		}
		for(IBond bond : templateCarbonRingBreakMatches){
			pk.removeBond(bond);
		}
	}
	
	/**
	 *  Opens pyrones in the carbon backbone while also recording the potential states of the carbons where bonds were broken 
	 * @param pk
	 * @param backboneCs list of backbone carbons, ordered
	 * @param ringContainingCurrentCarbon
	 * @param currentCarbon
	 * @param currentCarbonPosition the currentCarbons number in the backboneCs list
	 */
	private void oxygenInBackboneCheckForAndBreak(IAtomContainer pk, 
			List<IAtom> backboneCs, IRingSet ringContainingCurrentCarbon, IAtom currentCarbon, Integer currentCarbonPosition){

		IAtom nextAtom = backboneCs.get(currentCarbonPosition+1);
		IAtom next2Atom = null;
		try {next2Atom = backboneCs.get(currentCarbonPosition+2);}catch(IndexOutOfBoundsException e){return;}
		if (!nextAtom.getAtomTypeName().equals("O.sp3")){ //See if oxygen is next atom in backbone and if it is sp3 hybridised (connected to another carbon without anything else)
		return;
		}
		if (RingSetManipulator.getAtomCount(ringContainingCurrentCarbon) % 2 == 0){ //checks if ring size is even, if so we cannot tell which beta carbon had what starting unit, so make it variable
			pk.removeBond(currentCarbon, nextAtom);
			pk.removeBond(nextAtom, next2Atom);
			pk.addBond(pk.getAtomNumber(currentCarbon),pk.getAtomNumber(nextAtom), IBond.Order.DOUBLE);
			List<PolyKetideDomainEnums> domains = new ArrayList<PolyKetideDomainEnums>();
			domains.add(PolyKetideDomainEnums.HYDROXYL);
			domains.add(PolyKetideDomainEnums.KETONE);
			if(isLinear){ //For some linear polyketides this is also sometimes an epoxide, so must consider this as a possiblity
				domains.add(PolyKetideDomainEnums.DOUBLEBOND);
			}
			variableOxiState.put(currentCarbon, domains);
			variableOxiState.put(nextAtom, domains);
			
		}else { //if it is odd we know that the beta carbon had a ketone, so figure out which is the beta and give that the ketone after the break
			if (isBetaCarbon(backboneCs, currentCarbon)){
				pk.removeBond(nextAtom, next2Atom);
				pk.removeBond(currentCarbon, nextAtom);
				pk.addBond(pk.getAtomNumber(currentCarbon),pk.getAtomNumber(nextAtom), IBond.Order.DOUBLE);
				if(isLinear){ //For some linear polyketides this is also sometimes an epoxide, so must consider this as a possiblity
					List<PolyKetideDomainEnums> domains = new ArrayList<PolyKetideDomainEnums>();
					domains.add(PolyKetideDomainEnums.KETONE);
					domains.add(PolyKetideDomainEnums.DOUBLEBOND);
					variableOxiState.put(currentCarbon, domains);
				}
			}else{
				pk.removeBond(currentCarbon, nextAtom);
				pk.removeBond(nextAtom, next2Atom);
				pk.addBond(pk.getAtomNumber(nextAtom),pk.getAtomNumber(next2Atom), IBond.Order.DOUBLE);
				if(isLinear){ //For some linear polyketides this is also sometimes an epoxide, so must consider this as a possiblity
					List<PolyKetideDomainEnums> domains = new ArrayList<PolyKetideDomainEnums>();
					domains.add(PolyKetideDomainEnums.KETONE);
					domains.add(PolyKetideDomainEnums.DOUBLEBOND);
					variableOxiState.put(next2Atom, domains);
				}
			}
		}
	} 
	/**
	 * Checks if it is an epoxide rings at the current carbon in the backbone and breaks it
	 * @param pk
	 * @param backboneCs list of backbone carbons, ordered
	 * @param ringContainingCurrentCarbon
	 * @param currentCarbon
	 * @param currentCarbonPosition the currentCarbons number in the backboneCs list
	 */
	private void epoxideCheckForAndBreak(IAtomContainer pk, 
			List<IAtom> backboneCs,	IRingSet ringContainingCurrentCarbon, IAtom currentCarbon, Integer currentCarbonPosition){
		IAtom nextAtom = null;
		try {nextAtom = backboneCs.get(currentCarbonPosition+1);}catch(IndexOutOfBoundsException e){return;}
		
		if(RingSetManipulator.getAtomCount(ringContainingCurrentCarbon) == 3 
		&& ringContainingCurrentCarbon.contains(nextAtom)){ //checks if epoxide
			for (IAtom connectedAtom:pk.getConnectedAtomsList(currentCarbon)){	
				if (connectedAtom.getAtomTypeName().contains("O.")){ //confirms the third member is an oxygen

					removeBondsAndAtom(pk, connectedAtom);
					pk.removeBond(currentCarbon,nextAtom); //must remove the single bond in order to add the double bond in the next step
					pk.addBond(pk.getAtomNumber(currentCarbon),pk.getAtomNumber(nextAtom), IBond.Order.DOUBLE);
				}
			}
		}
	}

	private static void removeBondsAndAtom(IAtomContainer pk, 
			IAtom toBeRemoved) {
		for(IAtom connected: pk.getConnectedAtomsList(toBeRemoved)){
			pk.removeBond(connected, toBeRemoved);
		}
		pk.removeAtom(toBeRemoved);
	}
	/**
	 * @param backboneCs list of backbone carbons, ordered
	 * @param currentCarbon
	 * @return true if beta carbon, false if not
	 */
	private static boolean isBetaCarbon(List<IAtom> backboneCs, 
			IAtom currentCarbon) { //counts the carbons before the current carbon, if even it is a beta carbon.
		int counter = 0;
		for (IAtom backboneCarbon: backboneCs){
			counter++;
			if (backboneCarbon == currentCarbon){
				break;
			}
		}
		if (counter % 2 != 0){
			return true;
		}
		else{
			return false;
		}
	}
	public HashMap<IAtom, List<PolyKetideDomainEnums>> getVariableBetaOxidativeStates(){
		return variableOxiState; 
	}
}

