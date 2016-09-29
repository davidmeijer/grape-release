package ca.mcmaster.magarveylab.grape.pk.chem;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import ca.mcmaster.magarveylab.grape.enums.DomainEnums.*;
import ca.mcmaster.magarveylab.grape.pk.chem.BackboneRings;
import ca.mcmaster.magarveylab.grape.pk.modules.PKsubstrate;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.smsd.ring.HanserRingFinder;

/**
 * Finds and orders the carbon backbone (getBackbone) and then gets the oxidation state of each of the beta carbons (getBetaOxidationState).
 * Will also remove extra groups off of the alpha carbons and store them for later analysis
 * @author cDejong
 */

public class BackboneAnalyser {
	
		/**
		 * public method for use outside of this class
		 * @param cleavedMacrolide the macrolide to be checked
		 * @param startCarbon the first carbon in the backbone
		 * @param endCarbon the last carbon in the backbone
		 * @return wholesome, ordered, list of carbons in the backbone
		 */
	
	public static List<IAtom> getBackbone(IAtomContainer cleavedMacrolide, 
			IAtom startCarbon, IAtom endCarbon){
		List<IAtom> backboneCs = new ArrayList<IAtom>();
		backboneCs.addAll(PathTools.getShortestPath(cleavedMacrolide, startCarbon, endCarbon));
		return backboneCs;
	}	
	
	public static int getDomainCount(List<List<PolyKetideDomainEnums>> allDomains,
			PolyKetideDomainEnums domain) {
		Integer domainCount = 0;
		for(List<PolyKetideDomainEnums> domains : allDomains){
			for(PolyKetideDomainEnums singleDomain : domains){
				if(singleDomain == domain){
				domainCount = domainCount + (1/domains.size()); //to consider when there are several possibles, so not over counted
				}
			}
		}
		return(domainCount);
	}
	
	private IAtomContainer cleavedMacrolide; 
	private IAtom startCarbon; //first carbon in the backbone
	private IAtom endCarbon; //last carbon in the backbone
	private List<List<PolyKetideDomainEnums>> allDomains = new ArrayList<List<PolyKetideDomainEnums>>(); //the oxidation states of each beta carbon in the macrolide, ordered
	private List<PKsubstrate> loadingUnits = new ArrayList<PKsubstrate>(); //the alpha carbon fragments
	private BackboneRings ringRemover;	
	
	/**
	 * @param cleavedMacrolide the macrolide to be analyzed
	 * @param startCarbon first carbon in the backbone
	 * @param endCarbon last carbon in the backbone
	 */
	public BackboneAnalyser(IAtomContainer cleavedMacrolide, 
			IAtom startCarbon, IAtom endCarbon, boolean isLinear){
		this.cleavedMacrolide = cleavedMacrolide;
		this.startCarbon = startCarbon;
		this.endCarbon = endCarbon;
		ringRemover = new BackboneRings(cleavedMacrolide, startCarbon, endCarbon, isLinear);
		deriveBetaOxidationState();
		derivePolyketideAlphaUnits();
	}
	
	/**
	 * Finds the oxidation states of the beta carbons of the cleavedMacrolide and sets the variable betaCarbonsOxidationState
	 */
	private  void deriveBetaOxidationState(){
		
		List<IAtom> backboneCarbons = getBackbone(cleavedMacrolide, startCarbon, endCarbon);
		HashMap<IAtom, List<PolyKetideDomainEnums>> variableOxiState = ringRemover.getVariableBetaOxidativeStates();
		
		for (int carbon = 0; carbon < backboneCarbons.size(); carbon = 2 + carbon){ //cycle through he carbon backbone, querying the beta carbons
			IAtom currentC = backboneCarbons.get(carbon);
			List<PolyKetideDomainEnums> domains = new ArrayList<PolyKetideDomainEnums>();			
			
			//checks to see if this carbon has been edited (such as a ring opening event, and gives appropriate value
			if(variableOxiState.containsKey(currentC) == true){
				domains = variableOxiState.get(currentC);
			}
			else{
				for (IAtom connectedAtom:cleavedMacrolide.getConnectedAtomsList(currentC)){ //cycle through connected atoms,
					if(connectedAtom.getAtomicNumber() == 8
					&& cleavedMacrolide.getBond(currentC, connectedAtom).getOrder() == IBond.Order.valueOf("SINGLE") //checks for hydroxy group (single bond to O) attached to beta carbon
							){
						domains.add(PolyKetideDomainEnums.HYDROXYL);
						break;		
					}
					else if(connectedAtom.getAtomicNumber() == 8
					&& cleavedMacrolide.getBond(connectedAtom, currentC).getOrder() == IBond.Order.valueOf("DOUBLE")//checks for ketone group (double bond to O) attached to beta carbon
					){
						domains.add(PolyKetideDomainEnums.KETONE);
						break;
					}
					else if(connectedAtom.getAtomicNumber() == 6 
					&& cleavedMacrolide.getBond(connectedAtom, currentC).getOrder() == IBond.Order.valueOf("DOUBLE")//checks for a double bonded C to the beta carbon
					){
						domains.add(PolyKetideDomainEnums.DOUBLEBOND);
						break;
					}
				}
			}
			if(domains.size() < 1){
				domains.add(PolyKetideDomainEnums.SINGLEBOND); // state of carbon with only carbon bonds and two hydrogens (no groups or extra bonds attached to the beta carbon)
			}
			
			allDomains.add(domains);
		}
		//Remove the last oxidation state -> it is wrong with this logic because it is not taking into consideration possibility of two sepirate O bonds. It is always 3, so replace with it.
		allDomains.remove(allDomains.size()-1);
		List<PolyKetideDomainEnums> start = new ArrayList<PolyKetideDomainEnums>();
		start.add(PolyKetideDomainEnums.KETONE);
		start.add(null);
		allDomains.add(0, start);
	}
	
	/**
	 * Determines the alpha unit type. This allows determination of the type of loaded carbon unit.
	 */
	private void derivePolyketideAlphaUnits(){ //TODO: DO MORE TESTING. Update: Tested by Andrew Webber, working as expected, though more testing should be done on a diverse set of alpha subunits
		
		IRingSet RingAtoms = getAllRings();
		List<IAtom> backboneCarbons = getBackbone(cleavedMacrolide, startCarbon, endCarbon);
		
		for (int carbon = 1; carbon < backboneCarbons.size(); carbon = 2 + carbon){ //cycle through the carbon backbone, querying the alpha carbons
			
			IAtom currentC = backboneCarbons.get(carbon);
			PKsubstrate fragment = null;
			Integer sideChains = 0;
			
			for (IAtom connectedAtom:cleavedMacrolide.getConnectedAtomsList(currentC)){ //cycle through connected atoms,
				checkMethoxylMalonyl(backboneCarbons, connectedAtom, fragment);
				
				if(!backboneCarbons.contains(connectedAtom)
				&& connectedAtom.getAtomTypeName().contains("C.")
				){
					checkBenzoyl(RingAtoms, backboneCarbons, connectedAtom, fragment);
					sideChains ++;
				}
			}
			if(fragment == null
			&& sideChains > 0){
				if(sideChains == 1){ //ETHYLMALONYL METHYLMALONYL
					fragment = getSingleCarbonSideChainFragment(
							backboneCarbons,currentC);
				}else if (sideChains == 2){ //ISOBUTERYL AND 2_METHYLBUTERYL
					fragment = getDoubleCarbonSideChainFragment(
							backboneCarbons, currentC);
				}else{
					fragment = PKsubstrate.UNKNOWN;
				}
			}else{
				fragment = PKsubstrate.MALONYL;//if nothing else...
			}
			loadingUnits.add(fragment);
		}
		loadingUnits.add(0,PKsubstrate.UNKNOWN);
	}

	/**
	 * @param backboneCarbons list of backbone carbons, ordered
	 * @param currentC the carbon that is currently being analyzed
	 * @return returns the unit type if there are two carbons not from the backbone on the alpha carbon (4 carbons total)
	 */
	private PKsubstrate getDoubleCarbonSideChainFragment(
			List<IAtom> backboneCarbons, IAtom currentC) {
		PKsubstrate fragment;
		int connectedCarbonsCount = 0;
		for (IAtom connectedAtom:cleavedMacrolide.getConnectedAtomsList(currentC)){ //cycle through connected atoms,
			if(!backboneCarbons.contains(connectedAtom)
			&& connectedAtom.getAtomTypeName().contains("C.")	
			){
				for (IAtom connectedAtom2:cleavedMacrolide.getConnectedAtomsList(connectedAtom)){ //cycle through connected atoms,
					if(connectedAtom2.getAtomTypeName().contains("C.")){
						connectedCarbonsCount ++;
					}
				}			
			}
		}
		if (connectedCarbonsCount == 2){
			fragment = PKsubstrate.ISOBUTRYL;
		}
		else if(connectedCarbonsCount ==3){
			fragment = PKsubstrate.METHYLBUTERYL2;
		}
		else{
			fragment = PKsubstrate.UNKNOWN;
		}
		return fragment;
	}

	/**
	 * @param backboneCarbons list of backbone carbons, ordered
	 * @param currentC the carbon that is currently being analyzed
	 * @return returns the unit type if there is one carbon not from the backbone on the alpha carbon (3 carbons total)
	 */
	private PKsubstrate getSingleCarbonSideChainFragment(List<IAtom> backboneCarbons,
			IAtom currentC) {
		PKsubstrate fragment;
		int connectedCarbonsCount = 0;
		for (IAtom connectedAtom:cleavedMacrolide.getConnectedAtomsList(currentC)){ //cycle through connected atoms,
			if(!backboneCarbons.contains(connectedAtom)
			&& connectedAtom.getAtomTypeName().contains("C.")	
			){
				for (IAtom connectedAtom2:cleavedMacrolide.getConnectedAtomsList(connectedAtom)){ //cycle through connected atoms,
					if(connectedAtom2.getAtomTypeName().contains("C.")){
						connectedCarbonsCount ++;
					}
				}
				
			}
		}
		if (connectedCarbonsCount == 1){
			fragment = PKsubstrate.METHYLMALONYL;
		}
		else if(connectedCarbonsCount == 2){
			fragment = PKsubstrate.ETHYLMALONYL;
		}
		else{
			fragment = PKsubstrate.UNKNOWN;
		}
		return fragment;
	}

	/**
	 * @return all the rings in the molecule
	 * NOTE: May be better refactored to a general class, OR refactored to a parent object if this is done in other classes
	 */
	private IRingSet getAllRings(){
		HanserRingFinder ringSearch = new HanserRingFinder();
		
		IRingSet RingAtoms = null; // holds the atoms that are in rings
		
		try {
			RingAtoms = ringSearch.getRingSet(cleavedMacrolide); // get all atoms in rings
			
		} catch (CDKException e) {
			e.printStackTrace();			
		}
		return RingAtoms;
	}

	/**
	 * @param backboneCarbons list of backbone carbons, ordered
	 * @param connectedAtom atom connected to the alpha carbon
	 * @param fragment piece connected to alpha carbon
	 */
	private void checkMethoxylMalonyl(List<IAtom> backboneCarbons,
			IAtom connectedAtom, PKsubstrate fragment) {
		if(!backboneCarbons.contains(connectedAtom)
		&& connectedAtom.getAtomTypeName().contains("O.")
		){
			for (IAtom connectedToO: cleavedMacrolide.getConnectedAtomsList(connectedAtom)){
				if(connectedToO.getAtomTypeName().contains("C.")
				&&!backboneCarbons.contains(connectedToO)// need to check if this is also connected to a carbon that is not on the backbone
				){
					fragment = PKsubstrate.METHOXYLMALONYL;
				}
			}
		}
	}

	/**
	 * @param RingAtoms all ring atoms of molecule
	 * @param backboneCarbons list of backbone carbons, ordered
	 * @param connectedAtom atom connected to the alpha carbon
	 * @param fragment piece connected to alpha carbon
	 */
	private void checkBenzoyl(IRingSet RingAtoms, List<IAtom> backboneCarbons,
			IAtom connectedAtom, PKsubstrate fragment) {
		for (IAtom connectedToC: cleavedMacrolide.getConnectedAtomsList(connectedAtom)){
			if(connectedToC.getAtomTypeName().contains("C.")
			&&RingAtoms.contains(connectedToC)
			&&!backboneCarbons.contains(connectedToC)// need to check if this is also connected to a carbon that is not on the backbone
			){
				fragment = PKsubstrate.BENZOYL;
			}
		}
	}
	
	public List<List<PolyKetideDomainEnums>> getAllDomains(){
		return allDomains;
	}
	
	public List<PKsubstrate> getLoadingUnits(){
		return loadingUnits;
	}

}
