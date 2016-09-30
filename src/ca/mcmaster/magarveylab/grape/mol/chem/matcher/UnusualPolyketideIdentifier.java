package ca.mcmaster.magarveylab.grape.mol.chem.matcher;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;

import ca.mcmaster.magarveylab.grape.enums.DomainEnums.PolyKetideDomainEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.TailoringDomainEnums;
import ca.mcmaster.magarveylab.grape.mol.chem.Fragment;
import ca.mcmaster.magarveylab.grape.mol.chem.Fragment.FragmentType;
import ca.mcmaster.magarveylab.grape.pk.chem.BackboneAnalyser;
import ca.mcmaster.magarveylab.grape.pk.chem.PolyketideModulePredictor;
import ca.mcmaster.magarveylab.grape.pk.modules.PKsubstrate;
import ca.mcmaster.magarveylab.grape.util.ChemicalUtilities;
import ca.mcmaster.magarveylab.grape.util.io.SmilesIO;

public class UnusualPolyketideIdentifier {

	private IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
	private SubstrateMatcher sm;

	public UnusualPolyketideIdentifier(SubstrateMatcher sm){
		this.sm = sm;
	}
	
	/**
	 * Identify as a linear polyketide
	 * @param currentFragment
	 * @param macrolideType
	 */
	public void identifyAsLinearPK(Fragment currentFragment, int macrolideType) {
		// If there are is a lactone bond and an amino bond
		if(currentFragment.getLactoneHydroxylC() != null && currentFragment.getAminoCs().size() > 0) {
			//identifyAsLinearPK(currentFragment, currentFragment.getAminoC(), currentFragment.getLactoneHydroxylC());
			identifyAsLinearPK(currentFragment, null, currentFragment.getAminoCs().get(0), macrolideType);
		}
		// If there is one connected amino acid
		else if(currentFragment.getAminoCs().size() > 0 && currentFragment.getAminoNs().size() > 0) {
			identifyAsLinearPK(currentFragment, null, null, macrolideType);
			if(currentFragment.getFragmentType() != null) {
				currentFragment.setFragmentType(FragmentType.FA_OR_PK);
			}
		}
		// If this does not appear to be an amino acid
		else {
			identifyAsLinearPK(currentFragment, null, null, macrolideType);
		}
		// Check if there is a saturated carbon chain indicative of a possible fatty acid
		if(currentFragment.getFragmentType() == FragmentType.POLYKETIDE) {
			SMARTSQueryTool querytoolSingle = null;
			SMARTSQueryTool querytoolDouble = null;
			querytoolSingle = new SMARTSQueryTool("[#6;A][#6;A;D2][#6;A;D2][#6;A;D2][#6;A;D2][#6;A]", builder);
			querytoolDouble = new SMARTSQueryTool("[#6;A]=[#6;A;D2]\\[#6;A;D2]=[#6;A;D2]\\[#6;A;D2]=[#6;A]", builder);
			boolean matches = false;
			try {
				currentFragment.addImplicitHydrogens();
				matches = querytoolSingle.matches(currentFragment.getAtomContainer());
				if(!matches) matches = querytoolDouble.matches(currentFragment.getAtomContainer());
			} catch (CDKException e) {
				e.printStackTrace();
			}
			if(matches) {
				currentFragment.setFragmentType(FragmentType.FA_OR_PK);
			}
		}
		
		// Check if this is possibily a modified amino acid - nitrogen containing with at most 20 atoms and two or fewer PK units
		if(currentFragment.getAtomContainer().getAtomCount() <= 20) {
			boolean hasNitrogen = false;
			for(int i = 0; i < currentFragment.getAtomContainer().getAtomCount(); i++) {
				if(currentFragment.getAtomContainer().getAtom(i).getAtomicNumber() == 7) {
					hasNitrogen = true;
				}
			}
			if(hasNitrogen && currentFragment.getLoadingUnits() != null && currentFragment.getLoadingUnits().size() <= 2) {
				currentFragment.setFragmentType(FragmentType.FA_OR_PK);
			}
		}
	}


	/**
	 * General method for perfoming PK checks on a monomer fragment. Modifies the monomer fragment to update
	 * FragmentType if necessary and chemical changes corresponding to PK processing.
	 * Returns true if the Fragment is determined to be a PK.
	 * @param currentFragment
	 * @param start
	 * @param end
	 * @return
	 */
	public void identifyAsLinearPK(Fragment currentFragment, IAtom start, IAtom end, int macrolideType) {
		// If there are two connected amino acids
		if(end != null && start != null && ChemicalUtilities.hasCarbonPath(currentFragment.getAtomContainer(), start, end)) {
			try {
				PolyketideModulePredictor pkPredictor =
						new PolyketideModulePredictor(
								currentFragment.getAtomContainer(),
								start,
								end,
								macrolideType
								);
				if(pkPredictor.isPK()) {
					List<List<PolyKetideDomainEnums>> pkDomains = pkPredictor.getDomains();
					List<PKsubstrate> loadingUnits = pkPredictor.getLoadingUnits();
					//Modifications modifications = pkPredictor.getTailorsWithCount();
					currentFragment.setFragmentType(FragmentType.POLYKETIDE);
					currentFragment.setPkDomains(pkDomains);
					currentFragment.setLoadingUnits(loadingUnits);
					//currentFragment.setPkModifications(modifications);
				}
			} catch(Exception e) {
				//e.printStackTrace();
				System.err.println("Error in linear PK predictor - skipping a portion of PK predictions");
			}
		}
		// If there is an end
		else if(end != null) {
			PolyketideModulePredictor pkPredictor = null;
			try {
				pkPredictor = new PolyketideModulePredictor(
						currentFragment.getAtomContainer(),
						end,
						macrolideType);
				if(pkPredictor.isPK()) {
					currentFragment.setFragmentType(FragmentType.POLYKETIDE);
					List<List<PolyKetideDomainEnums>> pkDomains = pkPredictor.getDomains();
					List<PKsubstrate> loadingUnits = pkPredictor.getLoadingUnits();
					//Modifications modifications = pkPredictor.getTailorsWithCount();
					currentFragment.setPkDomains(pkDomains);
					currentFragment.setLoadingUnits(loadingUnits);
					//currentFragment.setPkModifications(modifications);
				}
			} catch (Exception e) {
				//e.printStackTrace();
				System.err.println("Error in linear PK predictor - skipping a portion of PK predictions");
			}
			
			
		}
		else {	
			PolyketideModulePredictor pkPredictor = null;
			try {
			pkPredictor = new PolyketideModulePredictor(
					currentFragment.getAtomContainer(),
					macrolideType
					);
			if(pkPredictor.isPK()) {
				currentFragment.setFragmentType(FragmentType.POLYKETIDE);
				List<List<PolyKetideDomainEnums>> pkDomains = pkPredictor.getDomains();
				List<PKsubstrate> loadingUnits = pkPredictor.getLoadingUnits();
				//Modifications modifications = pkPredictor.getTailorsWithCount();
				currentFragment.setPkDomains(pkDomains);
				currentFragment.setLoadingUnits(loadingUnits);
				//currentFragment.setPkModifications(modifications);
			}
			} catch(Exception e) {
				//e.printStackTrace();
				System.err.println("Error in linear PK predictor - skipping a portion of PK predictions");
			}
		}
		List<List<PolyKetideDomainEnums>> allDomains = currentFragment.getPkDomains();
		int numModules = allDomains.size();
		if(currentFragment.getFragmentType() != null ){
			if(currentFragment.getFragmentType().equals(FragmentType.POLYKETIDE) && numModules > 2){
				int numDHdomains = BackboneAnalyser.getDomainCount(allDomains,PolyKetideDomainEnums.DOUBLEBOND) + BackboneAnalyser.getDomainCount(allDomains,PolyKetideDomainEnums.SINGLEBOND);
				if(allDomains.get(allDomains.size() - 1).contains(PolyKetideDomainEnums.HYDROXYL))	numModules --;//get second last domain
				if(numModules - numDHdomains < 2){ // fully reduced other than the start and the potental beta hydroxyl then could be a fa_or_pk
					currentFragment.setFragmentType(FragmentType.FA_OR_PK);
				}
			}
		}
	}

	
	/**
	 * Take a monomer fragment and, if it is a ketoextended amino acid, return a two-valued array of monomerfragments
	 * corresponding to the amino acid and polyketide respectively
	 * @param currentFragment
	 * @return
	 */
	public Fragment[] getKetoextendedAAPieces(Fragment currentFragment, int macrolideType) {
		IAtomContainer mol = currentFragment.getAtomContainer();
		if(macrolideType > 0) {
			return null;
		}
		for(IAtom aminoN : mol.atoms()){
			if(aminoN.getAtomicNumber() != 7){
				continue;
			}
			boolean appropriateAminoN = false;
			if(currentFragment.getAminoNs().contains(aminoN)){
				appropriateAminoN = true;
				for(IAtom atom : mol.getConnectedAtomsList(aminoN)){
					if(atom.getAtomicNumber() != 7 &&
							atom.getAtomicNumber() != 6 && 
							atom.getAtomicNumber() != 1){
						appropriateAminoN = false;
						break;
					}
				}
			}
			//if not appropriate nitrogen, see if it is primary
			if(!appropriateAminoN){
				if(ChemicalUtilities.getConnectedAtomsCountNonHydrogen(mol, aminoN) == 1){
					appropriateAminoN = true;
				}
			}
			
			boolean inSmallRing = false;
			List<Set<IAtom>> smallestRings = ChemicalUtilities.getSmallestRings(mol);
			for(Set<IAtom> ring : smallestRings){
				if(ring.contains(aminoN) && ring.size() < 7){
					inSmallRing = true;
					break;
				}
			}
			
			if(inSmallRing){
				appropriateAminoN = true;
			}
			
			//if still no appropriate aminoN then skip this fragment
			if(!appropriateAminoN){
				continue;
			}
			// Create the array that will contain the C fragment and N fragment respectively
			Fragment[] pieces = new Fragment[2];
			// Check if this is a ketoextended amino acid
			IAtom possiblePKcarbon = null;
			List<IBond> candidateBondsToBreak = new ArrayList<IBond>();
			
			//check if a cyclized amino acid
			if(inSmallRing){
				IAtomContainer template = null; 
				try {
					template = SmilesIO.readSmilesTemplates("CC(=O)CN");
				} catch(Exception e) {
					e.printStackTrace();
				}
				IBond templateBondToBreak = template.getBond(template.getAtom(0),template.getAtom(1));
				candidateBondsToBreak.addAll(ChemicalUtilities.findMatchingBondsFromTemplate(template, templateBondToBreak, mol));
			}
			
			//check if linear regular amino acid
			if(candidateBondsToBreak.size() == 0){
				IAtomContainer template = null; 
				try {
					template = SmilesIO.readSmilesTemplates("CCC(=O)C(C)N");
				} catch(Exception e) {
					e.printStackTrace();
				}
				IBond templateBondToBreak = template.getBond(template.getAtom(1),template.getAtom(2));
				candidateBondsToBreak.addAll(ChemicalUtilities.findMatchingBondsFromTemplate(template, templateBondToBreak, mol));
			}
			for(IBond bond : candidateBondsToBreak){
				boolean ketoCarbon = false;
				for(IAtom atom : bond.atoms()){
					for(IAtom a : mol.getConnectedAtomsList(atom)){
						if(a.getAtomicNumber() == 8 && mol.getBond(a, atom).getOrder().equals(Order.DOUBLE)){
							ketoCarbon = true;
						}
					}
					if(!ketoCarbon){
						possiblePKcarbon = atom;
						break;
					}
				}
			}
			//if neither try an extension out
			if(candidateBondsToBreak.size() == 0){
				for(IAtom a : mol.getConnectedAtomsList(aminoN)) {
					if(a.getAtomicNumber() == 6
							&& ChemicalUtilities.getConnectedAtomsCountNonHydrogen(mol, a) > 1) {
						possiblePKcarbon = a;
					}
				}
				for(IBond b : currentFragment.getAtomContainer().getConnectedBondsList(possiblePKcarbon)) {
					if(b.contains(aminoN)) 
						continue;
					candidateBondsToBreak.add(b);
				}
				
			}
			for(IBond b : candidateBondsToBreak) {
				// Try breaking this bond, and see if this forms a linear polyketide
				mol.removeBond(b);
				IAtomContainerSet partitions = ConnectivityChecker.partitionIntoMolecules(mol);
				if(partitions.getAtomContainerCount() != 2) {
					currentFragment.getAtomContainer().addBond(b);
					continue;
				}
				
				IAtomContainer possiblePK = null;
				if(partitions.getAtomContainer(0).contains(aminoN)) {
					possiblePK = partitions.getAtomContainer(1);
				}
				else {
					possiblePK = partitions.getAtomContainer(0);
				}
				
				if(!possiblePK.contains(possiblePKcarbon)){
					if(possiblePK.contains(b.getAtom(0))){
						possiblePKcarbon = b.getAtom(0);
					}else{
						possiblePKcarbon = b.getAtom(1);
					}
				}
				
				//Extend the "end" by 1 as this underwent decarboxylation
				IAtom carbonAtom = new Atom("C");
				carbonAtom.setAtomTypeName("C.sp3");
				possiblePK.addAtom(carbonAtom);
				possiblePK.addBond(new Bond(carbonAtom, possiblePKcarbon, IBond.Order.SINGLE));
				//possiblePKcarbon = carbonAtom;
				Fragment aminoAcidFrag = null;
				Fragment pkFrag = null;
				IAtom aminoCarboxylicCarbon = null;
				IAtom startCarbon = null;
				List<Fragment> fragmentPartitions;
				boolean match = false;
				try{
					match = sm.identifyAsSmallPolyketide(possiblePK);
				}catch(Exception e){
					//System.err.println("Can't keto extend, not connected");
					currentFragment.getAtomContainer().addBond(b);
					return null;
				}
				if(match){
					fragmentPartitions = currentFragment.partitionIntoMonomerFragments();
					
					if(fragmentPartitions.get(0).getAtomContainer().contains(aminoN)) {
						aminoAcidFrag = fragmentPartitions.get(0);
						pkFrag = fragmentPartitions.get(1);
					}
					else {
						aminoAcidFrag = fragmentPartitions.get(1);
						pkFrag = fragmentPartitions.get(0);
					}
					if(b.getAtom(0).equals(possiblePKcarbon)){
						aminoCarboxylicCarbon = b.getAtom(1);
					}else{
						aminoCarboxylicCarbon = b.getAtom(0);
					}
					startCarbon = possiblePKcarbon;
				}else{
					// Do this as a check
					boolean isPK = false;
					PolyketideModulePredictor pkPredictor = null;
					try {
						pkPredictor = new PolyketideModulePredictor(
									possiblePK,
									0
									);	
						if(pkPredictor.isPK()) {
							isPK = true;
						}
					} catch (Exception e) {
						System.err.println("Error in linear PK predictor - skipping a portion of PK predictions");
						//e.printStackTrace();
						
						isPK = false;
					}
					mol.addBond(b);
					if(isPK == false) {
						continue;
					}
					//Try to find the 'carbon backbone start' for the potential PK piece
					IAtom pkCarboxylicCarbon = PolyketideModulePredictor.predictEndCarbon(possiblePK);
					if(pkCarboxylicCarbon == null){
						for(IAtom carbon : possiblePK.atoms()){
							if(PolyketideModulePredictor.isCarboxylicCarbon(possiblePK, carbon)){
								pkCarboxylicCarbon = carbon;
							}
						}
					}
					//Trace the potential carbon backbone of the potential PK portion, this is used to count the number of carbons in the backbone
					//See if it a beta or alpha amino acids
					List<IAtom> path = null;
					try{
						path = PolyketideModulePredictor.getCarbonOnlyPath(mol, pkCarboxylicCarbon, aminoN);
					}catch(Exception e){
						System.err.println("Couldn't get carbon only path");
						continue;
					}
					if(path.size() < 5) {
						return null;
					}
					if(path.size() % 2 != 0){
						currentFragment.getAtomContainer().removeBond(path.get(path.size()-3), path.get(path.size()-4));
						aminoCarboxylicCarbon = path.get(path.size()-3);
						startCarbon = path.get(path.size()-4);
					}else{ //out dated
						currentFragment.getAtomContainer().removeBond(path.get(path.size()-4), path.get(path.size()-5));
						aminoCarboxylicCarbon = path.get(path.size()-4);
						startCarbon = path.get(path.size()-5);
					}
					
					fragmentPartitions = currentFragment.partitionIntoMonomerFragments();
				}
				
				if(fragmentPartitions.get(0).getAtomContainer().contains(aminoN)) {
					aminoAcidFrag = fragmentPartitions.get(0);
					pkFrag = fragmentPartitions.get(1);
				}
				else {
					aminoAcidFrag = fragmentPartitions.get(1);
					pkFrag = fragmentPartitions.get(0);
				}
				
				//Forming carboyxylic acid with the aminoCarboxylicCarbon
				//check what the carbon is connected to
				int singleBondO = 0;
				int doubleBondO = 0;
				for(IAtom atom : aminoAcidFrag.getAtomContainer().getConnectedAtomsList(aminoCarboxylicCarbon)){
					if(atom.getAtomicNumber() == 1 || atom.getAtomTypeName() == null){
						continue;
					}
					if(atom.getAtomTypeName().equals("O.sp3")){
						for(IAtom connectedToOxygen : aminoAcidFrag.getAtomContainer().getConnectedAtomsList(atom)){
							if(!connectedToOxygen.equals(aminoCarboxylicCarbon)){
								if(aminoAcidFrag.getAtomContainer().getConnectedAtomsCount(connectedToOxygen) == 1 && connectedToOxygen.getAtomicNumber() == 6){
									aminoAcidFrag.getAtomContainer().removeAtomAndConnectedElectronContainers(connectedToOxygen);
									aminoAcidFrag.addTailoringDomain(TailoringDomainEnums.O_METHYLTRANSFERASE);
								}
							}
						}
							singleBondO ++;
					}
					if(atom.getAtomTypeName().equals("O.sp2")){
						doubleBondO ++;
					}
				}

				//add if not there
				if(singleBondO == 0){
					Atom firstO = new Atom("O");
					firstO.setAtomTypeName("O.sp3");
					aminoAcidFrag.getAtomContainer().addAtom(firstO);
					aminoAcidFrag.getAtomContainer().addBond(new Bond(aminoCarboxylicCarbon, firstO, IBond.Order.SINGLE));
				}
				
				//add if not there
				if(doubleBondO == 0){
					Atom secondO = new Atom("O");
					secondO.setAtomTypeName("O.sp2");
					aminoAcidFrag.getAtomContainer().addAtom(secondO);
					aminoAcidFrag.getAtomContainer().addBond(new Bond(aminoCarboxylicCarbon, secondO, IBond.Order.DOUBLE));
				}
				//Add all aminoNs to the new fragments
				for(IAtom aminoNtoAdd : currentFragment.getAminoNs()){
					if(pkFrag.getAtomContainer().contains(aminoNtoAdd)){
						if(!pkFrag.getAminoNs().contains(aminoNtoAdd)){
							pkFrag.addAminoN(aminoNtoAdd);
						}
					}else{
						if(!aminoAcidFrag.getAminoNs().contains(aminoNtoAdd)){
							aminoAcidFrag.addAminoN(aminoNtoAdd);
						}
					}
				}
				
				aminoAcidFrag.addAminoC(aminoCarboxylicCarbon);
				aminoAcidFrag.setAtomAfterCTerminus(startCarbon);
				pkFrag.setAtomAfterNTerminus(aminoCarboxylicCarbon);
				pkFrag.setAsKetoExtension();
				identifyAsLinearPK(pkFrag, 0);				
				pieces[0] = aminoAcidFrag;
				pieces[1] = pkFrag;
				return pieces;
			}
		}
		return null;
	}	
	
	
}
