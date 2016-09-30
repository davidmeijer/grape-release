package ca.mcmaster.magarveylab.grape.nrp.chem;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import ca.mcmaster.magarveylab.grape.enums.DomainEnums.TailoringDomainEnums;
import ca.mcmaster.magarveylab.grape.enums.MoleculeClasses.ChemicalSubType;
import ca.mcmaster.magarveylab.grape.enums.MoleculeClasses.ChemicalType;
import ca.mcmaster.magarveylab.grape.nrp.chem.Fragment.FragmentType;
import ca.mcmaster.magarveylab.grape.nrp.chem.matcher.UnusualPolyketideIdentifier;
import ca.mcmaster.magarveylab.grape.nrp.chem.matcher.SubstrateMatcher;
import ca.mcmaster.magarveylab.grape.pk.chem.PolyketideWalker;
import ca.mcmaster.magarveylab.grape.pk.chem.MacrolideChecker;


/**
 * NRP Predictor.
 * Given a nonribosomal peptide molecule, determine its monomeric constituents.
 * @author gmchen cDejong
 *
 */
public class MoleculePredictor {
	private SubstrateMatcher substrateMatcher;
	private UnusualPolyketideIdentifier unusualPolyketideIdentifier;
	private PolyketideWalker polyketideWalker;
	
	/**
	* Constructor for this GrapePredictor instance.
	*/
	public MoleculePredictor(String aminoAcidPath) {
		substrateMatcher = new SubstrateMatcher(aminoAcidPath);
		unusualPolyketideIdentifier = new UnusualPolyketideIdentifier(substrateMatcher);
		polyketideWalker = new PolyketideWalker();
	}
	
	public ChemicalAbstraction getChemicalAbstraction(IAtomContainer molecule, String name) throws CDKException {
		// Get macrolide type
		IAtomContainer currentMoleculeClone = null;
		try {
			currentMoleculeClone = molecule.clone();
		} catch (CloneNotSupportedException e1) {
			e1.printStackTrace();
		}
		
		// TODO change to enum
		// 0 for not macrolide, 1 for lactone, 2 for latam
		int macrolideType = MacrolideChecker.isMacrolide(currentMoleculeClone);
		
		ChemicalSubType PKType = ChemicalSubType.STANDARD;
		String scaffold = null;
		try {
			polyketideWalker.analyzeMolecule(molecule, macrolideType);
			PKType = polyketideWalker.getChemicalSubType();
			scaffold = polyketideWalker.getMatchName();
			
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		ChemicalAbstraction chemicalAbstraction = new ChemicalAbstraction(name);
		List<Fragment> monomerFragments = getMonomerFragments(molecule, macrolideType, PKType, chemicalAbstraction);
		for(Fragment m : monomerFragments){
			if(m.getFragmentType() == null) {
				m.setFragmentType(FragmentType.UNKNOWN_OTHER);
			}
		}
		chemicalAbstraction.setMonomerFragments(monomerFragments);
		Map<FragmentType, Integer> numOfType = new HashMap<FragmentType, Integer>();
		for(FragmentType type : FragmentType.values()){
			numOfType.put(type, 0);
		}
		for(Fragment m : monomerFragments) {
			numOfType.put(m.getFragmentType(), numOfType.get(m.getFragmentType()) + 1);
		}
		
		// Only allow FA_OR_PK pieces for compounds with at least two amino acids
		if(numOfType.get(FragmentType.AMINO_ACID) < 2) {
			for(Fragment m : chemicalAbstraction.getMonomerFragments()) {
				if(m.getFragmentType() == FragmentType.FA_OR_PK) {
					m.setFragmentType(FragmentType.POLYKETIDE);
				}
			}
		}
		
		// Set the type
		// If there contains a sulfur beta lactam, set as NRP
		if(PKType != ChemicalSubType.STANDARD) {
			chemicalAbstraction.setChemicalSubType(PKType);
			chemicalAbstraction.setChemicalType(ChemicalType.PK);
			if(scaffold != null){
				chemicalAbstraction.setChemicalScaffoldName(scaffold);
			}
		}
		else if (macrolideType != 0) {
			chemicalAbstraction.setChemicalSubType(ChemicalSubType.MACROLIDE);
			if (macrolideType == 1) {
				chemicalAbstraction.setChemicalType(ChemicalType.PK);
			}else {
				chemicalAbstraction.setChemicalType(ChemicalType.NRP_PK_HYBRID);
			}
		}
		else if(numOfType.get(FragmentType.POLYKETIDE) + numOfType.get(FragmentType.FA_OR_PK) > 0 && numOfType.get(FragmentType.AMINO_ACID) < 2) {
			chemicalAbstraction.setChemicalSubType(ChemicalSubType.TYPE_1_PK);
			chemicalAbstraction.setChemicalType(ChemicalType.PK);

		}
		else if(numOfType.get(FragmentType.POLYKETIDE) > 0 || numOfType.get(FragmentType.FA_OR_PK) > 0) { // this means at least three amino acids
			chemicalAbstraction.setChemicalType(ChemicalType.NRP_PK_HYBRID);
		}
		else if(numOfType.get(FragmentType.AMINO_ACID) > 0) {
			chemicalAbstraction.setChemicalType(ChemicalType.NRP);
		}
		if(chemicalAbstraction.getChemicalType() == null) {
			// Check if aminoglycoside
			boolean hasCyclohexane = false;
			if(numOfType.get(FragmentType.SUGAR) > 0) {
				for(Fragment m : chemicalAbstraction.getMonomerFragments()) {
					if(m.getFragmentType() == FragmentType.UNKNOWN_OTHER) {
						substrateMatcher.identifyAsAminoglycosideCyclohexane(m);
						if(m.getFragmentType() == FragmentType.CYCLOHEXANE) {
							hasCyclohexane = true;
						}
					}
				}
			}
			if(hasCyclohexane) {
				chemicalAbstraction.setChemicalType(ChemicalType.AMINOGLYCOSIDE);
			}
		}
		if(chemicalAbstraction.getChemicalType() == null) {
			chemicalAbstraction.setChemicalType(ChemicalType.UNKNOWN);
		}
		for(Fragment m : chemicalAbstraction.getMonomerFragments()) {
			if(m.getTailoringDomains().contains(TailoringDomainEnums.SULFUR_BETA_LACTAM)) {
				chemicalAbstraction.setChemicalType(ChemicalType.NRP);
			}
		}
		
		return chemicalAbstraction;
	}
	
	/**
	 * Get an ordered list of monomer fragmentsSs
	 * @return
	 * @throws CDKException 
	 */
	public List<Fragment> getMonomerFragments(IAtomContainer nrp, int macrolideType, ChemicalSubType PKType, ChemicalAbstraction chemicalAbstraction) throws CDKException {
		
		try {
			nrp = nrp.clone();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		
		MoleculeModifier moleculeModifier = new MoleculeModifier(nrp, chemicalAbstraction);
		moleculeModifier.performAllNrpModifications();
		
		List<Fragment> monomerFragments = moleculeModifier.getFragments();
		
		for(int i = 0; i < monomerFragments.size(); i++) {
			Fragment fragmentInConsideration = monomerFragments.get(i);
			int numNitrogens = 0;
			for(int j = 0; j < fragmentInConsideration.getAtomContainer().getAtomCount(); j++) {
				if(fragmentInConsideration.getAtomContainer().getAtom(j).getAtomicNumber() == 7) {
					numNitrogens++;
				}
			}
			if(numNitrogens > 0) {
				fragmentInConsideration.setNumNitrogens(numNitrogens);
			}
		}
		for(int i = 0; i < monomerFragments.size(); i++) {
			Fragment currentFragment = monomerFragments.get(i);
			AtomContainerManipulator.percieveAtomTypesAndConfigureUnsetProperties(currentFragment.getAtomContainer());
			currentFragment.addImplicitHydrogens();
			if(currentFragment.getFragmentType() == FragmentType.MULTIPLE_AMINO_ACID_PIECE) {
				currentFragment.setFragmentType(null);				
			}
			if(currentFragment.getFragmentType() == null) {
				substrateMatcher.identifyAsMultipleAminoAcids(currentFragment);
			}			
			if(currentFragment.getFragmentType() == null) {
				substrateMatcher.identifyAsPerfectAminoAcid(currentFragment);
			}
			if(currentFragment.getFragmentType() == null || currentFragment.getFragmentType() == FragmentType.SUGAR) {
				substrateMatcher.identifyAsSugar(currentFragment);
			}
			if(currentFragment.getFragmentType() == null){
				substrateMatcher.identifyAsSmallPolyketide(currentFragment);
			}
			if(currentFragment.getFragmentType() == null) {
				substrateMatcher.identifyAsStandardFattyAcid(currentFragment);
			}
			if(currentFragment.getFragmentType() == null) {
				substrateMatcher.identifyAsAcylAdenlyatingSubstrate(currentFragment);
			}
			if(currentFragment.getFragmentType() == null) {
				substrateMatcher.identifyAsAminoAcid(currentFragment);
			}
			if(currentFragment.getFragmentType() == null) {
				Fragment[] pieces = unusualPolyketideIdentifier.getKetoextendedAAPieces(currentFragment, macrolideType);
				if(pieces != null) {
					Fragment CFragment = pieces[0];
					Fragment NFragment = pieces[1];
					int indexToReplace = monomerFragments.indexOf(currentFragment);
					monomerFragments.remove(indexToReplace);
					monomerFragments.add(indexToReplace, CFragment);
					monomerFragments.add(indexToReplace + 1, NFragment);					
					i--; // check these fragments again since the pk piece could be keto extended in the other direction, and the amino acid portion may have other modifications
					continue;
				}
			}
			if(currentFragment.getFragmentType() == null && !PKType.equals(ChemicalSubType.ENEDYINE) && !PKType.equals(ChemicalSubType.TYPE_2) && !PKType.equals(ChemicalSubType.NON_TYPE_2_AROMATIC) && !PKType.equals(ChemicalSubType.TERPENE)){
				// Check if it is a possible polyketide.
				unusualPolyketideIdentifier.identifyAsLinearPK(currentFragment, macrolideType);
			}
			if(currentFragment.getFragmentType() == null && !PKType.equals(ChemicalSubType.ENEDYINE) && !PKType.equals(ChemicalSubType.TYPE_2) && !PKType.equals(ChemicalSubType.NON_TYPE_2_AROMATIC) && !PKType.equals(ChemicalSubType.TERPENE)) {
				//identifyAsFattyAcid(currentFragment);
			}
			if(currentFragment.getFragmentType() == null) {
				substrateMatcher.identifyAsKnownOther(currentFragment);
			}
			if(currentFragment.getFragmentType() == null || currentFragment.getFragmentType() == FragmentType.POLYKETIDE || currentFragment.getFragmentType() == FragmentType.FA_OR_PK){
				substrateMatcher.identifySubstrctureAsStarter(currentFragment);
			}
			if(currentFragment.getFragmentType() == null && currentFragment.getStarters().size() > 0){
				currentFragment.setFragmentType(FragmentType.KNOWN_OTHER);
			}
			
			if(currentFragment.getFragmentType() == null && !PKType.equals(ChemicalSubType.ENEDYINE) && !PKType.equals(ChemicalSubType.TYPE_2) && !PKType.equals(ChemicalSubType.NON_TYPE_2_AROMATIC) && !PKType.equals(ChemicalSubType.TERPENE)) {				
				if(currentFragment.getNumNitrogens() > 0 && currentFragment.getAminoCs().size() > 0 || currentFragment.getAminoNs().size() > 1) {
					currentFragment.setFragmentType(FragmentType.MULTIPLE_AMINO_ACID_PIECE);
				}
			}
			// If there are at least four carbons, two nitrogens, and two oxygens,
			// then we consider the possibility that this consists of multiple amino acids
			
			if(currentFragment.getFragmentType() == null) {
				int numCarbons = 0;
				int numNitrogens = 0;
				int numOxygens = 0;
				for(int j = 0; j < currentFragment.getAtomContainer().getAtomCount(); j++) {
					switch(currentFragment.getAtomContainer().getAtom(j).getAtomicNumber()) {
					case 6:
						numCarbons++;
						break;
					case 7:
						numNitrogens++;
						break;
					case 8:
						numOxygens++;
						break;
					}
				}
				if(numCarbons >= 4 && numNitrogens >= 2 && numOxygens >= 4) {
					currentFragment.setFragmentType(FragmentType.MULTIPLE_AMINO_ACID_PIECE);
				}
			}
			
			if(currentFragment.getFragmentType() == null) {
				currentFragment.setFragmentType(FragmentType.UNKNOWN_OTHER);
			}
		}
		return(monomerFragments);
	}

}
