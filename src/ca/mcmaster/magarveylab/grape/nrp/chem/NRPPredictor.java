package ca.mcmaster.magarveylab.grape.nrp.chem;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Scanner;
import java.util.Set;

import org.openscience.cdk.Atom;
import org.openscience.cdk.Bond;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;

import ca.mcmaster.magarveylab.grape.enums.AcylAdenylatingSubstrates;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums;
import ca.mcmaster.magarveylab.grape.enums.FattyAcidsEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.AminoAcidEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.MultipleAminoAcidEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.PolyKetideDomainEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.SmallPKunits;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.SugarModificationsEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.TailoringDomainEnums;
import ca.mcmaster.magarveylab.grape.enums.KnownOtherEnums;
import ca.mcmaster.magarveylab.grape.enums.MoleculeClasses.ChemicalSubType;
import ca.mcmaster.magarveylab.grape.enums.MoleculeClasses.ChemicalType;
import ca.mcmaster.magarveylab.grape.enums.SugarEnums;
import ca.mcmaster.magarveylab.grape.nrp.chem.Fragment.FragmentType;
import ca.mcmaster.magarveylab.grape.pk.chem.BackboneAnalyser;
import ca.mcmaster.magarveylab.grape.pk.chem.CatigorizeOtherPK;
import ca.mcmaster.magarveylab.grape.pk.chem.PolyketideModulePredictor;
import ca.mcmaster.magarveylab.grape.pk.chem.MacrolideChecker;
import ca.mcmaster.magarveylab.grape.pk.modules.PKsubstrate;
import ca.mcmaster.magarveylab.grape.util.ChemicalUtilities;
import ca.mcmaster.magarveylab.grape.util.io.SmilesIO;


/**
 * NRP Predictor.
 * Given a nonribosomal peptide molecule, determine its monomeric constituents.
 * @author gmchen cDejong
 *
 */
public class NRPPredictor {
	private ArrayList<IMolecule> aminoAcids;
	private ArrayList<AminoAcidEnums> aminoAcidEnums;
	private ArrayList<ArrayList<TailoringDomainEnums>> aminoAcidTailoringEnums;
	private Map<String, IMolecule> sugarsMap = SugarEnums.getAll();
	private ArrayList<IMolecule> multipleAminoAcidMols = MultipleAminoAcidEnums.allMols();
	private ArrayList<AminoAcidEnums[]> mutlipleAminoAcidConstituents = MultipleAminoAcidEnums.allAminoAcids();
	private IMolecule[] fattyAcids = FattyAcidsEnums.mols();
	private CatigorizeOtherPK nonStandardPKIdentifier;
	private boolean fungal = false;

	/**
	* Constructor for this GrapePredictor instance.
	*/
	public NRPPredictor(String aminoAcidPath) {
		aminoAcids = new ArrayList<IMolecule>();
		aminoAcidEnums = new ArrayList<AminoAcidEnums>();
		aminoAcidTailoringEnums = new ArrayList<ArrayList<TailoringDomainEnums>>();
		nonStandardPKIdentifier = new CatigorizeOtherPK();
		
		readAminoAcidInput(aminoAcidPath);

		//System.out.println("Finished reading input.");
	}
	
	public NRPPredictor(String aminoAcidPath, boolean fungal) {
		aminoAcids = new ArrayList<IMolecule>();
		aminoAcidEnums = new ArrayList<AminoAcidEnums>();
		aminoAcidTailoringEnums = new ArrayList<ArrayList<TailoringDomainEnums>>();
		nonStandardPKIdentifier = new CatigorizeOtherPK();
		this.fungal = fungal;
		
		readAminoAcidInput(aminoAcidPath);

		//System.out.println("Finished reading input.");
	}
	
	public ChemicalAbstraction getChemicalAbstraction(IMolecule nrp) {
		// Get macrolide type
		IMolecule currentMoleculeClone = null;
		try {
			currentMoleculeClone = nrp.clone();
		} catch (CloneNotSupportedException e1) {
			e1.printStackTrace();
		}
		
		int macrolideType = MacrolideChecker.isMacrolide(currentMoleculeClone); // 0 for not macrolide, 1 for lactone, 2 for latam
		
		ChemicalSubType PKType = ChemicalSubType.STANDARD;
		String scaffold = null;
		try {
			nonStandardPKIdentifier.analyzeMolecule(nrp, macrolideType);
			PKType = nonStandardPKIdentifier.getChemicalSubType();
			scaffold = nonStandardPKIdentifier.getMatchName();
			
		} catch (CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		ChemicalAbstraction chemicalAbstraction = new ChemicalAbstraction();
		List<Fragment> monomerFragments = getMonomerFragments(nrp, macrolideType, PKType, chemicalAbstraction);
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
		
		// If there are fewer than 3 amino acids, set any FA_OR_PK pieces as polyketide
		/*
		if(numOfType.get(FragmentType.AMINO_ACID) < 3) {
			for(Fragment m : chemicalAbstraction.getMonomerFragments()) {
				if(m.getMonomerType() == FragmentType.FA_OR_PK) {
					m.setMonomerType(FragmentType.POLYKETIDE);
				}
			}
		}
		*/
		
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
						identifyAsAminoglycosideCyclohexane(m);
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
	 */
	public List<Fragment> getMonomerFragments(IMolecule nrp, int macrolideType, ChemicalSubType PKType, ChemicalAbstraction chemicalAbstraction) {
		
		try {
			nrp = nrp.clone();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		
		NRPModifier nrpModifier = new NRPModifier(nrp, chemicalAbstraction, fungal);
		
		List<Fragment> monomerFragments = nrpModifier.getMonomerFragments();
		
		for(int i = 0; i < monomerFragments.size(); i++) {
			Fragment fragmentInConsideration = monomerFragments.get(i);
			int numNitrogens = 0;
			for(int j = 0; j < fragmentInConsideration.getMolecule().getAtomCount(); j++) {
				if(fragmentInConsideration.getMolecule().getAtom(j).getAtomicNumber() == 7) {
					numNitrogens++;
				}
			}
			if(numNitrogens > 0) {
				fragmentInConsideration.setNumNitrogens(numNitrogens);
			}
		}
		
		for(int i = 0; i < monomerFragments.size(); i++) {
			Fragment currentFragment = monomerFragments.get(i);
			if(currentFragment.getFragmentType() == FragmentType.MULTIPLE_AMINO_ACID_PIECE) {
				currentFragment.setFragmentType(null);				
			}
			if(currentFragment.getFragmentType() == null) {
				identifyAsMultipleAminoAcids(currentFragment);
			}			
			if(currentFragment.getFragmentType() == null) {
				identifyAsPerfectAminoAcid(currentFragment);
			}
			if(currentFragment.getFragmentType() == null || currentFragment.getFragmentType() == FragmentType.SUGAR) {
				identifyAsSugar(currentFragment);
			}
			if(currentFragment.getFragmentType() == null){
				identifyAsSmallPolyketide(currentFragment);
			}
			if(currentFragment.getFragmentType() == null) {
				identifyAsStandardFattyAcid(currentFragment);
			}
			if(currentFragment.getFragmentType() == null) {
				identifyAsAcylAdenlyatingSubstrate(currentFragment);
			}
			if(currentFragment.getFragmentType() == null) {
				Fragment[] pieces = getKetoextendedAAPieces(currentFragment, macrolideType);
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
			if(currentFragment.getFragmentType() == null) {
				identifyAsAminoAcid(currentFragment);
			}
			if(currentFragment.getFragmentType() == null && !PKType.equals(ChemicalSubType.ENEDYINE) && !PKType.equals(ChemicalSubType.TYPE_2) && !PKType.equals(ChemicalSubType.NON_TYPE_2_AROMATIC) && !PKType.equals(ChemicalSubType.TERPENE)){
				// Check if it is a possible polyketide.
				identifyAsLinearPK(currentFragment, macrolideType);	
			}
			if(currentFragment.getFragmentType() == null && !PKType.equals(ChemicalSubType.ENEDYINE) && !PKType.equals(ChemicalSubType.TYPE_2) && !PKType.equals(ChemicalSubType.NON_TYPE_2_AROMATIC) && !PKType.equals(ChemicalSubType.TERPENE)) {
				//identifyAsFattyAcid(currentFragment);
			}
			if(currentFragment.getFragmentType() == null) {
				identifyAsKnownOther(currentFragment);
			}
			if(currentFragment.getFragmentType() == null || currentFragment.getFragmentType() == FragmentType.POLYKETIDE || currentFragment.getFragmentType() == FragmentType.FA_OR_PK){
				identifySubstrctureAsStarter(currentFragment);
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
			/*
			if(currentFragment.getMonomerType() == null) {
				int numCarbons = 0;
				int numNitrogens = 0;
				int numOxygens = 0;
				for(int j = 0; j < currentFragment.getMolecule().getAtomCount(); j++) {
					switch(currentFragment.getMolecule().getAtom(j).getAtomicNumber()) {
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
					currentFragment.setMonomerType(FragmentType.MULTIPLE_AMINO_ACID_PIECE);
				}
			}
			*/
			if(currentFragment.getFragmentType() == null) {
				currentFragment.setFragmentType(FragmentType.UNKNOWN_OTHER);
			}
		}
		//for(int i = 0; i < monomerFragments.size(); i++) {
		//	System.out.println(monomerFragments.get(i).getMonomerType());
		//}
		return(monomerFragments);
	}
	
	private void identifyAsKnownOther(Fragment fragment) {
		for(KnownOtherEnums knownOther : KnownOtherEnums.values()){
			if(knownOther.getMol() != null){
				double score = ChemicalUtilities.getTanimotoScore(fragment.getMolecule(), knownOther.getMol());
				if(score > 0.9){
					fragment.setFragmentType(FragmentType.KNOWN_OTHER);
					fragment.setKnownOther(knownOther);
					fragment.setTanimotoScore(score);
					return;
				}
			}
		}
		
	}

	private void identifyAsAcylAdenlyatingSubstrate(Fragment fragment) {
		
		for(AcylAdenylatingSubstrates starter  : AcylAdenylatingSubstrates.values()){
			IMolecule[] mols = starter.mols();
			for(int i = 0; mols.length > i; i++){
				if(ChemicalUtilities.getTanimotoScore(fragment.getMolecule(), mols[i]) > 0.9){
					fragment.addStarter(starter);
					fragment.setTanimotoScore(ChemicalUtilities.getTanimotoScore(fragment.getMolecule(), mols[i]));
					fragment.setFragmentType(FragmentType.ACYL_ADENYLATING);
					return;
				}
			}
		}		
	}

	private void identifyAsSmallPolyketide(Fragment currentFragment) {
		
		for(SmallPKunits smallPK  : SmallPKunits.values()){
			if(ChemicalUtilities.getTanimotoScore(currentFragment.getMolecule(), smallPK.mol()) > 0.95){
				currentFragment.setFragmentType(FragmentType.POLYKETIDE);
				currentFragment.setTanimotoScore(ChemicalUtilities.getTanimotoScore(currentFragment.getMolecule(), smallPK.mol()));
				List<PolyKetideDomainEnums> pkDomain = new ArrayList<PolyKetideDomainEnums>();
				for(int i = 0; smallPK.pkDomains().length > i; i++){
					pkDomain.add(smallPK.pkDomains()[i]);
				}
				List<List<PolyKetideDomainEnums>> pkDomains = new ArrayList<List<PolyKetideDomainEnums>>();
				pkDomains.add(pkDomain);
				currentFragment.setPkDomains(pkDomains);
				
				List<PKsubstrate> pkSubstrate = new ArrayList<PKsubstrate>();
				for(int i = 0; smallPK.pkSubstrates().length > i; i++){
					pkSubstrate.add(smallPK.pkSubstrates()[i]);
				}
				currentFragment.setLoadingUnits(pkSubstrate);
				return;
			}
		}
	}

	private void identifySubstrctureAsStarter(Fragment fragment) {
		IMolecule mol = null;
		try {
			mol = SmilesIO.readSmiles(SmilesIO.generateSmiles(fragment.getMolecule()));
		} catch (IOException | CDKException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		for(AcylAdenylatingSubstrates starter  : AcylAdenylatingSubstrates.values()){
			if(!starter.canCheckForSubstructure()) continue;
			try {
				IMolecule[] mols = starter.mols();
				for(int i = 0; mols.length > i; i++){
					if(UniversalIsomorphismTester.isSubgraph(mol, mols[i])){
						fragment.addStarter(starter);
						return;
					}
				}
			} catch (CDKException e){}
		}
	}

	private void identifyAsMultipleAminoAcids(Fragment fragment) {
		ArrayList<Double> tanimotoScores = ChemicalUtilities.getTanimotoScores(fragment.getMolecule(), multipleAminoAcidMols);
		int highestScoreIndex = ChemicalUtilities.getHighestTanimotoScoreIndex(tanimotoScores);
		
		if(tanimotoScores.get(highestScoreIndex) > 0.95) {
			fragment.addAminoAcidDomains(mutlipleAminoAcidConstituents.get(highestScoreIndex));
			fragment.setFragmentType(FragmentType.MULTIPLE_AMINO_ACID_PIECE);
			fragment.setTanimotoScore(tanimotoScores.get(highestScoreIndex));
		}
	}
	
	public void identifyAsPerfectAminoAcid(Fragment fragment) {
		ArrayList<Double> aminoAcidTanimotoScores = ChemicalUtilities.getTanimotoScores(fragment.getMolecule(), aminoAcids);
		int highestScoreIndex = ChemicalUtilities.getHighestTanimotoScoreIndex(aminoAcidTanimotoScores);
		if(aminoAcidTanimotoScores.get(highestScoreIndex) == 1.0) {
			fragment.addAminoAcidDomain(aminoAcidEnums.get(highestScoreIndex));
			fragment.setFragmentType(FragmentType.AMINO_ACID);
			fragment.setTanimotoScore(aminoAcidTanimotoScores.get(highestScoreIndex));
		}
	}
	
	public void identifyAsAminoAcid(Fragment fragment) {
		ArrayList<Double> aminoAcidTanimotoScores = ChemicalUtilities.getTanimotoScores(fragment.getMolecule(), aminoAcids);
		int highestScoreIndex = ChemicalUtilities.getHighestTanimotoScoreIndex(aminoAcidTanimotoScores);
		if(aminoAcidTanimotoScores.get(highestScoreIndex) >=  0.9) {
			fragment.addAminoAcidDomain(aminoAcidEnums.get(highestScoreIndex));
			fragment.setFragmentType(FragmentType.AMINO_ACID);
			fragment.setTanimotoScore(aminoAcidTanimotoScores.get(highestScoreIndex));
		}
		// If this piece contains a thiazole, check if it contains a cysteine as a substructure
		if(fragment.getTailoringDomains().contains(TailoringDomainEnums.THIAZOLE) ||
				fragment.getTailoringDomains().contains(TailoringDomainEnums.SULFUR_BETA_LACTAM)) {
			try {
				IMolecule cysteineSubstructure = SmilesIO.readSmiles("C(C(C(=O)O)N)S");
				if(UniversalIsomorphismTester.isSubgraph(fragment.getMolecule(), cysteineSubstructure)) {
					fragment.setFragmentType(FragmentType.AMINO_ACID);
					fragment.addAminoAcidDomain(AminoAcidEnums.Cysteine);
					double tanimoto = ChemicalUtilities.getTanimotoScore(fragment.getMolecule(), cysteineSubstructure);
					fragment.setTanimotoScore(tanimoto);
				}
			} catch (IOException e)  {
				
			}
			catch (CDKException e) {
				//e.printStackTrace();
			}
		}
	}
	/**
	 * Check if a molecule matches an amino acid
	 * @param tanimotoScores
	 * @return
	 */
	
	public boolean isAminoAcid(IMolecule molecule) {
		ArrayList<Double> aminoAcidTanimotoScores = ChemicalUtilities.getTanimotoScores(molecule, aminoAcids);
		int highestScoreIndex = ChemicalUtilities.getHighestTanimotoScoreIndex(aminoAcidTanimotoScores);
		if(aminoAcidTanimotoScores.get(highestScoreIndex) >= 0.9) {
			return true;
		}
		return false;
	}
	
	/**
	 * First, check if a monomer fragment is a sugar if it has not already been identified.
	 * Then, if it is identified as a sugar, determine sugar modifications
	 * @param fragment
	 */
	public void identifyAsSugar(Fragment fragment) {
		
		
		String bestMatchName = null;
		double highestScore = 0;
		
		for(Entry<String, IMolecule> sugar : sugarsMap.entrySet()){
			double score = ChemicalUtilities.getTanimotoScore(fragment.getMolecule(), sugar.getValue());
			if(score > highestScore){
				bestMatchName = sugar.getKey();
				highestScore = score;
			}
		}
		
		
		if(highestScore > 0.8) {
			fragment.setFragmentType(FragmentType.SUGAR);
			fragment.addSugarName(bestMatchName);
			fragment.setTanimotoScore(highestScore);
		}
		
		
		//define sugar branches
		if(fragment.getFragmentType() == FragmentType.SUGAR) {
			// Identify sugar modifications
			// Check for nitrogens
			for(int i = 0; i < fragment.getMolecule().getAtomCount(); i++) {
				if(fragment.getMolecule().getAtom(i).getAtomicNumber() == 7) {
					fragment.addSugarModification(SugarModificationsEnums.AMINOTRANSFER);
				}
				// N Methylations
				SMARTSQueryTool querytool = null;
				try {
					querytool = new SMARTSQueryTool("N[C;D1]");
				} catch (CDKException e) {
					e.printStackTrace();
				}
				boolean hasNMethyl = false;
				try {
					hasNMethyl = querytool.matches(fragment.getMolecule());
				} catch (CDKException e) {
					e.printStackTrace();
				}
				if(hasNMethyl) {
					// Check
					fragment.addSugarModification(SugarModificationsEnums.N_METHYLATION);
				}
				
				// O Methylations
				querytool = null;
				try {
					querytool = new SMARTSQueryTool("O[C;D1]");
				} catch (CDKException e) {
					e.printStackTrace();
				}
				boolean hasOMethyl = false;
				try {
					hasOMethyl = querytool.matches(fragment.getMolecule());
				} catch (CDKException e) {
					e.printStackTrace();
				}
				if(hasOMethyl) {
					// Check
					fragment.addSugarModification(SugarModificationsEnums.O_METHYLATION);
				}
				
			}
		}
	}
	
	public void identifyAsFattyAcid(Fragment fragment) {
		// Check if this matches a known fatty acid
		if(fragment.getFragmentType() == null) {
			double highestFattyAcidScore = ChemicalUtilities.getHighestTanimotoScore(fragment.getMolecule(), Arrays.asList(fattyAcids));
			if(highestFattyAcidScore > 0.8) {
				fragment.setFragmentType(FragmentType.FATTY_ACID);
				fragment.setTanimotoScore(highestFattyAcidScore);
			}
		}
		
		/*
		// Check that the atoms are C, H, or O
		for(int i = 0; i < fragment.getMolecule().getAtomCount(); i++) {
			int atomicNumber = fragment.getMolecule().getAtom(i).getAtomicNumber();
			if(atomicNumber != 1 && atomicNumber != 6 && atomicNumber != 8) {
				return;
			}
		}
		
		SMARTSQueryTool querytool = null;
		try {
			querytool = new SMARTSQueryTool("C[C;D2][C;D2][C;D2][C;D2]C");
		} catch (CDKException e) {
			e.printStackTrace();
		}
		boolean matches = false;
		try {
			matches = querytool.matches(fragment.getMolecule());
		} catch (CDKException e) {
			e.printStackTrace();
		}
		if(matches) {
			// Check
			fragment.setMonomerType(FragmentType.FATTY_ACID);
		}
		*/
	}
	
	/**
	 * If this fragment appears to be a fatty acid or 3-OH fatty acid, change its FragmentType to fatty acid.
	 * A molecule will be identified as a fatty acid or subtructure of the first four carbons at the 
	 * carboxylic acid end of a normal carboxylic fatty acid or 3-OH fatty acid, contains only C, H, and O, and has no more than three oxygens.
	 * @param fragment
	 */
	public void identifyAsStandardFattyAcid(Fragment fragment) {
		// exit if the molecule contains anything other than a carbon, oxygen, hydrogen
		for(int i = 0; i < fragment.getMolecule().getAtomCount(); i++) {
			int atomicNumber = fragment.getMolecule().getAtom(i).getAtomicNumber();
			if(atomicNumber != 1 && atomicNumber != 6 && atomicNumber != 8) {
				return;
			}
		}
		IMolecule fattyAcid3OHTemplate = null;
		IMolecule fattyAcidTemplate = null;
		try {
			fattyAcidTemplate = SmilesIO.readSmiles("CC(O)=O");
			fattyAcid3OHTemplate = SmilesIO.readSmiles("CC(O)CC(O)=O");
		} catch (IOException | CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		boolean matchesNormalFattyAcid = false;
		boolean matches3OHFattyAcid = false;
		try {
			
			matchesNormalFattyAcid = UniversalIsomorphismTester.isSubgraph(fragment.getMolecule(), fattyAcidTemplate);
			matches3OHFattyAcid = UniversalIsomorphismTester.isSubgraph(fragment.getMolecule(), fattyAcid3OHTemplate);
		} catch (CDKException e) {
			return;
		}
		if(!matchesNormalFattyAcid && !matches3OHFattyAcid) {
			return;
		}
		// Count the number of carbons and oxygens, and quit if a non carbon, oxygen, or hydrogen is found
		int numOxygens = 0;
		int numCarbons = 0;
		//int numTertiaryCarbons = 0;
		for(int i = 0; i < fragment.getMolecule().getAtomCount(); i++) {
			int atomicNumber = fragment.getMolecule().getAtom(i).getAtomicNumber();
			if(atomicNumber == 8) {
				numOxygens++;
			}
			if(atomicNumber == 6) {
				numCarbons++;
				if(fragment.getMolecule().getConnectedAtomsCount(fragment.getMolecule().getAtom(i)) == 3) {
				//	numTertiaryCarbons++;
				}
			}
		}
		
		// Check that there are no carbon rings
		ArrayList<IAtom> atomsInCarbonRings = ChemicalUtilities.getAtomsInCarbonRings(fragment.getMolecule());
		if(!atomsInCarbonRings.isEmpty()) {
			return;
		}
		// Must be at least four carbons
		if(numCarbons < 4) {
			return;
		}
		
		// If there the 3-OH fatty acid substructure is not present, then at most two oxygens are permitted
		if(!matches3OHFattyAcid && numOxygens > 2) {
			return;
		}
		// Otherwise, at most three oxygens may be present
		if(numOxygens > 3) {
			return;
		}
		
		// Check that there are at least four carbons in a row
		SMARTSQueryTool querytool = null;
		try {
			querytool = new SMARTSQueryTool("C[C;D2][C;D2][C;D2][C;D2]C");
		} catch (CDKException e) {
			e.printStackTrace();
		}
		boolean matches = false;
		try {
			matches = querytool.matches(fragment.getMolecule());
		} catch (CDKException e) {
			e.printStackTrace();
		}
		if(!matches) {
			return;
		}
		
		// At this point, the monomer fragment has satisfied all conditions.
		fragment.setFragmentType(FragmentType.FATTY_ACID);
	}
	
	/**
	 * Identify as a linear polyketide
	 * @param currentFragment
	 */
	public void identifyAsLinearPK(Fragment currentFragment, int macrolideType) {
		// If there are is a lactone bond and an amino bond
		if(currentFragment.getLactoneHydroxylC() != null && currentFragment.getAminoCs().size() > 0) {
			//identifyAsLinearPK(currentFragment, currentFragment.getAminoC(), currentFragment.getLactoneHydroxylC());
			identifyAsLinearPK(currentFragment, null, currentFragment.getAminoCs().get(0), macrolideType);
			
			int pathLength = PathTools.getShortestPath(currentFragment.getMolecule(), currentFragment.getAminoCs().get(0), currentFragment.getLactoneHydroxylC()).size();
			if(currentFragment.getFragmentType() != null && pathLength == 3) {
				currentFragment.setFragmentType(FragmentType.FA_OR_PK);
			}
		}
		// If there is one connected amino acid
		else if(currentFragment.getAminoCs().size() > 0 && currentFragment.getAminoNs().size() > 0) {
			identifyAsLinearPK(currentFragment, null, null, macrolideType);
			if(currentFragment.getFragmentType() != null) {
				currentFragment.setFragmentType(FragmentType.FA_OR_PK);
			}
		}
		
		else if(currentFragment.getLactoneCarboxylC() != null) {
			identifyAsLinearPK(currentFragment, null, currentFragment.getLactoneCarboxylC(), macrolideType);
		}
		// If this does not appear to be an amino acid
		else {
			identifyAsLinearPK(currentFragment, null, null, macrolideType);
		}
		// Check if there is a saturated carbon chain indicative of a possible fatty acid
		if(currentFragment.getFragmentType() == FragmentType.POLYKETIDE) {
			SMARTSQueryTool querytoolSingle = null;
			SMARTSQueryTool querytoolDouble = null;
			try {
				querytoolSingle = new SMARTSQueryTool("[#6;A][#6;A;D2][#6;A;D2][#6;A;D2][#6;A;D2][#6;A]");
				querytoolDouble = new SMARTSQueryTool("[#6;A]=[#6;A;D2]\\[#6;A;D2]=[#6;A;D2]\\[#6;A;D2]=[#6;A]");
			} catch (CDKException e) {
				e.printStackTrace();
			}
			boolean matches = false;
			try {
				matches = querytoolSingle.matches(currentFragment.getMolecule());
				if(!matches) matches = querytoolDouble.matches(currentFragment.getMolecule());
			} catch (CDKException e) {
				e.printStackTrace();
			}
			if(matches) {
				currentFragment.setFragmentType(FragmentType.FA_OR_PK);
			}
		}
		
		// Check if this is possibily a modified amino acid - nitrogen containing with at most 20 atoms and two or fewer PK units
		if(currentFragment.getMolecule().getAtomCount() <= 20) {
			boolean hasNitrogen = false;
			for(int i = 0; i < currentFragment.getMolecule().getAtomCount(); i++) {
				if(currentFragment.getMolecule().getAtom(i).getAtomicNumber() == 7) {
					hasNitrogen = true;
				}
			}
			if(hasNitrogen && currentFragment.getLoadingUnits() != null && currentFragment.getLoadingUnits().size() <= 2) {
				currentFragment.setFragmentType(FragmentType.FA_OR_PK);
			}
		}
		// Check if a ketone has been converted to an anime for incorporation into a peptide bond
		//if(currentFragment.getAminoC() != null && currentFragment.getAminoN() != null) {
			//identifyAsAmineKetoneLinearPK(currentFragment);
		//}
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
		
		if(end != null && start != null && ChemicalUtilities.hasCarbonPath(currentFragment.getMolecule(), start, end)) {
			try {
				PolyketideModulePredictor pkPredictor =
						new PolyketideModulePredictor(
								currentFragment.getMolecule(),
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
				//System.out.println("Error in linear PK predictor - skipping");
			}
		}
		// If there is an end
		else if(end != null) {
			PolyketideModulePredictor pkPredictor = null;
			try {
				
				pkPredictor = new PolyketideModulePredictor(
						currentFragment.getMolecule(),
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
				//System .out.println("Error in linear PK predictor - skipping");
			}
			
			
		}
		else {	
			PolyketideModulePredictor pkPredictor = null;
			try {
			pkPredictor = new PolyketideModulePredictor(
					currentFragment.getMolecule(),
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
				e.printStackTrace();
				//System.out.println("Error in linear PK predictor - skipping");
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
	 * Check if this piece contains a hexane ring (and no other ring structures), consistent with the cyclohexane in an aminoglycoside
	 * @param currentFragment
	 */
	public void identifyAsAminoglycosideCyclohexane(Fragment currentFragment) {
		List<Set<IAtom>> smallRings = ChemicalUtilities.getSmallestRings(currentFragment.getMolecule());
		if(smallRings.size() != 1) {
			return;
		}
		// check if the one ring contains a non-carbon
		for(IAtom a : smallRings.get(0)) {
			if(a.getAtomicNumber() != 6) {
				return;
			}
			if(!a.getAtomTypeName().equals("C.sp3")) {
				return;
			}
		}
		currentFragment.setFragmentType(FragmentType.CYCLOHEXANE);
	}
	
	/**
	 * Take a monomer fragment and, if it is a ketoextended amino acid, return a two-valued array of monomerfragments
	 * corresponding to the amino acid and polyketide respectively
	 * @param currentFragment
	 * @return
	 */

	public Fragment[] getKetoextendedAAPieces(Fragment currentFragment, int macrolideType) {
		//if(currentFragment.getAminoNs().size() < 1 || macrolideType > 0 || PolyketideModulePredictor.getNumCarboxylicAcids(currentFragment.getMolecule()) > 1) return null;
		if(currentFragment.getAminoNs().size() < 1 || macrolideType > 0) return null;
		for(IAtom aminoN : currentFragment.getAminoNs()){
			boolean appropriateAminoN = true;
			for(IAtom atom : currentFragment.getMolecule().getConnectedAtomsList(aminoN)){
				if(atom.getAtomicNumber() != 6){
					appropriateAminoN = false;
					break;
				}
//				else if(currentFragment.getMolecule().getConnectedAtomsCount(atom) != 2){
//					appropriateAminoN = false;
//					break;
//				}
			}
			if(!appropriateAminoN) continue;
			
			// Create the array that will contain the C fragment and N fragment respectively
			Fragment[] pieces = new Fragment[2];
			// Check if this is a ketoextended amino acid
			IAtom adjacentCarbon = null;
			for(IAtom a : currentFragment.getMolecule().getConnectedAtomsList(aminoN)) {
				if(a.getAtomicNumber() == 6
						&& currentFragment.getMolecule().getConnectedAtomsCount(a) > 1) {
					adjacentCarbon = a;
				}
			}
			ArrayList<IBond> candidateBondsToBreak = new ArrayList<IBond>();
			for(IBond b : currentFragment.getMolecule().getConnectedBondsList(adjacentCarbon)) {
				if(b.contains(aminoN)) continue;
				candidateBondsToBreak.add(b);
			}
			for(IBond b : candidateBondsToBreak) {
				// Try breaking this bond, and see if this forms a linear polyketide
				
				currentFragment.getMolecule().removeBond(b);
				
				IAtom possibleEndCarbon = null;
				if(b.getAtom(0).equals(adjacentCarbon)) {
					possibleEndCarbon = b.getAtom(1);
				}
				else {
					possibleEndCarbon = b.getAtom(0);
				}
				
				IMoleculeSet partitions = ConnectivityChecker.partitionIntoMolecules(currentFragment.getMolecule());
				IMolecule possiblePK = null;
				if(partitions.getMoleculeCount() != 2) {
					currentFragment.getMolecule().addBond(b);
					continue;
				}
				
				if(partitions.getMolecule(0).contains(aminoN)) {
					possiblePK = partitions.getMolecule(1);
				}
				else {
					possiblePK = partitions.getMolecule(0);
				}
				
				IMolecule carbonToAdd = null;
				try {
					carbonToAdd = SmilesIO.readSmiles("C");} catch (Exception e1){}
				IAtom carbonAtom = carbonToAdd.getAtom(0);
				possiblePK.addAtom(carbonAtom);
				possiblePK.addBond(new Bond(carbonAtom, possibleEndCarbon, IBond.Order.SINGLE));
				possibleEndCarbon = carbonAtom;
				
				boolean isPK = false;
				// Do this as a check
				PolyketideModulePredictor pkPredictor = null;
				try {
					//Extend the "end" by 1 as this underwent decarboxcylation
					// Two cases: if there was something attached to the C terminus, or if this is a terminal piece
					pkPredictor = new PolyketideModulePredictor(
								possiblePK,
								0
								);	
					if(pkPredictor.isPK()) {
						isPK = true;
					}
				} catch (Exception e) {
					e.printStackTrace();
					isPK = false;
				}
				
				currentFragment.getMolecule().addBond(b);
				if(isPK == false) continue;
				//SmilesIO.drawMolecule(currentFragment.getMolecule(), "keto_before");
				
				//Try to find the 'carbon backbone start' for the potential PK piece
				IAtom pkCarboxylicCarbon = PolyketideModulePredictor.predictEndCarbon(currentFragment.getMolecule());
				for(IAtom carbon : currentFragment.getMolecule().atoms()){
					if(PolyketideModulePredictor.isCarboxylicCarbon(currentFragment.getMolecule(), carbon)){
						pkCarboxylicCarbon = carbon;
					}
				}
				
				//Trace the potential carbon backbone of the potential PK portion, this is used to count the number of carbons in the backbone
				//See if it a beta or alpha amino acids
				List<IAtom> path = null;
				try{
					path = PolyketideModulePredictor.getCarbonOnlyPath(currentFragment.getMolecule(), pkCarboxylicCarbon, aminoN);
				}catch(Exception e){
					continue;
				}
				if(path.size() < 5) return null;
				IAtom aminoCarboxylicCarbon = null;
				IAtom startCarbon = null;
				if(path.size() % 2 != 0){
					currentFragment.getMolecule().removeBond(path.get(path.size()-2), path.get(path.size()-3));
					aminoCarboxylicCarbon = path.get(path.size()-2);
					startCarbon = path.get(path.size()-3);
				}else{
					currentFragment.getMolecule().removeBond(path.get(path.size()-3), path.get(path.size()-4));
					aminoCarboxylicCarbon = path.get(path.size()-3);
					startCarbon = path.get(path.size()-4);
				}
				
				List<Fragment> fragmentPartitions = currentFragment.partitionIntoMonomerFragments();
				
				Fragment aminoAcidFrag = null;
				Fragment pkFrag = null;
				
				if(fragmentPartitions.get(0).getMolecule().contains(aminoN)) {
					aminoAcidFrag = fragmentPartitions.get(0);
					pkFrag = fragmentPartitions.get(1);
				}
				else {
					aminoAcidFrag = fragmentPartitions.get(1);
					pkFrag = fragmentPartitions.get(0);
				}
				
				//Add carboxylic acid to the amino acid piece
				Atom firstO = new Atom("O");
				firstO.setAtomTypeName("O.sp3");
				Atom secondO = new Atom("O");
				firstO.setAtomTypeName("O.sp2");
				Atom carboxylCarbon = new Atom("C");
				carboxylCarbon.setAtomTypeName("C.sp2");
				aminoAcidFrag.getMolecule().addAtom(firstO);
				aminoAcidFrag.getMolecule().addAtom(secondO);
				aminoAcidFrag.getMolecule().addAtom(carboxylCarbon);
				aminoAcidFrag.getMolecule().addBond(new Bond(carboxylCarbon, firstO, IBond.Order.SINGLE));
				aminoAcidFrag.getMolecule().addBond(new Bond(carboxylCarbon, secondO, IBond.Order.DOUBLE));
				aminoAcidFrag.getMolecule().addBond(new Bond(carboxylCarbon, aminoCarboxylicCarbon, IBond.Order.SINGLE));
				//PK piece does not have any atoms added to it
				
				//Add all aminoNs to the new fragments
				for(IAtom aminoNtoAdd : currentFragment.getAminoNs()){
					if(pkFrag.getMolecule().contains(aminoNtoAdd)){
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
				identifyAsLinearPK(pkFrag, 0);
				
				//SmilesIO.drawMolecule(aminoAcidFrag.getMolecule(), "cfrag");
				//SmilesIO.drawMolecule(pkFrag.getMolecule(), "nfrag");
				
				pieces[0] = aminoAcidFrag;
				pieces[1] = pkFrag;
				return pieces;
			}
		}
		return null;
	}	
	
	/**
	 * Read Amino Acid input. 
	 */
	private void readAminoAcidInput(String aminoAcidPath) {
		// Generate array of amino acids and NRP IMolecules from the file.
		
		// Read amino acid file
		Scanner aaScanner = null;
		try {	
			aaScanner = new Scanner(new InputStreamReader(getClass().getResourceAsStream(aminoAcidPath)));
		} catch (Exception e) {
			try {	
				aaScanner = new Scanner(new File(aminoAcidPath));
			} catch (Exception ex) {
				ex.printStackTrace();
			}
		}
		while(aaScanner.hasNextLine()) {
			String line = aaScanner.nextLine();
			String[] vals = line.split("\t");
			if(vals[0].startsWith("#")) {
				continue;
			}
			if(vals.length < 4) {
				System.err.println("Found line with fewer than four values: " + line);
				continue;
			}
			try {
				aminoAcids.add(SmilesIO.readSmiles(vals[2]));
				aminoAcidEnums.add(DomainEnums.getAminoAcidEnumFromAbbreviation(vals[3]));
				ArrayList<TailoringDomainEnums> tailoringDomains = new ArrayList<TailoringDomainEnums>();
				for(int i = 4; i < vals.length; i++) {
					// These are tailoring domains
					TailoringDomainEnums tailoringDomain = DomainEnums.getTailoringDomainFromAbbreviation(vals[i]);
					if(tailoringDomain == null) {
						System.out.println("Warning: In amino acid file, found unknown tailoring domain " + vals[i] + " - skipping");
					}
					else {
						tailoringDomains.add(tailoringDomain);
					}
				}
				if(tailoringDomains.size() > 0) {
					aminoAcidTailoringEnums.add(tailoringDomains);
				}
				else {
					aminoAcidTailoringEnums.add(null);
				}
				
			} catch (Exception e) {
				System.err.println("Found bad line: " + line);
				e.printStackTrace();
			}	
		}
	}
}
