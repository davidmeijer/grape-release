package ca.mcmaster.magarveylab.grape.nrp.chem.matcher;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;

import ca.mcmaster.magarveylab.grape.enums.AcylAdenylatingSubstrates;
import ca.mcmaster.magarveylab.grape.enums.FattyAcidsEnums;
import ca.mcmaster.magarveylab.grape.enums.KnownOtherEnums;
import ca.mcmaster.magarveylab.grape.enums.SugarEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.AminoAcidEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.MultipleAminoAcidEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.PolyKetideDomainEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.SmallPKunits;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.SugarModificationsEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.TailoringDomainEnums;
import ca.mcmaster.magarveylab.grape.nrp.chem.Fragment;
import ca.mcmaster.magarveylab.grape.nrp.chem.Fragment.FragmentType;
import ca.mcmaster.magarveylab.grape.pk.modules.PKsubstrate;
import ca.mcmaster.magarveylab.grape.util.ChemicalUtilities;
import ca.mcmaster.magarveylab.grape.util.io.ReadAminoAcidFile;
import ca.mcmaster.magarveylab.grape.util.io.SmilesIO;

public class SubstrateMatcher {
	
	IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
	UniversalIsomorphismTester uit = new UniversalIsomorphismTester();
	private List<AminoAcidEnums[]> mutlipleAminoAcidConstituents = MultipleAminoAcidEnums.allAminoAcids();
	private IAtomContainer[] fattyAcids = FattyAcidsEnums.mols();
	private List<IAtomContainer> aminoAcids;
	private List<String>	aminoAcidSpecificName;
	private List<AminoAcidEnums> aminoAcidEnums;
	private List<List<TailoringDomainEnums>> aminoAcidTailoringEnums;
	private Map<String, IAtomContainer> sugarsMap = SugarEnums.getAll();
	private List<IAtomContainer> multipleAminoAcidMols = MultipleAminoAcidEnums.allMols();
	
	public SubstrateMatcher(String aminoAcidPath){
		aminoAcids = new ArrayList<IAtomContainer>();
		aminoAcidEnums = new ArrayList<AminoAcidEnums>();
		aminoAcidSpecificName = new ArrayList<String>();
		aminoAcidTailoringEnums = new ArrayList<List<TailoringDomainEnums>>();
		
		ReadAminoAcidFile.read(aminoAcidPath, this);
	}
	
	public void addAminoAcid(IAtomContainer mol, String name, AminoAcidEnums aminoAcidEnum, List<TailoringDomainEnums> tailorings) {
		aminoAcids.add(mol);
		aminoAcidSpecificName.add(name);
		aminoAcidEnums.add(aminoAcidEnum);
		aminoAcidTailoringEnums.add(tailorings);
	}
	
	

	public void identifyAsMultipleAminoAcids(Fragment fragment) {
		ArrayList<Double> tanimotoScores = ChemicalUtilities.getTanimotoScores(fragment.getAtomContainer(), multipleAminoAcidMols);
		int highestScoreIndex = ChemicalUtilities.getHighestTanimotoScoreIndex(tanimotoScores);
		if(tanimotoScores.get(highestScoreIndex) > 0.90) {
			fragment.addAminoAcidDomains(mutlipleAminoAcidConstituents.get(highestScoreIndex));
			fragment.setFragmentType(FragmentType.MULTIPLE_AMINO_ACID_PIECE);
			fragment.setTanimotoScore(tanimotoScores.get(highestScoreIndex));
		}
	}
	
	public void identifyAsPerfectAminoAcid(Fragment fragment) {
		ArrayList<Double> aminoAcidTanimotoScores = ChemicalUtilities.getTanimotoScores(fragment.getAtomContainer(), aminoAcids);
		int highestScoreIndex = ChemicalUtilities.getHighestTanimotoScoreIndex(aminoAcidTanimotoScores);
		if(aminoAcidTanimotoScores.get(highestScoreIndex) == 1.0) {
			fragment.addAminoAcidDomain(aminoAcidEnums.get(highestScoreIndex));
			fragment.addSpecificName(aminoAcidSpecificName.get(highestScoreIndex));
			fragment.setFragmentType(FragmentType.AMINO_ACID);
			fragment.setTanimotoScore(aminoAcidTanimotoScores.get(highestScoreIndex));
		}
	}
	
	public void identifyAsAminoAcid(Fragment fragment) {
		ArrayList<Double> aminoAcidTanimotoScores = ChemicalUtilities.getTanimotoScores(fragment.getAtomContainer(), aminoAcids);
		int highestScoreIndex = ChemicalUtilities.getHighestTanimotoScoreIndex(aminoAcidTanimotoScores);
		if(aminoAcidTanimotoScores.get(highestScoreIndex) >=  0.9) {
			fragment.addAminoAcidDomain(aminoAcidEnums.get(highestScoreIndex));
			fragment.setFragmentType(FragmentType.AMINO_ACID);
			fragment.addSpecificName(aminoAcidSpecificName.get(highestScoreIndex));
			fragment.setTanimotoScore(aminoAcidTanimotoScores.get(highestScoreIndex));
		}
		// If this piece contains a thiazole, check if it contains a cysteine as a substructure
		if(fragment.getTailoringDomains().contains(TailoringDomainEnums.THIAZOLE) ||
				fragment.getTailoringDomains().contains(TailoringDomainEnums.SULFUR_BETA_LACTAM)) {
			try {
				IAtomContainer cysteineSubstructure = SmilesIO.readSmilesTemplates("C(C(C(=O)O)N)S");
				if(uit.isSubgraph(fragment.getAtomContainer(), cysteineSubstructure)) {
					fragment.setFragmentType(FragmentType.AMINO_ACID);
					fragment.addAminoAcidDomain(AminoAcidEnums.Cysteine);
					fragment.addSpecificName("Cys-t");
					double tanimoto = ChemicalUtilities.getTanimotoScore(fragment.getAtomContainer(), cysteineSubstructure);
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
	
	public boolean isAminoAcid(IAtomContainer molecule) {
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
		
		for(Entry<String, IAtomContainer> sugar : sugarsMap.entrySet()){
			double score = ChemicalUtilities.getTanimotoScore(fragment.getAtomContainer(), sugar.getValue());
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
			for(int i = 0; i < fragment.getAtomContainer().getAtomCount(); i++) {
				if(fragment.getAtomContainer().getAtom(i).getAtomicNumber() == 7) {
					fragment.addSugarModification(SugarModificationsEnums.AMINOTRANSFER);
				}
				// N Methylations
				SMARTSQueryTool querytool = null;
				querytool = new SMARTSQueryTool("N[C;D1]", builder);
				boolean hasNMethyl = false;
				try {
					hasNMethyl = querytool.matches(fragment.getAtomContainer());
				} catch (CDKException e) {
					e.printStackTrace();
				}
				if(hasNMethyl) {
					// Check
					fragment.addSugarModification(SugarModificationsEnums.N_METHYLATION);
				}
				
				// O Methylations
				querytool = null;
				querytool = new SMARTSQueryTool("O[C;D1]", builder);
				boolean hasOMethyl = false;
				try {
					hasOMethyl = querytool.matches(fragment.getAtomContainer());
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
			double highestFattyAcidScore = ChemicalUtilities.getHighestTanimotoScore(fragment.getAtomContainer(), Arrays.asList(fattyAcids));
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
		for(int i = 0; i < fragment.getAtomContainer().getAtomCount(); i++) {
			int atomicNumber = fragment.getAtomContainer().getAtom(i).getAtomicNumber();
			if(atomicNumber != 1 && atomicNumber != 6 && atomicNumber != 8) {
				return;
			}
		}
		IAtomContainer fattyAcid3OHTemplate = null;
		IAtomContainer fattyAcidTemplate = null;
		try {
			fattyAcidTemplate = SmilesIO.readSmilesTemplates("CC(O)=O");
			fattyAcid3OHTemplate = SmilesIO.readSmilesTemplates("CC(O)CC(O)=O");
		} catch (IOException | CDKException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		boolean matchesNormalFattyAcid = false;
		boolean matches3OHFattyAcid = false;
		try {
			
			matchesNormalFattyAcid = uit.isSubgraph(fragment.getAtomContainer(), fattyAcidTemplate);
			matches3OHFattyAcid = uit.isSubgraph(fragment.getAtomContainer(), fattyAcid3OHTemplate);
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
		for(int i = 0; i < fragment.getAtomContainer().getAtomCount(); i++) {
			int atomicNumber = fragment.getAtomContainer().getAtom(i).getAtomicNumber();
			if(atomicNumber == 8) {
				numOxygens++;
			}
			if(atomicNumber == 6) {
				numCarbons++;
			}
		}
		
		// Check that there are no carbon rings
		ArrayList<IAtom> atomsInCarbonRings = ChemicalUtilities.getAtomsInCarbonRings(fragment.getAtomContainer());
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
		querytool = new SMARTSQueryTool("C[C;D2][C;D2][C;D2][C;D2]C", builder);
		boolean matches = false;
		try {
			matches = querytool.matches(fragment.getAtomContainer());
		} catch (CDKException e) {
			e.printStackTrace();
		}
		if(!matches) {
			return;
		}
		
		// At this point, the monomer fragment has satisfied all conditions.
		fragment.setFragmentType(FragmentType.FATTY_ACID);
	}
	
	public void identifySubstrctureAsStarter(Fragment fragment) {
		IAtomContainer mol = null;
		try {
			mol = SmilesIO.readSmilesTemplates(SmilesIO.generateSmiles(fragment.getAtomContainer()));
		} catch (IOException | CDKException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		for(AcylAdenylatingSubstrates starter  : AcylAdenylatingSubstrates.values()){
			if(!starter.canCheckForSubstructure()) continue;
			try {
				IAtomContainer[] mols = starter.mols();
				for(int i = 0; mols.length > i; i++){
					if(uit.isSubgraph(mol, mols[i])){
						fragment.addStarter(starter);
						return;
					}
				}
			} catch (CDKException e){}
		}
	}
	

	public void identifyAsKnownOther(Fragment fragment) {
		for(KnownOtherEnums knownOther : KnownOtherEnums.values()){
			if(knownOther.getMol() != null){
				double score = ChemicalUtilities.getTanimotoScore(fragment.getAtomContainer(), knownOther.getMol());
				if(score > 0.9){
					fragment.setFragmentType(FragmentType.KNOWN_OTHER);
					fragment.setKnownOther(knownOther);
					fragment.setTanimotoScore(score);
					return;
				}
			}
		}
		
	}

	public void identifyAsAcylAdenlyatingSubstrate(Fragment fragment) {
		for(AcylAdenylatingSubstrates acyl  : AcylAdenylatingSubstrates.values()){
			IAtomContainer[] mols = acyl.mols();
			for(int i = 0; mols.length > i; i++){
				if(ChemicalUtilities.getTanimotoScore(fragment.getAtomContainer(), mols[i]) > 0.9){
					fragment.addStarter(acyl);
					fragment.setTanimotoScore(ChemicalUtilities.getTanimotoScore(fragment.getAtomContainer(), mols[i]));
					fragment.setFragmentType(FragmentType.ACYL_ADENYLATING);
					return;
				}
			}
		}		
	}

	public void identifyAsSmallPolyketide(Fragment currentFragment) {
		
		for(SmallPKunits smallPK  : SmallPKunits.values()){
			if(ChemicalUtilities.getTanimotoScore(currentFragment.getAtomContainer(), smallPK.mol()) > 0.95){
				currentFragment.setFragmentType(FragmentType.POLYKETIDE);
				currentFragment.setTanimotoScore(ChemicalUtilities.getTanimotoScore(currentFragment.getAtomContainer(), smallPK.mol()));
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
	
	
	public boolean identifyAsSmallPolyketide(IAtomContainer mol) {
		for(SmallPKunits smallPK  : SmallPKunits.values()){
			if(ChemicalUtilities.getTanimotoScore(mol, smallPK.mol()) > 0.95){
				return true;
			}
		}
		return false;
	}
	
	/**
	 * Check if this piece contains a hexane ring (and no other ring structures), consistent with the cyclohexane in an aminoglycoside
	 * @param currentFragment
	 */
	public void identifyAsAminoglycosideCyclohexane(Fragment currentFragment) {
		List<Set<IAtom>> smallRings = ChemicalUtilities.getSmallestRings(currentFragment.getAtomContainer());
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
	
}
