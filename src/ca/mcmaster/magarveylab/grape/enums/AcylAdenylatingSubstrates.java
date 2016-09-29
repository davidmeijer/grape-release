package ca.mcmaster.magarveylab.grape.enums;

import java.util.LinkedHashMap;
import java.util.Map;

import org.openscience.cdk.interfaces.IAtomContainer;


/**
 * Substrates of a distinct clade of acyl-adenylating enzymes revealed by phylogenetic analysis. 
 * Starter units for NRP and PK
 * @author dejong skinnider
 *
 */
public enum AcylAdenylatingSubstrates {
	
	// fatty acids
//	MYRISTATE("myristate.hmm", "Myristate", "C<sub>14</sub>", "CCCCCCCCCCCCCC(I)=O"), //matchs nonspecific fatty acid
//	LONG_CHAIN_FATTY_ACID("long_chain_fatty_acid.hmm", "Long-chain fatty acid", "C<sub>10</sub>", "O=C(I)CCCCCCCCC"), //matchs nonspecific fatty acid
//	SHORT_CHAIN_FATTY_ACID("short_chain_fatty_acid.hmm", "Short-chain fatty acid", "C<sub>6</sub>", "O=C(I)CCCCC"), //matchs nonspecific fatty acid
	_3_AMINONON_5_ENOIC_ACID("3-aminonon-5-enoic_acid.hmm", "3-aminonon-5-enoic acid", "C<sub>9</sub>", new String[]{"CCC/C=C/CC(N)CC(I)=O"}, true),

	// aromatic starters
	_2_3_DIHYDROXYBENZOIC_ACID("2_3_dihydroxybenzoic_acid.hmm", "2,3-dihydroxybenzoic acid", 
	"2,3DHB", new String[]{"O=C(C1=C(O)C(O)=CC=C1)I", "O=C(I)C1=CC=CC(O)=C1O"}, true),
	CINNAMIC_ACID("cinnamic_acid.hmm", "Cinnamic acid", "Cinn", new String[]{"IC(/C=C/C1=CC=CC=C1)=O", "IC(/C=C/C1=CC=CC=C1)=O"}, true),
	_3_AMINO_5_HYDROXYBENZOIC_ACID("AHBA.hmm", "3-amino-5-hydroxybenzoic acid", "AHBA", new String[]{"O=C(I)C1=CC(N)=CC(O)=C1", "O=C(C1=CC(O)=CC(N)=C1)I", "COC1=CC(C)=CC(N)=C1"}, true),
	_3_FORMAMIDO_5_HYDROXYBENZOIC_ACID("FHBA.hmm", "3-formamido-5-hydroxybenzoic acid", "FHBA", new String[]{"IC(C1=C(O)C(NC=O)=CC=C1)=O", "IC(C1=CC=CC(NC=O)=C1O)=O", "OC(C1=CC=CC(NC=O)=C1O)=O", "OC(=O)C1=C(O)C(NC=O)=CC=C1", "Oc1c(NC=O)cccc1C=O"}, true),
	_3_HYDROXYPICOLINIC_ACID("3-HPA.hmm", "3-hydroxypicolinic acid", "3-HPA", new String[]{"IC(C1=C(O)C=CC=N1)=O", "IC(C1=NC=CC=C1O)=O"}, true),
	_3_HYDROXYQUINALDIC_ACID("3-HQA.hmm", "3-hydroxyquinaldic acid", "3-HQA", new String[]{"IC(C1=NC2=CC=CC=C2C=C1O)=O", "IC(C1=NC(C=CC=C2)=C2C=C1O)=O", "IC(C(N=C(C=CC=C1)C1=C2)=C2O)=O"}, true),
	QUINOXALINE("quinoxaline.hmm", "Quinoxaline-2-carboxylic acid", "QX", new String[]{"IC(C1=NC2=CC=CC=C2N=C1)=O", "IC(C1=NC(C=CC=C2)=C2N=C1)=O", "IC(C1=CN=C2C(C=CC=C2)=N1)=O"}, true),
	PHENYLACETATE("phenylacetate.hmm", "Phenylacetate", "PAA", new String[]{"O=C(I)CC1=CC=CC=C1", "O=C(CC1=CC=CC=C1)I"}, true),
	PARA_HYDROXYBENZOIC_ACID("PHBA.hmm", "Para-hydroxybenzoic acid", "PHBA", new String[]{"O=C(C1=CC=C(C=C1)O)I", "O=C(I)C(C=C1)=CC=C1O"}, true),
	PARA_AMINOBENZOIC_ACID("PABA.hmm", "Para-aminobenzoic acid", "PABA", new String[]{"O=C(C1=CC=C(C=C1)N)I", "NC1=CC=C(C=C1)C(I)=O"}, true),
	
	// alpha keto/alpha hydroxy acids
	ALPHA_KETOISOCAPROATE("alpha-ketoisocaproate.hmm", "&alpha;-ketoisocaproate", "&alpha;kL", new String[]{"IC(C(CC(C)C)OF)=O"}, false),
	ALPHA_KETOISOVALERATE("alpha-ketoisovalerate.hmm", "&alpha;-ketoisovalerate", "&alpha;kV", new String[]{"IC(C(C(C)C)O)=O"}, false),
	PYRUVATE("pyruvate.hmm", "Pyruvate", "Pyr", new String[]{"IC(C(C)O)=O", "CC(=O)C(O)=O"}, false),
	_3_METHYL_2_OXOPENTANOIC_ACID("3-methyl-2-oxopentanoate.hmm", "3-methyl-2-oxopentanoate", "&alpha;kI", new String[]{"IC(C(OF)C(CC)C)=O"}, false),
	PHENYLPYRUVATE("phenylpyruvate.hmm", "Phenylpyruvate", "&alpha;kF", new String[]{"IC(C(OF)CC1=CC=CC=C1)=O", "FOC(CC1=CC=CC=C1)C(I)=O"}, true),
	_2_HYDROXY_2_METHYL_4_OXOPENTANOIC_ACID("", "2-hydroxy-2-methyl-4-oxopentanoic acid", "Hmop", new String[]{"OC(CC(C)=O)(C)C(O)=O"}, false),
	_2_OXOBUTANOIC_ACID("","2-oxobutanoic acid", "OxobA", new String[]{"O=C(CC)C(O)=O", "CCC(C(O)=O)=O"}, false),
	_2_HYDROXYBUTYRIC_ACID("", "2-hydroxybutyric acid", "3-HBA", new String[]{"CCC(O)C(=O)O"}, true),
	_2_OXOBUTARATE("", "2-oxobutarate", "2-OXO", new String[]{"CCC(=O)C([O-])=O"}, true),
	
	// small starters
	BETA_AMINOALANINAMIDE("beta_aminoalaninamide.hmm", "&beta;-aminoalaninamide", "&beta;-Aln", new String[]{"ICC(C(N)=O)N"}, true),
	DIHYDROXYCYCLOHEXANECARBOXYLIC_ACID("DHCHC.hmm", "Dihydroxycyclohexanecarboxylic acid", 
			"DHCHC", new String[]{"OC1C(CCC(C1)C(I)=O)O", "COC1CC(CI)CCC1O"}, true),
	CYCLOHEXANECARBOXYLATE("CHC.hmm", "Cyclohexanecarboxylic acid", "CHC", new String[]{"IC(C1CCCCC1)=O"}, true),
	ALKYLMALONYL_COA("alkylmalonyl_CoA.hmm", "Alkylmalonyl-CoA", "CoL", new String[]{"IC(CF)=O"}, false),
	AMINOLEVULINIC_ACID("aminolevulinic_acid.hmm", "5-aminolevulinic acid", "5-ALA", new String[]{"IC1=C(O)CCC1=O"}, true),
	_3_HYDROXYBUTANOIC_ACID("3-hydroxybutanoic_acid.hmm", "3-hydroxybutanoic acid", "OHBu", new String[]{"IC(CC(O)C)=O"}, false),
	SALICYLIC_ACID("salicylic_acid.hmm", "Salicylic acid", "Sal", new String[]{"IC(C1=C(O)C=CC=C1)=O", "OC1=CC=CC=C1C(I)=O"}, true),
	//Different structures used for GRAPE
	LACTATE("", "lactate", "lac", new String[]{"CC(C(O)=O)O", "C[C@@H](C(=O)O)O", "C1=CC=C(C=C1)CC(C(=O)O)O", "C1=CC(=CC=C1CC(C(=O)O)O)O", "O=C(O)C(C)O", "CC(O)C(=O)C(O)=O"}, false),
	VALERAIC_ACID("", "Valeraic acid", "VLA", new String[]{"CCCCC(O)=O", "CC(C)(O)C(O)C(O)=O", "C(C(=O)O)(CCCC1=CC=C(C=C1)O)N", "COC1=CC=C(C=C1)CCCC(C(=O)O)N", "C(CC(C(=O)O)O)(C)C", "CCC(C)CC(=O)O", "C1=CC=C(C=C1)CCCC(C(=O)O)N", "O=C(O)C(O)C(C)C"}, false),
	_3_HYDROXYLPATANOLIC_ACID("", "3-hydroxylpatanolic acid", "3-HLPA", new String[]{"OC(CC(O)CC)=O"}, false),
	BENZOIC_ACID("", "Benzoic acid", "BZA", new String[]{"OC1=CC=CC=C1C(O)=O", "OC(=O)C1=C(O)C=CC=C1"}, false),
	HYDROXY_3_METHYL_PENTANOIC_ACID("", "Hydroxy 3-methyl pentanoic acid", "HMP", new String[]{"OC(C(O)C(CC)C)=O", "C(C(C(CC)C)O[S](=O)(=O)O)(=O)O"}, false),
	PRENYLATION("", "Prenylation", "PRN", new String[]{"CC(C)(C)C=C", "C/C(C)=C\\CO"}, false),
	;
	
	private final String hmm;
	private final String name;
	private final String abbreviation;
	private final String smiles[];
	private final boolean uniqueSubstructure;
	
	private AcylAdenylatingSubstrates(final String hmm, final String name, final String abbreviation, 
			final String[] smiles, boolean canCheckForSubstructure) {
		this.hmm = hmm;
		this.smiles = smiles;
		this.name = name;
		this.abbreviation = abbreviation;
		this.uniqueSubstructure = canCheckForSubstructure;
	}

	public String hmm() { 
		return hmm; 
	}
	
	public String fullName() { 
		return name; 
	}
	
	public String abbreviation() { 
		return abbreviation; 
	}
	
	public String[] smiles() { 
		return smiles; 
	}
	
	public boolean canCheckForSubstructure(){
		return uniqueSubstructure;
	}
	
	public String[] realSmiles() {
		String[] fixedSmiles = new String[smiles.length];
		for(int i = 0; smiles.length > i; i++){
			String fixedSmile = smiles[i].replace("I", "");
			fixedSmiles[i] = fixedSmile.replace("F", "O"); 
		}
		return fixedSmiles;
	}
	
	public IAtomContainer[] mols(){
		IAtomContainer[] mols = new IAtomContainer[smiles.length];
		for(int i = 0; smiles.length > i; i++){
			mols[i] = MoleculeClasses.getMol(realSmiles()[i]);
		}
		return mols;
	}
	
	public static Map<String, IAtomContainer[]> getAll(){
		Map<String, IAtomContainer[]> all = new LinkedHashMap<String, IAtomContainer[]>(); 
		for(AcylAdenylatingSubstrates single : AcylAdenylatingSubstrates.values()){
			all.put(single.toString(), single.mols());
		}
		return all;
	}
	
	public static Map<String, IAtomContainer[]> getSubstructureSearchable(){
		Map<String, IAtomContainer[]> all = new LinkedHashMap<String, IAtomContainer[]>(); 
		for(AcylAdenylatingSubstrates single : AcylAdenylatingSubstrates.values()){
			if(single.canCheckForSubstructure()) all.put(single.toString(), single.mols());
		}
		return all;
	}

}

