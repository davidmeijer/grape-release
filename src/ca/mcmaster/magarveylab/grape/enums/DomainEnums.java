package ca.mcmaster.magarveylab.grape.enums;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.grape.pk.modules.PKsubstrate;
import ca.mcmaster.magarveylab.grape.util.io.SmilesIO;

/**
 * Contains all of the shared GRAPE enums. This file has been modified from the PRISM enum file.
 * @author gmchen
 */
public class DomainEnums {
	
		
	public enum SugarModificationsEnums {
		AMINOTRANSFER("AMINO"),
		O_METHYLATION("OMT"),
		N_METHYLATION("NMT"),
		ACETYLTRANSFER("ACETYL");
		private final String abbreviation;
		private SugarModificationsEnums(final String abbreviation) {
			this.abbreviation = abbreviation;
		}
		
		/**
		 * Get the abbreviation associated with this domain.
		 * @return	the abbreviation 
		 */
		public String getAbbreviation() {
			return(this.abbreviation);
		}
		}
	
	public enum PolyKetideDomainEnums {
		KETONE("KS-AT-T"),
		HYDROXYL("KS-AT-KR-T"),
		DOUBLEBOND("KS-AT-KR-DH-T"),
		SINGLEBOND("KS-AT-KR-DH-ER-T"),
		END("END"),
		BAD("not a proper module");
		
		private final String full; 
		
		private PolyKetideDomainEnums(String full){
			this.full = full;
		}
		
		public String getFull(){
			return full;
		}
	}
	
	/**
	 * Amino acid enums
	 * @author gmchen
	 */
	public enum AminoAcidEnums {
		Beta_Alanine("bAla"),
		Alanine("Ala"),
		BetaMethylPhenylalanine("MePhe"),
		Phenylalanine("Phe"),
		beta_Phenylalanine("bPhe"),
		Beta_HyxroxyPhenylalanine("OHPhe"),
		HydroxyLeucine("OHLeu"),
		Leucine("Leu"),
		Isoleucine("Ile"),
		Glycine("Gly"),
		Hydroxyphenylglycine("Hpg"),
		DiHydroxyphenylglycine("Dhpg"),
		Salicylic_Acid_or_2_3_Dihydroxy_Benzoic_Acid("Sal"),
		Hydroxyasparagine("OHAsn"),
		Asparagine("Asn"),
		Aspartic_Acid("Asp"),
		Hydroxyaspartic_Acid("OHAsp"),
		Beta_Methyl_Aspartic_Acid("MeAsp"),
		Benzoic_Acid("Bz"),
		Capreomycidine("Cap"),
		Citrulline("Cit"),
		Cysteine("Cys"),
		Glutamine("Gln"),
		Glutamic_Acid("Glu"),
		Histidine("His"),
		Methionine("Met"),
		Norvaline("Nva"),
		Isovaline("Iva"),
		Valine("Val"),
		Beta_hydroxy_Valine("OHVal"),
		HydroxyOrnithine("OHOrn"),
		Ornithine("Orn"),
		Tryptophan("Trp"),
		Threonine("Thr"),
		Butenyl_Methyl_Threonine("Bmt"),
		Arginine("Arg"),
		MethylProline("MePro"),
		Proline("Pro"),
		Serine("Ser"),
		Tyrosine("Tyr"),
		HydroxyTyrosine("OHTyr"),
		Lysine("Lys"),
		Beta_Lysine("bLys"),
		Adipic_Acid("Aad"),
		Lactate("Lac"),
		Lactic_Acid("Lac"),
		Propionyl("Hap"),
		Hydroxyisovalerate("Hiv"),
		Kynurenine("Kyn"),
		Aminobutyric_Acid("Abu"),
		Dehydro_Aminobutyric_Acid("Dhab"),
		Aminoisobutyric("Aib"),
		Coronamic("Cma"),
		Dolaproine("Dap"),
		Enduracididine("End"),
		Hydroxy_3_methylpentanoic_Acid("Hmp"),
		Diaminopropionic_Acid("Dpr"),
		Valeric_Acid("Vaa"),
		Diaminobutyric_Acid("Dab"),
		Pipecolic_Acid("Pip"),
		MethylGlutamate("MeGlu"),
		Pyruvate("Pyr"),
		Alpha_keto_isocaproate("akIL"),
		Alpha_keto_isovalerate("akIV"),
		Oxopentoate("Oxo"),
		Epoxi_OXODECANOIC_ACID("Aeo"),
		OH_quinaldic_acid("OHQA"),
		piperazic_acid("Piz"),
		quinoxaline_carboxylic_acid("QxCA"),
		_3_hydroxy_anthranilic_acid("3HA"),
		//fungal
		_4_hydroxyphenylpyruvate("OHPP"),
		DehydroalanOne("Dha"),
		HydroxyglutamOne("OHGln"),
		_2_hydroxy_3_methylpentanoOc_acOd("HMP"),
		Homo_serOne("Hse"),
		Homo_tyrosOne("Hty"),
		HydroxyhomotyrosOne_sulfate("OHHty"),
		Ondole_pyruvOc_acOd("OPA"),
		N_methoxytryptophan("OMeTrp"),		
		//ambiguous set
		Serine_or_Cysteine("SerCysA"),
		;
		
		private final String abbreviation;
		private AminoAcidEnums(final String abbreviation) {
			this.abbreviation = abbreviation;
		}
		/**
		 * Get the abbreviation associated with this domain.
		 * @return	the abbreviation 
		 */
		public String getAbbreviation() {
			return(this.abbreviation);
		}
	}
	public static AminoAcidEnums getAminoAcidEnumFromAbbreviation(String abbreviation) {
		for(AminoAcidEnums currentEnum : AminoAcidEnums.values()) {
			if(currentEnum.abbreviation.equals(abbreviation)) {
				return(currentEnum);
			}
		}
		System.out.println("Warning: amino acid domain abbreviation " + abbreviation + " does not exist - skipping.");
		return(null);
	}
	
	public enum MultipleAminoAcidEnums { //NCCCC(C(N)C(O)O)C1=CC=CC2=C1NC=C2C[C@H](N)C(O)=O
		SER_ASP_ASN("Cc1c(N)nc(nc1C(O)=O)C(CC(N)=O)NCC(N)C(N)=O", new AminoAcidEnums[]{AminoAcidEnums.Serine, AminoAcidEnums.Asparagine, AminoAcidEnums.Asparagine}),
		TRP_GLU("[N][C]([C](O[C]c1c(C([O])=O)n([O])c2[c][c][c]c([C][O])c12)[C]([O])C([O])=O)C([O])=O", new AminoAcidEnums[]{AminoAcidEnums.Tryptophan, AminoAcidEnums.Glutamic_Acid}),
		TRP_LYS("NCCCC(C(N)C(O)O)C1=CC=CC2=C1NC=C2C[C@H](N)C(O)=O", new AminoAcidEnums[]{AminoAcidEnums.Tryptophan, AminoAcidEnums.Lysine}),
		;
		
		private final String smiles;
		private final AminoAcidEnums[] aminoAcids;
		
		private MultipleAminoAcidEnums(final String smiles, final AminoAcidEnums[] aminoAcids){
			this.smiles = smiles;
			this.aminoAcids = aminoAcids;
		}
		public String smiles(){
			return smiles;
		}
		public AminoAcidEnums[] aminoAcids(){
			return aminoAcids;
		}
		public static Map<IAtomContainer, AminoAcidEnums[]> getAll() {
			Map<IAtomContainer, AminoAcidEnums[]> allMultiMap = new HashMap<IAtomContainer, AminoAcidEnums[]>();
			for(MultipleAminoAcidEnums single : MultipleAminoAcidEnums.values()){
				IAtomContainer mol = null;
				try{
					mol = SmilesIO.readSmilesTemplates(single.smiles());
				}catch(Exception e){}
				allMultiMap.put(mol, single.aminoAcids);
			}
			
			return allMultiMap;
		}
		public static ArrayList<IAtomContainer> allMols() {
			ArrayList<IAtomContainer> mols = new ArrayList<IAtomContainer>();
			for(MultipleAminoAcidEnums single : MultipleAminoAcidEnums.values()){
				IAtomContainer mol = null;
				try{mol = SmilesIO.readSmilesTemplates(single.smiles());}catch(Exception e){}
				mols.add(mol);
			}
			return mols;
		}
		public static ArrayList<AminoAcidEnums[]> allAminoAcids() {
			ArrayList<AminoAcidEnums[]> aminoAcids = new ArrayList<AminoAcidEnums[]>();
			for(MultipleAminoAcidEnums single : MultipleAminoAcidEnums.values()){
				aminoAcids.add(single.aminoAcids);
			}
			return aminoAcids;
		}
	}
	/**
	 * Enum containing a list of all domains used in a Prism genome search.
	 */
	public enum TailoringDomainEnums {
		// Tailoring enzymes
		GLYCOSYLTRANSFERASE("GTr"),
		HALOGENATION("Cl"),
		SULFOTRANSFERASE("ST"),
		TRYPTOPHAN_DIOXYGENASE("Kyn"),
		PROLINE_DEHYDROGENASE("DHO"),
		ISOPENICILLIN_N_SYNTHASE("IPNS"),
		ISOPENICILLIN_N_ACYLTRANSFERASE("IAT"),
		DEACETOXYCEPHALOSPORIN_C_SYNTHASE("DOACS"),
		P450("P450"),
		O_METHYLTRANSFERASE("OMT"),
		N_METHYLTRANSFERASE("NMT"),
		C_METHYLTRANSFERASE("CMT"),
		THIAZOLE("THZ"),
		OXAZOLE("OXZ"),
		SULFUR_BETA_LACTAM("SULFUR_BETA_LACTAM"),
		PRENYLATION("PRN"), 
		C_HYDROXYLATION("CH"), 
		LANTITHIONE_LINKAGE("LL"),
		PYRIDINE_CLYLYZATION("PC"),
		DURAMYCIN_LIKE_LINKAGE("DLL"),
		BOTTROMYCIN_LIKE_LINKAGE("BLL"), 
		THIOL_AMIDE("TA"),
		DECARBOXYLATION("DCA"),
		DEHYDRATASE("DEH"), 
		OXAZOLINE("OXZN"),
		THIAZOLINE("THZN"), 
		AVI_CYSTEINE_LINKAGE("ACL"), 
		LABIONIN_LINKAGE("LAL"),
		TRIFOLOTOXIN_REDUCTION("TR"), 
		THIOL_ESTHER("THE"), 
		LACTAM("LCTM"), 
		TERT_THIO_ETHER("TTE"),
		;
		
		private final String abbreviation;
		private TailoringDomainEnums(final String abbreviation) {
			this.abbreviation = abbreviation;
		}
		
		/**
		 * Get the abbreviation associated with this domain.
		 * @return	the abbreviation 
		 */
		public String getAbbreviation() {
			return(this.abbreviation);
		}
	}
	
	public static TailoringDomainEnums getTailoringDomainFromAbbreviation(String abbreviation) {
		for(TailoringDomainEnums currentEnum : TailoringDomainEnums.values()) {
			if(currentEnum.abbreviation.equals(abbreviation)) {
				return(currentEnum);
			}
		}
		System.out.println("Warning: tailoring domain abbreviation " + abbreviation + " does not exist - skipping.");
		return(null);
	}

	
		
	/**
	 * Enum for fatty acid tails added as starter units by fatty acyl-AMP ligase.
	 * @author skinnider
	 *
	 */
	public enum FattyAcidsEnum {
		_4_CARBONS("CCCC(I)=O", "C4"),
		_6_CARBONS("CCCCCC(I)=O", "C6"),
		_8_CARBONS("CCCCCCCC(I)=O", "C8"),
		_9_CARBONS("O=C(I)CCCCCCCC", "C9"),
		_10_CARBONS("O=C(I)CCCCCCCCC", "C10"),
		_11_CARBONS("O=C(I)CCCCCCCCCC", "C11"),
		_12_CARBONS("O=C(I)CCCCCCCCCCC", "C12"),
		_3_HYDROXYBUTANOIC_ACID("CC(O)CC(I)=O" ,"3HB"),
		PROLINE("IC(C1NCCC1)=O", "Pro");		
		private final String smiles;
		private final String abbreviation;
		private FattyAcidsEnum(final String smiles, final String abbreviation) {
			this.smiles = smiles;
			this.abbreviation = abbreviation;
		}
		
		public String smiles() {
			return smiles;
		}
		
		public String abbreviation() {
			return abbreviation;
		}
		
		public IAtomContainer mol() {
			return MoleculeClasses.getMol(smiles);
		}
	}
	
	/**
	 * Enum for matching small polyketide units that the polyketide predictor cannot determine due to size
	 * @author cDejong
	 *
	 */
	public enum SmallPKunits {
		ethylmalonate_K("O=C(O)CCC", new PKsubstrate[]{PKsubstrate.ETHYLMALONYL}, new PolyKetideDomainEnums[]{PolyKetideDomainEnums.KETONE}),
		_2_methoxyacetic_acid("OC(COC)=O", new PKsubstrate[]{PKsubstrate.MALONYL}, new PolyKetideDomainEnums[]{PolyKetideDomainEnums.HYDROXYL}),
		ethanol("CCO", new PKsubstrate[]{PKsubstrate.MALONYL}, new PolyKetideDomainEnums[]{PolyKetideDomainEnums.HYDROXYL}),
		_2_methylbut_2_enoic_acid("OC(/C(C)=C/C)=O", new PKsubstrate[]{PKsubstrate.METHOXYLMALONYL, PKsubstrate.MALONYL}, new PolyKetideDomainEnums[]{PolyKetideDomainEnums.KETONE, PolyKetideDomainEnums.DOUBLEBOND}), //TODO FIX
		_2_methylbutanoic_acid("CCC(C)C(O)=O", new PKsubstrate[]{PKsubstrate.MALONYL, PKsubstrate.METHYLMALONYL}, new PolyKetideDomainEnums[]{PolyKetideDomainEnums.KETONE, PolyKetideDomainEnums.KETONE}),
		isobutryic_acid("OC(C(C)C)=O", new PKsubstrate[]{PKsubstrate.ISOBUTRYL}, new PolyKetideDomainEnums[]{PolyKetideDomainEnums.KETONE}),
		_3_methoxyacrylic_acid("OC(/C=C/OC)=O", new PKsubstrate[]{PKsubstrate.MALONYL}, new PolyKetideDomainEnums[]{PolyKetideDomainEnums.KETONE}),
		acrylic_acid("OC(C=C)=O", new PKsubstrate[]{PKsubstrate.ETHYLMALONYL}, new PolyKetideDomainEnums[]{PolyKetideDomainEnums.KETONE}),
		methyl_malonate("CCC(O)=O", new PKsubstrate[]{PKsubstrate.METHYLMALONYL}, new PolyKetideDomainEnums[]{PolyKetideDomainEnums.KETONE}),
		methyl_malonate2("CC(O)=O", new PKsubstrate[]{PKsubstrate.METHYLMALONYL}, new PolyKetideDomainEnums[]{PolyKetideDomainEnums.KETONE}),
		methyl_malonate_hydroxyl("CCC(O)=O", new PKsubstrate[]{PKsubstrate.METHYLMALONYL}, new PolyKetideDomainEnums[]{PolyKetideDomainEnums.KETONE}),
		malonate("CC(O)C(O)=O", new PKsubstrate[]{PKsubstrate.MALONYL}, new PolyKetideDomainEnums[]{PolyKetideDomainEnums.KETONE}),
		methoxy_malonate("COC(C)C(O)=O", new PKsubstrate[]{PKsubstrate.METHOXYLMALONYL}, new PolyKetideDomainEnums[]{PolyKetideDomainEnums.HYDROXYL}), //OCC(O)=O
		hydroxy_malonate("OCC(O)=O", new PKsubstrate[]{PKsubstrate.MALONYL}, new PolyKetideDomainEnums[]{PolyKetideDomainEnums.HYDROXYL}),
		;
		private final String smiles;
		private final PKsubstrate[] pkSubstrate;
		private final PolyKetideDomainEnums[] pkDomain;
		private SmallPKunits(final String smiles, final PKsubstrate[] pkSubstrate, final PolyKetideDomainEnums[] pkDomain) {
			this.smiles = smiles;
			this.pkSubstrate = pkSubstrate;
			this.pkDomain = pkDomain;
		}
		
		public String smiles() {
			return smiles;
		}
		public PKsubstrate[] pkSubstrates(){
			return pkSubstrate;
		}
		public PolyKetideDomainEnums[] pkDomains(){
			return pkDomain;
		}
		public IAtomContainer mol() {
			return MoleculeClasses.getMol(smiles);
		}
	}
		
	/**
	 * Enum for C-starter peptide starter units.
	 * @author skinnider
	 *
	 */
	public enum CStarterSubstrate {
	CdaPS1			("O=C(C1C(O1)CCC)"),
	tioR_C1_start	("C(C1=NC2=C(C=C1O)C=CC=C2)=O"),
	antE			("C(C1=CC=CC(NC=O)=C1O)=O"),
	Skyllamycin		("C\\C=C\\C1=CC=CC(CCC=O)=C1"),
	Echinomycin		("C(C1=NC2=CC=CC=C2N=C1)=O"),
	WS9326			("C\\C=C\\C1=CC=CC(CCC=O)=C1"),
	act2_C1_start	("CC1=CC=C(C=O)C(N)=C1O"),
	bacil2_C1_start	("OC1=CC=CC(C=O)=C1O"),
	prist1_C1_start	("OC1=C(C=O)N=CC=C1"),
	cdaps1_C1_start	("CCCC1OC1C=O"),
	fk506			("CC1CC[C@@H](O)[C@H](O)C1"),
	_3_amino_5_hydroxybenzoic_acid	("OC1=CC(N)=CC(C(O)=O)=C1"); 
		
		private final String smiles;
		private CStarterSubstrate(final String smiles) {
			this.smiles = smiles;
		}
		
		public String smiles() {
			return smiles;
		}
		
		public static Map<CStarterSubstrate, IAtomContainer> all(){
			Map<CStarterSubstrate, IAtomContainer> all = new HashMap<CStarterSubstrate, IAtomContainer>();
			
			for(CStarterSubstrate substrate : CStarterSubstrate.values()){
				IAtomContainer mol = null;
				try{mol = SmilesIO.readSmilesTemplates(substrate.smiles());}catch(Exception e){}
				all.put(substrate, mol);
			}
			return all;
		}
		
		public static String[] names() {
		    return Arrays.toString(CStarterSubstrate.values()).replaceAll("^.|.$", "").split(", ");
		}
		
	}
	
}
