package ca.mcmaster.magarveylab.grape.enums;

import java.util.LinkedHashMap;
import java.util.Map;

import org.openscience.cdk.interfaces.IAtomContainer;

public enum SugarEnums{
	GLUCOSE("Glucose" ,"OC[C@H]1OC([C@@H]([C@H]([C@@H]1O)O)O)O"),
	N_ACETYLGLUCOSAMINE("N-acetylglucosamine" ,"O=C(C)N[C@@H]1[C@H]([C@@H]([C@@H](CO)OC1O)O)O"),
	MANNOSE("Mannose" ,"OC(O1)[C@@H](O)[C@@H](O)[C@H](O)[C@H]1CO"),
	GULOSE("Gulose" ,"O[C@H](C(O)O1)[C@@H](O)[C@H](O)[C@@H]1CO"),
	L_ACULOSE("L-aculose" ,"CC1O[C@H](C=CC1=O)O"),
	L_CINERULOSE_A("L-cinerulose A" ,"CC1O[C@H](CCC1=O)O"),
	L_RHODINOSE("L-rhodinose" ,"C[C@@H]1O[C@@H](CC[C@@H]1O)O"),
	REDNOSE("Rednose" ,"CC1OC(C(N)=CC1=O)O"),
	L_CINERULOSE_B("L-cinerulose B" ,"CC1O[C@H]([C@@H](O)CC1=O)O"),
	O_METHYL_L_AMICETOSE("O-methyl-L-amicetose" ,"COC1CC[C@@H](O[C@H]1C)O"),
	_4_O_METHYL_L_RHODINOSE("4-O-methyl-L-rhodinose" ,"CO[C@H]1CC[C@@H](O[C@H]1C)O"),
	L_DAUNOSAMINE("L-daunosamine" ,"C[C@@H]1OC(C[C@@H]([C@@H]1O)N)O"),
	L_RISTOSAMINE("L-ristosamine" ,"CC1OC(CC(C1O)N)O"),
	D_DIGITOXOSE("D-digitoxose" ,"CC1OC(CC(C1O)O)O"),
	L_DIGITOXOSE("L-digitoxose" ,"CC1OC(CC(C1O)O)O"),
	_2_DEOXY_L_FUCOSE("2-deoxy-L-fucose" ,"CC1OC(CC(C1O)O)O"),
	D_OLIVOSE("D-olivose" ,"CC1O[C@@H](C[C@H]([C@@H]1O)O)O"),
	D_OLIOSE("D-oliose" ,"CC1O[C@H](C[C@H]([C@H]1O)O)O"),
	_4_OXO_L_VANCOSAMINE("4-oxo-L-vancosamine" ,"C[C@@H]1OC(C[C@](N)(C1=O)C)O"),
	D_FOROSAMINE("D-forosamine" ,"CC1OC(CC[C@@H]1N(C)C)O"),
	L_ACTINOSAMINE("L-actinosamine" ,"COC1C(OC(CC1N)O)C"),
	L_VANCOSAMINE("L-vancosamine" ,"OC1O[C@H]([C@@H](O)[C@](C1)(N)C)C"),
	L_VICENISAMINE("L-vicenisamine" ,"CN[C@@H]1C(CC(OC1C)O)O"),
	D_CHALCOSE("D-chalcose" ,"CO[C@H]1C[C@H](OC([C@@H]1O)O)C"),
	D_MYCAROSE("D-mycarose" ,"CC1OC(C[C@](O)([C@H]1O)C)O"),
	L_OLEANDROSE("L-oleandrose" ,"CO[C@H]1C[C@H](O[C@H]([C@@H]1O)C)O"),
	OLIVOMOSE("Olivomose" ,"COC1C(CC(OC1C)O)O"),
	D_MYCOSAMINE("D-mycosamine" ,"C[C@H]1O[C@H]([C@H]([C@H]([C@@H]1O)N)O)O"),
	_4_DEOXY_4_THIO_D_DIGITOXOSE("4-deoxy-4-thio-D-digitoxose" ,"CC1O[C@H](C[C@H](C1S)O)O"),
	D_FUCOFURANOSE("D-fucofuranose" ,"C[C@H]([C@H]1O[C@H]([C@H](C1O)O)O)O"),
	D_FUCOSE("D-fucose" ,"CC1OC([C@@H]([C@H]([C@H]1O)O)O)O"),
	L_RHAMNOSE("L-rhamnose" ,"C[C@@H]1O[C@@H]([C@@H]([C@@H]([C@H]1O)O)O)O"),
	_4_N_ETHYL_4_AMINO_3_O_METHOXY_245_TRIDEOXYPENTOSE("4-N-ethyl-4-amino-3-O-methoxy-2,4,5-trideoxypentose" ,"CCN[C@H]1CO[C@H](C[C@H]1OC)O"),
	D_3_N_METHYL_4_O_METHYL_L_RISTOSAMINE("D-3-N-methyl-4-O-methyl-L-ristosamine" ,"CN[C@H]1CC(OC([C@@H]1OC)C)O"),
	NN_DIMETHYL_L_PYRROLOSAMINE("N,N-dimethyl-L-pyrrolosamine" ,"CC1OC(CC(C1N(C)C)O)O"),
	D_DESOSAMINE("D-desosamine" ,"C[C@@H]1C[C@H](N(C)C)[C@@H](O)[C@H](OO)O1"),
	L_MEGOSAMINE("L-megosamine" ,"C[C@@H]1OC(C[C@@H](N(C)C)[C@H]1O)O"),
	NOGALAMINE("Nogalamine" ,"OC1O[C@@H](C)[C@H](O)[C@@H](N(C)C)[C@@H]1O"),
	L_RHODOSAMINE("L-rhodosamine" ,"CC1O[C@H](CC(N(C)C)[C@H]1O)O"),
	D_ANGOLOSAMINE("D-angolosamine" ,"C[C@@H]1OC(CC(N(C)C)[C@H]1O)O"),
	KEDAROSAMINE("Kedarosamine" ,"OC1O[C@@H](C)[C@@H](N(C)C)[C@@H](O)C1"),
	L_NOVIOSE("L-noviose" ,"CC1(C)[C@H](O)[C@@H](O)[C@@H](O)C(O)O1"),
	L_CLADINOSE("L-cladinose" ,"C[C@H]1[C@H](O)[C@](C)(OC)CC(O)O1"),
	_2_N_METHYL_D_FUCOSAMINE("2'-N-methyl-D-fucosamine" ,"CN[C@H]1C(O[C@@H]([C@@H]([C@@H]1O)O)C)O"),
	D_DIGITALOSE("D-digitalose" ,"CO[C@H]1[C@H]([C@H](O[C@H]([C@@H]1O)O)C)O"),
	_3_O_METHYL_RHAMNOSE("3-O-methyl-rhamnose" ,"COC1[C@H]([C@H](OC([C@H]1O)O)C)O"),
	_2_O_METHYL_RHAMNOSE("2-O-methyl-rhamnose" ,"COC1[C@@H](OC([C@@H](C1O)O)C)O"),
	_6_DEOXY_3_C_METHYL_L_MANNOSE("6-deoxy-3-C-methyl-L-mannose" ,"C[C@@H]1OC([C@@H]([C@](O)([C@H]1O)C)O)O"),
	_46_DIDEOXY_4_HYDROXYLAMINO_D_GLUCOSE("4,6-dideoxy-4-hydroxylamino-D-glucose" ,"CC1OC(C(C(C1NO)O)O)O"),
	_3_NN_DIMETHYL_L_EREMOSAMINE("3-N,N-dimethyl-L-eremosamine" ,"OC1C[C@](C)(N(C)C)[C@@H](O)[C@H](C)O1"),
	CHROMOSE_4_O_ACETYL_BETA_D_OLIOSE("Chromose (4-O-acetyl-beta-D-oliose)" ,"CC1OC(CC(C1OC(C)=O)O)O"),
	_4_O_CARBAMOYL_D_OLIVOSE("4-O-carbamoyl-D-olivose" ,"CC1OC(CC(C1OC(N)=O)O)O"),
	D_RAVIDOSAMINE("D-ravidosamine" ,"C[C@H]1O[C@@H]([C@@H]([C@@H](N(C)C)[C@H]1O)O)O"),
	_3_NN_DIMETHYL_D_MYCOSAMINE("3-N,N-dimethyl-D-mycosamine" ,"C[C@H]1O[C@H]([C@@H]([C@@H](N(C)C)[C@@H]1O)O)O"),
	_23_O_DIMETHYL_L_RHAMNOSE("2,3-O-dimethyl-L-rhamnose" ,"COC1[C@@H](OC([C@@H](C1OC)O)C)O"),
	_24_O_DIMETHYL_L_RHAMNOSE("2,4-O-dimethyl-L-rhamnose" ,"CO[C@H]1C(O[C@H](C(C1O)OC)O)C"),
	_34_O_DIMETHYL_L_RHAMNOSE("3,4-O-dimethyl-L-rhamnose" ,"CO[C@H]1C(O[C@H](C(C1OC)O)O)C"),
	_2_THIOGLUCOSE("2-thioglucose" ,"OC1[C@H](S)[C@@H](O)[C@H](O)[C@@H](CO)O1"),
	OLIVOMYCOSE("Olivomycose" ,"C[C@@H]1O[C@H](C[C@](O)([C@H]1OC(C)=O)C)"),
	_4_NN_DIMETHYLAMINO_4_DEOXY_5_C_METHYL_L_RHAMNOSE("4-N,N-dimethylamino-4-deoxy-5-C-methyl-l-rhamnose" ,"CN(C1C(C(C(OC1(C)C)O)O)O)C"),
	_234_TRI_O_METHYLRHAMNOSE("2,3,4-tri-O-methylrhamnose" ,"CO[C@@H]1[C@H](O[C@@H]([C@H]([C@H]1OC)OC)O)C"),
	_4_O_ACETYL_L_ARCANOSE("4-O-acetyl-L-arcanose" ,"OC1C[C@@](C)(OC)[C@H](OC(C)=O)[C@H](C)O1"),
	_3_N_ACETYL_D_RAVIDOSAMINE("3-N-acetyl-D-ravidosamine" ,"C[C@H]1O[C@@H]([C@@H]([C@@H](N(C(C)=O)C)[C@H]1OC(C)=O)O)O"),
	_3_O_CARBAMOYL_L_NOVIOSE("3-O-carbamoyl-L-noviose" ,"CC1([C@@H]([C@H]([C@H](C(O1)O)O)OC(N)=O)O)C"),
	L_NOGALOSE("L-nogalose" ,"COC1C(OC(C(C1(OC)C)OC)O)C"),
	_4_O_ACETYL_D_RAVIDOSAMINE("4-O-acetyl-D-ravidosamine" ,"C[C@H]1O[C@@H]([C@@H]([C@@H](N(C)C)[C@H]1OC(C)=O)O)O"),
	_3_O_CARBAMOYL_4_O_METHYL_L_NOVIOSE("3-O-carbamoyl-4-O-methyl-L-noviose" ,"CC1(C)[C@H](OC)[C@@H](OC(N)=O)[C@@H](O)C(O)O1"),
	_3_N_ACETYL_4_O_ACETYL_D_RAVIDOSAMINE("3-N-acetyl-4-O-acetyl-D-ravidosamine" ,"C[C@H]1O[C@@H]([C@@H]([C@@H](N(C(C)=O)C)[C@H]1OC(C)=O)O)O"),
	_3_5_METHYL_2_PYRROLYLCARBONYL_4_O_METHYL_L_NOVIOSE("3-(5'-methyl-2'-pyrrolylcarbonyl-)4-O-methyl-L-noviose" ,"CO[C@@H]1[C@@H](C([C@@H](OC1(C)C)O)O)OC(C2=CC=C(N2)C)=O"),
	MADUROSE("Madurose" ,"O[C@@H]1[C@H](O)[C@@](O)(C)[C@H](N)CO1"),
	_4_N_METHYL_4_AMINO_3_O_METHOXY_245_TRIDEOXYPENTOSE("4-N-methyl-4-amino-3-O-methoxy-2,4,5-trideoxypentose" ,"O[C@@H]1C[C@H](OC)[C@@H](NC)CO1"),
	_2_DEOXYSTREPTAMINE("2-deoxystreptamine", "NC1CC(N)C(O)C(O)C1O"),
	GLUCOSAMINE("Glucosamine", "NC1C(O)OC(CO)C(O)C1O"),
	NEOSAMINE_C("Neosamine C", "NCC1OC(O)C(N)C(O)C1O"),
	ROBOSE("Ribose", "OCC1OC(O)C(O)C1O"),
	XYLOSE("Xylose", "OC1COC(O)C(O)C1O"),
	KANOSAMINE("Kanosamine", "NC1C(O)C(O)OC(CO)C1O"),
	_3_DEOXYGLUCOSAMINE("3-deoxyglucosamine", "NC1CC(O)C(CO)OC1O"),
	_2_DEOXYFORTAMINE("2-deoxyfortamine	", "CNC1C(O)C(O)C(N)CC1OC"),
	FORTAMINE("Fortamine", "CNC1C(O)C(O)C(N)C(O)C1OC"),
	_3_DEOXYNEOSAMINE_C("3-deoxyneosamine C	", "NCC1OC(O)C(N)CC1O"),
	STREPTIDINE("Streptidine","NC(N)=NC1C(O)C(O)C(O)C(N=C(N)N)C1O"),
	N_METHYL_GLUCOSAMINE("N-methyl glucosamine","CNC1C(O)OC(CO)C(O)C1O"),
	STREPTOSE("Streptose","CC1OC(O)C(O)C1(O)C=O"),
	DIHYDROSTREPTOSE("Dihydrostreptose","CC1OC(O)C(O)C1(O)CO"),
	BLUENSIDINE("Bluensidine","NC(=N)NC1C(O)C(O)C(O)C(NC(N)=O)C1O"),
	_2L_2_AMINO_2_DEOXY_4_5_O_METHYLENE_NEO_INOSITOL("2L-2-amino-2-deoxy-4,5-O-methylene-neo-inositol","NC1C(O)C2OCOC2C(O)C1O"),
	_5_DEHYDRO_L_FUCOFURANOSE("5-dehydro-L-fucofuranose","CC(=O)C1OC(O)C(O)C1O"),
	KASUGAMINE("Kasugamine","CC1OC(O)C(N)CC1NC(=N)C(O)=O"),
	INOSITOL("Inositol","OC1C(O)C(O)C(O)C(O)C1O"),
	ACTINAMINE("Actinamine","CNC1C(O)C(O)C(O)C(NC)C1O"),
	SPECTINOSE("Spectinose","CC1CC(=O)C(O)C(O)O1"),
	_3_4_DIDEOXY_PURPUROSAMINE("3,4-dideoxy purpurosamine","CC(N)C1CCC(N)C(O)O1"),
	;

	private final String name;
	private final String smiles;
	private SugarEnums(final String name, final String smiles) {
		this.name = name;
		this.smiles = smiles;
	}
	
	public static Map<String, IAtomContainer> getAll(){
		Map<String, IAtomContainer> all = new LinkedHashMap<String, IAtomContainer>(); 
		for(SugarEnums single : SugarEnums.values()){
			all.put(single.fullName(), single.mol());
		}
		return all;
	}
	
	public String fullName(){
		return name;
	}
	
	public String smiles(){
		return smiles;
	}
	
	public IAtomContainer mol(){
		IAtomContainer mol = MoleculeClasses.getMol(smiles);
		return mol;
	}
}
