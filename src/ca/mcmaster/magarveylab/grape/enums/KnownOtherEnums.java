package ca.mcmaster.magarveylab.grape.enums;

import java.io.IOException;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;

import ca.mcmaster.magarveylab.grape.util.io.SmilesIO;

public enum KnownOtherEnums{
	ChrA("ChrA", "azotobactins chromophore", "C1=C(C(=CC2=C1N4C3C(=C2)NC(N3CCC4C(O)=O)=O)O)O"),
	ChrD("ChrD", "5,6-dihydropyoverdin chromophore", "C1=C(C(=CC2=C1N3C(C(N)C2)NCCC3C(O)=O)O)O"),
	ChrI("ChrI", "isopyoverdin chromophore", "C1=C(C(=CC2=C1N3C(C(=C2)N)NC(CC3)C(O)=O)O)O"),
	ChrP("ChrP", "pyoverdin chromophore", "C1=C(C(=CC2=C1N3C(C(=C2)N)NCCC3C(O)=O)O)O"),
	ChrAct("ChrAct", "actinomycin chromophore", "CC1=C2C(=C(C=C1)C(=O)O)N=C3C(=C(C(=O)C(=C3O2)C)N)C(=O)O"),
	Pyr("Pyr", "pyrrolidone", "C1CC(=O)NC1"),
	Dpy("Dpy", "Dolapyrrolidone", "COC1=CC(=O)NC1CC2=CC=CC=C2"),
	OH_Pyr("OH-Pyr", "hydroxy pyrrolidone", "C1[C@@H](CNC1=O)O"),
	Pya("Pya", "pyruvate", "CC(=O)C(=O)O"),
	dPyr("dPyr", "dehydropyrrolidone", "C1(CC(C(=CC(=O)O)N1)N)=O"),
	NSpd("NSpd", "norspermidine", "C(CN)CNCCCN"),
	Spd("Spd", "spermidine", "C(CCNCCCN)CN"),
	GSpd("GSpd", "guanylspermidine", "C(CCNCCCN=C(N)N)CN"),
	Ist("Ist", "isostatine", "CCC(C)C(C(CC(=O)O)O)N"),
	Nst("Nst", "norstatine", "CC(C)CC(C(C(=O)O)O)N"),
	Sta("Sta", "statine", "CC(C)CC(C(CC(=O)O)O)N"),
	DMOG("DMOG", "DHP-methyloxazolinyl group", "C1(=NC(C(O1)C)C(O)=O)C2=CC=CC(=C2O)O"),
	Choi("Choi", "2-carboxy-6-hydroxyoctahydroindole", "C1C(NC2C1CCC(C2)O)C(=O)O"),
	Ibu("Ibu", "4-amino-2,2-dimethyl-3-oxopentanoic acid", "CC(C(C(C(O)=O)(C)C)=O)N"),
	Hpa("Hpa", "hydroxypicolinic acid", "C1=CC(=C(N=C1)C(=O)O)O"),
	Map("Map", "2-methyl-3-aminopentanoic acid", "C(C(C(C(=O)O)C)N)C"),
	MdCP("MdCP", "N-methyldichloropyrrole-2-carboxylic acid", "C1(=CC(=C([N]1C)Cl)Cl)C(=O)O"),
	Me_AOA("Me-AOA", "methyl-2-aminooctanoic acid", "C(CCC(C(=O)O)N)CCCC"),
	Daz("Daz", "2,6-diamino-7-hydroxyazelaic acid", "C(CC(C(CC(=O)O)O)N)CC(C(=O)O)N"),
	Doe("Doe", "Dolaphenine", "C1=CC=C(C=C1)CC(C2=NC=CS2)N"),
	MCP("MCP", "N-methylchloropyrrole", "CN1C=C(C=C1C(=O)O)Cl"),
	CO("CO", "oxo group", "C=O"),
	DHMDA("DHMDA", "8,10-Dimethyl-9-hydroxy-7-methoxytridecadienoic acid", "C(C(C(CCC)C)O)(C(C=CC=CCC(=O)O)OC)C"),
	Agdha("Agdha", "4-amino-7-guanidino-2,3-dihydroxyheptanoic acid", "C(C(C(C(O)=O)O)O)(CCCNC(=N)N)N"),
	Ahp("Ahp", "3-amino-6-hydroxy-2-piperidone", "C1C(N(C(C(C1)N)=O)C(C(C)O)C(=O)O)O"),
	Dov("Dov", "Dolavaline", "CC(C)C(N(C)C)C(=O)O"),
	Eta("Eta", "Ethanolamine", "C(CO)N"),
	DHPT("DHPT", "dihydroxyphenylthiazol group", "C1(=NC(CS1)C(O)=O)C2=C(C(=CC=C2)O)O"),
	Hpoe("Hpoe", "2-hydroxyphenyl-2-oxo-ethanoic acid", "C1=CC(=CC=C1C(=O)C(=O)O)O"),
	HseL("HseL", "homoserine lactone", "C1COC(=O)C1N"),
	Aca("Aca", "Anticapsin", "C1CC(=O)[C@H]2[C@@H]([C@H]1C[C@@H](C(=O)O)N)O2"),
	NMe_Lan("NMe-Lan", "N-Methyl-Lanthionine", "C(SCC(C(=O)O)N)C(C(O)=O)NO"),
	PTTA("PTTA", "4-propenoyl-2-tyrosylthiazole acid", "C(O)(=O)C=CC1=CSC(=N1)C(N)CC2=CC=C(C=C2)O"),
	Mab("Mab", "2-methyl-3-aminobutanoic acid", "CC(C(C)N)C(=O)O"),
	NMe_Dha("NMe-Dha", "N-Methyl-dehydroalanine", "CNC(=C)C(=O)O"),
	Pda("Pda", "pentanedioic acid", "C(CC(=O)O)CC(=O)O"),
	Pha("Pha", "phenylacetic acid", "C1=CC=C(C=C1)CC(=O)O"),
	PT("PT", "phosphinothricin", "CP(=O)(CCC(C(=O)O)N)O"),
	Put("Put", "putrescine", "C(CCN)CN"),
	COOH_Qui("COOH-Qui", "2-carboxyquinoxaline", "C1=CC=C2C(=C1)N=CC(=N2)C(=O)O"),
	Azd("Azd", "aziridine dicarboxylic acid", "C1(C(N1)C(=O)O)C(=O)O"),
	N_OH_Hta("N-OH-Hta", "N-Hydroxy-histamine", "C1(CCCN1)CCNO"),
	PALOA("PALOA", "propenoyl-alanyloxazole acid", "C(C=CC1=COC(=N1)C(N)C)(=O)O"),
	PAOA("PAOA", "propenoyl-2-aminobutanoyloxazole acid", "C(C=CC1=COC(=N1)C(N)CC)(=O)O"),
	PMST("PMST", "propenoyl-O-methylserinylthiazole acid", "C(C=CC1=CSC(=N1)C(N)COC)(=O)O"),
	DMAH("DMAH", "dialkylmaleic Anhydride", "O=C(C([C@@H](CC(O)=O)O)=C1C)OC1=O"),
	Ac("Ac", "acetate", "CC(O)=O"),
	IV("IV", "isovalerate", "CC(C)CC(O)=O"),
	Deu("Deu", "3’-deoxy-4'-enamine uridine", "O=C(N1)C=CN(C2O/C(CC2O)=C/N)C1=O"),
	Hita("Hita", "6-hydroxy-tetrahydroisoquinoline carboxylic acid", "CC1NC(C(O)=O)CC2=C1C=CC(O)=C2"),
	dHis("dHis", "decarboxy histidine", "ONCCC1=CNC=N1"),
	MeOH("MeOH", "methanol", "CO"),
	Cba("Cba", "carbamic acid", "OC(N)=O"),
	Iva("Iva", "isovaleric acid", "O=C(O)CC(C)C"),
	//siderophore
	OHput("OHput", "hydroxyputrescine", "NCC(O)CCNO"),
	DMHQ("DMHQ", "dimethyl_hydroquinolium", "NC1C[N+](C)(C)C2=C(C=C(O)C(O)=C2)C1"),
	Cad("Cad", "cadaverine", "NCCCCCN"),
	OH_Cad("OH-Cad", "hydroxycadaverine", "NCCCCCNO"),
	Acp("Acp", "azotobactin_chomophore", "O=C1NC2=CC3=CC(O)=C(O)C=C3[N+]4=C2N1CCC4C(O)=O"),
	HMPA("HMPA", "_5-hydroxy-3-methylpent-2-enoic_acid", "OCC/C(C)=C/C(O)=O"),
	Can("Can", "cinnamic_acid", "O=C(O)/C=C/C1=CC=CC=C1"),
	HArg("HArg", "hydroarginine", "O=C(O)C(N)CCCNC(N)N"),
	Dapp("Dapp", "1_3_diaminopropane", "NCCCN"),
	CA("CA", "citric_acid", "OC(CC(C(O)=O)(O)CC(O)=O)=O"),
	CGDHBA("CGDHBA", "C-glycosylated_dihydroxybenzoic_acid", "OCC1OC(C(O)C(O)C1O)C2=CC(O)=C(O)C(C(O)=O)=C2"),
	SA("SA", "succinic_acid", "OC(CCC(O)=O)=O"),
	;
	private final String abbreviation;
	private final String fullName;
	private final String smiles;
	private KnownOtherEnums(final String abbreviation, final String fullName, final String smiles) {
		this.abbreviation = abbreviation;
		this.fullName = fullName;
		this.smiles = smiles;
	}
	
	public String getAbbreviation() {
		return(this.abbreviation);
	}
	public String getfullName() {
		return(this.fullName);
	}
	public String getSmiles() {
		return(this.smiles);
	}
	public IAtomContainer getMol(){
		try {
			return(SmilesIO.readSmilesTemplates(this.smiles));
		} catch (IOException | CDKException e) {
			System.err.println("KnownOther: " + this.abbreviation + " Could not be converted");
		}
		return null;
	}

}
