package ca.mcmaster.magarveylab.grape.enums;

import java.util.LinkedHashMap;
import java.util.Map;
import org.openscience.cdk.interfaces.IMolecule;

import ca.mcmaster.magarveylab.grape.util.io.SmilesIO;

public class MoleculeClasses {
	
	public enum ChemicalType {
		NRP, NRP_PK_HYBRID, PK, AMINOGLYCOSIDE, UNKNOWN, BADSMILES
	}
	
	public enum ChemicalSubType{
		TYPE_1_PK, MACROLIDE, TYPE_2_PK, STANDARD, NON_TYPE_2_AROMATIC, TYPE_2, TERPENE, FUNGAL, ENEDYINE, CAROTENOID, NONE
	}
	
	public interface CompoundScaffold{}
	
	public enum AromaticFungal implements CompoundScaffold {
		Annularin("CCCC1CC(CC(=O)O1)OC"),
		Fungal1("C1CCC2CC3OCCCC3CC2C1"),
		Fungal2("C1CCC2CC3CCCCC3CC2C1"),
		Fungal5("C1CCC2C(C1)OC1CCCCC21"),
		Fungal6("CC1(C)COC2C1CC1CCCC3CCCC2C13"),
		Fungal7("CC1(C)COC2C1CC1CCOC3CCCC2C13"),
		Fungal8("CC1(C)COC2C1CC1COCC3CCCC2C13"),
		Fungal4("C1CCC2C(C1)COC1CCCCC21"),
		Lagopodin("CC1(C)CC(O)CC1C1=CC(=O)C=CC1=O"),
		Fungal3("C1CCC2CC3OC4CCCCC4CC3CC2C1"),
		Xanthone2("O=C1C2CCCCC2OC2C3C(CCC12)OC1OCCC31"),
		Xanthone("O=C1C2CCCCC2OC2CCCCC12"),
		Fuscin("CC1(C)CCC2C(CCC3COCCC23)O1"),
		Fuscinarin("CC1(C)CCC2C3COCC3CCC2O1"),
		Aporpinone("CC#C\\C=C1/OC(=O)C=C1"),
		Cyclopaidic_acid("OCC1C2C(O)OC(=O)C2C(O)CC1O"),
		Phenalenone1("CC(C)CCOC1CCC2CCCC3CCCC1C23"),
		Phenalenone4("CC(C)CCOC1CCC2CCOC3CCCC1C23"),
		Phenalenone5("CC(C)CCOC1CCC2COCC3CCCC1C23"),
		Phenalenone2("OC1CCC2CCCC3CCCC1C23"),
		Phenalenone3("OC1CCC2CCOC3CCCC1C23"),
		Phenalenone6("OC1CCC2COCC3CCCC1C23");
		
		private final String smiles;
		private AromaticFungal(final String smiles) {
			this.smiles = smiles;
		}
		
		public static Map<String, IMolecule> getAll(){
			Map<String, IMolecule> all = new LinkedHashMap<String, IMolecule>(); 
			
			for(AromaticFungal single : AromaticFungal.values()){
				all.put(single.toString(), getMol(single.smiles));
			}
			return all;
		}
	}
	
	public enum AromaticNotType2 implements CompoundScaffold {
		Rebeccamycin_agly("OC1NC(O)C2C1C1C(NC3C1CCCC3)C1C2C2C(N1)CCCC2"),
		Teradecomycin("CC1OCC2C1OC1CC3CCCCC3C2C1"),
		Bisindole_malleam("OC1NC(O)C(C2CNC3C2CCCC3)C1C1CNC2C1CCCC2"),
		Violacei("OC1NC(C\\C1C1/C(O)NC2C1CCCC2)C1CNC2C1CCCC2"),
		Bis_indole_quinon("N1CC(C2C1CCCC2)C1CCC(CC1)C1CNC2C1CCCC2"),
		Pre_violacei("OC1NC(CC1C1CNC2C1CCCC2)C1CNC2C1CCCC2"),
		Pulvinic_acid2("OC(O)C(C1OC(O)C(C1O)C1CCCCC1)C1CCCCC1"),
		Flavone("OC1CC(OC2CCCCC12)C1CCCCC1"),
		Pulvinic_acid1("OC1C(C(O)O\\C1C/C1CCCCC1)C1CCCCC1"),
		Bis_aryl_quinon("C1CCC(CC1)C1CCC(CC1)C1CCCCC1"),
		Quinolon("CCCCCC1C(O)C(O)C2C(N1)CCCC2"),
		Coumanin("NC1C(O)C2CCC(O)CC2OC1O"),
		Phenazin("C1CC2NC3CCCCC3NC2CC1"),
		Phenothazin("N1C2CCCCC2SC2C1CCCC2"),
		Phenoxazin("N1C2CCCCC2OC2C1CCCC2"),
		Carbolin("C1CC2C(CN1)NC1C2CCCC1"),
		Pyoluteorin_es1("OC(C1CCCN1)C1CCCCC1"),
		Pyoluteorin_es2("C(C1CCCN1)C1CCCCC1"),
		Chromophyrrolic_a("OC(O)C1CCC(N1)C(O)O"),
		Pyrrolintrin_fan("N1CCC(C1)C1CCCCC1"),
		Indolotryptolin("OC1C(O)C(O)NC1O"),
		//Maleimides1("OC1NC(O)CC1"),
		//Maleimides2("OC1NCCC1");
		;
		
		private final String smiles;
		private AromaticNotType2(final String smiles) {
			this.smiles = smiles;
		}
		
		public static Map<String, IMolecule> getAll(){
			Map<String, IMolecule> all = new LinkedHashMap<String, IMolecule>(); 
			
			for(AromaticNotType2 single : AromaticNotType2.values()){
				all.put(single.toString(), getMol(single.smiles));
			}
			return all;
		}
	}
	
	public enum AromaticType2 implements CompoundScaffold {
		FR_9001("CCC(C)C1OC2(OC1(C)C1C=CC3CC(CCC3C21C)C(O)O)C1(OC(C)O)C(O)OC(C)C1O"),
		Enterocin("OC(C1C2CC3CC(C1C(O)O3)C2C1CCCCO1)C1CCCCC1"),
		AZ154("C1CCC(CC1)C1CCCC2CC3CC4CCCCC4CC3CC12"),
		Fredricamycin4("C1CCC2CC3CC4(CCC3CC2C1)CCC1CC2CCNCC2CC1C4"),
		Fredericamycin3("C1C2CC3CCCCC3CC2CC11CCC2CC3CCNCC3CC2C1"),
		Fredericamycin2("C1C2CC3CCNCC3CC2CC11CCC2CC3CCCCC3CC2C1"),
		Fredericamycin1("C1C2CC3CCCCC3CC2CC11CC2CC3CCNCC3CC2C1"),
		Aureolic_acid1("CC(O)C(O)C(O)C(O)C1CCC2CC3CCCCC3CC2C1"),
		Lysolipi("C1CCC2OC3CC4CCC5CC6CCNCC6CC5C4CC3CC2C1"),
		Rubromycin4("C1C2CC3CCCCC3CC2OC11CC2CC3CCOCC3CC2O1"),
		Rubromycin3("C1C2CC3CCOCC3CC2OC11CCC2CC3CCCCC3CC2O1"),
		Rubromycin2("C1CCC2CC3OC4(CCC3CC2C1)CCC1CC2CCOCC2CC1O4"),
		Rubromycin1("C1C2CC3CCCCC3CC2OC11CCC2CC3CCOCC3CC2O1"),
		Chartari("C1CCC2C(C1)CC1COC3CCCC4COC2C1C34"),
		Resistomycin("C1CC2CCC3CCC4CCCC5CC(C1)C2C3C45"),
		Erdacimycin("C1C2CCC3CCCCC3C2C2CCC3CCCCC3C12"),
		Urdamycin("C1CCC2CC3C(CCC4CCCCCC34)CC2C1"),
		Pradimicin2("C1CCC2CC3CC4C(CCC5CCCCC45)CC3CC2C1"),
		Pradimicin1("C1CCC2OC3CC4CCC5CCCCC5C4CC3CC2C1"),
		Wailupemycin("C1CCC(OC1)C1CCCC2CCCCC12"),
		Angycyclines4("C1CCC2CC3C(CC2C1)NCC1CCCOC31"),
		Angucyclines3("C1CCC2CC3C(CCC4CCCOC34)CC2C1"),
		Angucyclines2("C1CCC2CC3C(CCC4CCCCC34)CC2C1"),
		Angucyclines1("C1CCC2CC3C(CC2C1)NCC1CCCCC31"),
		Brnzonaphthopyro("C1CCC2C(C1)CCC1C3CCCCC3COC21"),
		tetracyclines("C1CCC2CC3CC4CCCCC4CC3CC2C1"),
		Kinobscurin("C1C2CCCCC2C2CC3CCCCC3CC12"),
		Nogalonate("OC1C2CCCCC2C(O)C2CCCCC12"),
		Benziosochromanequinone("OC1C2CCCCC2CC2CCOCC12"),
		UT_X26("C1CCC2CC34CCCC(CCC3CC2C1)C4"),
		AB649("C1C2CCCCC2C2CCC3CCCCC3C12"),
		Juglomycin("OC1CCC(O)C2CCCCC12");

		
		private final String smiles;
		private AromaticType2(final String smiles) {
			this.smiles = smiles;
		}
		
		public static Map<String, IMolecule> getAll(){
			Map<String, IMolecule> all = new LinkedHashMap<String, IMolecule>(); 
			
			for(AromaticType2 single : AromaticType2.values()){
				all.put(single.toString(), getMol(single.smiles));
			}
			return all;
		}
	}
	
	public enum Terpenes implements CompoundScaffold {
		Mutilin("CCC1(C)CC[C@@H]2CCC[C@]3(CCC[C@@H]23)CC1"),
		Hopane("C1CC2CCC3C(CCC4C5CCCCC5CCC34)C2C1"),
		Sesquiterpene("C1CC2CCC3C(CCC4CCCCC34)C2C1"),
		Melledonal("CC1(C)CC2CCC3CCC3C2C1"),
		Abietadien("CC(C)C1CCC2C(CCC3C(C)(C)CCCC23C)C1"),
		Cyclocitrin("C1CC2CCC3C(CCC4CCCCC3C4)C2C1"),
		Cyathane("CC(C)C1CCC2CCC3CCCCCC3C12"),
		Pimaradiene("CC1(C)CCCC2(C)C3CCC(C)(CC3CCC12)CC"),
		Albaflavone("CC1(C)CC2CCCC22CCC1C2"),
		Acremostrictin("CC1CC2CC3(C)OC(O)C(C3C)C2C1"),
		Aszonapyrone("CC1(C)CCCC2(C)C1CCC1(C)C(CC3CCCOC3)CCCC21");

		private final String smiles;
		private Terpenes(final String smiles) {
			this.smiles = smiles;
		}
		
		public static Map<String, IMolecule> getAll(){
			Map<String, IMolecule> all = new LinkedHashMap<String, IMolecule>(); 
			
			for(Terpenes single : Terpenes.values()){
				all.put(single.toString(), getMol(single.smiles));
			}
			return all;
		}
	}
	
	public enum Carotenoid implements CompoundScaffold{
		Carotenoid
	}
	
	public enum Enedyine implements CompoundScaffold{
		NineMember(new String[]{"C1CC#CCCC#CC1", "C1CC#C\\C=C/C#CC1", "C1CC#C\\C=C/CC#C1"}),
		TenMember(new String[]{"C1CCC#C\\C=C/C#CC1"});

		private final String[] smiles;
		private Enedyine(final String[] smiles) {
			this.smiles = smiles;
		}
		
		public static Map<String, IMolecule[]> getAll(){
			Map<String, IMolecule[]> all = new LinkedHashMap<String, IMolecule[]>(); 
			
			for(Enedyine single : Enedyine.values()){
				all.put(single.toString(), getMols(single.smiles));
			}
			return all;
		}
	}
	
	public static IMolecule getMol(String smiles){
		IMolecule mol = null;
		try {
			mol = SmilesIO.readSmiles(smiles);
			mol = SmilesIO.readSmiles(SmilesIO.generateSmiles(mol)); //Regenerates the SMILES
		} catch (Exception e) {
			System.out.println("Bad SMILES string: " + smiles);
		}
		return mol;
	}
	
	public static IMolecule[] getMols(String[] smiles){
		IMolecule[] mols = new IMolecule[smiles.length];
		try {
			for(int i = 0; smiles.length > i; i++){
				IMolecule mol = SmilesIO.readSmiles(smiles[i]);
				mols[i] = (SmilesIO.readSmiles(SmilesIO.generateSmiles(mol))); //Regenerates the SMILES
			}
		} catch (Exception e) {
			System.out.println("Bad SMILES string: " + smiles);
		}
		return mols;
	}
}

