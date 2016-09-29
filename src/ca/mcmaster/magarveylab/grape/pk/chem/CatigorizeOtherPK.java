package ca.mcmaster.magarveylab.grape.pk.chem;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Map;
import java.util.Map.Entry;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;

import ca.mcmaster.magarveylab.grape.enums.MoleculeClasses.AromaticFungal;
import ca.mcmaster.magarveylab.grape.enums.MoleculeClasses.AromaticNotType2;
import ca.mcmaster.magarveylab.grape.enums.MoleculeClasses.Enedyine;
import ca.mcmaster.magarveylab.grape.enums.MoleculeClasses.ChemicalSubType;
import ca.mcmaster.magarveylab.grape.enums.MoleculeClasses.Terpenes;
import ca.mcmaster.magarveylab.grape.enums.MoleculeClasses.AromaticType2;
import ca.mcmaster.magarveylab.grape.util.io.SmilesIO;

public class CatigorizeOtherPK {
	
	private static Map<String, IMolecule> otherAromaticMolecules = AromaticNotType2.getAll();
	private static Map<String, IMolecule> type2Molecules = AromaticType2.getAll();
	private static Map<String, IMolecule> fungalAromaticMolecules = AromaticFungal.getAll();
	private static Map<String, IMolecule> terpeneMolecules = Terpenes.getAll();
	private static Map<String, IMolecule[]> enedyineMolecules = Enedyine.getAll();
	private ChemicalSubType type;
	private String match;
	
	public void analyzeMolecule(IMolecule mol, int macrolideType) throws CDKException{
		type = ChemicalSubType.STANDARD;
		match = null;
		IMolecule carotenoids = null;
		try{
			carotenoids = SmilesIO.readSmiles("CC(C)=C\\C=C\\C=C(/C)\\C=C\\C=C(C)C"); //specific terpene
		}catch(IOException e){}
	    if(UniversalIsomorphismTester.isSubgraph(mol, carotenoids)){
    		type = ChemicalSubType.TERPENE;
    		match = "Carotenoid";
    		return;
	    }
	    for(Entry<String, IMolecule[]> entry : enedyineMolecules.entrySet()){
	    	String name = entry.getKey();
	    	IMolecule[] scaffolds = entry.getValue();
	    	for(IMolecule scaffold : scaffolds){
		    	if(UniversalIsomorphismTester.isSubgraph(mol, scaffold)){
		    		type = ChemicalSubType.ENEDYINE;
		    		match = name;
		    		return;
		    	}
	    	}
	    }
	    if(macrolideType == 0){
	    	for(Entry<String, IMolecule> entry : fungalAromaticMolecules.entrySet()){
		    	String name = entry.getKey();
		    	IMolecule aro = entry.getValue();
		    	if(UniversalIsomorphismTester.isSubgraph(mol, aro)){
		    		type = ChemicalSubType.FUNGAL;
		    		match = name;
		    		return;
		    	}
		    }
			IMolecule molSingleBondsOnly = onlySingleBonds(mol);
		    for(Entry<String, IMolecule> entry : type2Molecules.entrySet()){
		    	String name = entry.getKey();
		    	IMolecule aro = entry.getValue();
		    	if(UniversalIsomorphismTester.isSubgraph(molSingleBondsOnly, aro)){
		    		type = ChemicalSubType.TYPE_2;
		    		match = name;
		    		return;
		    	}
		    }
		    for(Entry<String, IMolecule> entry : otherAromaticMolecules.entrySet()){
		    	String name = entry.getKey();
		    	IMolecule aro = entry.getValue();
		    	if(UniversalIsomorphismTester.isSubgraph(molSingleBondsOnly, aro)){
		    		type = ChemicalSubType.NON_TYPE_2_AROMATIC;
		    		match = name;
		    		return;
		    	}
		    }
		    for(Entry<String, IMolecule> entry : terpeneMolecules.entrySet()){
		    	String name = entry.getKey();
		    	IMolecule aro = entry.getValue();
		    	if(UniversalIsomorphismTester.isSubgraph(molSingleBondsOnly, aro)){
		    		type = ChemicalSubType.TERPENE;
		    		match = name;
		    		return;
		    	}
		    }
	    }
	}
	
	public ChemicalSubType getChemicalSubType(){
		return type;
	}
	public String getMatchName(){
		return match;
	}
	private static IMolecule onlySingleBonds(IMolecule mol) {
		IMolecule molClone = null;
		try{
			molClone = mol.clone();
		}catch(Exception e){}
		
		ArrayList<IBond> bonds = new ArrayList<IBond>();
		for(IBond bond : molClone.bonds()){
			bonds.add(bond);
		}
		for(IBond bond : bonds){
			ArrayList<IAtom> atoms = new ArrayList<IAtom>();
			for(IAtom atom : bond.atoms()){
					atoms.add(atom);
			}
			molClone.removeBond(bond);
			molClone.addBond(molClone.getAtomNumber(atoms.get(0)), molClone.getAtomNumber(atoms.get(1)), IBond.Order.SINGLE);
		}
		
		return molClone;
	}
}
