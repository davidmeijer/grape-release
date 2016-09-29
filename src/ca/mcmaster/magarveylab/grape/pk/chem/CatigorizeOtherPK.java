package ca.mcmaster.magarveylab.grape.pk.chem;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Map;
import java.util.Map.Entry;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.isomorphism.UniversalIsomorphismTester;

import ca.mcmaster.magarveylab.grape.enums.MoleculeClasses.AromaticFungal;
import ca.mcmaster.magarveylab.grape.enums.MoleculeClasses.AromaticNotType2;
import ca.mcmaster.magarveylab.grape.enums.MoleculeClasses.Enedyine;
import ca.mcmaster.magarveylab.grape.enums.MoleculeClasses.ChemicalSubType;
import ca.mcmaster.magarveylab.grape.enums.MoleculeClasses.Terpenes;
import ca.mcmaster.magarveylab.grape.enums.MoleculeClasses.AromaticType2;
import ca.mcmaster.magarveylab.grape.util.io.SmilesIO;

public class CatigorizeOtherPK {
	
	private static UniversalIsomorphismTester uit = new UniversalIsomorphismTester();
	private static Map<String, IAtomContainer> otherAromaticMolecules = AromaticNotType2.getAll();
	private static Map<String, IAtomContainer> type2Molecules = AromaticType2.getAll();
	private static Map<String, IAtomContainer> fungalAromaticMolecules = AromaticFungal.getAll();
	private static Map<String, IAtomContainer> terpeneMolecules = Terpenes.getAll();
	private static Map<String, IAtomContainer[]> enedyineMolecules = Enedyine.getAll();
	private ChemicalSubType type;
	private String match;
	
	public void analyzeMolecule(IAtomContainer mol, int macrolideType) throws CDKException{
		type = ChemicalSubType.STANDARD;
		match = null;
		IAtomContainer carotenoids = null;
		try{
			carotenoids = SmilesIO.readSmilesTemplates("CC(C)=C\\C=C\\C=C(/C)\\C=C\\C=C(C)C"); //specific terpene
		}catch(IOException e){}
	    if(uit.isSubgraph(mol, carotenoids)){
    		type = ChemicalSubType.TERPENE;
    		match = "Carotenoid";
    		return;
	    }
	    for(Entry<String, IAtomContainer[]> entry : enedyineMolecules.entrySet()){
	    	String name = entry.getKey();
	    	IAtomContainer[] scaffolds = entry.getValue();
	    	for(IAtomContainer scaffold : scaffolds){
		    	if(uit.isSubgraph(mol, scaffold)){
		    		type = ChemicalSubType.ENEDYINE;
		    		match = name;
		    		return;
		    	}
	    	}
	    }
	    if(macrolideType == 0){
	    	for(Entry<String, IAtomContainer> entry : fungalAromaticMolecules.entrySet()){
		    	String name = entry.getKey();
		    	IAtomContainer aro = entry.getValue();
		    	if(uit.isSubgraph(mol, aro)){
		    		type = ChemicalSubType.FUNGAL;
		    		match = name;
		    		return;
		    	}
		    }
			IAtomContainer molSingleBondsOnly = onlySingleBonds(mol);
		    for(Entry<String, IAtomContainer> entry : type2Molecules.entrySet()){
		    	String name = entry.getKey();
		    	IAtomContainer aro = entry.getValue();
		    	if(uit.isSubgraph(molSingleBondsOnly, aro)){
		    		type = ChemicalSubType.TYPE_2;
		    		match = name;
		    		return;
		    	}
		    }
		    for(Entry<String, IAtomContainer> entry : otherAromaticMolecules.entrySet()){
		    	String name = entry.getKey();
		    	IAtomContainer aro = entry.getValue();
		    	if(uit.isSubgraph(molSingleBondsOnly, aro)){
		    		type = ChemicalSubType.NON_TYPE_2_AROMATIC;
		    		match = name;
		    		return;
		    	}
		    }
		    for(Entry<String, IAtomContainer> entry : terpeneMolecules.entrySet()){
		    	String name = entry.getKey();
		    	IAtomContainer aro = entry.getValue();
		    	if(uit.isSubgraph(molSingleBondsOnly, aro)){
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
	private static IAtomContainer onlySingleBonds(IAtomContainer mol) {
		IAtomContainer molClone = null;
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
