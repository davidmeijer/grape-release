package ca.mcmaster.magarveylab.grape.util.io;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.codehaus.jackson.map.ObjectMapper;
import org.openscience.cdk.exception.CDKException;

import ca.mcmaster.magarveylab.grape.GrapeConfig;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.AminoAcidEnums;
import ca.mcmaster.magarveylab.grape.enums.DomainEnums.PolyKetideDomainEnums;
import ca.mcmaster.magarveylab.grape.mol.chem.ChemicalAbstraction;
import ca.mcmaster.magarveylab.grape.mol.chem.Fragment;
import ca.mcmaster.magarveylab.grape.mol.chem.Fragment.FragmentType;
import ca.mcmaster.magarveylab.grape.util.ChemicalUtilities;

public class JsonOutput {
	
	private String name;
	private String smiles;
	private String id;
	List<Fragment> monomerFragments;
	ChemicalAbstraction abstraction;
	Map<String,Object> data;
	
	public JsonOutput(String name, String smiles, String id, ChemicalAbstraction abstraction) {
		this.name = name;
		this.smiles = smiles;
		this.id = id;
		this.monomerFragments = abstraction.getMonomerFragments();
		this.abstraction = abstraction;
	}
	
	public void writeJson(File file) throws IOException {
		ObjectMapper mapper = new ObjectMapper();
		mapper.writeValue(file, data);
	}
	
	public String getJSONString() throws IOException{
		ObjectMapper mapper = new ObjectMapper();
		return mapper.writeValueAsString(data);
	}
	
	public String getName(){
		return name;
	}
		
	public void createJson() throws IOException, CDKException {
		List<ArrayList<Fragment>> orderedMonomerFragmentSequences = ChemicalUtilities.getOrderedMonomerFragmentSequences(monomerFragments); //TODO: this should be somewhere else
		
		
		DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Date date = new Date();
		Map<String,Object> grapeResults = new HashMap<String,Object>();	

		Map<String,Object> grapeRun = new HashMap<String,Object>();	
		GrapeConfig gc = GrapeConfig.getConfig();
		grapeRun.put("id", id);
		grapeRun.put("name", name);
		grapeRun.put("SMILES", smiles);
		grapeRun.put("version", gc.getVersion());
		grapeRun.put("date", dateFormat.format(date));
		grapeRun.put("chemical_type", abstraction.getChemicalType());
		grapeRun.put("chemical_subtype", abstraction.getChemicalSubType());
		grapeRun.put("chemical_scaffold", abstraction.getChemicalScaffoldName());
		if(abstraction.getErrorMessage() != null){
			grapeRun.put("Error", abstraction.getErrorMessage());
		}
		
		List<Object> frags = new ArrayList<Object>();	
		
		for(ArrayList<Fragment> orderedSequence : orderedMonomerFragmentSequences) {
			List<Object> subs = new ArrayList<Object>();	

			for(Fragment m : orderedSequence) {

				Map<String,Object> sub = new HashMap<String,Object>();	
				Map<String, Object> modifications = new HashMap<String, Object>();
				
				if(m.getFragmentType() == FragmentType.AMINO_ACID) {
					String abbreviation = new String();
					if(m.getAminoAcidDomains().size() == 0) {
						abbreviation = "unknown";
					}
					else {
						abbreviation = m.getAminoAcidDomains().get(0).getAbbreviation();
					}
					String exactName = new String();
					if(m.getSpecificNames().size() == 0) {
						exactName = "unknown";
					}
					else {
						exactName = m.getSpecificNames().get(0);
					}
					sub.put("name", abbreviation);
					sub.put("exact_name", exactName);
					String smiles = SmilesIO.generateSmiles(m.getAtomContainer());
					sub.put("known_start", m.hasKnownStart());
					sub.put("smiles", smiles);
					sub.put("tanimoto_score", m.getTanimotoScore());
					sub.put("num_amide_atoms", m.getAminoCs().size() + m.getAminoNs().size());
					
				}else if(m.getFragmentType() == FragmentType.POLYKETIDE || m.getFragmentType() == FragmentType.FA_OR_PK) {
					for(int i = 0; i < m.getPkDomains().size(); i++) {
						
						String abbreviation = m.getLoadingUnits().get(i).getAbbreviation();
						if (abbreviation.equals("Unk")) {
							abbreviation = "unknown";
						}
						List<String> unitStates = new ArrayList<String>();
						boolean state = false;
						//TODO this method should just return an empty list instead of null.
						if(m.getPkDomains().get(i) != null){
							for(PolyKetideDomainEnums pkDomainEnum : m.getPkDomains().get(i)){
								if(pkDomainEnum != null) unitStates.add(pkDomainEnum.name());
							}
							state = true;
						}
						m.getFragmentType();
						Map<String,Object> subpk = new HashMap<String,Object>();	

						subpk.put("name", abbreviation);
						subpk.put("type", m.getFragmentType());
						String smiles = SmilesIO.generateSmiles(m.getAtomContainer());
						subpk.put("smiles", smiles);

						if(state){
							subpk.put("unit_states", unitStates);
						}
						subs.add(subpk);
					}
				
				}else if (m.getFragmentType().equals(FragmentType.MULTIPLE_AMINO_ACID_PIECE)){
					for (AminoAcidEnums a : m.getAminoAcidDomains()) {
						Map<String,Object> subaa = new HashMap<String,Object>();	

						String abbreviation = a.getAbbreviation();
						String smiles = SmilesIO.generateSmiles(m.getAtomContainer());
						subaa.put("name", abbreviation);
						subaa.put("type", m.getFragmentType());
						subaa.put("smiles", smiles);
						subaa.put("tanimoto_score", m.getTanimotoScore());
						subs.add(subaa);
					}
					
				}else if(m.getFragmentType().equals(FragmentType.ACYL_ADENYLATING)){
					String abbreviation = new String();
					abbreviation = m.getStarters().get(0).abbreviation();
					
					sub.put("name", abbreviation);
					String smiles = SmilesIO.generateSmiles(m.getAtomContainer());
					sub.put("smiles", smiles);
					sub.put("tanimoto_score", m.getTanimotoScore());
					
				}else if(m.getKnownOther() != null){
					sub.put("name", m.getKnownOther().getAbbreviation());
					sub.put("smiles", SmilesIO.generateSmiles(m.getAtomContainer()));
					sub.put("tanimoto_score", m.getTanimotoScore());
					sub.put("type", m.getFragmentType());
					
				}else{
					if (m.getFragmentType().equals(FragmentType.SUGAR)) {
						
					
						
						if (m.getSugarModifications().size() > 0) {
							modifications.put("sugar_tailors",m.getSugarModifications());
						}
						
						sub.put("sugar_names", m.getSugarNames());
						sub.put("tanimoto_score", m.getTanimotoScore());
					}
					String smiles = SmilesIO.generateSmiles(m.getAtomContainer());
					sub.put("smiles", smiles);
				}
				
				sub.put("type", m.getFragmentType());
				
				if (m.getTailoringDomains() != null && m.getTailoringDomains().size() > 0) {
					modifications.put("tailors", m.getTailoringDomains());
				}
				if(!m.getFragmentType().equals(FragmentType.ACYL_ADENYLATING) && m.getStarters().size() > 0) {
					modifications.put("acyl_adenylating_units", m.getStarters());
				}
				
				if (modifications.size() > 0) {
					sub.put("modifications", modifications);
				}else{
					sub.put("modifications", null);
				}
				
				subs.add(sub);
			}
			if (subs.size() != 0) {
				frags.add(subs);
			}
		}
		grapeRun.put("fragments", frags);
		grapeResults.put("grape_results", grapeRun);
		data = grapeResults;
	}

	public static void writeJSON(String baseDir, String name, String smiles, String id, ChemicalAbstraction chemicalAbstraction) {
		JsonOutput jo = new JsonOutput(name, smiles, id, chemicalAbstraction);
		new File(baseDir + "json/").mkdirs();
		name = SmilesIO.getCleanFileName(name);
		try {
			jo.createJson();
		} catch (IOException | CDKException e) {
			e.printStackTrace();
		}
		try {
			if(!name.isEmpty()){
				jo.writeJson(new File(baseDir + "json/" + name + ".json"));
			}else if(!id.isEmpty()){
				jo.writeJson(new File(baseDir + "json/" + id + ".json"));
			}else{
				jo.writeJson(new File(baseDir + "json/output.json"));
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}	
}
