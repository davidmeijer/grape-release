
package ca.mcmaster.magarveylab.grape.nrp.chem;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Scanner;

import ca.mcmaster.magarveylab.grape.util.MonomerMatch;

/**
 * Get the domain representation of a monomer match
 * @author gmchen
 */
public class DomainMatcher 
{
	private HashMap<String, DomainUnit> domainMap;
	public DomainMatcher(String domainMapFilename) {
		initializeDomainMap(domainMapFilename);
	}
	
	public static boolean isMatch(ArrayList<ArrayList<String>> predictedDomains, ArrayList<String> actualDomains) {
		boolean withUnknowns = false;
		for(int i = 0; i < predictedDomains.size(); i++) {
			for(int j = 0; j < predictedDomains.get(i).size(); j++) {
				if("UNKNOWN".equals(predictedDomains.get(i).get(j))) {
					withUnknowns = true;
				}
			}
		}
		if(withUnknowns) {
			return isMatchWithUnknowns(actualDomains, predictedDomains, 0, new boolean[actualDomains.size()]);
		}
		return isMatchNoUnknowns(actualDomains, predictedDomains, 0, new boolean[predictedDomains.size()]);
	}
	
	public static boolean isMatchNoUnknowns(ArrayList<String> actualDomains, ArrayList<ArrayList<String>> predictedDomains, int currentActualDomainIndex, boolean[] predictedDomainsVisited) {
		if(actualDomains.size() != predictedDomains.size()) {
			return false;
		}
		
		// The exit condition of the recursion if it is a match
		if(currentActualDomainIndex >= actualDomains.size()) {
			return true;
		}
		
		for(int i = 0; i < predictedDomains.size(); i++) {
			if(predictedDomainsVisited[i] == false) {
				// check if the current actual domain matches one of the predicted domains of index i 
				for(int j = 0; j < predictedDomains.get(i).size(); j++) {
					if(actualDomains.get(currentActualDomainIndex).equals(predictedDomains.get(i).get(j))) {
						boolean[] newPredictedDomainsVisited = predictedDomainsVisited.clone();
						newPredictedDomainsVisited[i] = true;
						if( isMatchNoUnknowns(actualDomains, predictedDomains, currentActualDomainIndex + 1, newPredictedDomainsVisited) ) {
							return true;
						}
					}
				}
			}
		}
		
		return false;
	}
	/**
	 * Check if this is a match, possibly with unknowns
	 * @param actualDomains
	 * @param predictedDomains
	 * @param currentPredictedDomainIndex
	 * @param actualDomainsVisited
	 * @return
	 */
	public static boolean isMatchWithUnknowns(ArrayList<String> actualDomains, ArrayList<ArrayList<String>> predictedDomains, int currentPredictedDomainIndex, boolean[] actualDomainsVisited) {
		//All of the predicted domains must match to an actual domain. Since there are unknowns, there may be actual domains which are unmatched to a predicted domain.
		
		if(predictedDomains.size() > actualDomains.size()) {
			return false;
		}
		// The exit condition of the recursion if it is a match
		if(currentPredictedDomainIndex >= predictedDomains.size()) {
			return true;
		}
		
		for(int i = 0; i < actualDomains.size(); i++) {
			if(actualDomainsVisited[i] == false) {
				// check if the current actual domain matches one of the predicted domains of index currentPredictedDomainIndex
				for(int j = 0; j < predictedDomains.get(currentPredictedDomainIndex).size(); j++) {
					if(predictedDomains.get(currentPredictedDomainIndex).get(j).equals("UNKNOWN")) {
						continue;
					}
					if(actualDomains.get(i).equals(predictedDomains.get(currentPredictedDomainIndex).get(j))) {
						boolean[] newActualDomainsVisited = actualDomainsVisited.clone();
						newActualDomainsVisited[i] = true;
						if( isMatchWithUnknowns(actualDomains, predictedDomains, currentPredictedDomainIndex + 1, newActualDomainsVisited) ) {
							return true;
						}
					}
				}
			}
		}
		
		return false;
	}
	/**
	 * Create domain hashmap from file
	 * @param filename
	 */
	private void initializeDomainMap(String filename) {
		domainMap = new HashMap<String, DomainUnit>();
		Scanner scanner = null;
		try {	
			scanner = new Scanner(new File(filename));
		} catch (Exception e) {
			e.printStackTrace();
		}
		while(scanner.hasNextLine()) {
			String line = scanner.nextLine();
			String[] vals = line.split("\t");
			ArrayList<String> domainNames = new ArrayList<String>(Arrays.asList(Arrays.copyOfRange(vals, 2, vals.length)));
			domainMap.put(vals[0], new DomainUnit(domainNames));
		}
	}

	public ArrayList<DomainUnit> getDomains(MonomerMatch monomerMatch) {
		ArrayList<DomainUnit> domainUnits = new ArrayList<DomainUnit>();
		for(int i = 0; i < monomerMatch.getMonomerNames().size(); i++) {
			DomainUnit domainUnit = domainMap.get(monomerMatch.getMonomerNames().get(i));
			
			if(monomerMatch.getTanimotoScore() < 0.9) {
				// call the copy constructor
				domainUnit = new DomainUnit(domainUnit);
				domainUnit.getDomainNames().add("lowTanimoto");
			}
			
			boolean isNew = true;
			for(int j = 0; j < domainUnits.size(); j++) {
				if(domainUnits.get(j).equals(domainUnit)) {
					isNew = false;
					break;
				}
			}
			
			if(isNew) {
				domainUnits.add(domainUnit);
			}
		}
		return domainUnits;
	}
	/*
	private ArrayList<ArrayList<ArrayList<String>>> getDomains(ArrayList<ArrayList<String>> predictedMonomers) {
		ArrayList<ArrayList<ArrayList<String>>> domainsListListList = new ArrayList<ArrayList<ArrayList<String>>>();
		
		for(int i = 0; i < predictedMonomers.size(); i++) {
			ArrayList<ArrayList<String>> domainsListList = new ArrayList<ArrayList<String>>();
			for(int j = 0; j < predictedMonomers.get(i).size(); j++) {
				ArrayList<String> domainsList = domainMap.get(predictedMonomers.get(i).get(j));
				// If this domain list is the same as another domain list, then don't add it.
				boolean alreadyPresent = false;
				Collections.sort(domainsList);
				for(ArrayList<String> currentList:domainsListList) {
					
					if(domainsList.equals(currentList)) {
						alreadyPresent = true;
						break;
					}
				}
				// add to the domains.
				if(!alreadyPresent) {
					domainsListList.add(domainsList);
				}
			}
			domainsListListList.add(domainsListList);
		}
		
		return domainsListListList;
	}
	
	public ArrayList<ArrayList<ArrayList<String>>> getDomainsFromMatches(ArrayList<MonomerMatch> monomerMatches) {
		ArrayList<ArrayList<String>> predictedMonomers = new ArrayList<ArrayList<String>>();
		for(MonomerMatch m:monomerMatches) {
			predictedMonomers.add(m.monomerNames);
		}
		return getDomains(predictedMonomers);
	}
	*/
}
