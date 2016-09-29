package ca.mcmaster.magarveylab.grape.nrp.chem;

import java.util.ArrayList;

/**
 * For example, [C][A][T] is considered one DomainUnit.
 * This class gives a representation of a biosynthetic assembly line domain unit, which may consist of multiple actual assembly lines in series.
 *  * @author gmchen
 */
public class DomainUnit {

	private ArrayList<String> domainNames;
	
	public DomainUnit(ArrayList<String> domainNames) {
		this.domainNames = domainNames;
	}
	
	public DomainUnit(DomainUnit unitToCopy) {
		this.domainNames = new ArrayList<String>();
		for(String s:unitToCopy.getDomainNames()) {
			domainNames.add(s);
		}
	}

	public ArrayList<String> getDomainNames() {
		return domainNames;
	}

	public void setDomains(ArrayList<String> domainNames) {
		this.domainNames = domainNames;
	}
	
	@Override
	public boolean equals(Object other) {
		DomainUnit otherDomainUnit = (DomainUnit) other;
		
		if(domainNames.size() != otherDomainUnit.getDomainNames().size()) {
			return false;
		}
		
		for(int i = 0; i < domainNames.size(); i++) {
			if(!domainNames.get(i).equals(otherDomainUnit.getDomainNames().get(i))) {
				return false;
			}
		}
		
		return true;
	}
}
