package ca.mcmaster.magarveylab.grape.pk.modules;

public enum PKsubstrate {
	MALONYL("Mal"),
	METHYLMALONYL("MeM"),
	METHOXYLMALONYL("OMeMal"),
	BENZOYL("Bz"),
	ETHYLMALONYL("EtM"),
	ISOBUTRYL("IBu"),
	METHYLBUTERYL2("MBu"),
	UNKNOWN("Unk");
	
	private final String abbreviation;
	
	private PKsubstrate(String abbreviation){
		this.abbreviation = abbreviation;
	}
	
	public String getAbbreviation(){
		return abbreviation;
	}
	
}
