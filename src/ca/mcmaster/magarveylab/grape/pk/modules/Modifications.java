package ca.mcmaster.magarveylab.grape.pk.modules;

import java.util.HashMap;

public class Modifications {
	String macrolideType;
	Integer sugars = 0;
	Integer deoxySugars = 0;
	Integer epoxides = 0;
	Integer chlorines = 0;
	
	public Modifications(){	
	}
	public void setMacrolideType(String type){
		macrolideType = type;
	}
	public void addSugar(){
		sugars ++;
	}
	public void addDeoxySugar(){
		deoxySugars ++;
	}
	public void addEpoxide(){
		epoxides ++;
	}
	public void addChlorineCount(int numCl){
		chlorines = numCl;
	}
	public String getMacrolideType(){
		return(macrolideType);
	}
	public Integer getSugarCount(){
		return(sugars);
	}
	public Integer getDeoxySugarCount(){
		return(deoxySugars);
	}
	public Integer getEpoxideCount(){
		return(epoxides);
	}
	public Integer getChlorines(){
		return(chlorines);
	}
	public HashMap<String, String> getAllMods(){
		HashMap<String, String> allMods = new HashMap<String,String>();
		allMods.put("Macrolide Type:", macrolideType);
		if (chlorines > 0){
			allMods.put("Number of chlorine(s)", "" + chlorines); //toString better here?
		}
		if(sugars > 0){
			allMods.put("Number of sugar(s)", "" + sugars);
			}
		if(deoxySugars > 0){
			allMods.put("Number of deoxysugar(s)", "" + deoxySugars);
			}
		if(epoxides > 0){
			allMods.put("Number of epoxide(s)", "" + epoxides);
		}
		return(allMods);
	}
}
