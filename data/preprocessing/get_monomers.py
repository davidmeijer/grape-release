import re
import urllib2

f = open("norine_monomer_urls.txt", "r")
out = open("nrp_monomers.txt", "w")

urls = f.readlines()

monomer_data = []

counter = 1
for u in urls:
	print("Getting data from url number " + str(counter))
	response = urllib2.urlopen(u)
	text = response.read()
	# Regex that grabs the first thing between two tags after the text "complete name" in html
	complete_name_match = re.search("complete name[^>]*>([^<]*)<", text)
	complete_name = complete_name_match.group(1)
	complete_name = complete_name.strip()
	# Regex that grabs the short name
	short_name_match = re.search("<h2>([^<]*)</h2>", text)
	short_name = short_name_match.group(1)
	short_name = short_name.strip()
	
	# Regex for grabbing Smiles from html
	smiles_match = re.search("isomeric SMILES[^>]*>([^<]*)<", text)
	if(smiles_match == None):
		smiles_match = re.search("SMILES[^>]*>([^<]*)<", text)
	smiles = smiles_match.group(1)
	
	smiles = smiles.strip()
	
	monomer_data.append( (short_name, complete_name, smiles) )

	counter = counter + 1

for m in monomer_data:
	short_name, complete_name, smiles = m
	out.write(short_name + "\t" + complete_name + "\t" + smiles + "\n")

f.close()
out.close()
