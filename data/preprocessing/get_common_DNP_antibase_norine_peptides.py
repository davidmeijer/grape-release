import re

f_dnp = open("Combined_DNP_and_Antibase_SMILES_Database_-_Repeated_compounds_gone_obvious_non-peptides_removed_lipids.txt", "r")
f_norine = open("../norine_peptides.txt", "r")
f_out = open("dnp_antibase_norine_nrp_smiles_monomers.txt", "w")

dnp_lines = f_dnp.readlines()
norine_lines = f_norine.readlines()

dnp_data = [] # This will be a list of tuples, which each tuple being  (list of names, SMILES string)

firstline = True
for line in dnp_lines:
	if firstline == True:
		firstline = False
		continue
	vals = line.split("\t")
	names = vals[0].split(";")
	for i in range(len(names)):
		names[i] = names[i].strip()
	smiles = vals[2].strip()
	dnp_data.append( (names, smiles) )

norine_data = [] # This will contain a list of tuples, each tuple being (name, monomer_string)
for line in norine_lines:
	name = re.search("^(.*?)\t", line).group(1)
	monomer_text = re.search("^.*?\t(.*)", line).group(1)
	# Probably not necessary, but for good measure.
	name = name.strip()
	monomer_text = monomer_text.strip()
	norine_data.append( (name, monomer_text) )

for (dnp_names, smiles) in dnp_data:
	for (norine_name, monomer_text) in norine_data:
		for dnp_name in dnp_names:
			name1 = re.sub("[^0-1A-Z]", "", dnp_name.upper())
			name2 = re.sub("[^0-9A-Z]", "", norine_name.upper())
			if name1 == name2:
				f_out.write(norine_name + "\t" + smiles + "\t" + monomer_text + "\n")
				break

f_dnp.close()
f_norine.close()
f_out.close()
