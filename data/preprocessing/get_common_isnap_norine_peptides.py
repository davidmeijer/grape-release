import re

f_isnap = open("../nrp_from_isnap.txt", "r")
f_norine = open("../norine_peptides.txt", "r")
f_common = open("common_compounds.txt", "r")
f_out = open("isnap_norine_nrp.txt", "w")

isnap_lines = f_isnap.readlines()
norine_lines = f_norine.readlines()
common_lines = f_common.readlines()

for line in common_lines:
	line = line.strip()
	outline = ""
	norine_line = ""
	isnap_line = ""
	for l in norine_lines:		
		l = l.strip()
		if(line == re.sub("[^0-9A-Za-z]", "", l.split("\t")[0])):
			norine_line = l
			break
	for l in isnap_lines:
		l = l.strip()
		if(line == re.sub("[^0-9A-Za-z]", "", l.split("\t")[0])):
			isnap_line = l
			break
	# note that common_compounds may list compounds whose names are duplicated in one file (as opposed to being present in both files), so we'll skip this if it wasn't found in one file
	if norine_line == "" or isnap_line == "":
		continue
	outline = isnap_line
	outline = outline + "\t" + norine_line
	f_out.write(outline + "\n")

f_out.close()
f_isnap.close()
f_norine.close()
f_common.close()
