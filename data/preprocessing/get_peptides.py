import re
import urllib2

f = open("norine_peptide_urls.txt", "r")
out = open("norine_peptides.txt", "w")

urls = f.readlines()

monomer_data = []

counter = 1
for u in urls:
	print("Getting data from url number " + str(counter))
	response = urllib2.urlopen(u)
	text = response.read()
	
	# Regex that grabs the first thing between two tags after the text "complete name" in html
	name_match = re.search("<h3>(.*?)<", text, re.DOTALL)
	name = name_match.group(1)
	name = name.strip()
	
	# Regex for grabbing monomer short forms
	monomers = []
	text_match = re.search("Monomeric composition.*?<table>((.*?<table>.*?</table>.*?)*)?</table>", text, re.DOTALL)
	text2 = text_match.group(1);
	
	link_texts = re.findall("<a href.*?>.*?</a>", text2, re.DOTALL)
	
	for link_text in link_texts:
		link_match = re.match("<a href.*?>(.*?)</a>", link_text, re.DOTALL)
		mon_name = link_match.group(1).strip()
		monomers.append(mon_name)

	monomer_data.append( (name, monomers) )

	counter = counter + 1

for m in monomer_data:
	name, monomers = m
	out.write(name)
	for mon_name in monomers:
		out.write("\t" + mon_name)
	out.write("\n")

f.close()
out.close()
