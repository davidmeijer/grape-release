import re
import urllib2

f = open("NRP_DB.txt", "r")
out = open("nrp_from_isnap.txt", "w")

nrp_data = []

for line in f:
	vals = line.split(" = ")
	if len(vals) < 2:
		continue
	first, second = vals[0], vals[1]
	first = first.strip()
	second = second.strip()
	first = re.sub("^#define \$", "", first)
	first = re.sub("^;#define \$", "", first)
	out.write(first + "\t" + second + "\n")

f.close()
out.close()
