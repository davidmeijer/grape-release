

f = open("../norine_peptides.txt")

lines = f.readlines()

monomer_lists = []

for line in lines:
	vals = line.split("\t")
	monomer_list = []
	for i in range(len(vals)):
		if i == 0:
			continue
		monomer_list.append(vals[i].strip())
	monomer_list.sort()
	monomer_lists.append(monomer_list)

total = len(monomer_lists)

unique_lists = []
repeated_lists = []

for i in range(len(monomer_lists)):
	if not (monomer_lists[i] in unique_lists):
		unique_lists.append(monomer_lists[i])
	else:
		repeated_lists.append(monomer_lists[i])

num_not_uniquely_identifiable = 0
for i in range(len(monomer_lists)):
	if monomer_lists[i] in repeated_lists:
		num_not_uniquely_identifiable = num_not_uniquely_identifiable + 1

print "Num uniquely identifiable: " + str(total - num_not_uniquely_identifiable)
print "Number of nrps: " + str(total)
print "Number of unique monomer lists: " + str(len(unique_lists))
