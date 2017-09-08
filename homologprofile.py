import sys
argv = sys.argv

# Read from command line
if len(argv) != 4:
	print('Insufficient number of args')
	sys.exit()

database = argv[1]
family = argv[2]
matrix = argv[3]

# Read family file
f = open(family, 'r')
line = f.readline()
family = []

while line:
	line = line.rstrip()
	family.append(line)
	line = f.readline()
f.close()

# Read database
f = open(database, 'r')
database = f.readline()
database = database.rstrip()
f.close()

# Read input files and create necessary data structures
scoremat = [[0 for i in range(20)] for i in range(20)]
f = open(matrix, 'r')

# Remove header
line = f.readline()

#Set matrix indices to each residue using dictionary
line = line.split()    
aa = {}
count = 0
for i in line:
	aa[i] = count
	count += 1

# Start reading each row of matrix values
line = f.readline()
line = line.split()[1:len(line)]
for i in range(len(scoremat)):
	for j in range(len(line)):
		scoremat[i][j + (i)] =  int(line[j]) # offset the values that are blank in front of where actual values start
	for k in range(i):
		scoremat[i][k] = scoremat[k][i] # fill in symmetrically
	line = f.readline()
	line = line.split()[1:len(line)]

f.close()

# Create frequency matrix
freq = [[0 for i in range(len(family[1]))] for i in range(20)]

for i in range(len(family[1])):
	total = 0
	for y in family:
		total += 1
		x = y[i]
		freq[aa[x]][i] += 1
	for y in range(20):
		freq[y][i] = float(freq[y][i]) / float(total)

# Threshold score
threshold = []
count = 0
for sequence in family:
	seqscore = 0
	for i in range(len(sequence)):
		for j in range(20):
			seqscore += freq[j][i]*scoremat[j][aa[sequence[i]]]
			
	threshold.append(seqscore)

totthres = 0
for i in threshold:
	totthres += i

totthres = totthres / len(threshold)

# Search for homologs
front = 0
back = len(family[1])
window = database[front:back]

sequence = {}

while back <= len(database):
	seqscore = 0
	for i in range(len(window)):
		for j in range(20):
			seqscore += freq[j][i]*scoremat[j][aa[window[i]]]
	sequence[window] = seqscore

	back += 1
	front += 1
	window = database[front:back]

average = 0
for i in sequence:
	average += sequence[i]

average = float(average) / float(len(sequence))
hlogs = {}
tot = 0

for i in sequence:
	for j in sequence:
		if sequence[j] >= sequence[i]:
			tot += 1
	pvalue = float(tot) / float(len(sequence))

	if pvalue <= .05:
		hlogs[i] = sequence[i]
	pvalue = 0
	tot = 0

f = open('homologs.txt', 'w')
for i in hlogs:
	f.write(i)
	f.write(' , ')
	f.write(str(hlogs[i]))
	f.write('\n')
f.close()




