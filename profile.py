import sys
argv = sys.argv

#~~~~~~LINK A.A. SYMBOLS TO MATRIX POS (0, 1, 2,...,19)~~~~~~





#~~~~~~ CHECKS VALIDITY OF ARGUMENTS~~~~~~~~~~~~~~~~~~~~~~~
# 		Reads in filenames...parses information
print('checking input arguments...')
if len(argv) != 4:
	print 'insufficient arguments'
	print '\'python [FILENAME.py] [DATABASE] [FAMILY] [MATRIX]\''
	sys.exit()
#creation of database, family, 
db = argv[1]
fam = argv[2]
mat = argv[3]

#~~~~~~PARSE DATABASE~~~~~~~~~~#
f = open(db, 'r')
db = f.readline()
db = db.rstrip()
f.close()

#~~~~~~PARSE FAMILY~~~~~~~~~~~~#
f = open(fam, 'r')
line = f.readline()
fam = []
while line:
	line = line.rstrip()
	fam.append(line)
	line = f.readline()
f.close()


#~~~~~~~~PARSE MATRIX~~~~~~~~~~~~#
score = [[0 for x in range(20)] for x in range(20)] #20 amino acids!
f = open(mat, 'r')


line = f.readline()     #remove column labels [A, C, D,... ]
line = line.split()     #Take column labels and link A.A. symbol to pos in matrix
amino = {}
count = 0
for x in line:
	amino[x] = count
	count += 1


line = f.readline()

line = line.split()[1:len(line)] #remove row labels [A, C, D,... ]

#inputs values from mat.txt in score matrix
for y in range(len(score)): #for every row...
	for x in range(len(line)):
		score[y][x + (y)] =  int(line[x]) #fill in column
	for z in range(y):
		score[y][z] = score[z][y] #this here accounts for symmetric grid spaces left blank in .txt file...
	line = f.readline()
	line = line.split()[1:len(line)]

f.close()
#~~~~~~~~~~~~~~GENERATE FREQUENCY MATRIX~~~~~~~~~~~~~~~~~
freq = [[0 for x in range(len(fam[1]))] for x in range(20)]

for x in range(len(fam[1])):
	total = 0
	for y in fam:
		total += 1
		aa = y[x]
		freq[amino[aa]][x] += 1
	# ^^^ ~ count freq of a.a. for position x in alignment in each protein seq (y)
	# vvv ~ produces freqency value by dividing occurance of a.a. at pos x by total number of protein seqs 
	for y in range(20):
		freq[y][x] = float(freq[y][x]) / float(total)

# for x in range(20):
# 	print freq[x]


#~~~~~~~~~~~GENERATE TRESHOLD SCORE~~~~~~~~~~~~~~~~~~~~~
threshold = [] #store scores of each sequence in family
count = 0
for seq in fam:
	profScore = 0
	for aa in range(len(seq)):
		for x in range(20):
			profScore += freq[x][aa]*score[x][amino[seq[aa]]]
			
	threshold.append(profScore)
#generate average threshold score from the scores of each seq in family F
THRESHOLD = 0
for x in threshold:
	THRESHOLD += x

THRESHOLD = THRESHOLD / len(threshold)
print 'Average Profile Score of family F sequences: ', THRESHOLD







#~~~~~~~SLIDING FRAME: FIND SEQS IN DATABASE~~~~~~~~~~~~~~~
back = 0
front = len(fam[1])
window = db[back:front]

seq = {}

while front <= len(db):
	profScore = 0
	for aa in range(len(window)):
		for x in range(20):
			profScore += freq[x][aa]*score[x][amino[window[aa]]]
	seq[window] = profScore

	# if profScore >= THRESHOLD:
	# 	print "MATCH"
	#homologs.append(profScore)

	# if profScore > -20:
	# 	print profScore

	front += 1
	back += 1
	window = db[back:front]

average = 0
for x in seq:
	average += seq[x]

average = float(average) / float(len(seq))
print 'Average Profile Score of database Db sequences: ', average

f = open('scores.txt', 'w')
for x in seq:
	f.write(str(seq[x]))
	f.write('\n')
f.close()



homologs = {}
total = 0
for x in seq:
	for y in seq:
		if seq[y] >= seq[x]:
			total += 1
	pvalue = float(total) / float(len(seq))

	if pvalue <= .05:
		homologs[x] = seq[x]
	pvalue = 0
	total = 0

f = open('homologs.txt', 'w')
for x in homologs:
	f.write(x)
	f.write(' : ')
	f.write(str(homologs[x]))
	f.write('\n')
f.close()



