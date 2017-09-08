import sys, math, random
argv = sys.argv

codingfile = argv[1]
noncodingfile = argv[2]
humseqfile = argv[3]

# Parse Fastas
seqs = []

with open(codingfile) as fasta:
	for line in fasta:
		if line.startswith('>'):
			continue
		else:
			line = line.replace('N', random.choice('ACTG'))
			seqs.append(line.rstrip())

nonseq = ''

with open(noncodingfile) as fasta:
	for line in fasta:
		line = line.rstrip()
		line = line.replace(' ', '')

		if line.startswith('>'):
			continue

		if len(line) == 0:
			continue
		else:
			nonseq = nonseq + line

humanseq = ''

with open(humseqfile) as fasta:
	for line in fasta:
		if line.startswith('>'):
			continue
		else:
			humanseq = humanseq + line.rstrip()

# Create Markov Indices
mchain = {}

for seq in seqs:
	front = 0
	back = 6
	frame = seq[front:back]
	pentaframe = frame[0:5]

	while back <= len(seq):
		if pentaframe in mchain:
			mchain[pentaframe][frame] += 1
			mchain[pentaframe]['total'] += 1
		else:
			mchain[pentaframe] = {pentaframe+'A': 0, pentaframe+'C': 0, pentaframe+'T': 0, pentaframe+'G': 0, 'total':0}
			mchain[pentaframe][frame] += 1
			mchain[pentaframe]['total'] += 1

		front += 1
		back += 1
		frame = seq[front:back]
		pentaframe = frame[0:5]


mchainnon = {}

counter = 0
front = 0
back = 6
frame = seq[front:back]
pentaframe = frame[0:5]

while back <= len(nonseq):
	if pentaframe in mchainnon:
		mchainnon[pentaframe][frame] += 1
		mchainnon[pentaframe]['total'] += 1
	else:
		mchainnon[pentaframe] = {pentaframe+'A': 0, pentaframe+'C': 0, pentaframe+'T': 0, pentaframe+'G': 0, 'total':0}
		mchainnon[pentaframe][frame] += 1
		mchainnon[pentaframe]['total'] += 1

	front += 1
	back += 1
	frame = nonseq[front:back]
	pentaframe = frame[0:5]

# Select 100 random sequences

randseq = []
randnonseq = []

while len(randseq) < 100:
	randseq.append(random.choice(seqs))

front = 0
back = 100
rand = nonseq[front:back]

while len(randnonseq) < 100:
	randnonseq.append(rand)
	front += 100
	back += 100
	rand = nonseq[front:back]

# Coding Differential Calculations

diff = float(0)

f = open('codediff.txt', 'w')
for s in randseq:
	front = 0
	back = 6
	prexon = 0
	printron = 0
	frame = s[front:back]
	pentaframe = frame[0:5]

	while back <= len(s):
		if pentaframe in mchain:
			prexon += (float(mchain[pentaframe][frame])/float(mchain[pentaframe]['total']))

		if pentaframe in mchainnon:
			printron += (float(mchainnon[pentaframe][frame])/float(mchainnon[pentaframe]['total']))

		front += 1
		back += 1
		frame = s[front:back]
		pentaframe = frame[0:5]

	diff = prexon/printron
	codediff = math.log(diff, 2)

	f.write(str(codediff))
	f.write('\n')

for s in randnonseq:
	front = 0
	back = 6
	prexon = 0
	printron = 0
	frame = s[front:back]
	pentaframe = frame[0:5]

	while back <= len(s):
		if pentaframe in mchain:
			prexon += (float(mchain[pentaframe][frame])/float(mchain[pentaframe]['total']))

		if pentaframe in mchainnon:
			printron += (float(mchainnon[pentaframe][frame])/float(mchainnon[pentaframe]['total']))

		front += 1
		back += 1
		frame = s[front:back]
		pentaframe = frame[0:5]

	diff = prexon/printron
	codediff = math.log(diff, 2)
	f.write(str(codediff))
	f.write('\n')

f.close()

# GC Content Calculations

f = open('GC.txt', 'w')
for s in randseq:
	counter = 0
	GCcount = 0
	for base in s:
		if base == 'G' or base == 'C':
			counter += 1

	GCcount = float(counter)/float(len(s))
	f.write(str(GCcount))
	f.write('\n')

for s in randnonseq:
	counter = 0
	GCcount = 0
	for base in s:
		if base == 'G' or base == 'C':
			counter += 1

	GCcount = float(counter)/float(len(s))
	f.write(str(GCcount))
	f.write('\n')

f.close()





