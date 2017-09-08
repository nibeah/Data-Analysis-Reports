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


# Find 'translation-start', 'donor', 'acceptor' of original humanseq

atg = [i for i in range(len(humanseq)) if humanseq.startswith('ATG', i)]
gt = [i for i in range(len(humanseq)) if humanseq.startswith('GT', i)]
ag = [i for i in range(len(humanseq)) if humanseq.startswith('AG', i)]

f = open('ATG.txt', 'w')
for i in atg:
	f.write(str(i) + '\n')
f.close()

f = open('GT.txt', 'w')
for i in gt:
	f.write(str(i) + '\n')
f.close()

f = open('AG.txt', 'w')
for i in ag:
	f.write(str(i) + '\n')
f.close()

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

# Viterbi optimum parse

donor = [1 if humanseq.startswith('GT', i) else -1000 for i in range(len(humanseq))]
acceptor = [1 if humanseq.startswith('AG', i+1) else -1000 for i in range(len(humanseq))]
ve = [0] * len(humanseq)
vi = [0] * len(humanseq)

diff = float(0)
k = 0

for i in range(5, len(humanseq)):
	maxve = -1000

	if i < 150:
		k = 0
	else:
		k = i-150

	for j in range (k, i-4):

		if acceptor[j-1] == -1000:
			continue

		else:
			subseq = humanseq[j:i+1]
			front = 0
			back = 6
			prexon = 0
			printron = 0
			frame = subseq[front:back]
			pentaframe = frame[0:5]

			while back <= len(subseq):
				if pentaframe in mchain:
					prexon += (float(mchain[pentaframe][frame])/float(mchain[pentaframe]['total']))

				if pentaframe in mchainnon:
					printron += (float(mchainnon[pentaframe][frame])/float(mchainnon[pentaframe]['total']))

				front += 1
				back += 1
				frame = subseq[front:back]
				pentaframe = frame[0:5]

			diff = prexon/printron
			codediff = math.log(diff, 2)

			vescore = codediff + vi[j-1] + 1

			if vescore > maxve:
				maxve = vescore

	if maxve == -1000:
		ve[i] = 0
	else:
		ve[i] = maxve

	maxvi = -1000

	for j in range(k, i-4):

		if donor[j] == -1000:
			continue

		else:
			subseq = humanseq[j:i+1]
			front = 0
			back = 6
			prexon = 0
			printron = 0
			frame = subseq[front:back]
			pentaframe = frame[0:5]

			while back <= len(subseq):
				if pentaframe in mchain:
					prexon += (float(mchain[pentaframe][frame])/float(mchain[pentaframe]['total']))

				if pentaframe in mchainnon:
					printron += (float(mchainnon[pentaframe][frame])/float(mchainnon[pentaframe]['total']))

				front += 1
				back += 1
				frame = subseq[front:back]
				pentaframe = frame[0:5]

			diff = prexon/printron
			codediff = (math.log(diff, 2)) * -1

			viscore = codediff + ve[j-1] + 1

			if viscore > maxvi:
				maxvi = viscore

	if maxvi == -1000:
		vi[i] = 0
	else:
		vi[i] = maxvi

 	# print i,'     ', ve[i], '      ', vi[i] 

f = open('geneGFF.txt', 'w')
exons = ''
exonranges = []
firstexon = 0
lastexon = 0

for x in range(len(vi)):
	if vi[x] >= ve[x]:

		if ve[x-1] >= vi[x-1]:
			if lastexon >= 2:
				exonranges.append([firstexon, firstexon + lastexon])
			firstexon = 0
			lastexon = 0
		print 'I'

	else:
		exons = exons + humanseq[x]

		if vi[x-1] >= ve[x-1]:
			firstexon = x
		else:
			lastexon += 1
		print 'E'

for i in exonranges:
	f.write('X\t' + 'Assignment 5\t' + 'Exonic Region (Gene)\t' + str(i[0]) + '\t' + str(i[1]) + '\t' + '.\t' + '+\t' + '.' + '\n')
f.close()



