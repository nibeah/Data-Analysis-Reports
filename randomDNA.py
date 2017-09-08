import sys, random

#Problem 2
def randomDNA ():
	dnaoutput = open('dnaoutput.txt', 'r+')
	numseq = int(sys.argv[1])
	sizeseq = int(sys.argv[2])
	As = 0
	Cs = 0
	Gs = 0
	Ts = 0

	for i in range(0, numseq):
		dnaoutput.write(''.join(random.choice('CGTA') for _ in xrange(sizeseq)))
		dnaoutput.write('\n')

	with open('dnaoutput.txt') as f:
		for line in f:
			As = As + line.count('A')
			Cs = Cs + line.count('C')
			Gs = Gs + line.count('G')
			Ts = Ts + line.count('T')

	dnaoutput.write('Summary: ')
	#countnuc = dnaoutput.read()
	
	dnaoutput.write('# of A: ' + str(As))
	dnaoutput.write(' # of C: ' + str(Cs))
	dnaoutput.write(' # of G: ' + str(Gs))
	dnaoutput.write(' # of T: ' + str(Ts))



	return;

randomDNA()
