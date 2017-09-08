import sys
argv = sys.argv

family = argv[1]
threshold = float(argv[2])

with open(family) as f:
	line = f.readline()
	subset = []
	thresvec = []

	line = line.rstrip()
	subset.append(line)
	seqlen = len(line)

	for line in f:
		for i in xrange(len(subset)):
			s1 = line.rstrip()
			print s1
			s2 = subset[i]
			print s2
			diff = [i for i in xrange(len(s1)) if s1[i] == s2[i]]
			#print len(diff)/seqlen
			difference = (float(len(diff))/seqlen)*100
			print difference

			if (difference < threshold):
				thresvec.append(True)
				print 'ADDED'
			else:
				thresvec.append(False)

		if (thresvec.count(False) == 0):
			subset.append(s1)
			
		thresvec = []
			# elif (i == len(subset)):
			# 	subset.append(line)
			# 	print 'ADDED TO SUBSET'

	print subset

