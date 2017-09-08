import sys

#Problem 1

def locAL():
	in_arr = sys.argv
	locALoutput = open('locALoutput.txt', 'w')
	pair_array = []
	scores = [0, 0, 0]
	


	if '-m' not in in_arr or '-s' not in in_arr or '-d' not in in_arr:
		raise NameError('Error: inputs are incorrect')
	else:
		filename = in_arr[1]
		match = int(in_arr[in_arr.index('-m') + 1])
		mismatch = int(in_arr[in_arr.index('-s') + 1])
		indel = int(in_arr[in_arr.index('-d') + 1])

	with open(filename, 'r') as fasta:
		for line in fasta:
			if line.startswith('A') or line.startswith('C') or line.startswith('T') or line.startswith('G'):
				#pair_array.append(line)
				pair_array.append(line.replace('\r\n', ''))
				#pair_array.append(line.replace('\n', ''))

	len1 = len(pair_array[0])
	len2 = len(pair_array[1])
	
	almat = [[0 for x in range(len(pair_array[0])+1)] for x in range(len(pair_array[1])+1)]
	almatptr = [[0 for x in range(len(pair_array[0])+1)] for x in range(len(pair_array[1])+1)]

	#y values
	for num in range(1, len(pair_array[1])):
		almat[num][0] = indel + almat[num-1][0]

	#x values
	for num2 in range(1, len(pair_array[0])):
		almat[0][num2] = indel + almat[0][num2-1]

	for i in range(1, len(pair_array[1])+1):
		for j in range(1, len(pair_array[0])+1):
			if pair_array[1][i-1] == pair_array[0][j-1]:
				scores[0] = almat[i-1][j-1] + match
			else:
				scores[0] = almat[i-1][j-1] + mismatch

			scores[1] = almat[i][j-1] + indel
			scores[2] = almat[i-1][j] + indel

			almat[i][j] = max(scores)

			if scores.index(max(scores)) == 0:
				almatptr[i][j] = 'diagonal'
			elif scores.index(max(scores)) == 1:
				almatptr[i][j] = 'left'
			elif scores.index(max(scores)) == 2:
				almatptr[i][j] = 'right'

	# print almat
	# print almatptr

	i1 = len(pair_array[1])
	j1 = len(pair_array[0])
	seq1 = ''
	seq2 = ''
	allen = 0

	while i1 > 1 and j1 > 1:
		if almatptr[i1][j1] == 'diagonal':
			seq1 = pair_array[1][-1] + seq1
			seq2 = pair_array[0][-1] + seq2
			pair_array[1] = pair_array[1][:-1]
			pair_array[0] = pair_array[0][:-1]
			i1 = i1-1
			j1 = j1-1
			allen = allen + 1
		elif almatptr[i1][j1] == 'left':
			seq2 = pair_array[0][-1] + seq2
			pair_array[0] = pair_array[0][:-1]
			seq1 = '-' + seq1
			j1 = j1-1
			allen = allen + 1
		elif almatptr[i1][j1] == 'right':
			seq2 = '-' + seq2
			seq1 = pair_array[1][-1] + seq1
			pair_array[1] = pair_array[1][:-1]
			i1 = i1-1
			allen = allen + 1

	if j1 > 1:
		while j1 >= 1:
			seq2 = pair_array[0][-1] + seq2
			pair_array[0] = pair_array[0][:-1]
			seq1 = '-' + seq1
			j1 = j1-1
			allen = allen + 1
	elif i1 > 1:
		while i1 >= 1:
			seq2 = '-' + seq2
			seq1 = pair_array[1][-1] + seq1
			pair_array[1] = pair_array[1][:-1]
			i1 = i1-1
			allen = allen + 1

	seq1 = pair_array[1][-1] + seq1
	seq2 = pair_array[0][-1] + seq2
	allen = allen + 1

	locALoutput.write('Score ' + str(almat[len(pair_array[1])][len(pair_array[0])]) + '\n')
	locALoutput.write('Length ' + str(allen) + '\n')
	locALoutput.write('%s\n%s' % (seq1, seq2))


	return ;

locAL()





