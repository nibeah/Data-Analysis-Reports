import sys

# Haebin Liew A10639579

# This program takes an input of a dictionary file and database file to
# create an Aho Corasick trie and performs a search algorithm to search
# the database against the dictionary, resulting in the number of matches
# each keyword has to the dictionary.

if len(sys.argv) != 3:
	print ('Incorrect number of arguments')
	sys.exit()

file = open(sys.argv[1], 'r')
dictionary = []

results = dict()
line = file.readline()
longest = 0

while line:
	line = line.rstrip()

	if len(line) > longest:
		longest = len(line)

	dictionary.append(line)
	results.setdefault(line, 0)
	line = file.readline()

file.close()

root = dict()
curr = 0
stop = '_stop_'

for word in dictionary:
	curr = root

	for letter in word:
		curr = curr.setdefault(letter, {})

	curr[stop] = stop


dbfile = open(sys.argv[2], 'r')
db = dbfile.readline()
db.rstrip()
beg = 0
end = longest
frame = db[beg:end]
readframe = ''

while end != len(db):
	curr = root

	for letter in frame:
		if letter in curr:
			readframe = readframe + letter
			curr = curr[letter]

			if stop in curr:
				results[readframe] += 1

		else:
			frame = ''

	beg += 1
	end += 1
	frame = db[beg:end]
	readframe = ''

for string in results:
	print string, ', ', results[string]







