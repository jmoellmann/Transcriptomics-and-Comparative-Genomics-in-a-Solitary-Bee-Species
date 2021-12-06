import sys
import pysam

file_path = sys.argv[1]

samfile = pysam.AlignmentFile(file_path, "rb")
counts = [0 for i in range(10000)]

last_query_name = None
j = 0
count = 1

for (i, read) in enumerate(samfile.head(10000)):
	if read.query_name == last_query_name:
		count += 1
	else:
		j = i
		count = 1
	for k in range(j,i+1):
		counts[k] = count
	last_query_name = read.query_name
	
print(counts)
