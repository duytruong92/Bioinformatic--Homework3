#name: Duy Truong
#class: Bioinformatic

import random
from collections import OrderedDict



length = 400
f = open('NC_001133a.fna', 'r')
n = len(f.read())
#get the DNA sequences and get random fragment
stop= n-length
times = n/length*10
f.seek(0,0)
newline=f.read()
newline = newline.replace("\n","")
original = newline
array = [0 for i in range(0,times)]
for i in range (0, times):
	jump = random.randint(0,stop)
	array[i] = newline[jump:jump+400]
f.close()
#keep track the over lap
list(OrderedDict.fromkeys(array))
o = open("overlap1.fna","w")
overlap = [[0 for m in range(0,times)] for m in range(0,times)]

for m  in xrange(0,times):
	for n  in xrange(0,times):
		if(m==n):
			overlap[m][m]= 0
		for k  in xrange(0,length-1):
			check= array[m][length-k-1: length]
			check1= array[n][0:k+1]
			if(check1 == check):
				overlap[m][n]= len(check1)
			
		o.write(str(m)+", " +str(n)+", " +str(overlap[m][n])+"\n")
o.close()



newDNA = [0 for m in range(0,times)]
for m in range(0,times):
	newDNA[m] = array[m]

biggest_overlap= [0 for m in range(0,times)]
index_biggest =  [0 for m in range(0,times)]


for x in xrange(0,times):
	index = 0
	biggest_overlap[x] = max(overlap[x])
	while (overlap[x][index] != max(overlap[x])):
		index+=1
	index_biggest[x] = index
#merge the DNA gradment to get new DNA
temp = biggest_overlap
temp1 = index_biggest
for x in xrange(0,times):
	m=x
	while (temp [index_biggest[m]] != 0):
		newDNA[x] += array[index_biggest[m]][biggest_overlap[m]+1: length-1]
		temp [index_biggest[m]] =0
		m = index_biggest[m]


o = open("overlap_original.fna","w")
#compare with the original chromosomes				
for x in xrange(0, times):
	original_overlap = ""
	for y in xrange (0, len(newDNA[x])):
		n = len(newDNA[x])
		for z in xrange (0, len(original)- n):
			l = len (newDNA[x][y:len(newDNA[x])] )
			if( newDNA[x][y:len(newDNA[x])] == original[z:z+l] ):
				temp = newDNA[x][y:len(newDNA[x])]
			if(len(temp) > len (original_overlap)):
				original_overlap = temp
	o.write("newDNA number: "+str(x)+"\n" + "length: " + str(len(original_overlap)) + "\n" +"overlap with original: " + original_overlap + "\n")
o.close()


f = open('NC_001134a.fna', 'r')
n = len(f.read())
#get the DNA sequence and take random fragment
stop= n-length
times1 = n/length*10
f.seek(0,0)
newline=f.read()
newline = newline.replace("\n","")
original1 = newline
array1 = [0 for i in range(0,times1)]
for i in range (0, times1):
	jump = random.randint(0,stop)
	array1[i] = newline[jump:jump+400]
f.close()

array2= array+ array1
times2 = times + times1
o = open("overlap2.fna","w")
overlap1 = [[0 for m in range(0,times2)] for m in range(0,times2)]

for m  in xrange(0,times2):
	for n  in xrange(0,times2):
		if(m==n):
			overlap1[m][m]= 0
		for k  in xrange(0,length-1):
			check= array2[m][length-k-1: length]
			check1= array2[n][0:k+1]
			if(check1 == check):
				overlap1[m][n]= len(check1)
		o.write(str(m)+", " +str(n)+", " +str(overlap1[m][n])+"\n")
o.close()

newDNA1 = [0 for m in range(0,times2)]
for m in range(0,times2):
	newDNA1[m] = array2[m]
print newDNA1
biggest_overlap1= [0 for m in range(0,times2)]
index_biggest1 =  [0 for m in range(0,times2)]


for x in xrange(0,times2):
	index = 0
	biggest_overlap1[x] = max(overlap1[x])
	while (overlap1[x][index] != max(overlap1[x])):
		index+=1
	index_biggest1[x] = index

temp = biggest_overlap1
temp1 = index_biggest1
print temp
print temp1

#merge the DNA fragments to get new DNA
for x in xrange(0,times2):
	m=x
	while (temp [index_biggest1[m]] != 0):
		newDNA1[x] += array2[index_biggest1[m]][biggest_overlap1[m]+1: length-1]
		temp [index_biggest1[m]] =0
		m = index_biggest1[m]

#compare with the original chormosomes 1
o = open("overlap_original1__chromosome1.fna","w")
				
for x in xrange(0, times2):
	original_overlap = ""
	temp2 = ""
	for y in xrange (0, len(newDNA1[x])):
		n = len(newDNA1[x])
		for z in xrange (0, len(original1)- n):
			l = len (newDNA1[x][y:len(newDNA1[x])] )
			if( newDNA1[x][y:len(newDNA1[x])] == original[z:z+l] ):
				temp2 = newDNA1[x][y:len(newDNA1[x])]
			if(len(temp2) > len (original_overlap)):
				original_overlap = temp2
	print original_overlap
	o.write("newDNA number: "+str(x)+ "\n" + "length: " + str(len(original_overlap)) + "\n" +"overlap with original: " + str(original_overlap) + "\n")
o.close()
#compare with the original chormosomes 2
o = open("overlap_original1__chromosome2.fna","w")
				
for x in xrange(0, times2):
	original_overlap = ""
	temp2 = ""
	for y in xrange (0, len(newDNA1[x])):
		n = len(newDNA1[x])
		for z in xrange (0, len(original1)- n):
			l = len (newDNA1[x][y:len(newDNA1[x])] )
			if( newDNA1[x][y:len(newDNA1[x])] == original1[z:z+l] ):
				temp2 = newDNA1[x][y:len(newDNA1[x])]
			if(len(temp2) > len (original_overlap)):
				original_overlap = temp2
	print original_overlap
	o.write("newDNA number: "+str(x)+ "\n" + "length: " + str(len(original_overlap)) + "\n" +"overlap with original: " + str(original_overlap) + "\n")
o.close()