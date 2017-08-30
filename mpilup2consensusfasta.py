import sys


threshold=float(sys.argv[1])
header=sys.argv[2]

nucleotides=['A','C','T','G']
def basefun(ref,inp):
	string =''
	for thebase in inp:
		if thebase.isalpha() and thebase.upper() in nucleotides:
			string += thebase.upper()
		else:
			if thebase in ['.',',']:
				string += ref
	return string

newseq=''
for line in sys.stdin:
	col=line.split()
	refbase=col[2]
	depth=int(col[3])
	pileupcol=col[4]
	pileup=basefun(refbase,pileupcol)	
		
	alleles=list(set(pileup))
	allalleles=list(set(pileup+refbase))
	#find major base
	alleledict={}
	refcount=0
	for allele in alleles:
		alleledict[allele]=pileup.count(allele)

	counts=alleledict.values()
	if len(counts) <1:continue
	maxcount=max(counts)

	chosenallele=[]	
	for allele in alleledict.keys():
		if alleledict[allele] == maxcount:
			chosenallele.append(allele)
	if len(chosenallele) > 1:continue
	majorallele=chosenallele[0]
	majF=1.0*maxcount/depth
	
	if majF >= threshold:
		newallele=majorallele
		#newseq +=majorallele
	else:
		#newseq += 'N'
		newallele='N'
	newseq += newallele

	#print newallele,line, 
		
print '>'+header		
print newseq