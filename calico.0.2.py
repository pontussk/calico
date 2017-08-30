#pontus.skoglund@gmail.com


import sys
import math
from optparse import OptionParser
import subprocess

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("--notransitions", action="store_true", dest="notransitions",help="notransitions",default=False)
parser.add_option("--nodamage", action="store_true", dest="nodamage",help="nodamage",default=False)
parser.add_option("--mindepth", action="store", type="int", dest="mindepth",help="mindepth",default=3)
parser.add_option("--maxdepth", action="store", type="int", dest="maxdepth",help="maxdepth",default=1000)
parser.add_option("--verbose", action="store_true", dest="verbose",help="verbose",default=False)
parser.add_option("--noheader", action="store_true", dest="noheader",help="noheader",default=False)
parser.add_option("--indels", action="store_true", dest="indels",help="indels",default=True)
parser.add_option("--noindels", action="store_false", dest="indels",help="indels",default=True)
parser.add_option("--samtoolspath", action="store", type="string", dest="samtoolspath",help="samtoolspath",default='samtools')
parser.add_option("-m", "--requiremapq", action="store", type="int", dest="mapq",help="mapq",default=0)
parser.add_option("-q", "--requirebasequal", action="store", type="int", dest="qual",help="qual",default=0)
parser.add_option("--outputsites", action="store_true", dest="outputsites",help="outputsites",default=False)
parser.add_option("-r","--ref_file", action="store", type="string", dest="ref_file",help="ref_file",default=False)
parser.add_option("-b", "--bam", action="store", type="string", dest="bam",help="bam",default=False)
parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile",help="outfile",default=False)
parser.add_option("--distribution", action="store_true", dest="distribution",help="distribution",default=False)
(options, args) = parser.parse_args()



if options.outfile != False:
	sys.stdout = open(options.outfile, 'w')

def factorial(f_n): 
	if f_n < 2: return 1
	return reduce(lambda f_x, f_y: f_x*f_y, xrange(2, int(f_n)+1))

def prob(f_s, f_p, f_n):
	f_x = 1.0 - f_p
	f_a = f_n - f_s
	f_b = f_s + 1.0
	f_c = f_a + f_b - 1.0
	f_prob = 0
	for f_j in xrange(int(f_a), int(f_c) + 1):
		f_prob += factorial(f_c) / (factorial(f_j)*factorial(f_c-f_j)) * f_x**f_j * (1 - f_x)**(f_c-f_j)
	return f_prob


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


FNULL = open('/dev/null', 'w')
def mpileup_get (chrom,thebam):
	chrom=str(chrom)#'chr'+str(chrom)
	if options.addchr:
		chrom= 'chr'+chrom

	cmd_line=[options.samtoolspath, 'mpileup', '-q', str(options.mapq),'-Q', str(options.qual),'-f',options.ref_file, thebam] #,'-r',chrom

	if options.verbose: 
		for line in cmd_line: print line,

	if options.verbose: 
		outp_file=subprocess.Popen(cmd_line, stdout=subprocess.PIPE)
	else:
		outp_file=subprocess.Popen(cmd_line, stdout=subprocess.PIPE,stderr=FNULL)

	pileupline = outp_file.stdout.readlines()

	return pileupline



errorbase=0
nonerrorbase=0
endobase=0
contbase=0
diagnostic=0
diagnostic_list=[]

if options.bam:

	thepileup=mpileup_get(options.chromosome,options.bam)
else:
	thepileup=sys.stdin

for line in thepileup:
	col=line.split('\t')
	chromosome=col[0]
	position=col[1]
	refbase=col[2]
	depth=int(col[3])
	if depth < options.mindepth:continue
	if depth > options.maxdepth:continue

	if refbase not in nucleotides:continue
	pileupcol=col[4]
	pileup=basefun(refbase,pileupcol)
	
	#remove indels
	if options.indels ==False:
		#print hello
		if '+' in pileupcol or '-' in pileupcol or '*' in pileupcol:continue
	#print line,
	#print pileup
	alleles=list(set(pileup))
	allalleles=list(set(pileup+refbase))
	
	if options.notransitions:
		if 'C' in allalleles and 'T' in allalleles:continue
		if 'G' in allalleles and 'A' in allalleles:continue

		
	#print allalleles,alleles
	if len(allalleles) ==1:continue
	
	if options.distribution:
		if len(allalleles)>3:continue
		else:
			basefreq=1.0*pileup.count(refbase)/depth
			if basefreq==1.0:continue
			print depth,'\t',basefreq
			continue

	"""
	find major base
	"""
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
	majF=maxcount/depth
	#if majF <0.8:continue
	
	
	if majorallele == refbase:continue
	
	if options.nodamage:
		if 'C' == majorallele and 'T' == refbase:continue
		if 'G' == majorallele and 'A' == refbase:continue

	refcount=pileup.count(refbase)
	majcount=pileup.count(majorallele)
	contbase += refcount
	endobase += majcount

	
	diagnostic +=1
	if options.outputsites:
		diagnostic_list.append(str(position))
	if options.verbose:
		print position,'\t',refbase,'\t',pileup,'\t',contbase,'\t',contbase+endobase,'\t',1.0*contbase/(contbase+endobase),'\t',1.0*refcount/(refcount+majcount)
	
n=endobase+contbase
if n==0:
	freq='NA'
	confinterval='NA'
	confstring='NA\tNA'
	diagnostic_list=['NA']
	
else:
	freq=1.0*contbase/n
	freq=round(freq,3)
	if freq==0.0:
		for s_i in range(1,999):
			s_i=s_i/1000.0
			probability=prob(0,s_i,n) 
			#poisson=prob(0,s_i,n) 
			if probability <= 0.05:
				confstring=str(0.0)+'\t'+str(s_i)
				break
	else:			
		confinterval=1.96*math.sqrt((freq*(1.0-freq))/n)
		confstring=str(max(0,round(freq-confinterval,3)))+'\t'+str(min(1.0,round(freq+confinterval,3)))


if options.outputsites:
	confstring = confstring+'\t'+','.join(diagnostic_list)

if options.noheader ==False:
	print 'SITES','\t','MAJ','\t','MIN','\t','CONT','\t','CI_low','\t','CI_up',
	if options.outputsites:
		print '\t Informative_sites'
	print ''
print diagnostic,'\t',endobase,'\t',contbase,'\t',freq,'\t',confstring


