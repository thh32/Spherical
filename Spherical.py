#!/usr/bin/env python
from __future__ import division
from subprocess import call
import HTSeq
import random
import argparse
import numpy as np
import subprocess


parser = argparse.ArgumentParser() #simplifys the wording of using argparse as stated in the python tutorial
parser.add_argument("-i", type=str, action='store',  dest='input', help="Input the FASTA/Q file") # allows input of the forward read
parser.add_argument("-o", type=str, action='store',  dest='output', help="Name for the output files") # allows input of the forward read
parser.add_argument('-fasta', action='store_true', default=False, dest='fasta_switch', help='Input is fasta')
parser.add_argument('-fastq', action='store_true', default=False, dest='fastq_switch', help='Input is fastq')
parser.add_argument('-align', action='store', dest='alignmentrate', help='The alignment rate wanted before ending')
parser.add_argument('-iter', action='store', dest='iterations', help='Number of iterations before ending')
parser.add_argument('-m', action='store_true', default=False, dest='merge_switch', help='Merges all contig files into a singluar assembly.')
parser.add_argument('-bwa', action='store_true', default=False, dest='bwa_switch', help='Choose aligner BWA')
parser.add_argument('-bowtie2', action='store_true', default=False, dest='bowtie_switch', help='Choose aligner Bowtie2')
parser.add_argument('-velvet', action='store_true', default=False, dest='velvet_switch', help='Choose assembler Velvet')
parser.add_argument('-soapdenovo', action='store_true', default=False, dest='soapdenovo_switch', help='Choose assembler SOAPdenovo')
parser.add_argument('-k', action='store', dest='kmer', help='Enter Kmer size of choice')
parser.add_argument('-R', action='store', dest='RAM', help='Enter RAM available')

args = parser.parse_args()



# Place each of the input into a simple variable to call
INPUT = str(args.input)
OUTPUT = str(args.output)
iterations = int(args.iterations)
alignmentwanted = int(args.alignmentrate)
ksize = str(args.kmer)
RAM = int(args.RAM)


filetype = 0

if args.fasta_switch == True:
	inputfile = HTSeq.FastaReader( INPUT )
	filetype = 'fasta'
	newinput = INPUT
if args.fastq_switch == True:
	inputfile = HTSeq.FastqReader( INPUT, "phred")
	filetype = 'fastq'
	fastafile = open('Converted_input.fa','w')
	for read in inputfile:
		fastafile.write( ">" + read.name + '\n')
		fastfile.write(read.seq + '\n')
	fastafile.close()
	newinput = 'Converted_input.fa'





currentfile = newinput






# Provide an initial count of the raw reads so we can work out the alignment rate
totalreads = 0


for read in inputfile:
	totalreads +=1

currentiter = 0

unalignedfile = currentfile
currentalignrate = 0

while currentiter < iterations:
	if currentalignrate >= alignmentwanted:
		break
	else:
		currentiter +=1


		# Subsample from current file
		subsample = open('subsample.fa','w')
		amountneeded = RAM * 170000 
		count = 0
		if count < amountneeded:
			for read in HTSeq.FastaReader(currentfile):
				subsample.write('>' + read.name + '\n')
				subsample.write(read.seq + '\n')
				count +=1
		subsample.close()
		currentfile = 'subsample.fa'


		# Run assembly
		if args.velvet_switch == True:
			bashCommand = 'velveth out-dir ' + str(ksize) + ' -' + filetype + ' ' + currentfile
			notneeded = call(bashCommand, shell=True)
			bashCommand = 'velvetg out-dir -exp_cov auto'
			notneeded = call(bashCommand, shell=True)
			# Run velvet code
		elif args.soapdenovo_switch == True:
			# Run soapdenovo code
			print "Not yet available"


		# Run alignment
		if args.bowtie_switch == True:
			contigfilename = OUTPUT + '.' + str(currentiter)  
			bashCommand = 'mv out-dir/contigs.fa ' + contigfilename # Move file and rename so it isnt deleted
			notneeded = call(bashCommand, shell=True)
			bashCommand = 'bowtie2-build -f ' + currentfile + ' ' + 'Current_round_index' # Creates the bowtie index
			notneeded = call(bashCommand, shell=True)
			bashCommand = 'bowtie2 -f -N 1 --un Unaligned.fa -U ' + unalignedfile + '--al-gz /dev/null -x Current_round_index -S /dev/null' # Runs bowtie itself
			notneeded = call(bashCommand, shell=True)
			unalignedfile = 'Unaligned.fa'
			# Run bowtie code
		elif args.bwa_switch == True:
			# Run BWA code
			print "Not yet available"



		# Calculate current alignment rate
		unalignedreads = 0
		for reaf in HTSeq.FastaReader(unalignedfile):
			unalignedreads +=1
		print 'Unaligned reads; ', unalignedreads
		alignmentrate = (totalreads - unalignedreads) / totalreads * 100
		currentalignrate = alignmentrate
		print 'Percentage alignment; ', alignmentrate




# Merge output files

if args.merge_switch == True:
	#  Run code to merge all the contig files from each assembly into a single file where all reads have been renamed.
	bashCommand = 'cat ' + OUTPUT + '* > '  + OUTPUT + '.combined.fa'
	notneeded = call(bashCommand, shell=True)
	cfile =  OUTPUT + '.combined.fa'
	currentfile = HTSeq.FastaReader(cfile)
	listolengths = []
	for read in currentfile:
		listolengths.append(len(read.seq))
	print "Statistics for merged file."
	print 'Mean contig length; ', np.mean(listolengths)
	print '25 percentile; ' , np.percentile(listolengths,25)
	print '75 percentile; ', np.percentile(listolengths,75)
	print 'Standard deviation; ', np.std(listolengths)
	print "Final alignment rate; ", alignmentrate

# Provide final statistics on every round

for cround in range(0,currentiter):
	cround +=1
	cfile = OUTPUT + '.' + str(cround)  
	currentfile = HTSeq.FastaReader(cfile)
	print "For iteration ; ", cround
	listolengths = []
	for read in currentfile:
		listolengths.append(len(read.seq))
	print "Statistics for merged file."
	print 'Mean contig length; ', np.mean(listolengths)
	print '25 percentile; ' , np.percentile(listolengths,25)
	print '75 percentile; ', np.percentile(listolengths,75)
	print 'Standard deviation; ', np.std(listolengths)



