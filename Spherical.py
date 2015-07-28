#!/usr/bin/env python
from __future__ import division
from subprocess import call
import HTSeq
import random
import argparse
import numpy as np
import subprocess
import os
from datetime import datetime



parser = argparse.ArgumentParser() #simplifys the wording of using argparse as stated in the python tutorial


# Basic input output files
parser.add_argument("-i", type=str, action='store',  dest='input', help="Input the FASTA/Q file") # allows input of the forward read
parser.add_argument("-o", type=str, action='store',  dest='output', help="Name for the output files") # allows input of the forward read


# Input file type
parser.add_argument('-fasta', action='store_true', default=False, dest='fasta_switch', help='Input is fasta')
parser.add_argument('-fastq', action='store_true', default=False, dest='fastq_switch', help='Input is fastq')


# Limit assembly or not
parser.add_argument('-limit', action='store_true', default=False, dest='limit_switch', help='Limit velvet to produce contigs large than 300bp.')

# Choose assembler
parser.add_argument('-velvet', action='store_true', default=False, dest='velvet_switch', help='Choose assembler Velvet')
parser.add_argument('-soapdenovo', action='store_true', default=False, dest='soapdenovo_switch', help='Choose assembler SOAPdenovo')
parser.add_argument('-abyss', action='store_true', default=False, dest='abyss_switch', help='Choose assembler ABYSS')


# Choose aligner
parser.add_argument('-bwa', action='store_true', default=False, dest='bwa_switch', help='Choose aligner BWA')
parser.add_argument('-bowtie2', action='store_true', default=False, dest='bowtie_switch', help='Choose aligner Bowtie2')


# More indepth choices
parser.add_argument('-align', action='store', default= '70', dest='alignmentrate', help='The alignment rate wanted before ending, default is 70%.')
parser.add_argument('-iter', action='store', default= '5', dest='iterations', help='Number of iterations before ending, default is 5.')
parser.add_argument('-m', action='store_true', default=True, dest='merge_switch', help='Merges all contig files into a singluar assembly, default is true.')
parser.add_argument('-k', action='store', default= '31', dest='kmer', help='Enter Kmer size of choice, default is 31.')
parser.add_argument('-R', action='store', dest='RAM', help='Enter fraction of file to be used as sub-sample e.g. if -R 3 is used 1 third of the reads will be used in the sub-sample, no default')

args = parser.parse_args()



# Place each of the input into a simple variable to call
INPUT = str(args.input)
OUTPUT = str(args.output)
iterations = int(args.iterations)
alignmentwanted = int(args.alignmentrate)
ksize = str(args.kmer)
RAM = int(args.RAM)


limit = False
if args.limit_switch == True:
	limit = True


# Provide an initial count of the raw reads so we can work out the alignment rate
totalreads = 0



currentiter = 0

filetype = 'fasta'

currentalignrate = 0
if args.fasta_switch == True:
	newinput = INPUT
	filetype = 'fasta'
if args.fastq_switch == True:
	inputfile = HTSeq.FastqReader( INPUT, "phred")
	filetype = 'fasta'
	fastafile = open('Converted_input.fa','w')
	for read in inputfile:
		fastafile.write( ">" + read.name + '\n')
		fastfile.write(read.seq + '\n')
	fastafile.close()
	newinput = 'Converted_input.fa'


failed = False
failedfirst = False

unalignedfile = newinput
currentfile = newinput




print "Using file;" , currentfile

# Provide an initial count of the raw reads so we can work out the alignment rate
totalreads = 0



currentiter = 0

unalignedfile = currentfile


while currentiter < iterations:
	if currentalignrate >= alignmentwanted:
		break
	else:
		currentiter +=1
		print "Starting sub-sampling; " + str(datetime.now())

		# Subsample from current file
		subsample = open('subsample.fa','w')

		currentfasta = HTSeq.FastaReader(unalignedfile)
		for read in currentfasta:
			if random.random()%int(RAM) == 0:
				subsample.write('>' + read.name + '\n')
				subsample.write(read.seq + '\n')
			else:
				continue


		subsample.close()
		currentfile = 'subsample.fa'
		print "Subsampling completed."
		
		print "Ending sub-sampling; " + str(datetime.now())

		# Run assembly
		print "Starting Assembly; " + str(datetime.now())
		if args.velvet_switch == True:
			bashCommand = 'velveth out-dir ' + str(ksize) + ' -' + filetype + ' ' + currentfile + ' | cat >  Assembly_log'
			notneeded = call(bashCommand, shell=True)
			if limit == True:
				bashCommand = 'velvetg out-dir -exp_cov auto -min_contig_lgth 300 | cat > Assembly_log'
			else:
				bashCommand = 'velvetg out-dir -exp_cov auto | cat > Assembly_log'				
			notneeded = call(bashCommand, shell=True)
			# Run velvet code
		elif args.soapdenovo_switch == True:
			# Run soapdenovo code
			print "Not yet available"

		elif args.abyss_switch == True:
			# Run ABYSS code
			print "not yet available"
		print "Assembly complete."

		print "Ending Assembly; " + str(datetime.now())

		# Check if assembly produced anything




		# Run alignment
		print "Starting Alignment; " + str(datetime.now())
		if args.bowtie_switch == True:
			contigfilename = OUTPUT + '.' + str(currentiter)  
			bashCommand = 'mv out-dir/contigs.fa ' + contigfilename # Move file and rename so it isnt deleted
			notneeded = call(bashCommand, shell=True)
			if os.stat(contigfilename).st_size == 0:
				print "Assembly failed to produce any contigs."
				print "Spherical will now exit."
				failed = True
				bashCommand = 'rm ' + contigfilename
				notneeded = call(bashCommand, shell=True)
				if currentiter == 1:
					failedfirst = True
				break
			bashCommand = 'bowtie2-build -f ' + contigfilename + ' ' + 'Current_round_index | cat > Index_log ' # Creates the bowtie index
			notneeded = call(bashCommand, shell=True)
			bashCommand = 'bowtie2 -f -N 1 --un Unaligned.fa.' + str(currentiter) + ' -U ' + unalignedfile + ' --al /dev/null -x Current_round_index -S /dev/null | cat > Alignment_log ' # Runs bowtie itself
			notneeded = call(bashCommand, shell=True)
			unalignedfile = 'Unaligned.fa.' + str(currentiter)
			# Run bowtie code
		elif args.bwa_switch == True:
			# Run BWA code
			print "Not yet available"
		print "Alignment complete."
		print "Ending Alignment; " + str(datetime.now())



		# Calculate current alignment rate
		unalignedreads = 0
		for reaf in HTSeq.FastaReader(unalignedfile):
			unalignedreads +=1
		print 'Unaligned reads; ', unalignedreads
		alignmentrate = (totalreads - unalignedreads) / totalreads * 100
		currentalignrate = alignmentrate
		print 'Percentage alignment; ', alignmentrate






# Check if failed before output could be made or not

if failedfirst == False:
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
		print '\n'

	# Provide final statistics on every round

	if failed == True:
		currentiter = currentiter - 1

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
		print '\n'

else:
	print "Spherical failed with no successful assembly produced."
	print "Please try again with a different kmer size."



