#!/usr/bin/env python
from __future__ import division
from subprocess import call
import HTSeq
import random
import argparse
import numpy as np
import subprocess
import os
import sys
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
parser.add_argument('-abyss', action='store_true', default=False, dest='abyss_switch', help='Choose assembler ABYSS')
parser.add_argument('-megahit', action='store_true', default=False, dest='megahit_switch', help='Choose assembler Megahit')


# Choose aligner
parser.add_argument('-bwa', action='store_true', default=False, dest='bwa_switch', help='Choose aligner BWA')
parser.add_argument('-bowtie2', action='store_true', default=False, dest='bowtie_switch', help='Choose aligner Bowtie2')


# More indepth choices
parser.add_argument('-align', action='store', default= '70', dest='alignmentrate', help='The alignment rate wanted before ending, default is 70%.')
parser.add_argument('-iter', action='store', default= '5', dest='iterations', help='Number of iterations before ending, default is 5.')
parser.add_argument('-m', action='store_true', default=True, dest='merge_switch', help='Merges all contig files into a singluar assembly, default is true.')
parser.add_argument('-k', action='store', default= '31', dest='kmer', help='Enter Kmer size of choice, default is 31.')
parser.add_argument('-R', action='store', dest='RAM', help='Enter percentage of file to be used as sub-sample e.g. if -R 0.25 is used  25% of the reads will be used in the sub-sample, no default')
parser.add_argument("-x1", type=str, action='store',default= '',  dest='extra1', help="Allows additional options for assembly to be used in Velveth or ABYSS steps, used for options starting with -")
parser.add_argument("-x2", type=str, action='store',default= '',  dest='extra2', help="Allows additional options for assembly to be used in Velveth or ABYSS steps, used for options starting with --")
parser.add_argument("-u", type=str, action='store',default= ' ',  dest='bowtie_extra', help="Allows additional options for alignment to be used in Bowtie2")
parser.add_argument('-f', action='store_true', default=False, dest='scaffold_switch', help='Conducts a final assembly of the produced contigs.')


args = parser.parse_args()



# Place each of the input into a simple variable to call
INPUT = str(args.input)
OUTPUT = str(args.output)
EXTRA1 = str(args.extra1)
EXTRA2 = str(args.extra2)
BOWTIE_EXTRA = str(args.bowtie_extra)
iterations = int(args.iterations)
alignmentwanted = int(args.alignmentrate)
ksize = str(args.kmer)
RAM = float(args.RAM)


EXTRA = ''
temp = ''
if EXTRA1 != '':
	temp = EXTRA + '-' + EXTRA1
if EXTRA2 != '':
	EXTRA = temp + ' --' + EXTRA2
	

scaffold_switch = False
if  args.scaffold_switch == True:
	scaffold_switch = True

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
		fastafile.write(read.seq + '\n')
	fastafile.close()
	newinput = 'Converted_input.fa'


failed = False
failedfirst = False

unalignedfile = newinput
currentfile = newinput




print "Using file;" , currentfile
sys.stdout.flush()

# Provide an initial count of the raw reads so we can work out the alignment rate
totalreads = 0

for read in HTSeq.FastaReader(currentfile):
	totalreads += 1

currentiter = 0

unalignedfile = currentfile
RAM2 = RAM*100

while currentiter < iterations:
	if currentalignrate >= alignmentwanted:
		break
	else:
		currentiter +=1
		print "Starting sub-sampling; " + str(datetime.now())
		sys.stdout.flush()

		# Subsample from current file
		subsample = open('subsample.fa','w')

		currentfasta = HTSeq.FastaReader(unalignedfile)
		#print RAM
		for read in currentfasta:
			random_num = random.randint(0,100)
			if RAM2 >= random_num:
				subsample.write('>' + read.name + '\n')
				subsample.write(read.seq + '\n')
			else:
				continue


		subsample.close()
		currentfile = 'subsample.fa'
		print "Subsampling completed."
		sys.stdout.flush()
		
		print "Ending sub-sampling; " + str(datetime.now())
		sys.stdout.flush()

		# Run assembly
		print "Starting Assembly; " + str(datetime.now())
		if args.velvet_switch == True:
			bashCommand = 'velveth out-dir ' + str(ksize) + ' -' + filetype + ' ' + currentfile + ' ' + EXTRA + ' | cat >  Assembly_log'
			notneeded = call(bashCommand, shell=True)
			if limit == True:
				bashCommand = 'velvetg out-dir -exp_cov auto -min_contig_lgth 300 | cat > Assembly_log'
			else:
				bashCommand = 'velvetg out-dir -exp_cov auto | cat > Assembly_log'				
			notneeded = call(bashCommand, shell=True)
			# Run velvet code


		elif args.abyss_switch == True:
			# Run ABYSS code
			bashCommand = 'ABYSS -k' + str(ksize) +  ' ' + currentfile + ' -o temp_contigs.fa ' + ' ' + EXTRA + ' | cat >  Assembly_log'
			notneeded = call(bashCommand, shell=True)

		elif args.megahit_switch == True:
			# Run Megahit code
			contigfilename = OUTPUT + '.' + str(currentiter) 
			bashCommand = 'megahit ' + ' -r ' + currentfile + ' -o Iteration_' + str(currentiter) + ' ' + EXTRA + ' | cat >  Assembly_log'
			notneeded = call(bashCommand, shell=True)
			bashCommand = 'cp Iteration_' + str(currentiter) + '/final.contigs.fa ' + contigfilename
			notneeded = call(bashCommand, shell=True)
		print "Assembly complete."
		sys.stdout.flush()

		print "Ending Assembly; " + str(datetime.now())
		sys.stdout.flush()

		# Check if assembly produced anything




		# Run alignment
		print "Starting Alignment; " + str(datetime.now())
		sys.stdout.flush()

		if args.bowtie_switch == True:

			contigfilename = OUTPUT + '.' + str(currentiter)  
			if args.velvet_switch == True:
				bashCommand = 'mv out-dir/contigs.fa ' + contigfilename # Move file and rename so it isnt deleted
				notneeded = call(bashCommand, shell=True)
				if os.stat(contigfilename).st_size == 0:
					print "Assembly failed to produce any contigs."
					sys.stdout.flush()

					print "Spherical will now exit."
					sys.stdout.flush()

					failed = True
					bashCommand = 'rm ' + contigfilename
					notneeded = call(bashCommand, shell=True)
					if currentiter == 1:
						failedfirst = True
					sys.exit()
			elif args.abyss_switch == True:
				bashCommand = 'mv temp_contigs.fa ' + contigfilename # Move file and rename so it isnt deleted
				notneeded = call(bashCommand, shell=True)
				if os.stat(contigfilename).st_size == 0:
					print "Assembly failed to produce any contigs."
					sys.stdout.flush()

					print "Spherical will now exit."
					sys.stdout.flush()

					failed = True
					bashCommand = 'rm ' + contigfilename
					notneeded = call(bashCommand, shell=True)
					if currentiter == 1:
						failedfirst = True
					sys.exit()
			
			bashCommand = 'bowtie2-build -f ' + contigfilename + ' ' + 'Current_round_index | cat > Index_log ' # Creates the bowtie index
			notneeded = call(bashCommand, shell=True)
			bashCommand = 'bowtie2 -f -N 1 --un Unaligned.fa.' + str(currentiter) + ' -U ' + unalignedfile + ' ' +  BOWTIE_EXTRA + ' --al /dev/null -x Current_round_index -S /dev/null | cat > Alignment_log ' # Runs bowtie itself
			notneeded = call(bashCommand, shell=True)
			unalignedfile = 'Unaligned.fa.' + str(currentiter)
			# Run bowtie code
		elif args.bwa_switch == True:
			# Run BWA code
			print "Not yet available"
		print "Alignment complete."
		sys.stdout.flush()

		print "Ending Alignment; " + str(datetime.now())
		sys.stdout.flush()



		# Calculate current alignment rate
		unalignedreads = 0
		for reaf in HTSeq.FastaReader(unalignedfile):
			unalignedreads +=1
		print 'Unaligned reads; ', unalignedreads
		sys.stdout.flush()
		alignmentrate = (totalreads - unalignedreads) / totalreads * 100
		currentalignrate = alignmentrate
		print 'Percentage alignment; ', alignmentrate
		sys.stdout.flush()







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
		sys.stdout.flush()





	if args.scaffold_switch == True:
		print "Starting Scaffold step; " + str(datetime.now())
		if args.velvet_switch == True:
			bashCommand = 'velveth out-dir ' + str(ksize) + ' -' + filetype + ' ' + '*.combined.fa' + ' ' + EXTRA + ' | cat >  Assembly_log'
			notneeded = call(bashCommand, shell=True)
			if limit == True:
				bashCommand = 'velvetg out-dir -exp_cov auto -min_contig_lgth 300 | cat > Assembly_log'
			else:
				bashCommand = 'velvetg out-dir -exp_cov auto | cat > Assembly_log'				
			notneeded = call(bashCommand, shell=True)
			# Run velvet code


		elif args.abyss_switch == True:
			# Run ABYSS code
			bashCommand = 'ABYSS -k' + str(ksize) +  ' ' + currentfile + ' -o temp_contigs.fa ' + ' ' + EXTRA + ' | cat >  Assembly_log'
			notneeded = call(bashCommand, shell=True)

		elif args.megahit_switch == True:
			print "Scaffoling is not avaialble for Megahit" 

		print "Scaffold complete."
		sys.stdout.flush()

		print "Ending Scaffold; " + str(datetime.now())
		sys.stdout.flush()

		# Run alignment
		print "Starting Alignment; " + str(datetime.now())
		sys.stdout.flush()

		if args.bowtie_switch == True:

			contigfilename = OUTPUT + '.' + 'scaffold.fa'  
			if args.velvet_switch == True:
				bashCommand = 'mv out-dir/contigs.fa ' + contigfilename # Move file and rename so it isnt deleted
				notneeded = call(bashCommand, shell=True)
				if os.stat(contigfilename).st_size == 0:
					print "Scaffold failed to produce any contigs."
					sys.stdout.flush()

					print "Spherical will now exit."
					sys.stdout.flush()

					failed = True
					bashCommand = 'rm ' + contigfilename
					notneeded = call(bashCommand, shell=True)
					if currentiter == 1:
						failedfirst = True
					sys.exit()
			elif args.abyss_switch == True:
				bashCommand = 'mv temp_contigs.fa ' + contigfilename # Move file and rename so it isnt deleted
				notneeded = call(bashCommand, shell=True)
				if os.stat(contigfilename).st_size == 0:
					print "Assembly failed to produce any contigs."
					sys.stdout.flush()

					print "Spherical will now exit."
					sys.stdout.flush()

					failed = True
					bashCommand = 'rm ' + contigfilename
					notneeded = call(bashCommand, shell=True)
					if currentiter == 1:
						failedfirst = True
					sys.exit()


			elif args.abyss_switch == True:
				print "Scaffolding not available, no alignment conducted. Program will now likely crash due to bad programming. Sorry."
			bashCommand = 'bowtie2-build -f ' + contigfilename + ' ' + 'Current_round_index | cat > Index_log ' # Creates the bowtie index
			notneeded = call(bashCommand, shell=True)
			bashCommand = 'bowtie2 -f -N 1 --un Non_scaffolded_contigs.fa '   +  BOWTIE_EXTRA + ' --al Scaffolded_contigs.fa -x Current_round_index -S /dev/null | cat > Alignment_log ' # Runs bowtie itself
			notneeded = call(bashCommand, shell=True)


			print "Ending Scaffold alignment; " + str(datetime.now())
			sys.stdout.flush()





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
		sys.stdout.flush()


else:
	print "Spherical failed with no successful assembly produced."
	print "Please try again with a different kmer size."
	sys.stdout.flush()
