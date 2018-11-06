#!/usr/bin/env python3
# Python3 must be loaded prior to implementing this program.
#######################
## IMPORT STATEMENTS ##
#######################
import argparse
from argparse import RawTextHelpFormatter
import re
from subprocess import call
from subprocess import check_output
##########################
## FUNCTION DEFINITIONS ##
##########################
def parseSAM(lineSAM):
    """
    description: parseSAM takes 1 line from a SAM file and returns
    the UMI, chromosome #, POS, CIGAR, and the strand.
    INPUT: A SAM file read-line
    OUTPUT: UMI,CHROM,POS,CIGAR,STRAND
    """
    line = lineSAM.split("\t")
    header = line[0]
    umi = re.findall(".*:([A,C,G,T,N]{8}$)", header)
    barcode = umi[0]
    chromo = line[2]
    POS = line[3]
    CIGAR = line[5]
    flag = int(line[1])
    if((flag & 16) != 16):
        strand = "+"
    else:
        strand = "-"
    return barcode,chromo,POS,CIGAR,strand

def plusPOS(POS,CIGAR):
    """
    description: plusPOS checks if soft-clipping has occurred at the left-most mapped position on a "+" read.
    If soft-clipping has occurred, truePOS performs a calculation to determine the INT value for the original 
    left-most-mapped position and returns this corrected pos. If soft-clipping has not occurred, truePOS returns
    the original INT value of POS.
    INPUT: POS = A STR representing the left-most-mapped position of the read.
    CIGAR = A STR representing the CIGAR-string for the read.
    OUTPUT: pos = STRING: A soft-clip-corrected value for POS in reads where soft-clipping has occurred.
    POS = STRING: The original value for POS.
    """
    # Make POS an INT for calculations.
    POS = int(POS)
    # Search for left-most soft-clipping.
    sc = re.search("([0-9]+)S[0-9]+", CIGAR)
    # If it's found...
    if sc:
        # Set the difference to the INT value before the first S, for calculations.
        diff = int(sc.group(1))
        # Create the new corrected position, to check for PCR duplicates.
        pos = POS - diff
        # Return the corrected position.
        return str(pos)
    # If there's no soft-clipping...
    else:
        # Return the original position
        return str(POS)
    
def negPOS(POS,CIGAR):
    """
    description: negPOS checks if right-most-soft-clipping,deletions,Ns,or Ms have occurred in the CIGAR string.
    If right-most-soft-clipping,deletions,Ns,or Ms are present, negPOS adds the INT value to the original POS. 
    INPUT: POS = A STR representing the left-most-mapped position of the read.
    CIGAR = A STR representing the CIGAR-string for the read.
    OUTPUT: pos = A STR: A right-most value for POS in "-" strand reads.
    """
    # Make POS an INT for calculations.
    pos = int(POS)
    # Search for right-most soft-clipping.
    sc = re.search("([0-9]+)S$", CIGAR)
    # Search for mapped regions.
    mapped = re.findall("([0-9]+)M", CIGAR)
    # Search for deletions.
    dels = re.findall("([0-9]+)D", CIGAR)
    # Search for skipped regions.
    Ns = re.findall("([0-9]+)N", CIGAR)
    # If left-most soft-clipping is found...
    if sc:
        # Set the difference to the INT value before the first S, for calculations.
        diff = int(sc.group(1))
        # Create the new corrected position, to check for PCR duplicates.
        pos += diff
    elif len(dels) > 0:
        for n in dels:
            pos += int(n)
    elif len(Ns) > 0:
        for n in Ns:
            pos += int(n)
    for n in mapped:
        pos += int(n)
    return str(pos)
#####################
## ARGUMENT PARSER ##
#####################
# Create the argument parser.
parser = argparse.ArgumentParser(
        description=
        """
        INPUT:  A SAM file for deduplicating
                A file containing a list of UMIs
        OUTPUT: A deduplicated SAM file
                A text file containing basic STATS on deduplication
        """, formatter_class = RawTextHelpFormatter)
# Inform the parser of command line options.
# -h,--help is an argparse built-in function and will return the description w/ 
# options and option help messages seen below.
parser.add_argument('-f', '--file', help='Absolute SAM file path', 
                    type=str, required=True)
parser.add_argument('-p', '--paired', help='Designates file is paired-end',
                    action='store_true', required=False)
parser.add_argument('-u', '--umi', help='A file with a list of UMIs', 
                    type=str, required=False)
# Call the parser.
args = parser.parse_args()
# If -p/--paired is given on the command line, issue a technical-error message to st.out
if args.paired == True:
    msg = "\nTECHNICAL ERROR: This program does not currently handle paired-end reads.\n"
    print(msg)
    exit()
# If -u,--umi not specified, issue a technical-error message to st.out.
if args.umi == None:
    msg = """\nTECHNICAL ERROR: This program currently requires a file containing a list of UMIs.
    
    To specify your UMI filename use: -u,--umi <UMIfilename>
    \n"""
    print(msg)
    exit()
##################
## FILES PART 1 ##
##################
# Save the input file as fIN.
fIN = args.file
# Save the UMI list file as the variable umi.
umi = args.umi
# Split the SAMfile input at the .sam extention, store in a list.
filelist = fIN.split(".sam")
# Set up an input file, which is the output of the samtools sort command.
fin = filelist[0]+"_sorted.sam"
#######################
## UNIX COMMAND LINE ##
#######################
# Load requisite modules and run SAMtools sort on fIN.
check_output("module load samtools && module load python3/3.6.5 && \
            samtools sort -l 0 -m 4000M -o {} -T temp_sorted {}".format(fin,fIN),
            shell=True, timeout=None)   # ONLY ON TALAPAS
##################
## FILES PART 2 ##
##################
# Split the output file from samtools sort.
fOUTtemp = fin.split("_sorted.sam")
# Create an output file for this program.
fOUT = fOUTtemp[0]+"_deduped.sam"
# Create a statistical output text file.
stats = fOUTtemp[0]+"_STATS.txt"
################
## OPEN FILES ##
################
# Open all of the pertinent files.
with open(fin) as sam, open(umi) as u, open(fOUT,'w') as ded, open(stats, 'w') as pout:
    # Contain the unique molecular identifiers in a set.
    # Initiate an empty set for UMIs.
    UMIs = set()
    # Initiate an empty set for chromosomes "seen". 
    chromes = set()
    # Populate the set from the list of UMIs.
    for i in u:
        # Strip the UMI of invisible characters.
        i = i.strip()
        # Add the UMI to the set UMIs.
        UMIs.add(i)
    # Initiate Counters for the STATS output.
    # Count total duplicates.
    pcr_dup = 0
    # Count total reads.
    tot_rds = 0
    # Count unknown barcodes.
    unk_ids = 0
    # Iterate through each line of the input SAM file.
    for line in sam:
        # If the line is a header line...
        if re.match(r'^@[A-Z]{2}\s', line):
            # Write the header to fout.
            ded.write(line)
            # Continue iterating thru sam.
            continue
        # The line is not a header.
        else:
            # Increment the number of total reads.
            tot_rds += 1
            # Use parseSAM function to get the relevant SAM file fields to check for PCR duplicates.
            read_tup = parseSAM(line)
            # Store info in separate variables.
            this_UMI = read_tup[0]
            this_chr = read_tup[1]
            this_pos = read_tup[2]
            this_CIGAR = read_tup[3]
            this_strand = read_tup[4]
            # If current UMI is NOT IN UMIs...
            if this_UMI not in UMIs:
                # Increment the unknown umi counter.
                unk_ids += 1
                # Continue to the next line.
                continue
            # The UMI is in UMIs.
            else:
                # If this_chr is NOT in chromes...
                if this_chr not in chromes:
                    # Add the chromosome to the set, chromes.
                    chromes.add(this_chr)
                    # Reset a dictionary for + strands: key=UMI, value=set(POS1,POS2,POS3,...)
                    plusUMIdict = {}
                    # Reset a dictionary for - strands: key=UMI, value=set(POS1,POS2,POS3,...)
                    negUMIdict = {}
                    # If it's the positive strand.
                    if this_strand == "+":
                        # Find the true left-most mapped position.
                        new_pos = plusPOS(this_pos,this_CIGAR)
                        # Add the UMI and its set of positions to the plusUMIdict.
                        plusUMIdict[this_UMI] = set([new_pos])
                        # Write the line to ded.
                        ded.write(line)
                        # Continue to the next line.
                        continue
                    # It's the reverse strand...
                    else:
                        # Find the right-most read position.
                        new_pos = negPOS(this_pos,this_CIGAR)
                        # Add the UMI and its set of positions to the negUMIdict.
                        negUMIdict[this_UMI] = set([new_pos])
                        # Write the line to ded.
                        ded.write(line)
                        # Continue to the next line.
                        continue 
                # We're still on the same chromosome.
                else:
                    # If it's the positive strand.
                    if this_strand == "+":
                        # Find the true left-most mapped position.
                        new_pos = plusPOS(this_pos,this_CIGAR)
                        # Check if this_UMI is NOT in plusUMIdict...
                        if this_UMI not in plusUMIdict:
                            # Add the UMI as a KEY to the dictionary with the set of forward
                            # positions as the VALUE.
                            plusUMIdict[this_UMI] = set([new_pos])
                            # Write the line to ded.
                            ded.write(line)
                            # Continue to the next line.
                            continue
                        # We've already seen this UMI, we'll need to check for duplicates.
                        else:
                            # If the corrected position is in the set of positions for this UMI...
                            if new_pos in plusUMIdict[this_UMI]:
                                # Increment the PCR duplicate counter.
                                pcr_dup += 1
                                # Continue to the next line.
                                continue
                            # If the corrected position is not in the set of positions for this UMI...
                            else:
                                # Add the new position to this_UMI's set of positions, in
                                # plusUMIdict.
                                plusUMIdict[this_UMI].add(new_pos)
                                # Write the line to ded.
                                ded.write(line)
                                # Continue to the next line.
                                continue
                                
                    # It's the negative strand.
                    else:
                        # Find the right-most read position.
                        new_pos = negPOS(this_pos,this_CIGAR)
                        # Check if this_UMI is NOT in negUMIdict...
                        if this_UMI not in negUMIdict:
                            # Add the UMI as a KEY to the dictionary with the set of reverse 
                            # positions as the VALUE.
                            negUMIdict[this_UMI] = set([new_pos])
                            # Write the line to ded.
                            ded.write(line)
                            # Continue to the next line.
                            continue
                        # We've already seen this UMI, we'll need to check for duplicates.
                        else:
                            # If the corrected position is in this UMI's set of positions...
                            if new_pos in negUMIdict[this_UMI]:
                                # Increment the PCR duplicate counter.
                                pcr_dup += 1
                                # Continue to the next line.
                                continue
                            # The corrected position is not in this UMI's set of positions...
                            else:
                                # Add the corrected position to this UMI's set of positions.
                                negUMIdict[this_UMI].add(new_pos)
                                # Write the line to ded.
                                ded.write(line)
                                # Continue to the next line.
                                continue
    # Print output to the STATS file.
    print("Input SAM File:",fIN,file=pout)
    print("List of UMIs Used:",umi,file=pout)
    print("Deduplicated Output SAM File:",fOUT,file=pout)
    print("Total Reads:",tot_rds,file=pout)
    print("Total Duplicates:",pcr_dup,file=pout)
    print("Rate of PCR Duplication:",str(round(pcr_dup/tot_rds*100,2))+"%",file=pout)
    print("Total Unknown Barcodes:",unk_ids,file=pout)
    print("Rate of Unknown Barcodes:",str(round(unk_ids/tot_rds*100,2))+"%",file=pout)
# Create a template to remove the temporary samtools-sorted file.
rmcommand = "rm {}".format(fin)
call(rmcommand, shell=True)
