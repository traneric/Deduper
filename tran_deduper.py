#!/usr/bin/env python
import re
import argparse

def argument_parser():
    '''Return parser object.'''
    # Create parser object.
    parser = argparse.ArgumentParser(description='Given a SAM file of uniquely \
        mapped reads, this program removes all PCR duplicates (retain only a single \
        copy of each read). The SAM file must be sorted by chromosome number \
        prior to running the program. The program will output the unique reads, \
        PCR duplicate reads, and mismatched index reads in three files: \
        deduped.sam, duplicates.sam, and misindexed.sam')

    # Create flags for SAM file, paired-end option, and UMI file.
    parser.add_argument('-f', '--file', type=str, required=True, help='Absolute file path \
        for SAM file')
    parser.add_argument('-p', '--paired', action='store_true', help='Designates file is paired end')
    parser.add_argument('-u', '--umi', type=str, required=False, help='Designates file \
        containing the list of UMIs (unset if randomers instead of UMIs)')
    return parser.parse_args()

def is_unique_umi(umi, umi_list):
    '''This function returns True if the read has a unique UMI. Otherwise 
    it returns False.'''
    if umi in umi_list:
        return True
    else:
        return False

def is_reverse_complement(read):
    '''This function returns True if the read is a reverse complement. Otherwise
    it returns False.'''
    strand_integer = int(read.split("\t")[1])
    binaryNumber = '{0:05b}'.format(strand_integer)
    strand_bit_flag = binaryNumber[0]
    
    if strand_bit_flag == 1:
        return True
    else:
        return False

def is_soft_clipped(cigar_string):
    '''This function returns True if the read is soft clipped. Otherwise it
    returns False.'''
    if 'S' in cigar_string:
        return True
    else:
        return False

def is_matched(cigar_string):
    '''This function returns True if the read has matches. Otherwise it
    returns False.'''
    if 'M' in cigar_string:
        return True
    else:
        return False
    
def is_deleted(cigar_string):
    '''This function returns True if the read has deletions. Otherwise it
    returns False.'''
    if 'D' in cigar_string:
        return True
    else:
        return False

def is_skipped(cigar_string):
    '''This function returns True if the read is skipped. Otherwise it
    returns False.'''
    if 'N' in cigar_string:
        return True
    else:
        return False

def get_umi(read):
    '''This function accepts a read and returns its UMI value.'''
    umi = read.split("\t")[0][-8:]
    return umi

def get_cigar_string(read):
    '''This function returns a cigar string from a read'''
    return read.split("\t")[5]

def get_chromosome_number(read):
    '''This functions returns the read's chromosome number.'''
    chromosome_number = read.split("\t")[2]
    return chromosome_number

def get_start_position(read):
    '''This function returns the read's start position.'''
    start_position = int(read.split("\t")[3])
    return start_position

def adjust_start_position(read, start_position, reverse_complement):
    '''This function adjusts the start position of the read depending
    on its strandedness. Returns the adjusted start position.'''

    cigar_string = get_cigar_string(read)
    soft_clipped = is_soft_clipped(cigar_string)
    adjusted_start_position = start_position

    if reverse_complement == False and soft_clipped:
        # Subtract soft clipping from the forward read original position.
        soft_clipped_value = int((re.search("(\d+)S", cigar_string).group())[:-1])
        adjusted_start_position = start_position - soft_clipped_value
    elif reverse_complement == True:
        # The read is a reverse complement.
        # Add matches, deletions, skips. Substract end soft clipping.
        matched_value = 0
        deleted_value = 0
        skipped_value = 0
        soft_clipped_value = 0

        if is_matched(cigar_string):
            matched_value = sum(list(map(int, re.findall("(\d+)M", cigar_string))))
        if is_deleted(cigar_string):
            deleted_value = int((re.search("(\d+)D", cigar_string).group())[:-1])
        if is_skipped(cigar_string):
            skipped_value = int((re.search("(\d+)N", cigar_string).group())[:-1])
        if is_soft_clipped(cigar_string):
            # Search for right most soft clipping value.
            if cigar_string[-1] == 'S':
                soft_clipped_value = list(map(int, re.findall("(\d+)S", cigar_string)))[-1]
        adjusted_start_position = start_position + matched_value + deleted_value + skipped_value + soft_clipped_value
    return adjusted_start_position

def display_final_results(misindexed, duplicate, unique):
    '''This function displays the number of unique reads, duplicate reads, and
    misindexed reads.'''
    print("Misindexed reads: " + str(misindexed))
    print("Duplicate reads: " + str(duplicate))
    print("Unique reads: " + str(unique))

def close_files(samFile, outputSamFile, duplicateFile, misindexedFile):
    '''This function closes all files.'''
    samFile.close()
    outputSamFile.close()
    duplicateFile.close()
    misindexedFile.close()

def main():
    '''This is the main function.'''
    # Retrieve parser object.
    args = argument_parser()

    # Variables to store command-line inputs.
    absolute_file_path = args.file
    paired_end = args.paired
    umi_file_path = args.umi

    # Raise exception if the user specifies paired-end.
    if paired_end:
        raise Exception("This program does not support paired end functionality. The program will now exit. Please try again without the -p (paired-end) option.")

    # Generate a list of the 96 unique molecular indices.
    umi_list = open('STL96.txt').read().split('\n')

    # This list retains only a single copy of each read.
    non_pcr_duplicate_set = set()

    # Variables to track misindexed, duplicate, unique reads, and previous chromosome number.
    misindexed_read_count = 0
    duplicate_read_count = 0
    unique_read_count = 0
    previous_chromosome_number = "undeclared"

    # Open all the files.
    samFile = open(absolute_file_path, 'r')
    outputSamFile = open('deduped.sam', 'w')
    duplicateFile = open('duplicates.sam', 'w')
    misindexedFile = open('misindexed.sam', 'w')

    # Process the file.
    for line in samFile:
        read = line.strip()

        # Get UMI for the current read.
        umi = get_umi(read)

        # Continue to the next read if umi is not one of 96 unique UMI.
        if is_unique_umi(umi, umi_list) == False:
            misindexed_read_count += 1
            continue

        # Get variables to determine PCR duplicates.
        chromosome_number = get_chromosome_number(read)
        start_position = get_start_position(read)
        reverse_complement = is_reverse_complement(read)
        soft_clipped = is_soft_clipped(get_cigar_string(read))

        # Check if the tuple set needs to be purged.
        if previous_chromosome_number == "undeclared":
            previous_chromosome_number = chromosome_number
        elif chromosome_number != previous_chromosome_number:
            previous_chromosome_number = chromosome_number
            non_pcr_duplicate_set.clear()
        
        # Adjust the start position of the read.
        start_position = adjust_start_position(read, start_position, reverse_complement)

        # Store the current read information into a tuple.
        read_tuple = (umi, chromosome_number, start_position, reverse_complement)

        # Check for PCR duplicate and append tuple to list if the read is not a duplicate.
        if read_tuple not in non_pcr_duplicate_set:
            non_pcr_duplicate_set.add(read_tuple)
            outputSamFile.write(read + '\n')
            unique_read_count += 1
        else:
            duplicateFile.write(read + '\n')
            duplicate_read_count += 1

    close_files(samFile, outputSamFile, duplicateFile, misindexedFile)
    display_final_results(misindexed_read_count, duplicate_read_count, unique_read_count)

# Run program
main()
