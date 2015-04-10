# vim: fdm=indent
'''
author:     Fabio Zanini
date:       09/04/15
content:    Call SNPs with pysam from BAM files.
'''
# Modules



# Functions
def get_bamfile_mode(bamfilename):
    '''Guess what kind of file it is (BAM/SAM)'''
    ext = bamfilename.split('.')[-1]
    if ext.lower() == 'sam':
        read_mode = 'r'
    elif ext.lower() == 'bam':
        read_mode = 'rb'
    else:
        raise ValueError('BAM/SAM file not recognized')

    return read_mode


def estimate_genome_length(bamfilename):
    '''Estimate genome length from BAM/SAM file'''
    import pysam

    read_mode = get_bamfile_mode(bamfilename)
    length = 0
    with pysam.Samfile(bamfilename, read_mode) as bamfile:
        for read in bamfile:
            try:
                length = max(length, read.pos + 2 * abs(read.isize))
            except TypeError:
                continue

    return length


def get_allele_counts_read(read,
                           counts_out,
                           qual_min=30,
                           length=None,
                           VERBOSE=0,
                           paired_end=False,
                           phred_offset=33,
                           alpha=['A', 'C', 'G', 'T', '-', 'N']):
    '''Get allele counts from a single read

    Parameters:
       read (pysam.AlignedSegment): the read to analyze
       counts_out (ndarray, alphabet x sequence length): output data structure for counts
       qual_min (int): minimal phred value to accept
       length (int): genome length
       VERBOSE (int): verbosity level (0-4)
       paired_end (bool): True for paired end BAM files (illumina)
       phred_offset (int): 8-bit offset of phred scores
       alpha (list/ndarray): genome alphabet
    '''
    import numpy as np

    if read.is_unmapped or (paired_end and (read.isize == 0)):
        if VERBOSE >= 3:
            print('Read '+read.qname+': unmapped')
        return

    if paired_end and (not read.is_proper_pair):
        if VERBOSE >= 3:
            print('Read '+read.qname+': unpaired')
        return

    # Read CIGARs
    seq = np.fromstring(read.seq, 'S1')
    qual = np.fromstring(read.qual, np.int8) - phred_offset
    pos = read.pos

    # Iterate over CIGARs
    for ic, (block_type, block_len) in enumerate(read.cigar):

        # Check for pos: it should never exceed the length of the fragment
        if (length is not None) and (block_type in [0, 1, 2]) and (pos >= length):
            raise ValueError('Pos exceeded the length of the fragment')
    
        # Inline block
        if block_type == 0:
            seqb = seq[:block_len]
            qualb = qual[:block_len]
            # Increment counts
            for j, a in enumerate(alpha):
                posa = ((seqb == a) & (qualb >= qual_min)).nonzero()[0]
                if len(posa):
                    counts_out[j, pos + posa] += 1
    
            # Chop off this block
            if ic != len(read.cigar) - 1:
                seq = seq[block_len:]
                qual = qual[block_len:]
                pos += block_len
    
        # Deletion
        elif block_type == 2:
            # Increment gap counts
            counts_out[4, pos:pos + block_len] += 1
    
            # Chop off pos, but not sequence
            pos += block_len
    
        # Insertion
        # an insert @ pos 391 means that seq[:391] is BEFORE the insert,
        # THEN the insert, FINALLY comes seq[391:]
        # We do not track them here
        elif block_type == 1:   
            # Chop off seq, but not pos
            if ic != len(read.cigar) - 1:
                seq = seq[block_len:]
                qual = qual[block_len:]
    
        # Other types of cigar?
        else:
            raise ValueError('CIGAR type '+str(block_type)+' not recognized')


def get_allele_counts_from_file(bamfilename,
                                genome_length=None,
                                qual_min=30,
                                maxreads=-1,
                                VERBOSE=0,
                                paired_end=False,
                                phred_offset=33,
                                alpha=['A', 'C', 'G', 'T', '-', 'N']):
    '''Get the allele counts from a BAM file
    
    Parameters:
       bamfilename (str): path of the BAM/SAM file
       genome_length (int): genome length
       qual_min (int): minimal phred value to accept
       maxreads (int): number of reads to scan
       VERBOSE (int): verbosity level (0-4)
       paired_end (bool): True for paired end BAM files (illumina)
       phred_offset (int): 8-bit offset of phred scores
       alpha (list/ndarray): genome alphabet
    '''
    import numpy as np
    import pysam

    # Estimate genome length if necessary
    if genome_length is None:
        genome_length = estimate_genome_length(bamfilename)
        if VERBOSE >= 2:
            print('Genome length:', genome_length)

    # Make alpha a ndarray
    alpha = np.asarray(alpha, 'S1')

    # Prepare output structures
    counts = np.zeros((len(alpha), genome_length), int)

    # Open BAM/SAM file
    read_mode = get_bamfile_mode(bamfilename)
    with pysam.Samfile(bamfilename, read_mode) as bamfile:

        # Iterate over single reads
        for i, read in enumerate(bamfile):

            # Max number of reads
            if i == maxreads:
                if VERBOSE >= 2:
                    print('Max reads reached:', maxreads)
                break
        
            # Print output
            if (VERBOSE >= 3) and (not ((i +1) % 10000)):
                print(i+1)

            get_allele_counts_read(read, counts,
                                   length=genome_length,
                                   qual_min=qual_min,
                                   paired_end=paired_end,
                                   phred_offset=phred_offset,
                                   VERBOSE=VERBOSE,
                                   alpha=alpha)

    return counts



# Script
if __name__ == '__main__':

    # Parser for command-line args
    import argparse
    parser = argparse.ArgumentParser(description='Call SNPs from BAM/SAM file',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('bamfile',
                        help='path to the bamfile')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Number of reads to scan (for testing)')
    parser.add_argument('--quality', type=int, default=30,
                        help='Minimal phred quality to accept')
    parser.add_argument('--genome-length', type=int,
                        help='Genome length')
    parser.add_argument('--paired-end', action='store_true',
                        help='Set for paired-end sequencing (illumina)')


    args = parser.parse_args()
    bamfilename = args.bamfile
    VERBOSE = args.verbose
    maxreads = args.maxreads
    qual_min = args.quality
    genome_length = args.genome_length
    paired_end = args.paired_end
    
    counts = get_allele_counts_from_file(bamfilename,
                                         genome_length=genome_length,
                                         qual_min=qual_min,
                                         maxreads=maxreads,
                                         paired_end=paired_end,
                                         VERBOSE=VERBOSE)

