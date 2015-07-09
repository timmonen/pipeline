# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 14:42:50 2015

@author: timmonen
"""

#####################33
'''


Iterate over sam-file in read pairs
-------------------------------------
0.a  Get the reference positions of each (fwd and bwd) primer (only once) 
0.b Make dictionary of fragment and beginning and end positions

Preliminary clean-up of bad reads
Track overall number thrown out due to prelim cleaning

1. If reads unmapped, trash
    track number unmapped
2. If reads not proper pair, trash
    track number not proper pair
3. If insert size <300 or > 750, trash
    track number too small, number too big
3. See which fragment(s) insert belongs to: compare first and last position of 
    read to span of each fragment (if in array?)
    
    def is_slice_in_list(s,l):
    len_s = len(s) #so we don't recompute length of s on every iteration
    return any(s == l[i:len_s+i] for i in xrange(len(l) - len_s+1))
    
4. If insert in zero or more than one fragments, trash
    track how many thrown out due to 0 or more than 2 fragments
    keep how many in 3 fragments (should be none! sanity check)
    
 ---------------------------------------------------------------------   
We have exclusively mapped, paired, reasonably-sized inserts which belong to
exactly one fragment.

Trim reads in regions with primer overlap using pre-mapping information
--------------------------------------------------------------------------

    We know from pre-mapping the position of each read wrt to same reference:
    We know the length of each primer
    Say primer starts at position 2300 and ends at position 2320
    If read starts at position 2315, trim beginning to start at position 2021
    
    Since accept only insert sizes > 350, only need to worry about forward
    primers for forward reads, and backward primers for backward reads
    
    Reverse read will never reach into forward primer and vice versa
    
    If we need to accept insert sizes that are smaller, need to check forward
    primer agaist reverse read, and reverse read against forwards primer 
    because reads could start looping back
---------------------------------------------------------------------------

We already know which fragment each reads belongs to, so sort the trimmed reads
into 6 different bam/sam files with fragment name
'''
from pipeline.sample import Sample
import pysam 
import sys
import os
import subprocess as sp
import argparse
import numpy as np


def pair_generator(iterable):
    '''Generator for pairs in interleaved files, such as BAM files'''
    # Note: the last item is lost if odd
    it = iter(iterable)
    while True:
        try:
            p1 = it.next()
            p2 = it.next()
            yield (p1, p2)
        except StopIteration:
            raise    
            
        
def assign_to_fragment(reads, fragment_full, VERBOSE=0, **kwargs):
    '''Assign read pair to fragments'''
    i_fwd = reads[0].is_reverse
    # Insert coordinates
    ins_start = reads[i_fwd].pos
    ins_end = reads[i_fwd].pos + reads[i_fwd].isize

    # What fragments could the read pair come from?
    frags_pot = []
    for n_frag, (fr_start, fr_end) in enumerate(fragment_full.values()):
        if np.isscalar(fr_start):
            frag_ind = (ins_start >= fr_start) and (ins_end <= fr_end)
            if frag_ind:
                frags_pot.append(str(fragment_full.keys()[n_frag]))

    # If no fragments are compatible, it's cross-boundary (PCR crazyness)
    if len(frags_pot) == 0:
        pair_identity = 'Lost'

    # If it is only compatible with one primer, it's ok
    elif len(frags_pot) == 1:
        pair_identity = frags_pot[0]
        # If it is only compatible with one primer, it's ok
    elif len(frags_pot) > 1:
        pair_identity = "Ambiguous"
    if VERBOSE > 0:
        print fragment_full, ins_start, ins_end, pair_identity
    return pair_identity
      

def trim_primers(reads, fragment, fragments_trim, **kwargs):
    i_fwd = reads[0].is_reverse
    i_rev = not i_fwd
    tampered = "Intact"
    frag_pos = fragments_trim[fragment]
        # FWD primer
    if reads[i_fwd].pos < frag_pos[0]:
            tampered = "Trimmed"
            ref_pos = reads[i_fwd].pos
            read_pos = 0
            cigar = reads[i_fwd].cigar[::-1]
            for i, (bt, bl) in enumerate(reads[i_fwd].cigar):
                # Match
                if bt == 0:
                    if ref_pos + bl > frag_pos[0]:
                        cut = frag_pos[0] - ref_pos
                        cigar[-1] = (bt, bl - cut)
                        read_pos += cut
                        ref_pos = frag_pos[0]
                        break
                    cigar.pop(-1)
                    read_pos += bl
                    ref_pos += bl
                # Insertion
                elif bt == 1:
                    # Move up in read, nothing to move up in reference
                    cigar.pop(-1)
                    read_pos += bl
                # Deletion
                elif bt == 2:
                    # Starting with a deletion is not allowed
                    cigar.pop(-1)
                    # Move up in reference, nothing to move up in read
                    ref_pos += bl
                    if ref_pos > frag_pos[0]:
                        break
            cigar = cigar[::-1]
            # If you cut away everything, trash
            if not len(cigar):
                return "Trash"
                
           # Alter the read!!!
            seq = reads[i_fwd].seq
            qual = reads[i_fwd].qual
            reads[i_fwd].pos = ref_pos
            reads[i_fwd].seq = seq[read_pos:]
            reads[i_fwd].qual = qual[read_pos:]
            reads[i_fwd].cigar = cigar

        # REV primer
        # Get end position of read by summing up over cigars (only substitutions, dels, inserts)
    ref_pos = reads[i_rev].pos + sum(bl for (bt, bl) in reads[i_rev].cigar if bt in (0, 2))
    if ref_pos > frag_pos[1]:
        tampered = "Trimmed"
            # Get the last position in the read, trim backwards
        read_pos = reads[i_rev].rlen
        cigar = reads[i_rev].cigar
            # Go through all cigar statuses cutting them off as needed to match
        for i, (bt, bl) in enumerate(reads[i_rev].cigar[::-1]):
                # If substition
            if bt == 0:
                if ref_pos - bl < frag_pos[1]:
                    cut = ref_pos - frag_pos[1]
                    cigar[-1] = (bt, bl - cut)
                    read_pos -= ref_pos - frag_pos[1]
                    break
                cigar.pop(-1)
                read_pos -= bl
                ref_pos -= bl
                # If insertion
            elif bt == 1:
                cigar.pop(-1)
                read_pos -= bl
                # If deletion
            elif bt == 2:
                # Ending with a deletion is not allowed
                cigar.pop(-1)
                ref_pos -= bl
                if ref_pos < frag_pos[1]:
                    break

            # If you cut away everything, trash
            if not len(cigar):
                return "Trash"

            seq = reads[i_rev].seq
            qual = reads[i_rev].qual
            reads[i_rev].seq = seq[:read_pos]
            reads[i_rev].qual = qual[:read_pos]
            reads[i_rev].cigar = cigar

    # Fix mate pair
    if tampered == "Trimmed":
        reads[i_fwd].mpos = reads[i_rev].pos
        reads[i_rev].mpos = reads[i_fwd].pos
        isize = reads[i_rev].pos + sum(bl for bt, bl in reads[i_rev].cigar
                                       if bt in (0, 2)) - reads[i_fwd].pos
        reads[i_fwd].isize = isize
        reads[i_rev].isize = -isize

    return tampered


def trim_and_divide(sample, min_isize = 350, max_isize = 750, VERBOSE = 0,keep_trash = 1, **kwargs):
    

    fragments = sample.fragment_names
    fragment_full = sample.get_fragment_positions()
    fragments_trim = sample.get_fragment_trim_positions()    
        
    primers = sample.primer_dict
    primer_positions = sample.get_primer_positions()
        
    data_folder = Sample.get_data_foldername()
    input_filename = sample.get_pre_map_filename()['data']
        
    fragment_output_names = sample.get_fragment_output_names()
    trashed_read_names = sample.get_trashed_read_names()
    fn_outsum = sample.get_divide_summary_filenames()
  
    
    # Are reads ok?      
    unmapped = 0
    too_small = 0
    too_large = 0
    # Are fragments identified?
    frag_found = 0
    ambiguous = 0
    lost = 0
    # Do we trim primers or not?
    bad_primer = 0
    trimmed_primers = 0
    intact_primers = 0
    # How many of each type of fragment
    frag1 = 0
    frag2 = 0
    frag3 = 0
    frag4 = 0
    frag5 = 0
    frag6 = 0
    
    
    with pysam.Samfile(input_filename, 'r') as samfile:
        try:
            file_handles = [pysam.Samfile(ofn, 'w', template=samfile) 
                            for ofn in fragment_output_names[:len(fragments)]]
        
            if keep_trash == 1:
                trash_handles = [pysam.Samfile(otr, 'w', template = samfile)
                                for otr in trashed_read_names[:len(trashed_read_names)]]
    
            for irp, reads in enumerate(pair_generator(samfile)):
                trash = 0
                i_fwd = reads[0].is_reverse
    
                if (reads[0].is_unmapped) or (reads[1].is_unmapped) or (not reads[0].is_proper_pair) or (not reads[1].is_proper_pair):
                    unmapped += 1
                    if keep_trash == 1:
                        trash_handles[0].write(reads[0])
                        trash_handles[0].write(reads[1])
                    trash = 1
                    continue 
                        
                if (reads[i_fwd].isize < min_isize):
                    too_small += 1
                    if keep_trash == 1:
                        trash_handles[1].write(reads[0])
                        trash_handles[1].write(reads[1])
                    trash = 1
                    continue
                             
                if (reads[i_fwd].isize > max_isize):
                    too_large += 1   
                    if keep_trash == 1:
                        trash_handles[2].write(reads[0])
                        trash_handles[2].write(reads[1])
                    trash = 1
                    continue
                        
                if trash == 0:
                    fragment = assign_to_fragment(reads, fragment_full, VERBOSE = 0)
                    if fragment == "Ambiguous":
                        ambiguous += 1
                        if keep_trash == 1:
                            trash_handles[4].write(reads[0])
                            trash_handles[4].write(reads[1])
                    if fragment == "Lost":
                        lost += 1
                        if keep_trash == 1:
                            trash_handles[5].write(reads[0])
                            trash_handles[5].write(reads[1])
                    if fragment in fragment_full.keys():
                        frag_found += 1
                        n_frag = fragments.index(fragment)
                        primer_status = trim_primers(reads, fragment, fragments_trim)
                        if primer_status == "Trash":
                            bad_primer += 1
                            if keep_trash == 1:
                                trash_handles[3].write(reads[0])
                                trash_handles[3].write(reads[1])
                        if primer_status == "Trimmed":
                            trimmed_primers += 1
                            file_handles[n_frag].write(reads[0])
                            file_handles[n_frag].write(reads[1])
                        if primer_status == "Intact":
                            intact_primers += 1
                            file_handles[n_frag].write(reads[0])
                            file_handles[n_frag].write(reads[1])
                        if (primer_status == "Trimmed") or (primer_status == "Intact"):
                            if fragment == "F1":
                                frag1 += 1
                            if fragment == "F2":
                                frag2 += 1
                            if fragment == "F3":
                                frag3 += 1
                            if fragment == "F4":
                                frag4 += 1
                            if fragment == "F5":
                                frag5 += 1
                            if fragment == "F6":
                                frag6 += 1                            
        finally:
            for f in file_handles:
                f.close()
            if keep_trash == 1:
                for t in trash_handles:
                    t.close()
      
    summary = {'sample name': sample.name,
                   'sequencing run': sample.run,
                   'read pairs': irp + 1,
                   'unmapped reads': unmapped,
                   'too short inserts': too_small,
                   'too long inserts': too_large,
                   'ambigous fragments': ambiguous,
                   'non-assigned fragments': lost,
                   'matching fragments': frag_found,
                   'reads with untrimmed primers': intact_primers,
                   'reads with trimmed primers': trimmed_primers,
                   'reads lost due to trimming primers': bad_primer,
                   'F1': frag1,
                   'F2': frag2,
                   'F3': frag3,
                   'F4': frag4,
                   'F5': frag5,
                   'F6': frag6,
                  }    
    sample.write_json(summary, fn_outsum)
    
  # Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Divide into fragments and trip primers of reads')
    parser.add_argument('--sample', required=True,
                        help='MiSeq sample to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--min_isize', type=int, default = 350, help='Minimum insert size required')
    parser.add_argument('--max_isize', type=int, default = 750, help='Maximum insert size allowed')
    parser.add_argument('--keep_trash',type=int, default = 1, help = 'Keep trashed reads in separate files')
    args = parser.parse_args()

    sample = Sample(args.sample)
    sample.trim_and_divide(VERBOSE=args.verbose,min_isize = args.min_isize, max_isize = args.max_isize, keep_trash = args.keep_trash)
                    


