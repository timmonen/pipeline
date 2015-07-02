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
