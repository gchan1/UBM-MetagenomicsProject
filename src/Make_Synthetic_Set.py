#File: Make_Synthetic_Set.py
#Author: Grace Chandler
#Last edit: 18-07-14
#Description: This program takes in a set of known binding sites/motifs and reconstructs a snthetic data set from it

import time
import sys
import numpy as np
import os.path
import random

#This function takes in the lines from a file of aligned motifs and returns a list of the frequencies of each base. The order of the list is ATGC
def Get_Base_Frequencies(lines, base_frequencies):


#Train probabilities from strand
    a_count = 0
    t_count = 0
    c_count = 0
    g_count = 0
    base_count = 0



    #Iterate through each motif and count the total number of ATGC in the entire set of aligned motifs
    for motif in lines:
        #print motif 

        for base in motif:
    
            if base == 'T' or base == 't':
                t_count += 1
            elif base == 'G' or base == 'g':
                g_count += 1
            elif base == 'A' or base == 'a':
                a_count += 1
            elif base == 'C' or base == 'c':
                c_count += 1

            base_count += 1

    #print a_count
    #print t_count
    #print c_count
    #print g_count
    #print base_count

    base_count = float(base_count)
    
    #add the frequencies to the list
    base_frequencies.append(float(float(a_count)/base_count))
    base_frequencies.append(float(float(t_count)/base_count))
    base_frequencies.append(float(float(g_count)/base_count))
    base_frequencies.append(float(float(c_count)/base_count))

    print base_frequencies


def Logical(length, sites, lines, loci):
    
    #print 'length at beginning of logical ', length
    
    base_frq = []
    
    #Keep track of time to compare it to the Logical method
    start_time = time.time()

    Get_Base_Frequencies(lines, base_frq)

    print base_frq

    nucleotides = 'ATGC'
    selection_pool = []
    count = 0
    bs_len = len(lines[0])
    binding_sites = []
    
    
    #Make a selection pool to draw nucleotides from
    #This creates the means by which to obtain nucleotides with the proper probability
    for frequency in base_frq:
        frequency = int(frequency*1000)
        for i in range(0 , frequency):
            selection_pool.append(nucleotides[count])
        count += 1
    
    count = 0

    #This portion of the code makes binding sites based on the number of desired 
    for i in range(sites):
        temp_site = []
        while count < bs_len:
            index = random.randrange(len(selection_pool))
            #print 'index' , index
            nucleotide = selection_pool[index]
            temp_site.append(nucleotide)
            count += 1
        count = 0
        
        #make that puppy into a string for safe keeping
        temp_site = ''.join(temp_site)
        print temp_site
        binding_sites.append(temp_site)


    count = 0

    #Initializing a lot of lists
    #Maybe somehow find a better way to to this
    non_coding_nucleotides = []
    coding_nucleotides = []
    coding_region = []
    non_coding_region = [] 
    
    #Make a set of random coding frequency nucleotides to randomly select from
    while count < 4:
        for i in range(0,25):
            coding_nucleotides.append(nucleotides[count])
        
        if count < 2:
            for i in range(0,40):
                non_coding_nucleotides.append(nucleotides[count])
        if count > 1:
            for i in range(0,60):
                non_coding_nucleotides.append(nucleotides[count])
       
        count += 1
        
    count = 0   
    
   
    #assign 50 random nucleotides in the proportion for coding strand
    #use .25 for the coding strand
    for i in range(0,50):
        coding_region.append(coding_nucleotides[random.randrange(len(coding_nucleotides))])
    
  
    #how many binding sites should be assigned to the noncoding region
    leftovers = (length - (bs_len * sites))
   
   
    for i in range(0, leftovers):
        non_coding_region.append(non_coding_nucleotides[random.randrange(len(non_coding_nucleotides))])
        
    #This is the part of the show where Grace comes out and does a lot of number checking 
    print 'leftovers ', leftovers
    print 'length desired', length
    print 'num sites', sites
    
    
    print 'bs_len ', bs_len  
    
    print 'length of non_Coding ', len(non_coding_region)
    print 'length of Coding ', len(coding_region)  

    print 'total length should be 350'
    print 'total length is: ', (len(non_coding_region) + (sites * bs_len) + len(coding_region))
    
   #Now we insert the bs into the list of noncoding nucleotides!
    for site in binding_sites:
        locus = random.randrange(len(non_coding_region))
        
        #add this to the list of loci
        #also add the site for reference of what is at that loci
        loci.append(locus)
        loci.append(site)
        
        #add the binding site to the noncoding region
        non_coding_region.insert(locus, site)
        
    #add the 50 bp coding region to the end of the strand
    non_coding_region.extend(coding_region)
   
    set =  ''.join(non_coding_region)
            
    return set


    #calculate method runtime
    print "Time: " , time.time() - start_time , " seconds "

def Make_Synthetic_Set():
    
    #initialize the loci list
    loci = []

    total_start = time.time()

    outFile = open('Synth_Sets.txt', 'w')
    
    
    print "TF motif options: LexA, Fur, CodY "
    
    input_check = False
    while input_check != True:
        TF_input = raw_input("Please select TF motif: ")
        if TF_input in ['LexA', 'Fur', 'CodY']:
            input_check = True
        else:
            print "Incorrect input: please enter LexA, Fur or CodY"

        
    if input_check == True:
        

        file_path = os.path.join('..', 'data' , 'Aligned_motifs','%s.txt' %TF_input)
        f = open(file_path, "r")
        f.seek(0)
        lines = f.readlines()
        n = len(lines)
        
        strand_length = int(raw_input("Enter strand length: "))
        num_sites = int(raw_input("Enter number of sites to place in strand: "))
        num_sets = int(raw_input("Enter desired number of sets: "))
        set_count = 0
        while set_count < num_sets:


#print out which set you're using
            print "Using: %s" %TF_input

#determine number of sites to enrich strand with

            synthetic_set = Logical(strand_length, num_sites, lines, loci)
            print "length ", len(synthetic_set)
            print synthetic_set
            outFile.write(synthetic_set)
            outFile.write("\n")
                
            set_count += 1

#array of nucleotides

#randomly select one each time

#make biniding site a certain bp
            f.close()
#return the dataset
        outFile.close()
    
    print "Total time: " , (time.time() - total_start)

Make_Synthetic_Set()
