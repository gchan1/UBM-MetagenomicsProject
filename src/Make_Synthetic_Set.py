#File: Make_Synthetic_Set.py
#Author: Grace Chandler
#Last edit: 17-06-14
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


#This function utilizes the random.multinomial method to produce a synthetic dataset using the desired length and number of sites
def Multinomial(length, sites, lines):
    base_frequencies = []

    #Keep track of time to compare it to the Logical method
    start_time = time.time()

    base_frq = Get_Base_Frequencies(lines, base_frequencies)

    bs_len = len(lines[0])
    binding_sites = np.random.multinomial(bs_len, [base_frq]*4, size = sites)
    
    #shuffle the binding sites

    for site in binding_sites:
        #MAKE IT A LIST
        temp_list = []
        count = 0
        nucleotides = 'ATGC'

        for number  in site:
            print number 
            while count < number:
                temp_list.append(nucleotides[count])
            count += 1

        site = temp_list
        random.shuffle(site)
        site  = ''.join(site)

    #Create the other, non site bases
    #MAKE IT A LIST
    temp_list = []
    count = 0
    non_binding = np.random.multinomial(length - (bs_len * sites), [.25]*4, size = 1)
    for number in non_binding:
        while count < number:
            temp_list.append(nucleotides[count])
        count += 1

    non_binding = temp_list
    random.shuffle(non_binding)
    
            
    print binding_sites

    #Create the synthetic set

    for site in binding_sites:
        locus = random.range(len(non_binding))
        non_binding.insert(locus, site)


    print ''.join(non_binding)
    
    #calculate method runtime

    print "Time: " , time.time() - start_time , " seconds "


def Logical(length, sites, lines):
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
    
    for frequency in base_frq:
        frequency = int(frequency*1000)
        for i in range(0 , frequency):
            selection_pool.append(nucleotides[count])
        count += 1
    
    count = 0

    print 'length', len(selection_pool)

    for i in range(sites):
        temp_site = []
        while count < bs_len:
            index = random.randrange(len(selection_pool))
            #print 'index' , index
            nucleotide = selection_pool[index]
            temp_site.append(nucleotide)
            count += 1
        count = 0
        temp_site = ''.join(temp_site)
        print temp_site
        binding_sites.append(temp_site)


    count = 0
    random_nucleotides = []

    leftovers = int(int(length - (bs_len * sites)) / 4)
    print 'leftovers ', leftovers

    while count < 4:
        for i in range(0, leftovers):
            random_nucleotides.append(nucleotides[count])
        count += 1
        
    random.shuffle(random_nucleotides)

    for site in binding_sites:
        locus = random.randrange(len(random_nucleotides))
        random_nucleotides.insert(locus, site)


    set =  ''.join(random_nucleotides)
            
    return set


    #calculate method runtime
    print "Time: " , time.time() - start_time , " seconds "

def Make_Synthetic_Set():

    total_start = time.time()

    outFile = open('Synth_Sets.txt', 'w')
    
    
    print "TF motif options: LexA, Fur, CodY "
    
    input_check = False
    if input_check != True:
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
        prob_method = raw_input("Select method: enter M for multinomial, L for logistical ")
        num_sets = input("Enter desired number of sets: ")
        set_count = 0
        while set_count < num_sets:


#print out which set you're using
            print "Using: %s" %TF_input

#determine number of sites to enrich strand with

            if prob_method == 'M':
                synthetic_set = Multinomial(strand_length, num_sites, lines)

            elif prob_method == 'L':
                synthetic_set = Logical(strand_length, num_sites, lines)
                print "length ", len(synthetic_set)
                print synthetic_set
                outFile.write(synthetic_set)
                outFile.write("\n")
                
            set_count += 1
#Multinomial probability set


#array of nucleotides

#randomly select one each time

#make biniding site a certain bp
            f.close()
#return the dataset
        outFile.close()
    
    print "Total time: " , (time.time() - total_start)

Make_Synthetic_Set()
