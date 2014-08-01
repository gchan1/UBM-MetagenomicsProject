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

    freq_matrix = []

#Train probabilities from strand
    a_count = 0
    t_count = 0
    c_count = 0
    g_count = 0
    base_count = 0


    bs_len = len(lines[0])
    #Iterate through each motif and count the total number of ATGC in the entire set of aligned motifs
    
    #base_count = float(base_count)
    
    #add the frequencies to the list
   

   # print return frq_matrix

def transpose(arrays):
    return zip(*arrays)

def Logical(length, sites, lines, loci_list, loci, Synth_Set_Stats):
    
    #print 'length at beginning of logical ', length
    
    #base_frq = []
    
    #Keep track of time to compare it to the Logical method
    start_time = time.time()

    temp = []
    for line in lines:
        line = line.strip()
        line = list(line)
        if len(line) == 43:
            temp.append(line)
    
    lines = temp
        
            
            
    print 'len of lines: ', len(lines)
    print 'len of one line: ', len(lines[0])
    #Get_Base_Frequencies(lines, base_frq)
    cols = zip(*lines)
    
    #print base_frq

    nucleotides = 'ATGC'
    selection_pool = []
    count = 0
    bs_len = len(lines[0])
    binding_sites = []
    
    
    #Make a selection pool to draw nucleotides from
    #This creates the means by which to obtain nucleotides with the proper probability
   
    count = 0

    print cols
    print 'len: cols' , len(cols)

    #This portion of the code makes binding sites based on the number of desired 
    for i in range(sites):
        temp = "".join(random.choice(col) for col in cols)
        #print temp
        print 'added a BS'
        print 'len temp', len(temp)
        binding_sites.append(temp)
        
        #make that puppy into a string for safe keeping
        #temp_site = ''.join(temp_site)
       # print temp_site
       # binding_sites.append(temp_site)


    count = 0

    #print 'binding_sites =  ' , binding_sites

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
    print 'leftovers: ' , leftovers   
   
    for i in range(0, leftovers):
        non_coding_region.append(non_coding_nucleotides[random.randrange(len(non_coding_nucleotides))])
        
    #This is the part of the show where Grace comes out and does a lot of number checking 
    #print 'leftovers ', leftovers
    #print 'length desired', length
    #print 'num sites', sites
    
    
    print 'bs_len ', bs_len  
    
    print 'length of non_Coding ', len(non_coding_region)
    print 'length of Coding ', len(coding_region)  

    #print 'total length should be 350'
    print 'total length is: ', (len(non_coding_region) + (sites * bs_len) + len(coding_region))
    
    
    
    
   
   #Now we insert the bs into the list of noncoding nucleotides!
    for i in range(sites):

        site = binding_sites[i]
        
        locus = loci_list[i]
        
        #add this to the list of loci
        #also add the site for reference of what is at that loci
        loci.append(locus)
        loci.append(site)
        
        #add the binding site to the noncoding region
        non_coding_region.insert(locus, site)
        print "inserted site"

    print "scary bs len" , len(binding_sites[0])
    #print "noncoding region: " , non_coding_region
    print "coding region: " , coding_region
    #add the 50 bp coding region to the end of the strand
    non_coding_region.extend(coding_region)
   
    set =  ''.join(non_coding_region)

    #print 'set: ' , set
    print 'LEN OF ONE SET', len(set)
    return set


    #calculate method runtime
    print "Time: " , time.time() - start_time , " seconds "

def Make_Synthetic_Set():
    
    loci = []

    total_start = time.time()

    outFile = open('Synth_Sets.txt', 'w')
    #BS_File = open('Binding_Sites.txt', 'w')
    Synth_Set_Stats = open('Synth_Set_Stats.txt', 'w')
    
    print "TF motif options: LexA, Fur, CodY "
    
    input_check = False
    while input_check != True:
        TF_input = raw_input("Please select TF motif: ")
        if TF_input in ['LexA', 'Fur', 'CodY']:
            input_check = True
        else:
            print "Incorrect input: please enter LexA, Fur or CodY"

        
    if input_check == True:
        
        Synth_Set_Stats.write('TF motif')
        Synth_Set_Stats.write('\n')
        Synth_Set_Stats.write(TF_input)
        Synth_Set_Stats.write('\n')
    
    
        file_path = os.path.join('..', 'data' , 'Aligned_motifs','%s.txt' %TF_input)
        f = open(file_path, "r")
        f.seek(0)
        lines = f.readlines()
        n = len(lines)
        bs_len = len(lines[0])
        
        strand_length = int(raw_input("Enter strand length: "))
        num_sites = int(raw_input("Enter number of sites to place in strand: "))
        num_sets = int(raw_input("Enter desired number of sets: "))
        
        loci_list = []
        for i in range(num_sites):
            #you will fix this and make it better later
            loci_list.append(int(raw_input("Enter loci for binding sites. Keep in mind the length of each site ")))
       
        
        Synth_Set_Stats.write('Number of binding sites:')
        Synth_Set_Stats.write('\n')
        Synth_Set_Stats.write(str(num_sites))
        Synth_Set_Stats.write('\n')
        Synth_Set_Stats.write('Binding site length: ')
        Synth_Set_Stats.write('\n')
        Synth_Set_Stats.write(str(bs_len))
        Synth_Set_Stats.write('\n')
        
        Synth_Set_Stats.write('Loci:')
        Synth_Set_Stats.write('\n')
    
        count = 0
        
        for item in loci_list:
            
            #this accounts for the length of the bs changing the locus
            item += (int(bs_len)*count)
            Synth_Set_Stats.write(str(item))
            Synth_Set_Stats.write(', ')
            count += 1
        
        Synth_Set_Stats.write('\n')
        
        Synth_Set_Stats.write('Strand Length')
        Synth_Set_Stats.write('\n')
        Synth_Set_Stats.write(str(strand_length))
        Synth_Set_Stats.write('\n')
        
        
        set_count = 0
        
        
        while set_count < num_sets:


#print out which set you're using
            print "Using: %s" %TF_input

#determine number of sites to enrich strand with

            synthetic_set = Logical(strand_length, num_sites, lines, loci_list, loci, Synth_Set_Stats)
            print "length ", len(synthetic_set)
            #print synthetic_set
            #print 'LENGTH OF SYN SET', len(synthetic_set)
            outFile.write("< binding site loci: ")
                
            
           
            
            for locus in loci:
                try:
                    outFile.write(locus)
                    #if loci.index(locus) != 0 or loci.index(locus)%2 != 0:
                        #BS_File.write(locus) 
                        #BS_File.write("\n")
                except TypeError:
                    #print 'TYPE ERROR OCCURED. REVIEW NEEDED'
                    outFile.write(str(locus))
                outFile.write(" ")
            outFile.write("\n") 
            outFile.write(synthetic_set)
            outFile.write("\n")
                
            set_count += 1

#array of nucleotides

#randomly select one each time

#make biniding site a certain bp
            f.close()
#return the dataset
        outFile.close()
        #BS_File.close()
    
    print "Total time: " , (time.time() - total_start)
    print "Set count check: " , set_count

Make_Synthetic_Set()
