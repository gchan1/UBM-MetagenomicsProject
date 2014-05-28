
#File: gcContentHist.py
#Author: Grace Chandler
#Last edit: 04-23-14
#Description:
#This program creates a histogram based on gc content in the noncoding regions of the data. The procedure includes finding the GC base concentration in the non coding regions of the metagnomic sample (provided) which are 300 bp in length. This will then be displayed using a histogram to compare GC content within the given sample. 
#Files imported:
#Displays: Histogram
import sys
import re
import numpy as np
import os.path
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

def main():

    #This list will hold gc contentrations for each noncoding region
    #to be used for histogram later
    gcContents = []

    #make variable for the GC to use for search
    search = "GC"

    #create the path for the folder 
    data_folder = os.path.join('..', 'data' , 'LexA_Binding_and_Upstream_Regions', 'Pruned')
    
    for file in os.listdir(data_folder):
        #Open the file
        f  = open( os.path.join( data_folder , file), 'r+')
        f.seek(0)



        #Parse the file
        
            #skip lines
            #lines where data begins in the files
        start = 1
        
        #basecount will count how many bases we iterate through
        base_count =0

        #make temp string
        lines = f.readlines()
        
        n = len(lines)

        #Create a list to hold the 300 bp regions
        #upstream_regions = [lines[i:i+3] for i in range(0,n,5)] 


   #iterate through each 300 bp region in the provided data
        #for upstream_region in upstream_regions:
        #make a count to use for the gc concentration
        #the count needs to be re initialized for each 300 bp section
        gc_count = 0
        
   #iterate through each bp in the nonCoding region of 300 bp
        for line in lines:
            #if the bp is g or c, increment the count
            if line[0] != "<":
                for base in line:
                    if base in search or  base in search.lower():
                        gc_count += 1
                    base_count += 1
        #conentration will equal the number of g and c divided by total search width 300
        gcFrequency = float(float(gc_count) / float(base_count))
        
        #add this concentration to the array
        gcContents.append(gcFrequency)
        f.close()
        
#    print gcContents

    print np.std(gcContents), "-Standard Deviation"
    print np.mean(gcContents), "-Mean"

    #Make the histogram
    plt.hist(gcContents, bins = 100)
    plt.ylabel('Frequency')
    plt.xlabel('GC Content')
    plt.show()


main()
