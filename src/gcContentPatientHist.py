#File: gcContentHist.py
#Author: Grace Chandler
#Last edit: 05-28-14
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


    #make variable for the GC to use for search
    search = "GC"


   #keep track of number of files 
    file_count = 0

    #create the path for the folder 
    data_folder = os.path.join('..', 'data' , 'LexA_Binding_and_Upstream_Regions', 'Pruned')
    
    
    #Make an outfile to store the data
    outFile = open('data.txt', 'w')
    statFile = open('stat.txt', 'w')

    for file in os.listdir(data_folder):
        #Open the file
        f  = open( os.path.join( data_folder , file), 'r+')
        f.seek(0)

        
        gcContents = []
        region_count = 0


        #Parse the file
        
                
        #make temp string
        lines = f.readlines()
         
        n = len(lines)

        #Parse the file
        for line in lines:
            #if the bp is g or c, increment the count                                                                                                                                
            
            #Lines that start with this character do not contain DNA bases
            if line != "\n" and line != "\r\n" and line[0] != ">":
                
                #increment region count
                region_count += 1

                #set gc count to 0
                gc_count = 0 
                
                #set base count to 0

                base_count = 0

                #iterate through the bases
                for base in line:

                    #if the base is g or c then add to gc count
                    if base in search or  base in search.lower():
                        gc_count += 1

                    #iterate through all the bases
                    base_count += 1

        #conentration will equal the number of g and c divided by total search width 300                                                                                             
                gcFrequency = float(float(gc_count) / float(base_count))
                
                if gcFrequency < 0.1:
                    print line
        #add this concentration to the array                                                                                                                                         
                gcContents.append(gcFrequency)
        
                #Make the data to put in outfile
            #    for content in gcContents:
                    
        patientData = [len(filter(lambda x: (i/100.0) <= x < (i/100.0)+.01, gcContents)) for i in range(0,101)]

        
        patientDataLine = ",".join(map(str, patientData))
        outFile.write(patientDataLine)
        outFile.write("\n")
                     
        

        file_count += 1

        print "File: %d" %file_count
        statFile.write("File%d \n" %file_count)

    #print the number of regions for tracking reasons
        print "Number of regions: ", region_count
        statFile.write("Number of regions: %d \n" %region_count)
    #print the standard deviation
        print np.std(gcContents), "-Standard Deviation"
        statFile.write(" %f -Standard Deviation \n" %np.std(gcContents))

    #print the mean
        print np.mean(gcContents), "-Mean"
        statFile.write(" %f -Mean \n" %np.mean(gcContents))

    #Make the histogram
        plt.hist(gcContents, bins = 50 )
    
    #y axis
        plt.ylabel('Frequency')
    
    #x axis
        plt.xlabel('GC Content')
    
    #display the plot
        plt.savefig('Patient%d' %file_count)

        plt.close()

        region_count = 0
        f.close()

    outFile.close()
    statFile.close()
main()
