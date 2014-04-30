
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
    folder = os.path.join('/', 'Users', 'gracekathrynchandler', 'ErillLab', 'data', 'LexA_Binding_and_Upstream_Regions',)
    
    for file in os.listdir(folder):
        #Open the file
        f  = open( file, 'r+')
        
        #Create a list to hold the 300 bp regions
        bpPromoterRegions = []

        #Parse the file
        for line in f:
            
            #skip lines
            #lines where data begins in the files
            start = 1
            end = 4
            #make temp string
            tempString = f.readlines()[1:4]
            bpPromoterRegions.add(tempString)
            #Increment the start and end for the next section
            start += 5
            end += 5
   #iterate through each 300 bp region in the provided data
        for region in metagenomeData:
        #make a count to use for the gc concentration
        #the count needs to be re initialized for each 300 bp section
            count = 0
        #iterate through each bp in the nonCoding region of 300 bp
            for i in nonCoding:
            #if the bp is g or c, increment the count
                if i in search:
                    count += 1
        #conentration will equal the number of g and c divided by total search width 300
                    count = (count / 300.0)
        #add this concentration to the array
                    gcContents.add(count)


#make the gc content histogram
# make an x coordinate set
    #Put the GC content list in order
    gcContents = np.arrange(gcContents)
    x = dict()
    for content in gcContents:
        if content not in x:
            #value is a count of how many ( will be used as y values for historgram)
            x[content] = 1
        else:
            x[content] += 1

    plt.plot(x.keys(), x.values())
    


main()
