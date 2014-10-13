#File: UBM_Function_Library.py
#Author: Grace Chandler
#Last edit: 28-05-14
#Project: Metagenomics
#Description:
#This file contains functions to be reused amongst my other projects




#This function assigns locations to be linked to values in a list, for when a dictionary can't be used
#A situation where this would be needed would be for sorting.
def assign_location(array):
    new_array = []
    for i in range(len(array)):
        for item in array:
                temp_array = [(i+1), item]
                new_array.append(temp_array)
                
    return new_array
        
#This function is a slightly altered quicksort made to sort PSSM scores while maintaining their locus 
def pssm_quick_sort(array):
    
    #To store values during sort
    less_than = []
    greater_than = []
    equal_to = []
    
    if len(array) > 1:
        
        #this is the point the other  values will be compared to
        pivot_point = array[0]
    
        for item in array:
            
            #Only the second number is the sorting value
            x = array[1]
            if x < pivot_point:
                less_than.append(item)
            if x > pivot_point:
                greater_than.append(item)
            if x == pivot_point:
                equal_to.append(item)
                
        return sort(less_than) + equal + sort(greater_than)
    
    else: 
        
        return array
    
#This function is used to plot the histograms                                                                                                                                   
def Plot_Histogram(data_set, patient_num, num_bins, x_label, y_label):

    #Make the histogram                                                                                                                                                           
    plt.hist(data_set, bins = num_bins )

    #y axis                                                                                                                                                                        
    plt.ylabel(y_label)

    #x axis                                                                                                                                                  
    plt.xlabel(x_label)
    
    #display the plot                                                                                                                                                               
    plt.savefig('Patient%d' %patient_num)
    
    #Close the plot to avoid overlap                                                                                                                                               
    plt.close()



#This function sorts the data into counts and returns a set of the counts                                                                                                            
#It then writes to the outifle, making a CSV file                                                                                                                                    
def Get_Counts(data_set, outFile):

        #Make the data to put in outfile                                                                                                                                             
        #This function sorts the gcContents into percentile groups, showing a count of contents for each percentage in range 1-100                                                   
        #Patient data is a list of all the gcContent counts per percentage                                                                                                           
        patientData = [len(filter(lambda x: (i/100.0) <= x < (i/100.0)+.01, data_set)) for i in range(0,101)]


        #To make a CSV file, join the different counts with a comma, first converting them to strings                                                                                
        #Write the line and then a space to the outfile                                                                                                                              
        patientDataLine = ",".join(map(str, patientData))
        outFile.write(patientDataLine)
        outFile.write("\n")


#This function gets the gc content for a line of DNA sequence
def Get_GC_Content(line, data_set, gc_count, base_count):

        #make variable for the GC to use for search                                                                                                                                  
        search = "GC"

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
        data_set.append(gcFrequency)


#This function is the Markov Monster
#heavy edits are needed, I should probably not put this here
#Good point grace, put it here later