#File: UBM_Function_Library.py
#Author: Grace Chandler
#Last edit: 28-05-14
#Project: Metagenomics
#Description:
#This file contains functions to be reused amongst my other projects


#Function: Plot
#Input: set of data, x label y label
#Output: mean, standard deviation, histogram
def Plot( dataSet,  xlabel,  ylabel){


    #print the standard deviation                                                                                                                                                 
    print np.std(dataSet), "-Standard Deviation"

    #print the mean                                                                                                                                                                
    print np.mean(dataSet), "-Mean"

    #Make the histogram                                                                                                                                                             
    plt.hist(dataSet, bins = 100)

    #y axis                                                                                                                                                                         
    plt.ylabel(ylabel)

    #x axis                                                                                                                                                                         
    plt.xlabel(xlabel)

    #display the plot                                                                                                                                                               
    plt.show()

}
