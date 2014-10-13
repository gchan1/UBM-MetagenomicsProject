#File: PSSM_Calls.py
#Author: Grace Chandler
#Last edit: 07-31-2014
#Description: This file contains the function to make PSSM calls

import os
import sys

#This function reads in a file of PSSM scores produced by the file PSSM_Markov_1.py and makes them usable. 
def Get_PSSM_Sets():
    
    sets = []
    
    InFile = open('PSSM_Scores.txt', 'r')
    lines = InFile.readlines()
    for line in lines:
        line = line.strip()
        #print line
        line = line.split()
        #print line
        for item in line:
            item = float(item)
        sets.append(line)
        
    return sets

#This function makes a master list of all the scores
#This is useful for determining an overall mean and standard deviation for all scores
def Make_List_All(pssm_sets):
    
    total_set = []
    
    for set in pssm_sets:
        for score in set:
            total_set.append(float(score))
    
    return total_set
 
 
#this function checks if there is a binding site at an index
#this is based off of info from the ROC_Data.txt file created by the PSSM_Markov_1 file           
def Check_If_Binding_Site(index, bs_len, loci_list):
    
    #determine if the PSSM score's index lies within the range of the actual binding site
     
    #this section is for if you want to hard code loci   
    #Make the binding site ranges
    #binding_site_regions = [[10, (10+bs_len)], [(40+bs_len), (40+(2*(bs_len)))]]
    
    
    #Check the start and end regions
    for region in loci_list:
        start = int(region)
        end = start+bs_len
        if index in range(start,end):
            return True       
        else:
            return False
        
    #return false
    
#Making site calls
#Input: array of pssm scores (total), array of pssm scores(strand), threshold
#Output: technically void, this will mutate the true pos/false negs for the ROC curve
#This calls the check true positive function
#The threshold is a variable so this can be run many times for the ROC curve
def Make_Call(pos_neg, PSSM_Scores, threshold, bs_len, loci_list):
    
    #check to see if the score is above the threshold
    for score in PSSM_Scores:
       
        #where is the score in the strand
        index = PSSM_Scores.index(score)
        
        
        #fixing the negatives because they weren't being recognized
        if '-' in str(score):
            score = score.strip('-')
            score = -float(score)
       
        score = float(score)
        
        if score < threshold:
            #print 'less than 0'
            if Check_If_Binding_Site(index, bs_len, loci_list):
                pos_neg['false_neg'] += 1
            else:
                pos_neg['true_neg'] += 1
        
        #if its above the threshold, its a score
        elif score >= threshold:
            if Check_If_Binding_Site(index, bs_len, loci_list):
                pos_neg['true_pos'] += 1
            else:
                pos_neg['false_pos'] += 1
       
                    
    print pos_neg
    return pos_neg  
 
#this function plots an ROC curve for the PSSM data, showing how effective our method was at different thresholds
def Plot_ROC(true_pos_rates, false_pos_rates, thresholds):
    
    import numpy as np
    import matplotlib.pyplot as plt

    t = np.arange(0., 1., 0.01)

      
    extrax = [0, false_pos_rates[-1]]
    extray = [0, true_pos_rates[-1]]
    
    #plt.plot(extrax, extray) 
    y = true_pos_rates
    x = false_pos_rates
    
    #figure out a way to do this computationally
    z = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L']

    fig, ax = plt.subplots()
    plt.plot(false_pos_rates, true_pos_rates, '-b', false_pos_rates, true_pos_rates, 'ro')
    ax.plot(x, y)
    plt.plot(t, t, 'r--')
    
    for X, Y, Z in zip(x, y, z):
        # Annotate the points 5 _points_ above and to the left of the vertex
        ax.annotate('{}'.format(Z), xy=(X,Y), xytext=(-20, 20), ha='right',  textcoords='offset points')
       
    plt.ylabel('True Positive Rate') 
    plt.xlabel('False Positive Rate')
    #plt.title('ROC Curve for Site Identification at Different Thresholds')
    plt.savefig('ROC_Curve')
    plt.show()
    
    
#This function returns the true positive and false positive rates 
#These are used to plot the ROC curve
def Get_Rates(pos_neg):
    
    TPR = float(pos_neg['true_pos'])/ float((pos_neg['true_pos'] + pos_neg['false_neg']))
    FPR = float(pos_neg['false_pos'])/ float((pos_neg['false_pos'] + pos_neg['true_neg']))
    
    return TPR, FPR


def main():
     
     
    pos_neg_sets = []
    #getting info from the synthetic set file to be used
    Synth_Set_Stats = open('Synth_Set_Stats.txt', 'r')
    Threshold_Data = open('Threshold_Data.txt', 'w')
    lines = Synth_Set_Stats.readlines()
    For_Excel = open("For_Excel.txt", 'w')
    
    bs_len = int(lines[5])
    loci_list = lines[7].strip()
    loci_list = lines[7].strip(',')
    loci_list = loci_list.split(',') 
     
    for item in loci_list:
        try:
            item = int(item)
        except ValueError:
            loci_list.remove(item)
            
    #Establishing the dictionary to hold data to make TPR and FPR            
    pos_neg = dict()
    pos_neg['true_pos'] = 0
    pos_neg['false_pos'] = 0
    pos_neg['true_neg'] = 0
    pos_neg['false_neg'] = 0    
    
    keys = pos_neg.keys()
     
    pssm_sets = Get_PSSM_Sets()
    total_set = Make_List_All(pssm_sets)

   
    import numpy as np  
    total_set = np.array(total_set)
    
      
    mean = np.mean(total_set)
    standard_dev = np.std(total_set)
    
    thresholds = [ mean, (mean-standard_dev), (mean - (2*standard_dev)), (mean - (3* standard_dev)),  (mean+standard_dev), (mean + (2*standard_dev)), (mean + (3* standard_dev)), (mean + (4* standard_dev)), (mean + (5* standard_dev)), (mean + (6* standard_dev)),(mean + (7* standard_dev)),(mean + (8* standard_dev))]
   
    thresholds = sorted(thresholds)
    
    print thresholds
    true_pos_rates = []
    false_pos_rates = []
    
    count = 0 
    
    #This for excel file has all the data necessary to make the ROC curve so you can do it in other programs
    #maybe make it CSV?
    For_Excel.write("Thresholds \n")
    for threshold in thresholds:
        
        For_Excel.write(str(threshold))
        For_Excel.write('\n')
        count += 1
        print 'threshold #: ', count
        
        for set in pssm_sets:
            pos_neg.keys()
            pos_neg = Make_Call(pos_neg, set, threshold, bs_len, loci_list)
            pos_neg_sets.append(pos_neg)
        
        TPR, FPR = Get_Rates(pos_neg)
        true_pos_rates.append(TPR)
        false_pos_rates.append(FPR)
            
        Threshold_Data.write('Threshold: ')
        Threshold_Data.write(str(threshold))
        Threshold_Data.write('\n')
        
        for key in keys:
             Threshold_Data.write(key)
             Threshold_Data.write(': ')
             Threshold_Data.write(str(pos_neg[key]))
             Threshold_Data.write('\n')
    
    For_Excel.write('TPR \n')
    for item in true_pos_rates:
        For_Excel.write(str(item))
        For_Excel.write('\n')
        
    For_Excel.write('FPR \n')
    for item in false_pos_rates:
        For_Excel.write(str(item))
        For_Excel.write('\n')
        
    Plot_ROC(true_pos_rates, false_pos_rates, thresholds)    
    print 'Finished'
     
    Synth_Set_Stats.close()
    Threshold_Data.close()
    For_Excel.close()
main()           