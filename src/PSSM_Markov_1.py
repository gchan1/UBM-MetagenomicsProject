__author__ = 'jobott'
#Name:
#Inputs:
#Outputs:
#Description:

#Importing the math function (needed to perfom the log)
import math
import os.path
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

    
#Function to determine the PSFM, takes an input of the putative binding sites in the form of a list.
#Note: all of the binding sites are assumed to be same length
#Inputs: List of binding sites
#Outputs: PSFM
def PSFM(binding_sites):

   
    a_bases = []
    c_bases = []
    t_bases = []
    g_bases = []

    #Definging some other important variables
    total_sites = float(len(binding_sites))
    site_length = int(len(binding_sites[0]))
    window = site_length

    #Making the PSFM have the number of columns as the length of the first binding site
    a_bases.extend([float(0)]*site_length)
    c_bases.extend([float(0)]*site_length)
    t_bases.extend([float(0)]*site_length)
    g_bases.extend([float(0)]*site_length)

    #Viewing each site individually
    for site in binding_sites:

        #Creating the count aspect of the PSFM
        for base in range(len(site)):
            pos = base + 1
            if site[base] == 'A':
                a_bases[base] += 1
            if site[base] == 'C':
                c_bases[base] += 1
            if site[base] == 'T':
                t_bases[base] += 1
            if site[base] == 'G':
                g_bases[base] += 1

    #Making the counts into frequencies
    a_freq = [float(x/total_sites) for x in a_bases]
    c_freq = [float(x/total_sites) for x in c_bases]
    t_freq = [float(x/total_sites) for x in t_bases]
    g_freq = [float(x/total_sites) for x in g_bases]


    #Printing the PSFM
    print 'Obtained PSFM:'
    print 'A:', a_freq
    print 'C:', c_freq
    print 'T:', t_freq
    print 'G:', g_freq


    #Create a list to store the PSFM in (a list with the 4 lists of the frequencies of a,c,t,g)
    PSFM = [a_freq,c_freq,t_freq,g_freq]

    #Returning the PSFM
    return PSFM


#Applying the PSFM to determine motif-based probabilities when applying a sliding window
#Inputs: PSFM; sequences which we want to apply the PSFM to (in one long string)
#Outputs: A list of P(subsequence|motif) for each subsequence fromed by a sliding window
def P_motif(PSFM,sequences):

    #Creating a list to store the probability of each subsequece belonging to the TF motif
    list_of_motif_probs = []

    #Determining the window length based off of the length of the PSFM motif
    window = int(len(PSFM[0]))

    #Iterating the indexes of all of the initial bases for each case that the sliding window encounters
    for initial in range(len(sequences)-window+1):
        #Creating a subsequence that starts at the given initial index and is as long as the pre-determined window
        subsequence = sequences[initial:initial+window]
        #print 'subsequence', subsequence
        #Setting the initial value for the motif probability of the sliding window at 1
        motif_prob = float(1)
        #Iterating over the indexes of the bases in the subsequence
        for i in range(len(subsequence)):
          
            
            #Determines what the value of the base is at that position in the subsequence (a,c,t,g)
            if subsequence[i] == 'A':
                #Multiply the value of the motif_prob by the value of the frequency of that base at the
                #position of that base within the PSFM
                motif_prob *= float(PSFM[0][i])
            elif subsequence[i] == 'C':
                motif_prob *= float(PSFM[1][i])
            elif subsequence[i] == 'T':
                motif_prob *= float(PSFM[2][i])
            elif subsequence[i] == 'G':
                motif_prob *= float(PSFM[3][i])
        print'motif prob', motif_prob
        #Appending the list that stores the motif probabilities by the value motif_prob
        list_of_motif_probs.append(motif_prob)
        
        
        #Re-setting the value of motif_prob to 1 to apply to the next subsequence
        motif_prob = 1

    #Iterating over the indexes of the values in list_of_motif_prob
    for start in range(len(list_of_motif_probs)):
        #Printing out the values of motif_prob with appropriate labeling
        #General format:
        #'P_motif(initial base - final base in window) = appropriate value of P_motif at the
        #                                                given point of the sliding window
        print 'P_motif(bases %s'%(start+1), '- %s)='%(start+window),  '%s'%list_of_motif_probs[start]

    #Having the function return the list containing the motif probabilities
    return list_of_motif_probs






#Function used to find the background probabilities when using a Markov-1 model
#Inputs: sequence that is going to be used for both the background generation and to apply the sliding window on
#Outputs: A list of P(subsequence|background) for each subsequence fromed by a sliding window
def background_probs(sequence, binding_sites):

    #A way where the user is asked what length of a sliding window that they want to use
    #Also makes the window length an int
    
    input_check = False
    window = int(len(binding_sites[0]))
            
    #Might want to make this a function like def initialize_counts()
    #Making a dictionary to hold the counts of all possible base transitions
   
    counts  = {}
    counts["AA_count"] = 0
    counts["AC_count"] = 0
    counts["AT_count"] = 0
    counts["AG_count"] = 0
    counts["CA_count"] = 0
    counts["CC_count"] = 0
    counts["CT_count"] = 0
    counts["CG_count"] = 0
    counts["TA_count"] = 0
    counts["TC_count"] = 0
    counts["TT_count"] = 0
    counts["TG_count"] = 0
    counts["GA_count"] = 0
    counts["GC_count"] = 0
    counts["GT_count"] = 0
    counts["GG_count"] = 0
    
    keys = counts.keys()
    values = counts.values()

    #Defining a variable for sequence length
    sequence_length = len(sequence)


    print "sequence", sequence
    
    #Calculating the counts of the transitions in the sequence inputted
    #Iterates over the all the bases in a sequence except the last one
    #Not the last one becasue it does not have a base after it,
    #hence no transition between bases that start at that base
    for i in range(sequence_length-1):
        #The first set of if/elif statements asks whether the
        # intitial base is A, T, C or G
        if sequence[i] == 'A':
            if sequence[i+1] == 'A':
                counts["AA_count"] += 1
            elif sequence[i+1] == 'C':
                counts["AC_count"] += 1
            elif sequence[i+1] == 'T':
                counts["AT_count"] += 1
            elif sequence[i+1] == 'G':
                counts["AG_count"] += 1
        elif sequence[i] == 'C':
            if sequence[i+1] == 'A':
                counts["CA_count"] += 1
            elif sequence[i+1] == 'C':
                counts["CC_count"] += 1
            elif sequence[i+1] == 'T':
                counts["CT_count"] += 1
            elif sequence[i+1] == 'G':
                counts["CG_count"] += 1
        elif sequence[i] == 'T':
            if sequence[i+1] == 'A':
                counts["TA_count"] += 1
            elif sequence[i+1] == 'C':
                counts["TC_count"] += 1
            elif sequence[i+1] == 'T':
                counts["TT_count"] += 1
            elif sequence[i+1] == 'G':
                counts["TG_count"] += 1
        elif sequence[i] == 'G':
            if sequence[i+1] == 'A':
                counts["GA_count"] += 1
            elif sequence[i+1] == 'C':
                counts["GC_count"] += 1
            elif sequence[i+1] == 'T':
                counts["GT_count"] += 1
            elif sequence[i+1] == 'G':
                counts["GG_count"] += 1

    #Defining a variable for the total number of transitions observed
    transition_total = 0
    #Iterates over all the counts in the counts dictionary and adds them to the value
    #of the total number of transitions
    
    for key in counts:
        transition_total += counts[key]
    
        
    #This is a check to make sure that the total number of transitions observed
    #was correct
    #If it is, then the total number of transitions observed in our Markov-1 assembly
    #should be equal to the length of the sequence minus one
    if transition_total != sequence_length - 1:
        print 'Incorrect number of transitions observed'

    #Printing out information regarding the counts of each transition
    for key in keys:
        print key,
        print '='
        print counts[key]
        
        
    #Making the counts a float
    for key in keys:
        counts[key] = float(counts[key])

    #Printing out the sequence length
    print ' '
    print 'The length of the sequence was %s' %sequence_length

    for key in keys:
        print key, counts[key]
        
    #Defining variables that are the values of counts with the same first base (A,C,T,G)
    AX_count = float(counts['AA_count'] + counts['AC_count'] + counts['AT_count'] + counts['AG_count'])
    CX_count = float(counts['CA_count'] + counts['CC_count'] + counts['CT_count'] + counts['CG_count'])
    TX_count = float(counts['TA_count'] + counts['TC_count'] + counts['TT_count'] + counts['TG_count'])
    GX_count = float(counts['GA_count'] + counts['GC_count'] + counts['GT_count'] + counts['GG_count'])

    #Makingg a dictionary to hold all the values of the background probabilities
    #Note: prob_AC = P(C|A), the conditional probability of getting C as the second base
    #given that A is the first base; prob_TG = P(G|T); ect.
    background_probs = {
        'prob_AA' : counts['AA_count']/AX_count,
        'prob_AC' : counts['AC_count']/AX_count,
        'prob_AT' : counts['AT_count']/AX_count,
        'prob_AG' : counts['AG_count']/AX_count,
        'prob_CA' : counts['CA_count']/CX_count,
        'prob_CC' : counts['CC_count']/CX_count,
        'prob_CT' : counts['CT_count']/CX_count,
        'prob_CG' : counts['CG_count']/CX_count,
        'prob_TA' : counts['TA_count']/TX_count,
        'prob_TC' : counts['TC_count']/TX_count,
        'prob_TT' : counts['TT_count']/TX_count,
        'prob_TG' : counts['TG_count']/TX_count,
        'prob_GA' : counts['GA_count']/GX_count,
        'prob_GC' : counts['GC_count']/GX_count,
        'prob_GT' : counts['GT_count']/GX_count,
        'prob_GG' : counts['GG_count']/GX_count,
    }

    keys = background_probs.keys()
    values = background_probs.values()
    
    #Printing out the conditional probabilities
    print ' '
    print 'Conditional probabilities:'
    for key in keys:
        print key,
        print '=',
        print background_probs[key]

    #Overall section scanning an input sequence with a sliding window and
    #applying the Markov-1 background probabilities

    #Defining a list to store the background probabilities for each window observed
    list_of_back_probs = []

    #Iterating over all the bases that would be the start of a window
    #Hence the use of sequence_length - window becasue any of the bases within
    #a window length of the end would not have a full window length of window
    #of bases following it
    for i in range(sequence_length-window):
        #Defines a subsequence wich is the length of the window and starts at the first base
        subseq = sequence[i:i+window+1]
        #Defining a variable to hold the value of the Markov-1 background
        #probability for the given window
        m1_prob = 1
        #Only goes through this if it is the first term of the sequence:
        #Reverses the order of the sequence and then applies the markov-1 background
        #probability model to it. (This is explained in the comments on the next section)
        if i == 0:
            #Inverting the sequence order
            subseq = subseq[::-1]
            for i in range(len(subseq)-1):
                if subseq[i] == 'A':
                    if subseq[i+1] == 'A':
                        m1_prob *= background_probs['prob_AA']
                    elif subseq[i+1] == 'C':
                        m1_prob *= background_probs['prob_AC']
                    elif subseq[i+1] == 'T':
                        m1_prob *= background_probs['prob_AT']
                    elif subseq[i+1] == 'G':
                        m1_prob *= background_probs['prob_AG']
                elif subseq[i] == 'C':
                    if subseq[i+1] == 'A':
                        m1_prob *= background_probs['prob_CA']
                    elif subseq[i+1] == 'C':
                        m1_prob *= background_probs['prob_CC']
                    elif subseq[i+1] == 'T':
                        m1_prob *= background_probs['prob_CT']
                    elif subseq[i+1] == 'G':
                        m1_prob *= background_probs['prob_CG']
                elif subseq[i] == 'T':
                    if subseq[i+1] == 'A':
                        m1_prob *= background_probs['prob_TA']
                    elif subseq[i+1] == 'C':
                        m1_prob *= background_probs['prob_TC']
                    elif subseq[i+1] == 'T':
                        m1_prob *= background_probs['prob_TT']
                    elif subseq[i+1] == 'G':
                        m1_prob *= background_probs['prob_TG']
                elif subseq[i] == 'G':
                    if subseq[i+1] == 'A':
                        m1_prob *= background_probs['prob_GA']
                    elif subseq[i+1] == 'C':
                        m1_prob *= background_probs['prob_GC']
                    elif subseq[i+1] == 'T':
                        m1_prob *= background_probs['prob_GT']
                    elif subseq[i+1] == 'G':
                        m1_prob *= background_probs['prob_GG']
            
            
            list_of_back_probs.append(m1_prob)
            print subseq
            #Re-inverting the sequence order so that it is the same as before
            subseq = subseq[::-1]
            #Re-setting m1_prob to 1 for the markov-1 model to be applied again
            m1_prob = 1

        #Iterating over the bases in the formed subsequence
        for i in range(len(subseq)-1):
            #Set of if/elif statements that sorts based on the indexed base
            if subseq[i] == 'A':
                #Set of if/elif statements that sorts further based on the base
                #that follows the indexed base
                if subseq[i+1] == 'A':
                    #Once sorted, multiplying the current value of the m1_prob for
                    #the subsequence by the given conditional probability of the second
                    #base occuring given that the first base occured
                    m1_prob *= background_probs['prob_AA']
                elif subseq[i+1] == 'C':
                    m1_prob *= background_probs['prob_AC']
                elif subseq[i+1] == 'T':
                    m1_prob *= background_probs['prob_AT']
                elif subseq[i+1] == 'G':
                    m1_prob *= background_probs['prob_AG']
            elif subseq[i] == 'C':
                if subseq[i+1] == 'A':
                    m1_prob *= background_probs['prob_CA']
                elif subseq[i+1] == 'C':
                    m1_prob *= background_probs['prob_CC']
                elif subseq[i+1] == 'T':
                    m1_prob *= background_probs['prob_CT']
                elif subseq[i+1] == 'G':
                    m1_prob *= background_probs['prob_CG']
            elif subseq[i] == 'T':
                if subseq[i+1] == 'A':
                    m1_prob *= background_probs['prob_TA']
                elif subseq[i+1] == 'C':
                    m1_prob *= background_probs['prob_TC']
                elif subseq[i+1] == 'T':
                    m1_prob *= background_probs['prob_TT']
                elif subseq[i+1] == 'G':
                    m1_prob *= background_probs['prob_TG']
            elif subseq[i] == 'G':
                if subseq[i+1] == 'A':
                    m1_prob *= background_probs['prob_GA']
                elif subseq[i+1] == 'C':
                    m1_prob *= background_probs['prob_GC']
                elif subseq[i+1] == 'T':
                    m1_prob *= background_probs['prob_GT']
                elif subseq[i+1] == 'G':
                    m1_prob *= background_probs['prob_GG']

        #After iterating through the subsequence and obtaining its m1_prob,
        #apending that to our list of background probabilities
        list_of_back_probs.append(m1_prob)
        #Re-setting m1_prob to 1 so that the next subsequence's m1_prob can be calculated
        m1_prob = 1

    #Printing the background probabilities for each subsequence
    print ' '
    for start in range(len(list_of_back_probs)):
        print 'P_back(bases %s'%(start+1), '- %s)='%(start+window),  '%s'%list_of_back_probs[start]

    return list_of_back_probs

#to make histograms
def Plot_Histogram(data_set, file_name, num_bins, x_label, y_label):

    #Make the histogram                                                                                                                                                           
    plt.hist(data_set, bins = num_bins )

    #y axis                                                                                                                                                                        
    plt.ylabel(y_label)

    #x axis                                                                                                                                                  
    plt.xlabel(x_label)
    
    #display the plot         
    
    #plt.show()
                                                                                                                                                          
    plt.savefig('Sequence: %d' %file_name)
    
    #Close the plot to avoid overlap                                                                                                                                               
    plt.close()


#This function assigns locations to be linked to values in a list, for when a dictionary can't be used
#A situation where this would be needed would be for sorting.
def assign_location(array):
    new_array = []
    for i in range(len(array)):
    
        temp_array = [(i+1), array[i]]
        new_array.append(temp_array)
                
    return new_array
        
#This function is a slightly altered quicksort made to sort PSSM scores while maintaining their locus 
def pssm_quick_sort(array):
    
    
    print array
    
    #To store values during sort
    less_than = []
    greater_than = []
    equal_to = []
    
    
    
    if len(array) > 1:
        
        pivot_point = array[0][1]
        
        #this is the point the other  values will be compared to
    
        for item in array:
            
            
            #Only the second number is the sorting value
            x = item[1]
            if x < pivot_point:
                print 'add to less than'
                less_than.append(item)
            if x > pivot_point:
                print 'add to greater than'
                greater_than.append(item)
            if x == pivot_point:
                print 'add to equal to'
                equal_to.append(item)
                
        return pssm_quick_sort(less_than) + equal_to + pssm_quick_sort(greater_than)
    
    else: 
        return array

#Forming the PSSM
#Inputs: sequence, list of binding sites
#Outputs: PSSM scores for each subsequence in the sliding window
def PSSM(binding_sites,sequence, number):

    window = len(binding_sites[0])
    
    #Forming the motif probabilities
    list_of_motif_probs = P_motif(PSFM(binding_sites),sequence)

    #Forming the background probabilities
    list_of_back_probs = background_probs(sequence, binding_sites)

   
    
    #Creating a list to store the ratio of P(subsequence|motif)/P(subsequence|background)
    prob_ratios = []

    #Iterating over the index of each subsequence initial base
    for initial in range(len(list_of_motif_probs)):
         
        print "length of background probs", len(list_of_back_probs)
        print "length of motif probs", len(list_of_motif_probs)
         
        #Determining each initial index's P(subsequence|motif)/P(subsequence|background)
        #and appending it to prob_ratios
        prob_ratios.append(list_of_motif_probs[initial]/list_of_back_probs[initial])

    #Creating a list to store the PSSM scores in
    PSSM_scores = []
    
    
    #Iterating over the ratios in prob_ratios
    for ratio in prob_ratios:
        #Appending the PSSM list with the PSSM score, determined by taking the log base 2 for the ratio
        PSSM_scores.append(math.log(ratio,2))

    #Iterating over the indexes of the values in list_of_motif_prob
    for start in range(len(PSSM_scores)):
        #Printing out the PSSM values with appropriate labeling
        #General format:
        #'PSSM_score(initial base - final base in window) = appropriate PSSM value at the
        #                                                   given point of the sliding window
       
       
        print 'PSSM_score(bases %s'%(start+1), '- %s)='%(start+window),  '%s'%PSSM_scores[start]


    #file_name = raw_input("Please enter a name for the Histogram file: ")
    Plot_Histogram(PSSM_scores, number, 50 , 'Score', 'Counts')

    PSSM_scores = assign_location(PSSM_scores)
    PSSM_scores = pssm_quick_sort(PSSM_scores)
    return  PSSM_scores

    

def main():

    Sequence_Data = open('Synth_Sets.txt', 'r')
    Binding_Sites = open('Binding_Sites.txt', 'r')
    PSSM_Scores = open('PSSM_Scores.txt', 'w')
    Sequence_Data.seek(0)
    Binding_Sites.seek(0)
    Sequence_Data_lines = Sequence_Data.readlines()
    Binding_Sites_lines = Binding_Sites.readlines()

    
    sequences = []
    binding_sites = []
    
    for line in Sequence_Data_lines:
        if line != "\n" and line != "\r\n" and line[0] != "<":
            sequences.append(line.strip())
    print sequences
    for line in Binding_Sites_lines:
        if line != "\n" and line != "\r\n" and line[0] != "<":
            binding_sites.append(line.strip())
    
    #do the pssm
    for i in range(len(sequences)):
        print'sequence i ', sequences[i]
        PSSM_Score_List = PSSM(binding_sites, sequences[i] , i)
        PSSM_Scores.write("Sequence #%d"  %i)
        
        for item in PSSM_Score_List:
            PSSM_Scores.write(str(item))
       
        PSSM_Scores.write("\n")
    print 'Finished'
    #bsites = ['ACTGACTG', 'CTGACTGA' ,'TGACTGAC', 'GACTGACT', 'ACCTGAAT', 'ACCTGAAT', 'ACCCGATT','AACTGTAT']
    #x = 'AAATAAATCGAGCTACATAGAATATCTGTTCACCCTCGGGGAGCGTGGGGTGTAC' 
    
    #PSSM(bsites, x, 2)
main()
