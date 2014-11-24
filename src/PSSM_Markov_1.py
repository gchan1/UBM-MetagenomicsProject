#Name: Grace Chandler and Jacob O'Bott
#Description: This program is a PSSM calling script using, for this time, synthetically created data sets. It makes pssm binding site calls
#and based off of these calls, it determines a true pos/false neg rate and an ROC curve. This is an incredibly rough copy of the final 
#product, but should work


#Things that need to be done

#training this off of the sequences without binding sites
#problems: motif probs aren't updating. WHY?
#Tasks: make both the markov1 and markov0 options
#it should do both because the entire point is for comparison

import time
import math
import os.path
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
    
#Function to determine the PSFM, takes an input of the putative binding sites in the form of a list.
#Note: all of the binding sites are assumed to be same length
#Inputs: List of binding sites
#Outputs: PSFM
def PSFM(binding_sites):

    print 'entered psfm'
    print binding_sites
    start = time.time()
    
    a_bases = []
    c_bases = []
    t_bases = []
    g_bases = []

    #Definging some other important variables
    total_sites = float(len(binding_sites))
    site_length = int(len(binding_sites[0]))
    window = site_length

    print 'total sites ', total_sites
    print 'site length ' , site_length


    #Making the PSFM have the number of columns as the length of the first binding site
    a_bases.extend([float(1)]*site_length)
    c_bases.extend([float(1)]*site_length)
    t_bases.extend([float(1)]*site_length)
    g_bases.extend([float(1)]*site_length)

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
                t_bases[base]  += 1
            if site[base] == 'G':
                g_bases[base] += 1

    #Making the counts into frequencies
    a_freq = [x/(total_sites+4) for x in a_bases]
    c_freq = [x/(total_sites+4) for x in c_bases]
    t_freq = [x/(total_sites+4) for x in t_bases]
    g_freq = [x/(total_sites+4) for x in g_bases]


    #Create a list to store the PSFM in (a list with the 4 lists of the frequencies of a,c,t,g)
    PSFM = [a_freq,c_freq,t_freq,g_freq]


    print 'time: ' ,( time.time()- start)
    #Returning the PSFM
    print 'psfm: ', PSFM
    return PSFM


#Applying the PSFM to determine motif-based probabilities when applying a sliding window
#Inputs: PSFM; sequences which we want to apply the PSFM to (in one long string)
#Outputs: A list of P(subsequence|motif) for each subsequence from by a sliding window
def P_motif(PSFM,sequences):


    print PSFM
    print sequences

    #Creating a list to store the probability of each subsequece belonging to the TF motif
    list_of_motif_probs = []

    #Determining the window length based off of the length of the PSFM motif
    window = int(len(PSFM[0]))
    print 'window: ', window

    #Iterating the indexes of all of the initial bases for each case that the sliding window encounters
    for initial in range(len(sequences)-window+1):

        #Creating a subsequence that starts at the given initial index and is as long as the pre-determined window
        subsequence = sequences[initial:initial+window]
        
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
        
        #Appending the list that stores the motif probabilities by the value motif_prob
       
        #Added the motif prob
        
        list_of_motif_probs.append(motif_prob)
        #print 'added the motif prob'
        
        #Re-setting the value of motif_prob to 1 to apply to the next subsequence
        motif_prob = float(1)

    #Iterating over the indexes of the values in list_of_motif_prob
   # for start in range(len(list_of_motif_probs)):
        #Printing out the values of motif_prob with appropriate labeling
        #General format:
        #'P_motif(initial base - final base in window) = appropriate value of P_motif at the
        #                                                given point of the sliding window
        #print 'P_motif(bases %s'%(start+1), '- %s)='%(start+window),  '%s'%list_of_motif_probs[start]

    #Having the function return the list containing the motif probabilities

    print list_of_motif_probs
    return list_of_motif_probs






#Function used to find the background probabilities when using a Markov-1 model
#Inputs: sequence that is going to be used for both the background generation and to apply the sliding window on
#Outputs: A list of P(subsequence|background) for each subsequence fromed by a sliding window
def background_probs(sequence, sequence_without_sites, binding_sites):

    #A way where the user is asked what length of a sliding window that they want to use
    #Also makes the window length an int
    
    input_check = False
    window = int(len(binding_sites[0]))
            
    #Might want to make this a function like def initialize_counts()
    #Making a dictionary to hold the counts of all possible base transitions
   
    counts  = {}
    counts["AA"] = 0
    counts["AC"] = 0
    counts["AT"] = 0
    counts["AG"] = 0
    counts["CA"] = 0
    counts["CC"] = 0
    counts["CT"] = 0
    counts["CG"] = 0
    counts["TA"] = 0
    counts["TC"] = 0
    counts["TT"] = 0
    counts["TG"] = 0
    counts["GA"] = 0
    counts["GC"] = 0
    counts["GT"] = 0
    counts["GG"] = 0
    
    keys = counts.keys()
    values = counts.values()

    #Defining a variable for sequence length
    sequence_length = len(sequence)

    print 'sequence length', sequence_length
    print 'without sites length', len(sequence_without_sites)
    
    #Calculating the counts of the transitions in the sequence inputted
    #Iterates over the all the bases in a sequence except the last one
    #Not the last one becasue it does not have a base after it,
    #hence no transition between bases that start at that base
    
    #CHANGE THIS SEQUENCE TO BE THE SEQUENCE WITHOUT THE BINDING SITES
    for i in range(sequence_length-1):
        #The first set of if/elif statements asks whether the
        # intitial base is A, T, C or G
        keyNucleotides = sequence_without_sites[i]+sequence_without_sites[i+1]
        print keyNucleotides
        counts[keyNucleotides] += 1;

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
    #if transition_total != sequence_length - 1:
        #print 'Incorrect number of transitions observed'

    #Printing out information regarding the counts of each transition
    for key in keys:
        print key,
        print '='
        print counts[key]
        
        
    #Making the counts a float
    for key in keys:
        counts[key] = float(counts[key])

    #Printing out the sequence length
    #print ' '
    #print 'The length of the sequence was %s' %sequence_length

    for key in keys:
        print key, counts[key]
        
    #Defining variables that are the values of counts with the same first base (A,C,T,G)
    A_count = float(counts['AA'] + counts['AC'] + counts['AT'] + counts['AG'])
    C_count = float(counts['CA'] + counts['CC'] + counts['CT'] + counts['CG'])
    T_count = float(counts['TA'] + counts['TC'] + counts['TT'] + counts['TG'])
    G_count = float(counts['GA'] + counts['GC'] + counts['GT'] + counts['GG'])

    
    #Makingg a dictionary to hold all the values of the background probabilities
    #Note: prob_AC = P(C|A), the conditional probability of getting C as the second base
    #given that A is the first base; prob_TG = P(G|T); ect.
   
    Nucleotides = ['A', 'C','T','G']
    
    background_probs = {}
    for nucleotide1 in Nucleotides:
        for nucleotide2 in Nucleotides:
            background_probs[nucleotide1+nucleotide2] = (counts[nucleotide1+nucleotide2]/ (locals()[(nucleotide1+ '_count')]))
    
    keys = background_probs.keys()
    values = background_probs.values()
    

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
                keyNucleotides = sequence[i]+sequence[i+1]
                m1_prob *= background_probs[keyNucleotides]
            
            list_of_back_probs.append(m1_prob)
            print 'added back prob'
            print subseq
            #Re-inverting the sequence order so that it is the same as before
            subseq = subseq[::-1]
            #Re-setting m1_prob to 1 for the markov-1 model to be applied again
            m1_prob = 1

        m1_prob = 1
        
        #Iterating over the bases in the formed subsequence
        for i in range(len(subseq)-1):
            #Set of if/elif statements that sorts based on the indexed base
            #sequence i represents the 1st base, i+1 the second base
            keyNucleotides = sequence[i]+sequence[i+1]
            m1_prob *= background_probs[keyNucleotides]
            

        #After iterating through the subsequence and obtaining its m1_prob,
        #apending that to our list of background probabilities
        list_of_back_probs.append(m1_prob)
        print 'added background prob'
        #Re-setting m1_prob to 1 so that the next subsequence's m1_prob can be calculated
        m1_prob = 1

    #Printing the background probabilities for each subsequence
    print ' '
    for start in range(len(list_of_back_probs)):
        print 'P_back(bases %s'%(start+1), '- %s)='%(start+window),  '%s'%list_of_back_probs[start]

    return list_of_back_probs

#Markov 0 model
#alternative ot the other background probs function
#for comparison purposes mainly
def markov_0_probs(sequence, sequence_without_sites, binding_sites):

    #A way where the user is asked what length of a sliding window that they want to use
    #Also makes the window length an int
    
    input_check = False
    window = int(len(binding_sites[0]))
            
    #Might want to make this a function like def initialize_counts()
    #Making a dictionary to hold the counts of all possible base transitions
   
    counts  = {}
    counts["A"] = 0
    counts["C"] = 0
    counts["T"] = 0
    counts["G"] = 0
    
    
    keys = counts.keys()
    values = counts.values()

    #Defining a variable for sequence length
    sequence_length = len(sequence)

    print 'sequence length', sequence_length
    print 'without sites length', len(sequence_without_sites)
    
    #Calculating the counts of the transitions in the sequence inputted
    #Iterates over the all the bases in a sequence except the last one
    #Not the last one becasue it does not have a base after it,
    #hence no transition between bases that start at that base
    
    #CHANGE THIS SEQUENCE TO BE THE SEQUENCE WITHOUT THE BINDING SITES
    for i in range(sequence_length-1):
        #The first set of if/elif statements asks whether the
        # intitial base is A, T, C or G
        keyNucleotides = sequence_without_sites[i]
        print keyNucleotides
        counts[keyNucleotides] += 1;

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
    #if transition_total != sequence_length - 1:
        #print 'Incorrect number of transitions observed'

    #Printing out information regarding the counts of each transition
    for key in keys:
        print key,
        print '='
        print counts[key]
        
        
    #Making the counts a float
    for key in keys:
        counts[key] = float(counts[key])

    #Printing out the sequence length
    #print ' '
    #print 'The length of the sequence was %s' %sequence_length

    for key in keys:
        print key, counts[key]
        
    #Defining variables that are the values of counts with the same first base (A,C,T,G)
    A_count = float(counts['A'])
    C_count = float(counts['C'])
    T_count = float(counts['T'])
    G_count = float(counts['G'])

    
    #Makingg a dictionary to hold all the values of the background probabilities
    #Note: prob_AC = P(C|A), the conditional probability of getting C as the second base
    #given that A is the first base; prob_TG = P(G|T); ect.
   
    Nucleotides = ['A', 'C','T','G']
    
    background_probs = {}
    for nucleotide1 in Nucleotides:
        background_probs[nucleotide1] = (counts[nucleotide1]/ (locals()[(nucleotide1+ '_count')]))
    
    keys = background_probs.keys()
    values = background_probs.values()
    

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
        m0_prob = 1
        
        #Only goes through this if it is the first term of the sequence:
        #Reverses the order of the sequence and then applies the markov-1 background
        #probability model to it. (This is explained in the comments on the next section)
        if i == 0:
            #Inverting the sequence order
            subseq = subseq[::-1]
            for i in range(len(subseq)-1):
                keyNucleotides = sequence[i]
                m0_prob *= background_probs[keyNucleotides]
            
            list_of_back_probs.append(m0_prob)
            print 'added back prob'
            print subseq
            #Re-inverting the sequence order so that it is the same as before
            subseq = subseq[::-1]
            #Re-setting m1_prob to 1 for the markov-1 model to be applied again
            m0_prob = 1

        m0_prob = 1
        
        #Iterating over the bases in the formed subsequence
        for i in range(len(subseq)-1):
            #Set of if/elif statements that sorts based on the indexed base
            #sequence i represents the 1st base, i+1 the second base
            keyNucleotides = sequence[i]
            m0_prob *= background_probs[keyNucleotides]
            

        #After iterating through the subsequence and obtaining its m1_prob,
        #apending that to our list of background probabilities
        list_of_back_probs.append(m0_prob)
        print 'added background prob'
        #Re-setting m1_prob to 1 so that the next subsequence's m1_prob can be calculated
        m0_prob = 1

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

        
#Forming the PSSM
#Inputs: sequence, list of binding sites
#Outputs: PSSM scores for each subsequence in the sliding window
def PSSM_Markov_1(binding_sites,sequence, number, dna_without_sites):

    print binding_sites
    window = len(binding_sites[0])
    
    #Forming the motif probabilities
    list_of_motif_probs = P_motif(PSFM(binding_sites),sequence)
    print list_of_motif_probs
    #Forming the background probabilities
    list_of_back_probs = background_probs(sequence, dna_without_sites,binding_sites)
    print list_of_back_probs
    
    #Creating a list to store the ratio of P(subsequence|motif)/P(subsequence|background)
    prob_ratios = []


    print "length of background probs", len(list_of_back_probs)
    print "length of motif probs", len(list_of_motif_probs)
    #Iterating over the index of each subsequence initial base
    for initial in range(len(list_of_motif_probs)):
        #Determining each initial index's P(subsequence|motif)/P(subsequence|background)
        #and appending it to prob_ratios
        print 'motif prob', float(list_of_motif_probs[initial])
        print 'back prob', float(list_of_back_probs[initial])
        tempRatio = float(list_of_motif_probs[initial])/float(list_of_back_probs[initial])
        prob_ratios.append(tempRatio)

    #Creating a list to store the PSSM scores in
    PSSM_scores = []
    
    #print prob_ratios
    #Iterating over the ratios in prob_ratios
    for ratio in prob_ratios:
        #Appending the PSSM list with the PSSM score, determined by taking the log base 2 for the ratio
        print "ratio", ratio
        PSSM_scores.append(math.log(ratio,2))

    #Iterating over the indexes of the values in list_of_motif_prob
    for start in range(len(PSSM_scores)):
        #Printing out the PSSM values with appropriate labeling
        #General format:
        #'PSSM_score(initial base - final base in window) = appropriate PSSM value at the
#                                                   given point of the sliding window
       
       
        print 'PSSM_score(bases %s'%(start+1), '- %s)='%(start+window),  '%s'%PSSM_scores[start]


    #file_name = raw_input("Please enter a name for the Histogram file: ")
    #Plot_Histogram(PSSM_scores, number, 50 , 'Score', 'Counts')

    
    #PSSM_scores = assign_location(PSSM_scores)
    #PSSM_scores = pssm_quick_sort(PSSM_scores)
    return  PSSM_scores

#This is the pssm markov 0 version to use for background model comparison
def PSSM_Markov_0(binding_sites,sequence, number, dna_without_sites):

    print binding_sites
    window = len(binding_sites[0])
    
    #Forming the motif probabilities
    list_of_motif_probs = P_motif(PSFM(binding_sites),sequence)
    print list_of_motif_probs
    #Forming the background probabilities
    list_of_back_probs = markov_0_probs(sequence, dna_without_sites,binding_sites)
    print list_of_back_probs
    
    #Creating a list to store the ratio of P(subsequence|motif)/P(subsequence|background)
    prob_ratios = []


    print "length of background probs", len(list_of_back_probs)
    print "length of motif probs", len(list_of_motif_probs)
    #Iterating over the index of each subsequence initial base
    for initial in range(len(list_of_motif_probs)):
        #Determining each initial index's P(subsequence|motif)/P(subsequence|background)
        #and appending it to prob_ratios
        print 'motif prob', float(list_of_motif_probs[initial])
        print 'back prob', float(list_of_back_probs[initial])
        tempRatio = float(list_of_motif_probs[initial])/float(list_of_back_probs[initial])
        prob_ratios.append(tempRatio)

    #Creating a list to store the PSSM scores in
    PSSM_scores = []
    
    #print prob_ratios
    #Iterating over the ratios in prob_ratios
    for ratio in prob_ratios:
        #Appending the PSSM list with the PSSM score, determined by taking the log base 2 for the ratio
        print "ratio", ratio
        PSSM_scores.append(math.log(ratio,2))

    #Iterating over the indexes of the values in list_of_motif_prob
    for start in range(len(PSSM_scores)):
        #Printing out the PSSM values with appropriate labeling
        #General format:
        #'PSSM_score(initial base - final base in window) = appropriate PSSM value at the
#                                                   given point of the sliding window
       
       
        print 'PSSM_score(bases %s'%(start+1), '- %s)='%(start+window),  '%s'%PSSM_scores[start]


    #file_name = raw_input("Please enter a name for the Histogram file: ")
    #Plot_Histogram(PSSM_scores, number, 50 , 'Score', 'Counts')

    
    #PSSM_scores = assign_location(PSSM_scores)
    #PSSM_scores = pssm_quick_sort(PSSM_scores)
    return  PSSM_scores

def main():

    #open all the files
    Sequence_Data = open('Synth_Sets.txt', 'r')
    file_path = os.path.join('..', 'data' , 'Aligned_motifs','Fur.txt' )
    Binding_Sites = open(file_path, "r")
    Binding_Sites.seek(0)
    PSSM_Scores_0 = open('PSSM_Scores_0.txt', 'w')
    PSSM_Scores_1 = open('PSSM_Scores_1.txt', 'w')
    Synth_Sans_Sites = open('noSites.txt', 'r')
    
    #make sure all the files are properly opened and ready to write
    Sequence_Data.seek(0)
    Binding_Sites.seek(0)
    Synth_Sans_Sites.seek(0)
    
    #get all the lines from these infiles
    Sequence_Data_lines = Sequence_Data.readlines()
    Binding_Sites_lines = Binding_Sites.readlines()
    Synth_Sans_Sites_lines = Synth_Sans_Sites.readlines()
    
    
    print "synthsetsanssites", Synth_Sans_Sites_lines
    
 
    #initialize arrays for sequences ,binding sites, and synth dna without sites
    sequences = []
    binding_sites = []
    
    
    
    #strip all the whitespace from all the lines
    for line in Sequence_Data_lines:
        if line != "\n" and line != "\r\n" and line[0] != "<":
            sequences.append(line.strip())
        
    for line in Binding_Sites_lines:
        if line != "\n" and line != "\r\n" and line[0] != "<":
            binding_sites.append(line.strip())
    
    no_sites = ''.join(Synth_Sans_Sites_lines).strip()
   

    #do the pssm
    for i in range(len(sequences)):
        print'sequence i ', sequences[i]
        
        #change this so we pass in the no_sites
        PSSM_Scores_Markov_0 = PSSM_Markov_0(binding_sites, sequences[i] , i, no_sites)
        PSSM_Scores_Markov_1 = PSSM_Markov_1(binding_sites, sequences[i] , i, no_sites)
        
        #write the scores to the appropriate file
        for item in PSSM_Scores_Markov_0:
            PSSM_Scores_0.write(str(item))
            PSSM_Scores_0.write(' ')
            
        for item in PSSM_Scores_Markov_1:
            PSSM_Scores_1.write(str(item))
            PSSM_Scores_1.write(' ')
        
        #go to the next line because we are finished with one sequence
        PSSM_Scores_0.write("\n")
        PSSM_Scores_1.write("\n")
            
        
    PSSM_Scores_0.close()
    PSSM_Scores_1.close()
    Binding_Sites.close()
    Sequence_Data.close()
       
        
    print 'Finished'
    


main()



