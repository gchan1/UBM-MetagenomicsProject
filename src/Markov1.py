#Author: Jacob O'bott
#File: Markov_1.py
#Last edit: 06-27-14
#Authors: Jacob O'Bott, Grace Chandler
#Description:
#Input sequences and they will be scanned and given Markov-1 background probabilities
#Need to split into multiple functions

#Function used to find the background probabilities when a Markov-1 model
def background_probs(sequence):

    #A way where the user is asked what length of a sliding window that they want to use
    #Also makes the window length an int
    
    input_check = False
    while input_check == False:
        window = int(raw_input("Enter window size: "))
        #This is a way to make sure that the window inputted is appropriate
        if window > 0:
            input_check = True
        else:
            print 'Invalid sliding window: please enter a number greater than zero'
            
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

    #Calculating the counts of the transitions in the sequence inputted
    #Iterates over the all the bases in a sequence except the last one
    #Not the last one becasue it does not have a base after it,
    #hence no transition between bases that start at that base
    for i in range(sequence_length-1):
        #The first set of if/elif statements asks whether the
        # intitial base is A, T, C or G
        if sequence[i] == 'A':
            #The second set of if/elif statements is nested inside
            #each of the initial if/elif statements
            #It catagorizes what the base following the initial base is
            #From the information of the base and the base following it,
            #the computer then determiines what transition is observed
            if sequence[i+1] == 'A':
                #The following then adds one to the value of the count of the
                #appropriate base transition occuring
                counts['AA_count'] += 1
            elif sequence[i+1] == 'C':
                counts['AC_count'] += 1
            elif sequence[i+1] == 'T':
                counts['AT_count'] += 1
            elif sequence[i+1] == 'G':
                counts['AG_count'] += 1
        elif sequence[i] == 'C':
            if sequence[i+1] == 'A':
                counts['CA_count'] += 1
            elif sequence[i+1] == 'C':
                counts['CC_count'] += 1
            elif sequence[i+1] == 'T':
                counts['CT_count'] += 1
            elif sequence[i+1] == 'G':
                counts['CG_count'] += 1
        elif sequence[i] == 'T':
            if sequence[i+1] == 'A':
                counts['TA_count'] += 1
            elif sequence[i+1] == 'C':
                counts['TC_count'] += 1
            elif sequence[i+1] == 'T':
                counts['TT_count'] += 1
            elif sequence[i+1] == 'G':
                counts['TG_count'] += 1
        elif sequence[i] == 'G':
            if sequence[i+1] == 'A':
                counts['GA_count'] += 1
            elif sequence[i+1] == 'C':
                counts['GC_count'] += 1
            elif sequence[i+1] == 'T':
                counts['GT_count'] += 1
            elif sequence[i+1] == 'G':
                counts['GG_count'] += 1

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
    for key in counts:
        print key,
        print '='
        print counts[key]


    #Making the counts a float
    for key in counts:
        counts[key] = float(counts[key])

    #Printing out the sequence length
    print ' '
    print 'The length of the sequence was %s' %sequence_length

    #Defining variables that are the values of counts with the same first base (A,C,T,G)
    AX_count = counts['AA_count'] + counts['AC_count'] + counts['AT_count'] + counts['AG_count']
    CX_count = counts['CA_count'] + counts['CC_count'] + counts['CT_count'] + counts['CG_count']
    TX_count = counts['TA_count'] + counts['TC_count'] + counts['TT_count'] + counts['TG_count']
    GX_count = counts['GA_count'] + counts['GC_count'] + counts['GT_count'] + counts['GG_count']

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


x = 'AAGTAATCGAGCTACATAGAATATCTGTTCACCCTCGGGGAGCGTGGGGTGTAC'
list_of_back_probs = background_probs(x)
print list_of_back_probs

