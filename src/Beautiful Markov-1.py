__author__ = 'jobott'

#making a comment to test out git 
def site_length(binding_sites):
    window = int(len(binding_sites[0]))
    return window


def PSFM(binding_sites):

    #Defining variables for the total number of putative binding sites(float) and the length of the first site
    total_sites = float(len(binding_sites))
    site_length = int(len(binding_sites[0]))

    #Defining lists to store the counts for each base (the rows of the PSFM) and
    #making the PSFM have the number of columns as the length of the binding sites
    #Starts the value off at zero for each to create pseudo-counts
    a_bases = []; a_bases.extend([float(1)]*site_length)
    c_bases = []; c_bases.extend([float(1)]*site_length)
    t_bases = []; t_bases.extend([float(1)]*site_length)
    g_bases = []; g_bases.extend([float(1)]*site_length)

    #Iterating over the binding sites (viewing each site individually)
    for site in binding_sites:

        #Iterating over the indexes for the bases within each binding site
        for base in range(len(site)):

            #Determines what the value of the base is at that position in the binding site (a,c,t,g)
            if site[base] == 'A':
                #Adding one to the value of the count stored at the same position of the base in the *_bases list
                a_bases[base] += 1
            elif site[base] == 'C':
                c_bases[base] += 1
            elif site[base] == 'T':
                t_bases[base] += 1
            elif site[base] == 'G':
                g_bases[base] += 1

    #Defining new lists which hold the frequencies that are determined from the counts
    #General format:
    #frequency_list = the value of the count of that base at that specific location divided
    #                 by the total number of binding sites overall; this value is stored as a float
    a_freq = [float(x/total_sites) for x in a_bases]
    c_freq = [float(x/total_sites) for x in c_bases]
    t_freq = [float(x/total_sites) for x in t_bases]
    g_freq = [float(x/total_sites) for x in g_bases]

    #Printing the PSFM
    #Printing the title
    print 'Obtanied PSFM:'
    #Printing each set of frequencies with their appropriate labels for their bases
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
    list_of_log_of_motif_probs = []

    #Determining the window length based off of the length of the PSFM motif
    window = int(len(PSFM[0]))

    #Iterating the indexes of all of the initial bases for each case that the sliding window encounters
    for initial in range(len(sequences)-window+1):
        #Creating a subsequence that starts at the given initial index and is as long as the pre-determined window
        subsequence = sequences[initial:initial+window]
        #Setting the initial value for the motif probability of the sliding window at 1
        motif_prob = 0
        #Importing math to use the log function below
        import math
        #Iterating over the indexes of the bases in the subsequence
        for base in range(len(subsequence)):
            #Determines what the value of the base is at that position in the subsequence (a,c,t,g)
            if subsequence[base] == 'A':
                #Multiply the value of the motif_prob by the value of the frequency of that base at the
                #position of that base within the PSFM
                motif_prob += math.log(PSFM[0][base],2)
            elif subsequence[base] == 'C':
                motif_prob += math.log(PSFM[1][base],2)
            elif subsequence[base] == 'T':
                motif_prob += math.log(PSFM[2][base],2)
            elif subsequence[base] == 'G':
                motif_prob += math.log(PSFM[3][base],2)
        #Appending the list that stores the motif probabilities by the value motif_prob
        list_of_log_of_motif_probs.append(motif_prob)
        #Re-setting the value of motif_prob to 1 to apply to the next subsequence
        motif_prob = 0

    #Iterating over the indexes of the values in list_of_motif_prob
    for start in range(len(list_of_log_of_motif_probs)):
        #Printing out the values of motif_prob with appropriate labeling
        #General format:
        #'P_motif(initial base - final base in window) = appropriate value of P_motif at the
        #                                                given point of the sliding window
        print 'log[P_motif(bases %s'%(start+1), '- %s)]='%(start+window),  '%s'%list_of_log_of_motif_probs[start]

    #Having the function return the list containing the motif probabilities
    return list_of_log_of_motif_probs



#Function used to find the background probabilities when using a Markov-1 model
#Inputs: sequence that is going to be used for both the background generation and to apply the sliding window on
#Outputs: A list of P(subsequence|background) for each subsequence fromed by a sliding window
def background_probs(binding_sites,sequence):

    #A way where the user is asked what length of a sliding window that they want to use
    #Also makes the window length an int
    window = int(len(binding_sites[0]))

    #Might want to make this a function like def initialize_counts()
    #Making a dictionary to hold the counts of all possible base transitions
    counts = {}
    bases = ['A','C','T','G']

    for base_1 in bases:
        for base_2 in bases:
            key = (base_1,base_2)
            #Using zero to start off the counts
            counts[key] = 0

    #Defining a variable for sequence length
    sequence_length = len(sequence)

    #Calculating the counts of the transitions in the sequence inputted
    #Iterates over the all the bases in a sequence except the last one
    #Not the last one becasue it does not have a base after it,
    #hence no transition between bases that start at that base
    for i in range(sequence_length-1):
        counts[(sequence[i],sequence[i+1])] += 1

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
    AX_count = counts[('A','A')] + counts[('A','C')] + counts[('A','T')] + counts[('A','G')]
    CX_count = counts[('C','A')] + counts[('C','C')] + counts[('C','T')] + counts[('C','G')]
    TX_count = counts[('T','A')] + counts[('T','C')] + counts[('T','T')] + counts[('T','G')]
    GX_count = counts[('G','A')] + counts[('G','C')] + counts[('G','T')] + counts[('G','G')]

    #Makingg a dictionary to hold all the values of the background probabilities
    #Note: prob_AC = P(C|A), the conditional probability of getting C as the second base
    #given that A is the first base; prob_TG = P(G|T); ect.
    background_probs = {}
    for base_1 in bases:
        for base_2 in bases:
            key = (base_1,base_2)
            if base_1 == 'A':
                background_probs[key] = counts[key]/AX_count
            if base_1 == 'C':
                background_probs[key] = counts[key]/CX_count
            if base_1 == 'T':
                background_probs[key] = counts[key]/TX_count
            if base_1 == 'G':
                background_probs[key] = counts[key]/GX_count

    #Printing out the conditional probabilities
    print ' '
    print 'Conditional probabilities:'
    for key in background_probs:
        print key,
        print '=',
        print background_probs[key]

    #Overall section scanning an input sequence with a sliding window and
    #applying the Markov-1 background probabilities

    #Defining a list to store the background probabilities for each window observed
    list_of_log_of_back_probs = []

    #Importimg the math function to take the log below
    import math

    #Iterating over all the bases that would be the start of a window
    #Hence the use of sequence_length - window becasue any of the bases within
    #a window length of the end would not have a full window length of window
    #of bases following it
    for i in range(sequence_length-window):
        #Defines a subsequence wich is the length of the window and starts at the first base
        subseq = sequence[i:i+window+1]
        #Defining a variable to hold the value of the Markov-1 background
        #probability for the given window
        m1_prob = 0
        #Only goes through this if it is the first term of the sequence:
        #Reverses the order of the sequence and then applies the markov-1 background
        #probability model to it. (This is explained in the comments on the next section)
        if i == 0:
            #Inverting the sequence order
            subseq = subseq[::-1]
            for i in range(len(subseq)-1):
                m1_prob += math.log(background_probs[(subseq[i],subseq[i+1])],2)
            list_of_log_of_back_probs.append(m1_prob)
            print subseq
            #Re-inverting the sequence order so that it is the same as before
            subseq = subseq[::-1]
            #Re-setting m1_prob to 1 for the markov-1 model to be applied again
            m1_prob = 0

        #Iterating over the bases in the formed subsequence
        for i in range(len(subseq)-1):
                m1_prob += math.log(background_probs[(subseq[i],subseq[i+1])],2)
        #After iterating through the subsequence and obtaining its m1_prob,
        #apending that to our list of background probabilities
        list_of_log_of_back_probs.append(m1_prob)
        #Re-setting m1_prob to 1 so that the next subsequence's m1_prob can be calculated
        m1_prob = 0

    #Printing the background probabilities for each subsequence
    print ' '
    for start in range(len(list_of_log_of_back_probs)):
        print 'log[P_back(bases %s'%(start+1), '- %s)]='%(start+window),  '%s'%list_of_log_of_back_probs[start]

    #Having the function return the list containing the Markov-1 background probabilities
    return list_of_log_of_back_probs



#Forming the PSSM
#Inputs: sequence, list of binding sites
#Outputs: PSSM scores for each subsequence in the sliding window
def PSSM(binding_sites,sequence):

    #Defining the length of the window
    window = site_length(binding_sites)

    #Forming the motif probabilities
    list_of_log_of_motif_probs = P_motif(PSFM(binding_sites),sequence)

    #Forming the background probabilities
    list_of_log_of_back_probs = background_probs(binding_sites,sequence)

    #defining variable for the number of subsequences used
    number_of_subseq = len(list_of_log_of_motif_probs)

    #Creating a list to store the PSSM scores in
    PSSM_scores = []
    #Iterating over the ratios in prob_ratios
    for initial in range(number_of_subseq):
        #Appending the PSSM list with the PSSM score, determined by taking the log base 2 for the ratio
        PSSM_scores.append(list_of_log_of_motif_probs[initial]-list_of_log_of_back_probs[initial])

    #Iterating over the indexes of the values in list_of_motif_prob
    for start in range(number_of_subseq):
        #Printing out the PSSM values with appropriate labeling
        #General format:
        #'PSSM_score(initial base - final base in window) = appropriate PSSM value at the
        #                                                   given point of the sliding window
        print 'PSSM_score(bases %s'%(start+1), '- %s)='%(start+window),  '%s'%PSSM_scores[start]

    #Returning the PSSM scores list
    return PSSM_scores



#Test example
bsites = ['ACTGACTG','CTGACTGA','TGACTGAC','GACTGACT','ACCTGAAT','ACCTGAAT','ACCCGATT','AACTGTAT']
x = 'AAGTAAATCGAGCTACATAGAATATCTGTTCACCCTCGGGGAGCGTGGGGTGTAC'
PSSM(bsites,x)