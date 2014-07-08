#__author__ = 'jobott'
#Name:
#Inputs:
#Outputs:
#Description:


#Function to determine the PSFM
def PSFM(binding_sites,sequence):
#Note: the binding_sites variable here is a list of the binding sites, all with the same length

    #Defining lists to store counts for the rows of the PSFM
    a_bases = []
    c_bases = []
    t_bases = []
    g_bases = []

    #Definging some other important variables
    total_sites = float(len(binding_sites))
    site_length = float(len(binding_sites[0]))
    window = int(site_length)

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
            if site[base] == 'a':
                a_bases[base] += 1
            if site[base] == 'c':
                c_bases[base] += 1
            if site[base] == 't':
                t_bases[base] += 1
            if site[base] == 'g':
                g_bases[base] += 1

    #Making the counts into frequencies
    a_freq = [float(float(x)/total_sites) for x in a_bases]
    c_freq = [float(float(x)/total_sites) for x in c_bases]
    t_freq = [float(float(x)/total_sites) for x in t_bases]
    g_freq = [float(float(x)/total_sites) for x in g_bases]


    #Printing the PSFM
    print 'Obtanied PSFM:'
    print a_freq
    print c_freq
    print t_freq
    print g_freq

    #!!!!!Markov-1 Background forming!!!!!
    aa_count = 0
    ac_count = 0
    at_count = 0
    ag_count = 0
    ca_count = 0
    cc_count = 0
    ct_count = 0
    cg_count = 0
    ta_count = 0
    tc_count = 0
    tt_count = 0
    tg_count = 0
    ga_count = 0
    gc_count = 0
    gt_count = 0
    gg_count = 0

    #Defining some other variables
    errors_count = 0
    seq_length = len(sequence)

    #Calculating the counts of the transitions in the sequence inputted
    for i in range(seq_length-1):
        if sequence[i] == 'a' and sequence[i+1] == 'a':
            aa_count += 1
        elif sequence[i] == 'a' and sequence[i+1] == 'c':
            ac_count += 1
        elif sequence[i] == 'a' and sequence[i+1] == 't':
            at_count += 1
        elif sequence[i] == 'a' and sequence[i+1] == 'g':
            ag_count += 1
        elif sequence[i] == 'c' and sequence[i+1] == 'a':
            ca_count += 1
        elif sequence[i] == 'c' and sequence[i+1] == 'c':
            cc_count += 1
        elif sequence[i] == 'c' and sequence[i+1] == 't':
            ct_count += 1
        elif sequence[i] == 'c' and sequence[i+1] == 'g':
            cg_count += 1
        elif sequence[i] == 't' and sequence[i+1] == 'a':
            ta_count += 1
        elif sequence[i] == 't' and sequence[i+1] == 'c':
            tc_count += 1
        elif sequence[i] == 't' and sequence[i+1] == 't':
            tt_count += 1
        elif sequence[i] == 't' and sequence[i+1] == 'g':
            tg_count += 1
        elif sequence[i] == 'g' and sequence[i+1] == 'a':
            ga_count += 1
        elif sequence[i] == 'g' and sequence[i+1] == 'c':
            gc_count += 1
        elif sequence[i] == 'g' and sequence[i+1] == 't':
            gt_count += 1
        elif sequence[i] == 'g' and sequence[i+1] == 'g':
            gg_count += 1
        else:
            errors_count += 1


    #Variable with the total number of transitions observed
    trans_total = aa_count+ac_count+at_count+ag_count+ca_count+cc_count+ct_count+cg_count+ta_count+tc_count+tt_count+tg_count+ga_count+gc_count+gt_count+gg_count+errors_count

    aa_count = round(aa_count,7)
    ac_count = round(ac_count,7)
    at_count = round(at_count,7)
    ag_count = round(ag_count,7)
    ca_count = round(ca_count,7)
    cc_count = round(cc_count,7)
    ct_count = round(ct_count,7)
    cg_count = round(cg_count,7)
    ta_count = round(ta_count,7)
    tc_count = round(tc_count,7)
    tt_count = round(tt_count,7)
    tg_count = round(tg_count,7)
    ga_count = round(ga_count,7)
    gc_count = round(gc_count,7)
    gt_count = round(gt_count,7)
    gg_count = round(gg_count,7)

    if aa_count+ac_count+at_count+ag_count != 0:
        prob_aa = float(aa_count/(aa_count+ac_count+at_count+ag_count))
    if aa_count+ac_count+at_count+ag_count != 0:
        prob_ac = float(ac_count/(aa_count+ac_count+at_count+ag_count))
    if aa_count+ac_count+at_count+ag_count != 0:
        prob_at = float(at_count/(aa_count+ac_count+at_count+ag_count))
    if aa_count+ac_count+at_count+ag_count != 0:
        prob_ag = float(ag_count/(aa_count+ac_count+at_count+ag_count))
    if ca_count+cc_count+ct_count+cg_count != 0:
        prob_ca = float(ca_count/(ca_count+cc_count+ct_count+cg_count))
    if ca_count+cc_count+ct_count+cg_count != 0:
        prob_cc = float(cc_count/(ca_count+cc_count+ct_count+cg_count))
    if ca_count+cc_count+ct_count+cg_count != 0:
        prob_ct = float(ct_count/(ca_count+cc_count+ct_count+cg_count))
    if ca_count+cc_count+ct_count+cg_count != 0:
        prob_cg = float(cg_count/(ca_count+cc_count+ct_count+cg_count))
    if ta_count+tc_count+tt_count+tg_count != 0:
        prob_ta = float(ta_count/(ta_count+tc_count+tt_count+tg_count))
    if ta_count+tc_count+tt_count+tg_count != 0:
        prob_tc = float(tc_count/(ta_count+tc_count+tt_count+tg_count))
    if ta_count+tc_count+tt_count+tg_count != 0:
        prob_tt = float(tt_count/(ta_count+tc_count+tt_count+tg_count))
    if ta_count+tc_count+tt_count+tg_count != 0:
        prob_tg = float(tg_count/(ta_count+tc_count+tt_count+tg_count))
    if ga_count+gc_count+gt_count+gg_count != 0:
        prob_ga = float(ga_count/(ga_count+gc_count+gt_count+gg_count))
    if ga_count+gc_count+gt_count+gg_count != 0:
        prob_gc = float(gc_count/(ga_count+gc_count+gt_count+gg_count))
    if ga_count+gc_count+gt_count+gg_count != 0:
        prob_gt = float(gt_count/(ga_count+gc_count+gt_count+gg_count))
    if ga_count+gc_count+gt_count+gg_count != 0:
        prob_gg = float(gg_count/(ga_count+gc_count+gt_count+gg_count))

    print ' '
    print 'Conditional probabilities:'
    if aa_count+ac_count+at_count+ag_count != 0:
       print 'P(A|A) = %s' %prob_aa
    if aa_count+ac_count+at_count+ag_count != 0:
        print 'P(C|A) = %s' %prob_ac
    if aa_count+ac_count+at_count+ag_count != 0:
        print 'P(T|A) = %s' %prob_at
    if aa_count+ac_count+at_count+ag_count != 0:
        print 'P(G|A) = %s' %prob_ag
    if ca_count+cc_count+ct_count+cg_count != 0:
        print 'P(A|C) = %s' %prob_ca
    if ca_count+cc_count+ct_count+cg_count != 0:
        print 'P(C|C) = %s' %prob_cc
    if ca_count+cc_count+ct_count+cg_count != 0:
        print 'P(T|C) = %s' %prob_ct
    if ca_count+cc_count+ct_count+cg_count != 0:
        print 'P(G|C) = %s' %prob_cg
    if ta_count+tc_count+tt_count+tg_count != 0:
        print 'P(A|T) = %s' %prob_ta
    if ta_count+tc_count+tt_count+tg_count != 0:
        print 'P(C|T) = %s' %prob_tc
    if ta_count+tc_count+tt_count+tg_count != 0:
        print 'P(T|T) = %s' %prob_tt
    if ta_count+tc_count+tt_count+tg_count != 0:
        print 'P(G|T) = %s' %prob_tg
    if ga_count+gc_count+gt_count+gg_count != 0:
        print 'P(A|G) = %s' %prob_ga
    if ga_count+gc_count+gt_count+gg_count != 0:
        print 'P(C|G) = %s' %prob_gc
    if ga_count+gc_count+gt_count+gg_count != 0:
        print 'P(T|G) = %s' %prob_gt
    if ga_count+gc_count+gt_count+gg_count != 0:
        print 'P(G|G) = %s' %prob_gg

    count = 0
    list_of_back_probs = []
    errors_markov = 0

    for i in range(len(sequence)-window):
        subseq = sequence[i:i+window+1]
        '''subseqs.append(subseq)
        print subseqs'''
        m1_prob = 1

        if i == 0:
            subseq = subseq[::-1]
            for i in range(len(subseq)-1):
                if subseq[i] == 'a' and subseq[i+1] == 'a':
                    m1_prob *= prob_aa
                elif subseq[i] == 'a' and subseq[i+1] == 'c':
                    m1_prob *= prob_ac
                elif subseq[i] == 'a' and subseq[i+1] == 't':
                    m1_prob *= prob_at
                elif subseq[i] == 'a' and subseq[i+1] == 'g':
                    m1_prob *= prob_ag
                elif subseq[i] == 'c' and subseq[i+1] == 'a':
                    m1_prob *= prob_ca
                elif subseq[i] == 'c' and subseq[i+1] == 'c':
                    m1_prob *= prob_cc
                elif subseq[i] == 'c' and subseq[i+1] == 't':
                    m1_prob *= prob_ct
                elif subseq[i] == 'c' and subseq[i+1] == 'g':
                    m1_prob *= prob_cg
                elif subseq[i] == 't' and subseq[i+1] == 'a':
                    m1_prob *= prob_ta
                elif subseq[i] == 't' and subseq[i+1] == 'c':
                    m1_prob *= prob_tc
                elif subseq[i] == 't' and subseq[i+1] == 't':
                    m1_prob *= prob_tt
                elif subseq[i] == 't' and subseq[i+1] == 'g':
                    m1_prob *= prob_tg
                elif subseq[i] == 'g' and subseq[i+1] == 'a':
                    m1_prob *= prob_ga
                elif subseq[i] == 'g' and subseq[i+1] == 'c':
                    m1_prob *= prob_gc
                elif subseq[i] == 'g' and subseq[i+1] == 't':
                    m1_prob *= prob_gt
                elif subseq[i] == 'g' and subseq[i+1] == 'g':
                    m1_prob *= prob_gg
                else:
                    errors_markov += 1
            list_of_back_probs.append(m1_prob)
            print subseq
            subseq = subseq[::-1]
            m1_prob = 1
            count += 1



        for i in range(len(subseq)-1):
            if subseq[i] == 'a' and subseq[i+1] == 'a':
                m1_prob *= prob_aa
            elif subseq[i] == 'a' and subseq[i+1] == 'c':
                m1_prob *= prob_ac
            elif subseq[i] == 'a' and subseq[i+1] == 't':
                m1_prob *= prob_at
            elif subseq[i] == 'a' and subseq[i+1] == 'g':
                m1_prob *= prob_ag
            elif subseq[i] == 'c' and subseq[i+1] == 'a':
                m1_prob *= prob_ca
            elif subseq[i] == 'c' and subseq[i+1] == 'c':
                m1_prob *= prob_cc
            elif subseq[i] == 'c' and subseq[i+1] == 't':
                m1_prob *= prob_ct
            elif subseq[i] == 'c' and subseq[i+1] == 'g':
                m1_prob *= prob_cg
            elif subseq[i] == 't' and subseq[i+1] == 'a':
                m1_prob *= prob_ta
            elif subseq[i] == 't' and subseq[i+1] == 'c':
                m1_prob *= prob_tc
            elif subseq[i] == 't' and subseq[i+1] == 't':
                m1_prob *= prob_tt
            elif subseq[i] == 't' and subseq[i+1] == 'g':
                m1_prob *= prob_tg
            elif subseq[i] == 'g' and subseq[i+1] == 'a':
                m1_prob *= prob_ga
            elif subseq[i] == 'g' and subseq[i+1] == 'c':
                m1_prob *= prob_gc
            elif subseq[i] == 'g' and subseq[i+1] == 't':
                m1_prob *= prob_gt
            elif subseq[i] == 'g' and subseq[i+1] == 'g':
                m1_prob *= prob_gg
            else:
                errors_markov += 1
        list_of_back_probs.append(m1_prob)
        print subseq
        m1_prob = 1
        count += 1

    for start in range(len(list_of_back_probs)):
        print 'P_back(bases %s'%(start+1), '- %s)='%(start+window),  '%s'%list_of_back_probs[start]





    #!!!!!Start of calculating P_motif!!!!!!
    list_of_motif_probs = []

    for i in range(len(sequence)-window):
        subseq = sequence[i:i+window]
        motif_prob = 1
        if i == 0:
            subseq = subseq[::-1]
            for i in range(len(subseq)):
                if subseq[i] == 'a':
                    motif_prob *= a_freq[i]
                elif subseq[1] == 'c':
                    motif_prob *= c_freq[i]
                elif subseq[i] == 't':
                    motif_prob *= t_freq[i]
                elif subseq[i] == 'g':
                    motif_prob *= g_freq[i]
            list_of_motif_probs.append(motif_prob)
            subseq = subseq[::-1]
            motif_prob = 1

        for i in range(len(subseq)):
            if subseq[i] == 'a':
                motif_prob *= a_freq[i]
            elif subseq[1] == 'c':
                motif_prob *= c_freq[i]
            elif subseq[i] == 't':
                motif_prob *= t_freq[i]
            elif subseq[i] == 'g':
                motif_prob *= g_freq[i]
        list_of_motif_probs.append(motif_prob)
        motif_prob = 1

    for start in range(len(list_of_motif_probs)):
        print 'P_motif(bases %s'%(start+1), '- %s)='%(start+window),  '%s'%list_of_motif_probs[start]




    #!!!!!Combining the P_back and P_motif and forming the scores!!!!!
    prob_ratios = []
    for i in range(len(list_of_motif_probs)):
        prob_ratios.append(list_of_motif_probs[i]/list_of_back_probs[i])

    PSSM_scores = []
    import math
    for i in prob_ratios:
        PSSM_scores.append(math.log(i,2))

    for start in range(len(PSSM_scores)):
        print 'PSSM_score(bases %s'%(start+1), '- %s)='%(start+window),  '%s'%PSSM_scores[start]











bsites = ['actgactg','ctgactga','tgactgac','gactgact','acctgaat','acctgaat','acccgatt','aactgtat']
x = 'aagtaaatcgagctacatagaatatctgttcaccctcggggagcgtggggtgtac'
PSFM(bsites,x)



