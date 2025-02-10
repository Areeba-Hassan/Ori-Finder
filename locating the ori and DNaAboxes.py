#checking frequencies of C's in half starnds of bac genome
#Ori is where the minimum number of C is met because that's where the reverse half strand transitions to forward half strand
def patcount(symbol, genome):  # counts the occurences of a given symbol in the genome
    count = 0
    for i in range(len(genome) - len(symbol) + 1):
        if genome[i:i + len(symbol)] == symbol:
            count = count + 1
    return count

def symbolarray(genome, symbol):
    #generates a frequency map for a given symbol in the genome
    array = {}
    n = len(genome)
    ext_genome = genome + genome [0:(n//2)]
    #initializing first index as the number of Cs in first position of the sliding window
    array[0] = patcount(symbol, genome[0:n//2])
    #sliding the window
    for i in range(1,n - (n//2)):
    #assigning array[i-1] to array[i]
        array[i] = array[i-1]
        #incrementing or decrementing array[i] by 1 based on the new elemnt inside the sliding window
        if ext_genome[i-1] == symbol:
            array[i] = array[i]-1
        if ext_genome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

#finding the minimum skew
#minimum skew is where G-C is minimum in the genome
#C-G decreseas along the reverse half strand and increases along the positive half strand
#so it is minimum at the transition point, which is where the ORI is supposedly located

def skew(genome):
    #computes G-C along the genome
    skew = [0]
    n = len(genome)
    for i in range (n):
        if genome[i] == "C":
            skew.append(skew[i]-1)
        elif genome[i] == "G":
            skew.append(skew[i]+1)
        else:
            skew.append(skew[i])
    return skew

def min_skew(genome):
    #finds position with minimum skew, which are likely to be ori
    pos = []
    raw_skew = skew(genome)
    m = min(raw_skew)
    for i in range(len(raw_skew)):
        if raw_skew[i] == m:
            pos.append(i)
    return pos

#generating a frequency map that contains all k-mers in a sequence and their number of occurences
def freqmap (seq, k): #where k is the size of k-mer or pattern
    freq = {}  # initializing empty dictionary
    n = len(seq)
    for i in range(n - k + 1):
        pattern = seq[i:i + k]
        freq[pattern] = freq.get(pattern, 0) + 1
    return freq

def getfreqword(seq, k):
    #fetching the most frequent kmer
    #by passing ori as seq, and probable lenght of DNaA box as k, we'll get candidate DNAaboxes
    freq = freqmap(seq, k)
    max_freq = max(freq.values())
    return [key for key in freq if freq[key] == max_freq]

def rev (seq):
    r = ''
    for char in seq:
        r = char + r
    return r

def comp (seq):
    c = ''
    for char in seq:
        if char == "T" or char == "t":
            c = c + "A"
        elif char == "A" or char == "a":
            c = c + "T"
        elif char == "C" or char == "c":
            c = c + "G"
        elif char == "G" or char == "g":
            c = c + "C"
    return c

def revcomp (seq):
    #finding reverse complement of a sequence
    pat = comp(rev(seq))
    return pat

def ham_dist (p,q):
    #finding the hamming distance
    #since biologically, DNaA can bind even to imperfect sequences
    #the hamming distance defines the amount of mismatches below which an imperfect sequence will be considered a DNaA box candidate
    dist = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            dist += 1
    return dist

def app_pat_match (text, pat, ham):
    # fetching the DNaA boxes that satisfy the hamming distance threshold
    pos = []
    for i in range (len(text)-len(pat)+1):
        if ham_dist (pat, text[i:i+len(pat)]) <= ham:
            pos.append(i)
    return pos

def find_ori_and_dnaboxes(genome, k , ham):
    print ("Finding Ori candidates...")
    ori_positions = min_skew(genome)
    if not ori_positions:
        raise ValueError("No minimum skew positions found. Check input sequence.")
    ori_start = ori_positions[0]
    ori_end = min(ori_start + 500, len(genome)) #assuming the window size to be 500
    ori_region = genome[ori_start:ori_end]

    print(f"Ori region located at {ori_start} to {ori_end}")

    print("Identifying DNaAboxes...")
    candidate_dnaboxes = getfreqword(ori_region,k)

    print("Identifying reverse complements...")
    dnabox_pairs = [(box, revcomp(box)) for box in candidate_dnaboxes]

    print("Finding matches in Ori...")
    confirmed_dnaboxes = {}
    for box, rev_box in dnabox_pairs:
        positions = app_pat_match(ori_region, box, ham)
        rev_positions = app_pat_match(ori_region, rev_box, ham)

        if positions or rev_positions:
            confirmed_dnaboxes[box] = positions
            confirmed_dnaboxes[rev_box] = rev_positions

    print("Results:")
    print(f"Ori: {ori_start} to {ori_end}")
    print(f"Candidate DNaA boxes: {confirmed_dnaboxes}")
    return ori_start, ori_end, confirmed_dnaboxes

#Test variables
genome_sequence = "ACGTACGTGCGCGATATGCGCGATATGCATCGATCGTGCATGCATG" #genomic sequence
kmer_length = 9  #Expected DNaA box length
hamming_threshold = 2  #mismatches allowed

find_ori_and_dnaboxes(genome_sequence, kmer_length, hamming_threshold)
