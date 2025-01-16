#generating a frequency map that contains all k-mers in a sequence and their numbers of occurence
def freqmap (seq, k): #where k is the size of k-mer or pattern
    freq = {} #initializing empty dictionary
    n = len(seq)
    for i in range (n-k+1): #make sure the range is error free
        pattern = seq[i:i+k]
        freq[pattern] = 0 #initializing the occurence of each k-mer with 0
    for i in range (n-k+1):
        pattern = seq[i:i+k]
        freq[pattern]=freq[pattern]+1 #incrementing the occurence of every k-mer by 1 if it is detected in the sequence
    return freq

def getfreqword (seq, k):
    words = []
    freq = freqmap(seq, k) #calling freqmap as a subroutine
    m = max(freq.values()) #finding the maximum of occurences
    for key in freq:
        if freq[key] == m: #finding the k-mer with the highest occerence
            words.append(key)
    return words

print(getfreqword("ATCGCGGATAGGCGTAGGATAGCTAGACTGAT", 3))
#in the above function, pass the ori sequence of a bacterium as seq, and the expected length of its DNaA box as k
#and we will have the most frequently occuring patterns as candidate DNaA boxes

#now, finding the reverse complement of a pattern/sequence
def revcomp (seq):
    pat = comp(rev(seq))
    return pat

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

print(revcomp("ATCGCG")) '''pass the most frequently occuring k-mers found using getfreqwords()
and check if any of them are reverse complements of one another
The most repeatedly occuring k-mer and its reverse-coplemnetary seq would now be probable DNaA boxes
however, to make sure that they really are a special characteristic of the ori, and not just patterns which are
common all across the genome, we will check to see how many times they occur in the entire bacterial genome'''

#writing function to match a pattern to the genome
def patternmatch (pat, genome):
    positions = []
    pat_len = len(pat)
    gen_len = len(genome)
    for i in range (gen_len - pat_len + 1):
        if genome[i:i+pat_len] == pat:
            positions.append(i) #storing positions of pattern in the genome
    return positions

print (patternmatch("ATCGCG", "ATCGCGGATAGGCGTAGGATAGCTAGACTGAT"))
#in the above function, pass the most frequently occuring pattern and the genomic sequence
#if the pattern occuring frequently in the ori is not so much frequent in the rest of the genome
#we can further confirm that it might be the DNaA box sequence