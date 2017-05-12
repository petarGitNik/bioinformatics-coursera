#!/usr/bin/python

# example how to run this program in shell:
# ./replication.py < datasets/PatternCount.txt
# ./replication.py [output_file_name] < datasets/[input_file_name]

# Problem: Count the number of occurences of a Pattern within a Text. This is a
#          sliding window problem, so it should count overlapping occurences.
# Input: String Text and Pattern
# Output: Count(Text, Pattern)
def PatternCount(Text, Pattern):
    count = 0
    n = len(Text)
    k = len(Pattern)
    # or for i in range(n - k): ? there is no +1 in original algorithm
    # it does not count the last pattern! So it wont work! Look at the
    # debug datasets file.
    for i in range(n-k+1):
        if Text[i:i+k] == Pattern:
            count += 1
    return count

# Problem: Find the most frequent k-mers in a string
# Input: A string Text and an integer k
# Output: All most frequent k-mers in Text
def CountDict(Text, k):
    """
    The index in the count array corresponds to the position of the i-th pattern
    in a text.
    """
    Count = {}
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        Count[i] = PatternCount(Text, Pattern)
    return Count

def RemoveDuplicates(Words):
    """
    Remove duplicates from the Words array. Return array where every word is
    unique.
    """
    UniqueWords = []
    for word in Words:
        if word not in UniqueWords:
            UniqueWords.append(word)
    return UniqueWords

# Charging station - The Frequency Array:
def PatternToNumber(Pattern):
    code_book = { 'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3 }
    coded_pattern = 0
    for exponent, letter in enumerate(Pattern[::-1]):
        coded_pattern += code_book[letter] * (4 ** exponent)
    return coded_pattern

def NumberToPattern(index, k):
    """
    Convert a number from a base-10 to a string of characters from a 4-letter
    alphabet. That is, from a base-10 number to a base-4 'number'. If a resulting
    number from division is less than k signs long, add padding (with zeroes).
    For example if NumberToPattern([some_number], 3) = G, add two 'A's in front
    of it i.e. AAG is the resulting pattern.
    """
    decode_book = { 0 : 'A', 1 : 'C', 2 : 'G', 3 : 'T' }
    residues = []
    decoded_pattern = ''
    while index >= 4:
        index_temp = index / 4
        residues.append(index % 4)
        index = index_temp
    # Add the final value resulting from division/conversion
    residues.append(index)
    while len(residues) < k:
        residues.append(0)
    for residue in residues[::-1]:
        decoded_pattern += decode_book[residue]
    return decoded_pattern

def ComputingFrequencies(Text, k):
    FrequencyArray = {}
    for i in range(4**k):
        FrequencyArray[i] = 0
    # Alternatively: [0] * (4**k) = [0, ..., 0] where len([given_list]) = 4**k
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        j = PatternToNumber(Pattern)
        FrequencyArray[j] += 1
    return FrequencyArray

def FasterFrequentWords(Text, k):
    """
    Find frequent words using a frequency array. This algorithm is faster for
    smaller values of k.
    """
    pass

def FrequentWords(Text, k):
    """
    Take every k-mer from a text, and count the number of its occurences in the
    text. Find the maximum occuring kmers and, according to their indices, form
    an array of frequent patterns. Remove any duplicates and return the array.
    """
    # Why does it fails the debug dataset 2? See: FrequentWordsProblem.pdf
    FrequentPatterns = []
    Count = CountDict(Text, k)
    maxCount = max(Count.values())
    for i in Count:
        if Count[i] == maxCount:
            FrequentPatterns.append(Text[i:i+k])
    FrequentPatternsWithoutDuplicates = RemoveDuplicates(FrequentPatterns)
    return FrequentPatternsWithoutDuplicates

def FrequentWordsWithDuplicates(Text, k):
    FrequentPatterns = []
    Count = CountDict(Text, k)
    maxCount = max(Count.values())
    for i in Count:
        if Count[i] == maxCount:
            FrequentPatterns.append(Text[i:i+k])
    return FrequentPatterns

# Problem: Find a reverse complement of a DNA string
# Input: A DNA string named 'Pattern'
# Output: 'Pattern_rc' i.e. the reverse complement of 'Pattern'
def ReverseStrand(Pattern):
    return Pattern[::-1]

def ReverseComplement(Pattern):
    return ReverseStrand(Pattern).replace('A', 't').replace('T', 'a').replace('G', 'c').replace('C', 'g').upper()

# Problem: Find all occurences of a pattern in a string
# Input: Strings Pattern and Genome
# Output: All starting positions in Genome Where Pattern appears as a substring
def PatternMatching(Pattern, Genome):
    k = len(Pattern)
    n = len(Genome)
    position = []
    for i in range(n-k+1):
        if Genome[i:i+k] == Pattern:
            position.append(i)
    return position

"""
The following function does not search for a specific clump of a pattern, but
rather all patterns that form clumps in genome. It is hope that one of those
clumps will reveal the location of 'ori'.
Question: Wich of these clumps would reveal 'ori'? Or, will *any* of these
          clumps reveal 'ori'?
Answer: We don't know.
Remark: We know that DnaA boxes are usually 9 nucleoties long.
Following Q and A are from the course Authors:
Question: How does DnaA know which DnaA box to bind to?
Answer: DnaA does not necessarily bind to just one DnaA box. In fact, it may
        bind to all of them.

A k-mer clump is a called such if it appears many times within a short interval
of genome. A k-mer 'Pattern' forms an (L,t)-clump inside a 'Genome' if there is
an interval of 'Genome' of length 'L' in which k-mer appears at least 't' times.

Remark: It can be solved in several ways, by using:
        - FrequentWords
        - FasterFrequentWords i.e. using frequency array
        - i.e. finding frequent words by sorting
"""
# Problem: Find patterns forming clumps in a string
# Input: A string Genome, and integers k, L, and t
# Output: All distinct k-mers forming (L,t)-clumps in Genome
def Clump(Genome, k, L, t):
    clump = []
    for i in range(len(Genome) - L + 1):
        clumpCandidates = FrequentWordsWithDuplicates(Genome[i:i+L], k)
        # ispitaj mi ovaj clump
        # napravi mi recnik svih k-mera koje sam dobio
        potentialClumps = {}
        for pattern in RemoveDuplicates(clumpCandidates):
            potentialClumps[pattern] = 0
        # prebroj im njihovu ucestanost
        for pattern in clumpCandidates:
            potentialClumps[pattern] += 1
        # ispitaj da li je neki od njih najmanje t puta tu, i zapamti
        for pattern in potentialClumps:
            if potentialClumps[pattern] >= t:
                clump.append(pattern)
    return RemoveDuplicates(clump)

    # pogledaj prozor i:i+L
    # pronadji u njemu sve najcesce k-mere
    # pogledaj da li se neki od njih sadrzi vise ili jednako od 't' puta
    # ako se nalazi, zapamti ga
    # ponovi proces

# --- Program testing --- #
import sys
lines = sys.stdin.read().splitlines()

with open('./outputs/' + sys.argv[1], 'w') as f:
    # --- PatternCount --- #
    #f.write(str(PatternCount(lines[0], lines[1])) + '\n')

    # --- FrequentWords --- #
    #f.write(' '.join(FrequentWords(lines[0], int(lines[1]))) + '\n')

    # --- ComputingFrequencies --- #
    #listOfFrequencies = ComputingFrequencies(lines[0], int(lines[1])).values()
    #f.write(' '.join([str(num) for num in listOfFrequencies]) + '\n')

    # --- ReverseComplement --- #
    #f.write(ReverseComplement(lines[0]) + '\n')

    # --- Pattern Matching --- #
    #listOfIndices = PatternMatching(lines[0], lines[1])
    #f.write(' '.join([str(num) for num in listOfIndices]) + '\n')

    # --- Clump Finding Problem --- #
    k, L, t = [int(num) for num in lines[1].split()]
    f.write(' '.join(Clump(lines[0], k, L, t)) + '\n')
