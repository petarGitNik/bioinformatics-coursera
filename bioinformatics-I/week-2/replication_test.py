#!/usr/bin/python

# example how to run this program in shell:
# ./replication.py < datasets/PatternCount.txt
# ./replication.py [output_file_name] < datasets/[input_file_name]

# --- Week 1 --- #

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

# --- Week 2 --- #
from collections import OrderedDict

# Problem: Calculate the skew array i.e. the difference betwee 'C' and 'G'
# Input: Genome
# Output: An array of values which represents the difference between 'C' and 'G'
def Skew(Genome):
    skew = { 0 : 0 }
    for index, nucleotide in enumerate(Genome):
        if nucleotide == 'G':
            skew[index + 1] = skew[index] + 1
        elif nucleotide == 'C':
            skew[index + 1] = skew[index] - 1
        else:
            skew[index + 1] = skew[index]
    return OrderedDict(sorted(skew.items()))

# Problem: Find a position in a genome where the skew diagram attains a minimum
# Input: A DNA string 'Genome'
# Output: All integer(s) i minimizing Skew_i(Genome) among all values of 'i'
#         (from 0 to |Genome|)
def MinSkew(Genome):
    skew_array = Skew(Genome)
    minimum = min(skew_array.values())
    minimum_indices = []
    for position, count in skew_array.iteritems():
        if count == minimum:
            minimum_indices.append(position)
    return minimum_indices

# Problem: Compute the Hamming distance between two strings
# Input: Two strings of equal length
# Output: The Hamming distance between these strings
def HammingDistance(p, q):
    distance = 0
    for idx, nucleotide in enumerate(p):
        if nucleotide != q[idx]:
            distance += 1
    return distance

# Problem: Find all approximate occurrences of a pattern in a string
# Input: String 'Pattern' and 'Text' along with an integer 'd'
# Output: All string positions where 'Pattern' appears as a substring of 'Text'
#         with at most 'd' mismatches
def ApproximatePatternMatching(Pattern, Text, d):
    positions = []
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Pattern, Text[i:i+len(Pattern)]) <= d:
            positions.append(i)
    return positions

# ApproximatePatternCount
# Problem: Find the repetition number of Pattern in a Text with at most d
#          mismatches
# Input: 'Pattern', 'Text', and 'd'
# Output: Number of patterns
def CountWithMismatches(Pattern, Text, d):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Pattern, Text[i:i+len(Pattern)]) <= d:
            count += 1
    return count

# Problem: Find the most frequent k-mers with mismatches in a string
# Input: A string 'Text' as well as integers 'k' abd 'd'
#        Assume that k .lte. 12, and d .lte. 3
# Output: All most frequent k-mers with up to 'd' mismatches in 'Text'
def ComputingFrequenciesWithMismatches(Text, k, d):
    FrequencyArray = {}
    for i in range(4**k):
        FrequencyArray[i] = 0
    # Alternatively: [0] * (4**k) = [0, ..., 0] where len([given_list]) = 4**k
    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        Neighborhood = Neighbor(Pattern, d)
        for ApproximatePattern in Neighborhood:
            FrequencyArray[j] += 1
    return FrequencyArray

def FasterFrequentWordsWithMismatches():
    pass

def ImmediateNeighbors(Pattern):
    Neighborhood = [Pattern]
    for i in range(0,len(Pattern)):
        symbol = Pattern[i]
        listOfNucleotides = ['A', 'C', 'G', 'T']
        listOfNucleotides.remove(symbol)
        for nucleotide in listOfNucleotides:
            PatternPrime = list(Pattern)
            PatternPrime[i] = nucleotide
            Neighbor = ''.join(PatternPrime)
            Neighborhood.append(Neighbor)
    return Neighborhood

# ========
# INFESTED
# ========
# THERE'S AN ERROR HERE, IT RETURNS MORE DIIFF STRINGS THAN THE DEFINED DISTANCE
def IterativeNeighbors(Pattern, d):
    Neighborhood = [Pattern]
    for j in range(0,d):
        for NeighPattern in Neighborhood:
        #NeighPattern = Neighborhood[0]
            Neighborhood += ImmediateNeighbors(NeighPattern)
            #[Neighborhood.append(num) for num in ImmediateNeighbors(NeighPattern)]
            Neighborhood = RemoveDuplicates(Neighborhood)
    return Neighborhood
# ========
# INFESTED
# ========

def Neighbors(Pattern, d):
    if d == 0: # ako je hamingova razdaljina nula, samo mi vrati obrazac
        return [Pattern]
    if len(Pattern) == 1:
        # ako je obrazac samo jedan nukleotid, vrati mi sve posto su to sve
        # varijacije na dati obrazac
        return ['A', 'C', 'G', 'T']
    # u normalnom radu, skrati mi pattern (suffix) i onda za svakog suseda iz
    # sufiksa, ako je hamingova razdaljina manja od 'd' izmedju sufiksa ulaznog
    # obrasca i teksta iz skupa sufiksa, onda dodaj kao prefiks na taj teks sve
    # nukleotide. Ako nije to slucaj, onda samo dodaj natrag prvo slovo obrasca
    # na tekst. I naravno, vrati mi rezultat
    Neighborhood = []
    SuffixNeighbors = Neighbors(Suffix(Pattern), d)
    for Text in SuffixNeighbors:
        if HammingDistance(Suffix(Pattern), Text) < d:
            for nucleotide in 'ACGT':
                Neighborhood.append(nucleotide + Text)
        else:
            Neighborhood.append(FirstSymbol(Pattern) + Text)
    return Neighborhood

def FirstSymbol(Pattern):
    return Pattern[0]

def Suffix(Pattern):
    return Pattern[1:]
