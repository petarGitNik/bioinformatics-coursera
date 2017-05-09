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

# --- Program testing --- #
import sys
#lines = sys.stdin.read().splitlines()

#with open('./outputs/' + sys.argv[1], 'w') as f:
    # --- PatternCount --- #
    #f.write(str(PatternCount(lines[0], lines[1])) + '\n')
    # --- FrequentWords --- #
#    f.write(' '.join(FrequentWords(lines[0], int(lines[1]))) + '\n')
