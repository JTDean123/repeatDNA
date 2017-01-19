# repeatDNA
A program that identifies repetitive DNA sequences

### Objective
This objective of this program is to idenitfy [repetitive DNA sequences](https://en.wikipedia.org/wiki/Repeated_sequence_(DNA)) of specified length that occur in a user supplied DNA sequence.  A repetitive DNA sequence is defined as a sequence that occurs more than once.  

###  Development
I found that searching for repetitive sequences in short (100 character) inputs was easily solved with brute force by searching the entire input sequence for a string of specified length and iterating then querying the input sequence by one character at a time.  The following example and pseudocode illustrates this example:

Repeat length: 2  
Input:  aattcgaa

First candidate repeat:  aa

Iterate over length of input:  
aa == at  false  
aa == tt  false  
aa == tc  false  
aa == cg  false  
aa = ga  false  
aa = aa  TRUE  

This method identifies 'aa' as a repeat sequence, occuring twice.  Next, 'at' is queried in the same way.  This algorithm has asymptotic complexity of O(n^2), where n is the length of the input sequence. As such, this algorithm was adequeate for short input sequences but was found to slow to a screaching halt when larger (10^8) sequence inputs were provided.

This algorithm was substantially improved with hash tables and linked lists.  To accomplish this, the input sequence was disassembled into pieces with length equal to the repeat length being evaluated.  When building the hash table if a sequence was found to already exist a new linked node was added.  Evaluating the number of repeats then becomes as simple as summing the number of nodes connected to each bucket.  This modified algorithm reduces the asymptotic complexity to O(n), and genome scale inputs can now be analyzed.

###  Usage
Example usage of this program follows:

First, repeatDNA.c is compiled.

```
$ gcc -o repeatDNA repeatDNA.c
```

Next, an input sequence 'gDNA.txt' is queried for repeats between 6 and 10 characters in length.

```
$ ./repeatDNA gDNA.txt 6 10
```

Terminal output:

total 6 nt repeat occurences: 222  
total 7 nt repeat occurences: 72  
total 8 nt repeat occurences: 20  
total 9 nt repeat occurences: 4  
total 10 nt repeat occurences: 0  
success

The output is saved to 'results.txt' in the directory.
