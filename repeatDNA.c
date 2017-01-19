/*
this program finds repetetive DNA sequences.  the input to
the program is a DNA file, the smallest length repeat to search for,
and the largest length repeat to search for.  the search lengths must be
>2 and <10000.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <stdbool.h>

//these variables will be used to build a hash table of repeat DNA sequences
typedef struct node
{
    char srepeat[10000 + 1];
    struct node* next;
} node;

unsigned long hashfxn(char* str);
int size;

//define the number of buckets to use in the backbone of the hash table
const int buckets = 50000;
const int rbuckets = 50000;

//make the hashtable a global variable
struct node* hashtable[buckets];
struct node* rhashtable[buckets];

//prototype to make nodes for the hash table
struct node* make_node();

//prototype to test if repeat has already been identified
int check(char* repeat);

//prototype to test if repeat has already been identified
bool check_rhash(char* repeat);



//let the good times roll.....
int main(int argc, char* argv[])
{
    //ensure that the user has submitted correct input
    if (argc != 4)
    {
        printf("argument error \ninput format: ./gRepeat gDNA_sequence.txt lower upper\n");
        return 0;
    }

    int small = atoi(argv[2]);
    int large = atoi(argv[3]);

    if (large > 10000)
    {
        printf("argument error: large repeat length must be <= 2000\n");
        return 0;
    }

    if (small < 2)
    {
        printf("argument error: small repeat length must be >1\n");
        return 0;
    }

    if (small >= large)
    {
        printf("argument error:  small repeat must be < large repeat\n");
        return 0;
    }

    //open the gDNA sequence file that the user has specified
    char* gDNAtxt = argv[1];
    FILE* gDNAin = fopen(gDNAtxt, "r");

    if (gDNAin == NULL)
    {
        printf("could not open the input file\n");
    }

    //determine file size and assign a pointer that will point to a string containing the gDNA sequence
    fseek(gDNAin, 0, SEEK_END);
    size = ftell(gDNAin);
    fseek(gDNAin, 0 ,SEEK_SET);
    char* temp_DNA = malloc(sizeof(char)*size+1);
    char* gDNA = calloc(size+1, sizeof(char));

    if (temp_DNA == NULL)
   {
        printf("memory allocation error\n");
        return 0;
    }

    if (gDNA == NULL)
    {
        printf("memory allocation error\n");
        return 0;
    }

    //convert DNA sequence to lowercase and ensure that the input only contains a, t, c, and g
    int dna_count = 0;

    for (int i = 0; i < size; i++)
    {
        *(temp_DNA+i) = tolower(fgetc(gDNAin));

        //ensure that the position is a, t, c, or g
        if (*(temp_DNA+i) == 'a' || *(temp_DNA+i) == 't' || *(temp_DNA+i) == 'c' || *(temp_DNA+i) == 'g')
        {
            *(gDNA+dna_count) = *(temp_DNA+i);
            dna_count++;
        }
    }

    free(temp_DNA);

    //now that all of the arguments and input have been handled, initialize hash tables
    for (int i = 0; i < buckets; i++)
    {
        hashtable[i] = make_node();
        rhashtable[i] = make_node();
    }

    //open files to output the data to
    FILE* output = NULL;
    output = fopen("results.txt", "w");

    if (output == NULL)
    {
        printf("error opening output file\n");
    }

    //repeat sequences will go into a hash table
    unsigned long hash;
    unsigned long rhash;

    //find repetetive sequences of different lengths
    for (int x = small; x <= large; x++)
    {
        //allocate memory for the nucleotide (potential) repeat
        char* repeat = calloc(x+1, sizeof(char));

        if (repeat == NULL)
        {
            printf("memory allocation error\n");
            return 0;
        }

        //build the hash table
        for (int i = 0; i < dna_count - x - 2; i++)
        {
            strncpy(repeat, gDNA+i, x);

            //put the repeat into a hash table
            hash = hashfxn(repeat);

            //make a temp node that will be inserted into the list
            node* insert = make_node();

            //insert into the hash table
            strcpy(insert->srepeat, repeat);
            insert->next = hashtable[hash];
            hashtable[hash] = insert;
        }

        //count the total number of repeats for a given size
        int repeat_count = 0;

        //see how many times each item occurs in the hash table.  if > 1, add it to a new hash table consisting of repeats only
        for (int i = 0; i < dna_count - x - 2; i++)
        {
            strncpy(repeat, gDNA+i, x);

            int counter = 0;
            counter = check(repeat);
            if (counter > 1)
            {
                repeat_count = repeat_count + counter;

                //put the repeat into a hash table
                rhash = hashfxn(repeat);

                //make a temp node that will be inserted into the list
                node* rinsert = make_node();

                //insert into the hash table
                strcpy(rinsert->srepeat, repeat);
                rinsert->next = rhashtable[rhash];
                rhashtable[rhash] = rinsert;

                //check to see if repeat already exists.  if this is the first instance, write it to the output
                if (check_rhash(repeat) == false)
                {
                    //format data for output and write it to output
                    char* outdata = malloc((strlen(repeat)+sizeof(counter)+40)*sizeof(char));
                    sprintf(outdata, "sequence: %s occurences: %i\n", repeat, counter);
                    fputs(outdata, output);
                    free(outdata);
                }
            }
        }

        fprintf(output, "\ntotal %i nt repeat occurences: %i\n\n", x, repeat_count);
        printf("total %i nt repeat occurences: %i\n", x, repeat_count);
        //printf(".");
        free(repeat);
    }

    //close all of the open files
    fclose(gDNAin);
    fclose(output);

    //free pointers
    free(gDNA);

    //let the user know that the program ran successfully
    printf("success\n");

    return 0;
}


//the hash function used here (djb2) is from the following website:  http://www.cse.yorku.ca/~oz/hash.html
unsigned long hashfxn(char* str)
{
    unsigned long hash = 5381;

    while (*str != '\0')
    {
        hash = ((hash << 5) + hash) + *str;
        str++;
    }
    hash = hash % buckets;

    return hash;
}

//a function to initialize a new node
struct node* make_node()
{
    struct node* temp = malloc(sizeof(node));
    temp->next = NULL;
    return temp;
}

//check to see if the repeat exists in the hash table
int check(char* repeat)
{
    node* checker = NULL;

    //get the hashfxn index
    unsigned long hash = hashfxn(repeat);

    //check if repeat is in the hash table
    int counter = 0;
    for (checker = hashtable[hash]; checker->next != NULL; checker = checker->next)
    {
        if (strcmp(checker->srepeat, repeat) == 0)
        {
            counter++;
        }
    }

    return counter;
}

//check to see if the repeat exists in the hash table
bool check_rhash(char* repeat)
{
    node* checker = NULL;

    //get the hashfxn index
    unsigned int hash = hashfxn(repeat);

    //check if repeat is already in the hash table
    int counter = 0;
    for (checker = rhashtable[hash]; checker->next != NULL; checker = checker->next)
    {
        if (strcmp(checker->srepeat, repeat) == 0)
        {
            counter++;
            if (counter == 2)
            {
                return true;
            }
        }
    }

    return false;
}
