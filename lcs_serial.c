#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#define ALPHABET 26

int* calculate_P(char *b, int b_size){
    int i, j;
    const int span = b_size + 1;
    int *P = (int *) calloc(ALPHABET * (b_size + 1), sizeof(int));

    for (i = 0; i < ALPHABET; i++) {
        int block = i * span;
        char letter1 = 'A' + i;
        for (j = 1; j < b_size + 1; j++) {
            if (b[j - 1] == letter1) {
                P[block + j] = j;
            } else {
                P[block + j] = P[block + j - 1];
            }
        }
    }
    return P;
}

int* calculate_S(char *a, int a_size, char *b, int b_size, int* P){
    const int span = b_size + 1;
    int i, j;
    int *S = (int*) calloc((a_size + 1) * (b_size + 1), sizeof(int));
    if (S == NULL) {
        printf("Error! Not enough memory!\n");
        exit(-1);
    }
    int t, s;
    for (i = 1; i < a_size + 1; i++) {
        int block = i * span;
        for (j = 1; j < b_size + 1; j++) {
            int c = a[i - 1] - 'A';
            int p = P[c * span + j];

            if((0 - p) < 0) t = 1;
            else t = 0;

            if((0 - (S[block - span + j] - t * S[block - span + p - 1])) < 0) s = 1;
            else s = 0;

            S[block + j] = S[block - span + j] + t * (s ^ 1);
        }
    }
    return S;
}

char* calculate_LCS(char *a, int a_size, char *b, int b_size, int* S){
    int i = a_size;
    int j = b_size;
    const int span = b_size + 1;
    int block = i * span;
    int last_cell = S[block + j];
    char *LCS = (char*) calloc(last_cell + 1, sizeof(char));
    LCS[last_cell] = '\0';

    while(last_cell > 0) {
        block = i * span;
        if (a[i - 1] == b[j - 1]) {
            LCS[last_cell - 1] = a[i - 1];
            i--; j--; last_cell--;
        }
        else if (S[block + j - 1] < S[block + j - span]) i--;
        else j--;
    }
    return LCS;
}

int main(int argc, char const *argv[]) {

    FILE *file;
    char *a, *b;
    int a_size = 0, b_size = 0;
    double start_time, stop_time;

    file = fopen(argv[1], "r");
    if(file == NULL){
        printf("%s\n", "File cannot be opened!");
        exit(-1);
    }
    else if(file != NULL) {
        fscanf(file, "%d %d", &a_size, &b_size);
        a = (char *) calloc(a_size + 1, sizeof(char));
        b = (char *) calloc(b_size + 1, sizeof(char));
        fscanf(file, "%s %s", a, b);
    }
    fclose(file);

    start_time=MPI_Wtime();

    int i, j;
    int *P = calculate_P(b, b_size);
    int *S = calculate_S(a, a_size, b, b_size, P);

    const int span = b_size + 1;
    int length = S[a_size * span + b_size];
    char *LCS = calculate_LCS(a,a_size,b,b_size,S);

    stop_time=MPI_Wtime();

    printf("Length of LCS: %d\n", length);
    printf("LCS:");
    for(i = 0; i < length; i++) printf("%c", LCS[i]);
    printf("\nExecution time: %lf\n", stop_time - start_time);

    free(S);
    free(P);
    free(a);
    free(b);
    free(LCS);

    return 0;
}