#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define ALPHABET 26

int Block_lower_bound(int taskid, int numtasks, int len) {
    return taskid*len/numtasks;
}
int Block_upper_bound(int taskid, int numtasks, int len) {
    return Block_lower_bound(taskid + 1, numtasks, len) - 1;
}
int Block_length(int taskid, int numtasks, int len) {
    return Block_upper_bound(taskid, numtasks, len) - Block_lower_bound(taskid, numtasks, len) + 1;
}
int Block_of_task(int low,int numtasks,int high) {
    return (numtasks * (low + 1) - 1) / high;
}
int Width_of_block(int len) {
    if (len >= 100) return len / 50;
    else return len / 4;
}

int Num_task_blocks(int len, int bw) {
    float x = len /(bw * 1.0);
    int x_low = x;
    if (x - x_low > 0.5) x_low++;
    return x_low;
}

int Length_last_block(int len, int bw) {
    float x = len /(bw * 1.0);
    int x_low = x;
    if (x - x_low > 0.5) return len - (bw * x_low);
    else return bw + (len - (bw * x_low));
}

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

int* calculate_S(char *a, int a_size, char *b, int b_size, int* P, int first, int taskid, int numtasks, int size, MPI_Status status){

    int block_lower_bound = Block_lower_bound(taskid, numtasks, a_size + 1);
    int block_upper_bound = Block_upper_bound(taskid, numtasks, a_size + 1);
    int block_length = Block_length(taskid, numtasks, a_size + 1);
    int width_of_block = Width_of_block(b_size + 1);
    int num_task_blocks = Num_task_blocks(b_size + 1, width_of_block);
    int length_last_block = Length_last_block(b_size + 1, width_of_block);

    int *S = (int*) calloc(size * (b_size + 1), sizeof(int));
    if (S == NULL) {
        printf("Error! Not enough memory!\n");
        exit(-1);
    }

    int k,j;

    for(k = 0; k < num_task_blocks; k++) {
        int col = width_of_block;
        int cur_col = k * col;

        if (k + 1 == num_task_blocks) col = length_last_block;
        if (taskid != first && block_length > 0) MPI_Recv(S + cur_col, col, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

        int lb,i;

        if(block_lower_bound == 0) lb = 1;
        else lb = block_lower_bound;

        const int span = b_size + 1;
        int t, s;

        for (i = 1; lb < block_upper_bound + 1 && i < a_size + 1; i++, lb++) {
            int block = i * span;
            for (j = 0; j < col; j++) {
                if (k != 0 || j != 0) {
                    int col1 = cur_col + j;
                    int c = a[lb - 1] - 'A';
                    int p = P[c * span + col1];

                    if((0 - p) < 0) t = 1;
                    else t = 0;

                    if (p == 0){
                        if((0 - (S[block - span + col1] - t * S[block - span])) < 0) s = 1;
                        else s = 0;
                    }else{
                        if((0 - (S[block - span + col1] - t * S[block - span + p - 1])) < 0) s = 1;
                        else s = 0;
                    }

                    S[block + col1] = S[block - span + col1] + t * (s ^ 1);
                }
            }
        }
        if (block_length > 0 && block_upper_bound < a_size) {
            int n = Block_of_task(block_upper_bound + 1, numtasks, a_size + 1);
            if (taskid == first) MPI_Send(&S[block_upper_bound * span + cur_col], col, MPI_INT, n, 0, MPI_COMM_WORLD);
            else MPI_Send(&S[block_length * span + cur_col], col, MPI_INT, n, 0, MPI_COMM_WORLD);
        }
    }
    return S;
}
double start_time = 0;

void calculate_LCS(int size, int b_size, int a_size, int taskid, int numtasks, int last, MPI_Status status, int first){
    int i = size - 1;
    int j = b_size;
    char xi, yj;
    int last_cell;

    int *lcs;
    const int span = b_size + 1;
    int block = i * span;

    int length;
    if (taskid != last) {

        MPI_Recv(&length, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        lcs = (int*) calloc(length + 1, sizeof(int));
        MPI_Recv(lcs, length + 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);

        j = lcs[0];
        last_cell = S[block + j];

    } else {
        last_cell = S[block + j];
        length = last_cell;
        lcs = (int*) calloc(length + 1, sizeof(int));
    }

    int block_lower_bound = Block_lower_bound(taskid, numtasks, a_size + 1);

    while(i > 0 && last_cell > 0) {
        int block = i * span;
        if (taskid == first) xi = a[(block_lower_bound + i) - 1];
        else xi = a[(block_lower_bound + i) - 2];

        yj = b[j - 1];
        if (xi == yj) {
            lcs[last_cell] = xi;
            last_cell--;
            i--;
            j--;
        } else if (S[block + j - span] > S[block + j - 1]) {
            i--;
        } else {
            j--;
        }
    }
    double stop_time;
    if (taskid != first) {
        lcs[0] = j;
        int pr = Block_of_task(block_lower_bound - 1, numtasks, a_size + 1);
        MPI_Send(&length, 1, MPI_INT, pr , 0, MPI_COMM_WORLD);
        MPI_Send(lcs, length + 1, MPI_INT, pr , 0, MPI_COMM_WORLD);

    } else {
        stop_time = MPI_Wtime();
        printf("Execution time:: %lf\n",stop_time-start_time);
        printf("%d\n", length);
        for (i = 1; i < length + 1; i++) {
            printf("%c", lcs[i]);
        }
        printf("\n");

    }

    free(lcs);
}

int main(int argc, char const *argv[]) {
    int taskid, numtasks;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Status status;

    FILE *file;
    char *a, *b;
    int a_size = 0, b_size = 0;

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

    start_time = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);
    int *P = calculate_P(b, b_size);
    MPI_Barrier(MPI_COMM_WORLD);

    int block_lower_bound = Block_lower_bound(taskid, numtasks, a_size + 1);
    int block_upper_bound = Block_upper_bound(taskid, numtasks, a_size + 1);
    int block_length = Block_length(taskid, numtasks, a_size + 1);
    int first = Block_of_task(0, numtasks, a_size + 1);
    int last = Block_of_task(a_size, numtasks, a_size + 1);
    int size;

    if (taskid == first) size = block_length;
    else size = block_length + 1;

    //Calculating S table
    int *S = calculate_S(a,a_size,b,b_size,P,first,taskid,numtasks,size,status);

    calculate_LCS(size,b_size,a_size,taskid,numtasks,last,&status, first);

    free(S);
    free(P);
    free(a);
    free(b);


    MPI_Finalize();

    return 0;



}