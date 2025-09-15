// chkenc.c - check that encoding matrix is invertible

#include <string.h>
#include <stdint.h>
#include <stdio.h>

#define MAX_CHECK 63 /* Size is limited by using uint64_t to represent subsets */
#define M_MAX     0x20
#define K_MAX     0x10
#define ROWS      M_MAX
#define COLS      K_MAX

#define SYSTEMATIC 0

// only set one of these
#define CAUCHY 0
#define RS 0
#define VND 1

// if 0 generate standard Vandermonde matrix
#define SYSTEMATIC 0

static unsigned char abExp[512];        /* antilog table */
static unsigned char abLog[256];        /* log table */
static unsigned char abInv[256];        /* inverse table */

static unsigned char gf_mul(unsigned char b0, unsigned char b1)
{
        if(b0 == 0 || b1 == 0)
            return 0;
        return(abExp[(int)abLog[b0]+(int)abLog[b1]]);
}

static unsigned char gf_pow(unsigned char b0, unsigned int i1)
{
    if(i1 == 0)
        return 1;
    if(b0 == 0)
        return 0;
    return(abExp[(abLog[b0]*i1)%(255u)]);
}

static unsigned char gf_inv(unsigned char b0)
{
        if(b0 == 0){
            printf("divide by zero\n");
            return(0);
        }
        return abInv[b0];
}

static unsigned char gf_div(unsigned char b0, unsigned char b1)
{
        return gf_mul(b0, gf_inv(b1));
}

// low level gf multiply
static unsigned char gf_mul0(unsigned char b0, unsigned char b1)
{
int i;
int product;
        product = 0;
        for(i = 0; i < 8; i++){
            product <<= 1;
            if(product & 0x100u){
                product ^= 0x11du;}
            if(b0 & 0x80u){
                product ^= b1;}
            b0 <<= 1;}
        return((unsigned char)product);
}

static void initgf(void)
{
int i;
unsigned char t;

        t = 1;                              /* init abExp[] */
        for(i = 0; i < 512; i++){
            abExp[i] = t;
            t = gf_mul0(t, 2);
        }

        abLog[0] = 0xffu;                   /* init abLog[] */
        for(i = 0; i < 255; i++)
            abLog[abExp[i]] = i;

        abInv[0] = 0xffu;                   /* init abInv[] */
        for(i = 1; i < 256; i++)
            abInv[i] = gf_div(1,i);
}

int
gf_invert_matrix(unsigned char *in_mat, unsigned char *out_mat, const int n)
{
        int i, j, k;
        unsigned char temp;

        // Set out_mat[] to the identity matrix
        for (i = 0; i < n * n; i++) // memset(out_mat, 0, n*n)
                out_mat[i] = 0;

        for (i = 0; i < n; i++)
                out_mat[i * n + i] = 1;

        // Inverse
        for (i = 0; i < n; i++) {
                // Check for 0 in pivot element
                if (in_mat[i * n + i] == 0) {
                        // Find a row with non-zero in current column and swap
                        for (j = i + 1; j < n; j++)
                                if (in_mat[j * n + i])
                                        break;

                        if (j == n) // Couldn't find means it's singular
                                return -1;

                        for (k = 0; k < n; k++) { // Swap rows i,j
                                temp = in_mat[i * n + k];
                                in_mat[i * n + k] = in_mat[j * n + k];
                                in_mat[j * n + k] = temp;

                                temp = out_mat[i * n + k];
                                out_mat[i * n + k] = out_mat[j * n + k];
                                out_mat[j * n + k] = temp;
                        }
                }

                temp = gf_inv(in_mat[i * n + i]); // 1/pivot
                for (j = 0; j < n; j++) {         // Scale row i by 1/pivot
                        in_mat[i * n + j] = gf_mul(in_mat[i * n + j], temp);
                        out_mat[i * n + j] = gf_mul(out_mat[i * n + j], temp);
                }

                for (j = 0; j < n; j++) {
                        if (j == i)
                                continue;

                        temp = in_mat[j * n + i];
                        for (k = 0; k < n; k++) {
                                out_mat[j * n + k] ^= gf_mul(temp, out_mat[i * n + k]);
                                in_mat[j * n + k] ^= gf_mul(temp, in_mat[i * n + k]);
                        }
                }
        }
        return 0;
}

void
gf_gen_rs_matrix(unsigned char *a, int m, int k)
{
        int i, j;
        unsigned char p, gen = 1;

        memset(a, 0, k * m);
        for (i = 0; i < k; i++)
                a[k * i + i] = 1;

        for (i = k; i < m; i++) {
                p = 1;
                for (j = 0; j < k; j++) {
                        a[k * i + j] = p;
                        p = gf_mul(p, gen);
                }
                gen = gf_mul(gen, 2);
        }
}

void
gf_gen_cauchy1_matrix(unsigned char *a, int m, int k)
{
        int i, j;
        unsigned char *p;

        // Identity matrix in high position
        memset(a, 0, k * m);
        for (i = 0; i < k; i++)
                a[k * i + i] = 1;

        // For the rest choose 1/(i + j) | i != j
        p = &a[k * k];
        for (i = k; i < m; i++)
                for (j = 0; j < k; j++)
                        *p++ = gf_inv(i ^ j);
}

void
gf_gen_vnd_matrix(unsigned char *a, int m, int k)
{
        int i, j;
        unsigned char g;
#if SYSTEMATIC
        int n;
        unsigned char d;
#endif        
        memset(a, 0, k * m);
        // generate simple Vandermonde matrix
        g = 0;
        for (i = 0; i < m; i++) {
                for (j = 0; j < k; j++) {
                        a[k*i+j] = gf_pow(g, j);
                }
                g = g+1;
        }
#if SYSTEMATIC
        // modified gaussian jordan inversion of k by k sub-matrix of a
        // and update row k to m-1
        for (i = 0; i < k; i++) {
                p = a[k*i+i];                   /* p = pivot */
                d = gf_inv(p);                  /* d = 1/p */
                a[k*i+i] = 1;                   /* pivot = 1 */
                for(n = 0; n < k; n++)          /* divide row by p */
                        a[k*i+n] = gf_mul(a[k*i+n], d);
                for(j = 0; j < m; j++){         /* update other rows */
                        if(j == i)
                            continue;
                        p = a[k*j+i];           /* p = pivot */
                        a[k*j+i] = 0;           /* pivot = 0 */
                        for(n = 0; n < k; n++)
                                a[k*j+n] ^= gf_mul(a[k*i+n], p);
                }
        }
        // set a[...k] = identity matrix
        memset(a, 0, k*k);
        for(i = 0; i < k; i++)
                a[k*i+i] = 1;
#endif
}

static inline uint64_t
min(const uint64_t a, const uint64_t b)
{
        if (a <= b)
                return a;
        else
                return b;
}

void
gen_sub_matrix(unsigned char *out_matrix, const uint64_t dim, unsigned char *in_matrix,
               const uint64_t rows, const uint64_t cols, const uint64_t row_indicator,
               const uint64_t col_indicator)
{
        uint64_t i, j, r, s;

        for (i = 0, r = 0; i < rows; i++) {
                if (!(row_indicator & ((uint64_t) 1 << i)))
                        continue;

                for (j = 0, s = 0; j < cols; j++) {
                        if (!(col_indicator & ((uint64_t) 1 << j)))
                                continue;
                        out_matrix[dim * r + s] = in_matrix[cols * i + j];
                        s++;
                }
                r++;
        }
}

/* Gosper's Hack */
uint64_t
next_subset(uint64_t *subset, uint64_t element_count, uint64_t subsize)
{
        uint64_t tmp1 = *subset & (0-*subset);
        uint64_t tmp2 = *subset + tmp1;
        *subset = (((*subset ^ tmp2) >> 2) / tmp1) | tmp2;
        if (*subset & (((uint64_t) 1 << element_count))) {
                /* Overflow on last subset */
                *subset = ((uint64_t) 1 << subsize) - 1;
                return 1;
        }

        return 0;
}

int
are_submatrices_singular(unsigned char *vmatrix, const uint64_t rows, const uint64_t cols)
{
        unsigned char matrix[COLS * COLS];
        unsigned char invert_matrix[COLS * COLS];
        uint64_t subsize;

        /* Check all square subsize x subsize submatrices of the rows x cols
         * vmatrix for singularity*/
        for (subsize = 1; subsize <= min(rows, cols); subsize++) {
                const uint64_t subset_init = (1ULL << subsize) - 1ULL;
                uint64_t col_indicator = subset_init;
                do {
                        uint64_t row_indicator = subset_init;
                        do {
                                gen_sub_matrix(matrix, subsize, vmatrix, rows, cols, row_indicator,
                                               col_indicator);
                                if (gf_invert_matrix(matrix, invert_matrix, (int) subsize))
                                        return 1;

                        } while (next_subset(&row_indicator, rows, subsize) == 0);
                } while (next_subset(&col_indicator, cols, subsize) == 0);
        }

        return 0;
}

int
main(int argc, char **argv)
{
        unsigned char vmatrix[(ROWS + COLS) * COLS];
        uint64_t rows, cols;

        if (K_MAX > MAX_CHECK) {
                printf("K_MAX too large for this test\n");
                return 0;
        }
        if (M_MAX > MAX_CHECK) {
                printf("M_MAX too large for this test\n");
                return 0;
        }
#if 0
        if (M_MAX > K_MAX) {
                printf("M_MAX must be <= K_MAX");
                return 0;
        }
#endif
        initgf();
                                                  
#if CAUCHY
        printf("Checking gf_gen_caucy1_matrix for k <= %d and m <= %d.\n", K_MAX, M_MAX);
#elif RS
        printf("Checking gf_gen_rs_matrix for k <= %d and m <= %d.\n", K_MAX, M_MAX);
#else
        printf("Checking gf_gen_vnd_matrix for k <= %d and m <= %d.\n", K_MAX, M_MAX);
#endif
        printf("gen_rs_matrix creates erasure codes for:\n");

        for (cols = 1; cols <= K_MAX; cols++) {
                for (rows = 1; rows <= M_MAX - cols; rows++) {
#if CAUCHY
                        gf_gen_cauchy1_matrix(vmatrix, (int)(rows + cols), (int)cols);
#elif RS
                        gf_gen_rs_matrix(vmatrix, (int)(rows + cols), (int)cols);
#else
                        gf_gen_vnd_matrix(vmatrix, (int)(rows + cols), (int)cols);
#endif
                        /* Verify the Vandermonde portion of vmatrix contains no
                         * singular submatrix */
                        if (are_submatrices_singular(&vmatrix[cols * cols], rows, cols))
                                break;
                }
                printf("   k = %2u, m <= %2u\n", (unsigned) cols, (unsigned) (rows + cols - 1));
        }
        return 0;
}
