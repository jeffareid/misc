void
gf_gen_bch_matrix(unsigned char *a, int m, int k)
{
        unsigned char g[64];
        int i, j, n;
        unsigned char p;
        unsigned char d;
        // generate bch polynomial
        memset(g, 0, m-k+1);
        g[0] = 1;
        p = 1;
        for(j = 0; j < m-k; j++){
                for(i = j; i >= 0; i--){
                    g[i+1] ^= gf_mul(g[i], p);}
                p = gf_mul(p, 2);
        }
        // generate bch non-systematic matrix
        memset(a, 0, k * m);
        for(i = 0; i < k; i++){
                n = 0;
                for(j = 0; j < m; j++){
                        if(n != (m-k+1))
                                a[(i+j)*k+i] = g[n++];
                }
        }
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
}
