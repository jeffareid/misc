void
gf_gen_bch_matrix(unsigned char *a, int m, int k)
{
        unsigned char g[256];
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
        // gaussian reduction (column swap not needed)
        for (i = 0; i < k; i++) {               /* for all columns */
                p = a[k*i+i];                   /* p = pivot */
                d = gf_inv(p);                  /* d = 1/p */
                for(j = 0; j < m; j++)          /* divide column by p */
                        a[k*j+i] = gf_mul(a[k*j+i], d);
                for(n = 0; n < k; n++){         /* update other columns */
                        if(n == i)
                            continue;
                        p = a[k*i+n];
                        for(j = 0; j < m; j++)
                            a[k*j+n] ^= gf_mul(p, a[k*j+i]);
                }
        }
}
