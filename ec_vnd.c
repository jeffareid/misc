void
gf_gen_vnd_matrix(unsigned char *a, int m, int k)
{
        int i, j, n;
        unsigned char p, g;
        unsigned char d;
        memset(a, 0, k * m);
        // generate Vandermonde matrix
        a[0] = 1;
        g = 1;
        for (i = 1; i < m; i++) {
                p = 1;
                for (j = 0; j < k; j++) {
                        a[k*i+j] = p;
                        p = gf_mul(p, g);
                }
                g = gf_mul(g, 2);
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
