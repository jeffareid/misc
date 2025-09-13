void
gf_gen_vnd_matrix(unsigned char *a, unsigned char *b, int m, int k)
{
        int i, j, n;
        unsigned char p, g;
        unsigned char d;
        memset(a, 0, k * m);
        // generate Vandermonde matrix
        g = 1;
        for (i = 0; i < m; i++) {
                p = 1;
                for (j = 0; j < k; j++) {
                        a[k*i+j] = p;
                        p = gf_mul(p, g);
                }
                g = gf_mul(g, 2);
        }
        // modified gaussian jordan inversion of k by k sub-matrix of a
        for (i = 0; i < k; i++) {
                p = a[k*i+i];                   /* p = pivot */
                d = gf_inv(p);                  /* d = 1/p */
                a[k*i+i] = 1;                   /* pivot = 1 */
                for(n = 0; n < k; n++)          /* divide row by p */
                        a[k*i+n] = gf_mul(a[k*i+n], d);
                for(j = 0; j < k; j++){         /* update other rows */
                        if(j == i)
                            continue;
                        p = a[k*j+i];           /* p = pivot */
                        a[k*j+i] = 0;           /* pivot = 0 */
                        for(n = 0; n < k; n++)
                                a[k*j+n] ^= gf_mul(a[k*i+n], p);
                }
        }
        // a[k...] = b[k...] = mpy by inversed k by k sub-matrix of a
        for(i = k; i < m; i++){
                for(j = 0; j < k; j++){
                        b[k*i+j] = 0;
                        for(n = 0; n < k; n++)
                                b[k*i+j] ^= gf_mul(a[k*j+n], a[k*i+n]);
                }
        }
        memcpy(&a[k*k], &b[k*k], m*(m-k));
        // set a[...k] = identity matrix
        memset(a, 0, k*k);
        for(i = 0; i < k; i++)
                a[k*i+i] = 1;
}
