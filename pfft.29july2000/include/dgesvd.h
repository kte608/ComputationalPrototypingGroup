/* ==== CLAPACK routines ==== */
int dgesvd_(char *jobu, char *jobvt, long int  *m, long int  *n, 
	    double *a, long int  *lda, double *s, 
	    double *u, long int  *ldu, double *vt, long int  *ldvt, 
	    double *work, long int *lwork, long int  *info);
