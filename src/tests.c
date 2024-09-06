#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>


/* ***********************************************************************
 *                                                                       *
 *                          Declaration of functions                     *
 *                                                                       *
 * ********************************************************************* */



/* Functions coming from the package ade4 */

void vecpermut (double *A, int *num, double *B);
double alea (void);
void aleapermutvec (double *a);
void trirapideintswap (int *v, int i, int j);
void trirapideint (int *x , int *num, int gauche, int droite);
void sqrvec (double *v1);
void getpermutation (int *numero, int repet);
void prodmatABC (double **a, double **b, double **c);
void prodmatAtAB (double **a, double **b);
void prodmatAtBC (double **a, double **b, double **c);
void prodmatAAtB (double **a, double **b);
void prodmatAtBrandomC (double **a, double **b, double **c, int*permut);
void taballoc (double ***tab, int l1, int c1);
void vecalloc (double **vec, int n);
void vecintalloc (int **vec, int n);
void freetab (double **tab);
void freevec (double *vec);
void freeintvec (int *vec);
void matcentrage (double **A, double *poili, char *typ);
void matmodiffc (double **tab, double *poili);
void matmodifcp (double **tab, double *poili);
void matmodifcs (double **tab, double *poili);
void matmodifcn (double **tab, double *poili);
void matmodifcm (double **tab, double *poili);
void DiagobgComp (int n0, double **w, double *d, int *rang);





/* Functions from the package adehabitat */
void rpath(double **xp, double *rcx, double *rcy, double **asc, 
	   double **tabdist, double *dt, 
	   double *angles, double *xc, double *yc,
	   double *cs, int r);
void randpath(double *xpr, double *rcrx, double *rcry, double *ascr, 
	      double *xcr, double *ycr, double *csr,
	      double *tabdistr, double *dtr, double *anglesr, 
	      int *nlasc, int *ncasc, int *nltdr, int *nlocsr);
void dtmp(double x1, double x2, double y1, double y2, 
	  double *di);
void fptt(double *x, double *y, double *t, int pos, double radius, 
	  double *fptto, int nlo);
void fipati(double *x, double *y, double *t, 
	    int nlo, int nra, double *rad, 
	    double **fpt);
void fipatir(double *xr, double *yr, double *tr, 
	     int *nlocs, double *radius, int *nrad, 
             double *fptr);
void perclu(double **map, int nr, int nc, double *x, double *y,
	    int nmax, int *nreel, double *pm);
void perclur(double *mapr, int *nrm, int *ncm, double *probamr,
	     double *xr, double *yr, int *nmaxr, int *nreel);
void resolpol(double a, double b, double c, double *x1, double *x2, 
	      int *warn);
void discretraj(double *x, double *y, double *dat, double *xn, 
		double *yn, int n, int nn, double *datn, 
		double u, int *neff);
void discretrajr(double *xr, double *yr, double *datr, double *xnr, 
		 double *ynr, int *nr, int *nnr, double *datnr, 
		 double *xdeb, double *ydeb, double *ur, double *dat0, 
		 int *neff);
void permutR2n(double *xyr, int *nro, int *nrepr, 
	       double *R2nr, double *dtr, double *dtsimr);
void runsltr(int *xr, int *nr, double *res, int *nrepr);
void testindepangl (double *sim, double *ang, int *nang, int *debut, 
                    int *fin, int *ndeb, int *ni);
void testindepdist (double *sim, double *di, int *ndi, int *debut, 
		    int *fin, int *ndeb, int *ni);
void prepquart (double *dtur, int *nur, double *dtrr, double *resr, int *nrr, 
		int *ncr, double *rur);
void optcut (double **Pid, double **mk, int *maxk);
void optcutr (double *Pidr, double *mkr, int *maxk, int *lr, int *pr, 
	      double *baye);
void partraj(double **Pid, int *maxk, double **Mk, double **Mkd, 
	     double **res);
void partrajr(double *Pidr, double *curmar, int *curmodr, int *curlocr, 
	      int *lr, int *Dr, int *Kmr);
void acfdist (double *sim, double *di, int *ndi, double *diper, int *nbobs, 
	      int *indxnona, int *nsim, int *lag);
void acfangl (double *sim, double *ang, int *nang, double *angper, 
	      int *nbobs, int *indxnona, int *nsim, int *lag);

/*********************************************************************
 *********************************************************************
 *********                                                       *****
 *********               The sources of ADE-4                    *****
 *********               --------------------                    *****
 *********************************************************************
 *********************************************************************
 */



/**************************/
double alea (void)
{
    double w;
    GetRNGstate();
    w = unif_rand();
    PutRNGstate();
    return (w);
}

/*************************/
void aleapermutvec (double *a)
{
    /* Randomly permutes the elements of a vector a
       Manly p. 42 The vector is modified
       from Knuth 1981 p. 139 */
    int lig, i,j, k;
    double z;
    
    lig = a[0];
    for (i=1; i<=lig-1; i++) {
	j=lig-i+1;
	k = (int) (j*alea()+1);
	/* k = (int) (j*genrand()+1); */
	if (k>j) k=j;
	z = a[j];
	a[j]=a[k];
	a[k] = z;
    }
}


/*******************/	
void vecpermut (double *A, int *num, double *B)
{
/*---------------------------------------
 * A is a vector with n elements
 * B is a vector with n elements
 * num is a random permutation of the n first integers
 * B contains in output the permuted elements of A
 * ---------------------------------------*/
    
    int lig, lig1, lig2, i, k;
    
    lig = A[0];
    lig1 = B[0];
    lig2 = num[0];
    
    
    if ( (lig!=lig1) || (lig!=lig2) ) {
	/* err_message ("Illegal parameters (vecpermut)");
	   closelisting(); */
    }
    
    for (i=1; i<=lig; i++) {
	k=num[i];
	B[i] = A[k];
    }
}

/********* Centring accrding to row weights poili **********/	
void matcentrage (double **A, double *poili, char *typ)
{
    
    if (strcmp (typ,"nc") == 0) {
	return;
    } else if (strcmp (typ,"cm") == 0) {
	matmodifcm (A, poili);
	return;
    } else if (strcmp (typ,"cn") == 0) {
	matmodifcn (A, poili);
	return;
    } else if (strcmp (typ,"cp") == 0) {
	matmodifcp (A, poili);
	return;
    } else if (strcmp (typ,"cs") == 0) {
	matmodifcs (A, poili);
	return;
    } else if (strcmp (typ,"fc") == 0) {
	matmodiffc (A, poili);
	return;
    } else if (strcmp (typ,"fl") == 0) {
	matmodifcm (A, poili);
	return;
    }
}

/*********************/
void matmodifcm (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a complete disjonctive table with n rows and m columns
 * poili is a vector with n components
 * The process returns tab centred by column
 * with weighting poili (sum=1)
 * centring type multple correspondances
 --------------------------------------------------*/
{
    double		poid;
    int 			i, j, l1, m1;
    double		*poimoda;
    double		x, z;
    
    l1 = tab[0][0];
    m1 = tab[1][0];
    vecalloc(&poimoda, m1);
    
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=m1;j++) {
	    poimoda[j] = poimoda[j] + tab[i][j] * poid;
	}
    }
    
    for (j=1;j<=m1;j++) {
	x = poimoda[j];
	if (x==0) {
	    for (i=1;i<=l1;i++) tab[i][j] = 0;
	} else {
	    
	    for (i=1;i<=l1;i++) {
		z = tab[i][j]/x - 1.0;
		tab[i][j] = z;
	    }
	}
    }
    freevec (poimoda);
}

/*********************************************************/
void matmodifcn (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table n rows and p columns
 * poili is a vector with n components
 * the function returns tab normed by column
 * with the weighting poili (sum=1)
 --------------------------------------------------*/
{
    double		poid, x, z, y, v2;
    int 			i, j, l1, c1;
    double		*moy, *var;
    
    l1 = tab[0][0];
    c1 = tab[1][0];
    
    vecalloc(&moy, c1);
    vecalloc(&var, c1);
    
    
/*--------------------------------------------------
 * centred and normed table
 --------------------------------------------------*/
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=c1;j++) {
	    moy[j] = moy[j] + tab[i][j] * poid;
	}
    }
    
    for (i=1;i<=l1;i++) {
	poid=poili[i];
	for (j=1;j<=c1;j++) {
	    x = tab[i][j] - moy[j];
	    var[j] = var[j] + poid * x * x;
	}
    }
    
    for (j=1;j<=c1;j++) {
	v2 = var[j];
	if (v2<=0) v2 = 1;
	v2 = sqrt(v2);
	var[j] = v2;
    }
    
    for (i=1;i<=c1;i++) {
	x = moy[i];
	y = var[i];
	for (j=1;j<=l1;j++) {
	    z = tab[j][i] - x;
	    z = z / y;
	    tab[j][i] = z;
	}
    }
    
    freevec(moy);
    freevec(var);
    
}

/*********************************************************/
void matmodifcs (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table n rows, p columns
 * poili is a vector with n components
 * The function returns tab standardised by column
 * for the weighting poili (sum=1)
 --------------------------------------------------*/
{
	double		x,poid, z, y, v2;
	int 			i, j, l1, c1;
	double		*var;
	
	l1 = tab[0][0];
	c1 = tab[1][0];
	vecalloc(&var, c1);
	

/*--------------------------------------------------
 * calculation of the standardised table
 --------------------------------------------------*/
	
	for (i=1;i<=l1;i++) {
	    poid=poili[i];
	    for (j=1;j<=c1;j++) {
		x = tab[i][j];
		var[j] = var[j] + poid * x * x;
	    }
	}
	
	for (j=1;j<=c1;j++) {
	    v2 = var[j];
	    if (v2<=0) v2 = 1;
	    v2 = sqrt(v2);
	    var[j] = v2;
	}
	
	for (i=1;i<=c1;i++) {
	    y = var[i];
	    for (j=1;j<=l1;j++) {
		z = tab[j][i];
		z = z / y;
		tab[j][i] = z;
	    }
	}
	freevec(var);
}


/**********/
void matmodifcp (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table with n rows and p colonnes
 * poili is a vector with n components
 * The function returns tab centred by column
 * for the weighting poili (sum=1)
 --------------------------------------------------*/
{
    double		poid;
    int 			i, j, l1, c1;
    double		*moy, x, z;
    
    l1 = tab[0][0];
    c1 = tab[1][0];
    vecalloc(&moy, c1);
    
    
/*--------------------------------------------------
 * Centred table
 --------------------------------------------------*/
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=c1;j++) {
	    moy[j] = moy[j] + tab[i][j] * poid;
	}
    }
    
    
    for (i=1;i<=c1;i++) {
	x = moy[i];
	for (j=1;j<=l1;j++) {
	    z = tab[j][i] - x;
	    tab[j][i] = z;
	}
    }
    freevec(moy);
}

/*********************/
void matmodiffc (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table with n rows and m columns
 * of number >=0
 * poili is a vector with n components
 * The function returns tab doubly centred
 * for the weighting poili (sum=1)
 * centring type simple correspondance analysis
 --------------------------------------------------*/
{
    double		poid;
    int 			i, j, l1, m1;
    double		*poimoda;
    double		x, z;
    
    l1 = tab[0][0];
    m1 = tab[1][0];
    vecalloc(&poimoda, m1);
    
    
    for (i=1;i<=l1;i++) {
	x = 0;
	for (j=1;j<=m1;j++) {
	    x = x + tab[i][j];
	}
	if (x!=0) {
	    for (j=1;j<=m1;j++) {
		tab[i][j] = tab[i][j]/x;
	    }
	}	
    }
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=m1;j++) {
	    poimoda[j] = poimoda[j] + tab[i][j] * poid;
	}
    }
    
    for (j=1;j<=m1;j++) {
	x = poimoda[j];
	if (x==0) {
	    /* err_message("column has a nul weight (matmodiffc)"); */
	}
	
	for (i=1;i<=l1;i++) {
	    z = tab[i][j]/x - 1.0;
	    tab[i][j] = z;
	}
    }
    freevec (poimoda);
}









/*****************/
void getpermutation (int *numero, int repet)
/*----------------------
 * affects a random permutation of the first n integers
 * in an integer vector of length n
 * First vecintalloc is needed
 * *numero is a vector of integer
 * repet is an integer which can take any arbitrary value
 * used in the seed of the pseudo-random number generation process
 * if it is increased in repeated calls (e.g. simulation), it is ensured that
 * two calls returns different results (seed=clock+repet)
 ------------------------*/
{
    int i, n;
    int *alea;
    
    n=numero[0];
    vecintalloc (&alea,n);
    
    /*-------------
     * numbering in numero
     -----------*/
    for (i=1;i<=n;i++) {
	numero[i]=i;
    }
    
    /*-------------
     * affects random numbers in alea
     ----------------*/
    for (i=1;i<=n;i++) {
	GetRNGstate();
	alea[i] = ((int) (1e8)*(unif_rand()));
	PutRNGstate();
    }
    
    trirapideint (alea , numero, 1, n);
    freeintvec (alea);
}

/*****************************************/
/* Sorting: used in getpermutation */

void trirapideint (int *x , int *num, int gauche, int droite)
{
    int j, dernier, milieu, t;
    
    if ( (droite-gauche)<=0) return;
    
    milieu = (gauche+droite)/2;
    trirapideintswap (x, gauche, milieu);
    trirapideintswap (num, gauche, milieu);
    
    t=x[gauche];
    dernier=gauche;
    for (j = gauche+1; j<=droite; j++) {
	if (x[j] < t) {
	    dernier = dernier + 1;
	    trirapideintswap (x, dernier, j);	
	    trirapideintswap (num, dernier, j);
	}
    }
    trirapideintswap (x, gauche, dernier);
    trirapideintswap (num, gauche, dernier);
    
    trirapideint (x, num, gauche, dernier-1);
    trirapideint (x, num, dernier+1, droite);
    
}

/**************************************/
/* Sorting: used in trirapideint */

void trirapideintswap (int *v, int i, int j)
{
    int provi;
    
    provi=v[i];
    v[i]=v[j];
    v[j]=provi;
}

/***********************************************************************/
void sqrvec (double *v1)
/*--------------------------------------------------
 * Square root of the elements of a vector
 --------------------------------------------------*/
{
    int i, c1;
    double v2;
    
    c1 = v1[0];
    
    for (i=1;i<=c1;i++) {
	v2 = v1[i];
	/* if (v2 < 0.0) err_message("Error: Square root of negative number (sqrvec)"); */
	v2 = sqrt(v2);
	v1[i] = v2;
    }
}

/***********************************************************************/
void DiagobgComp (int n0, double **w, double *d, int *rang)
/*--------------------------------------------------
 * Eigenstructure of a matrix. See
 * T. FOUCART Analyse factorielle de tableaux multiples,
 * Masson, Paris 1984,185p., p. 62. D'apr?s VPROP et TRIDI,
 * de LEBART et coll.
 --------------------------------------------------*/
{
    double			*s;
    double			a, b, c, x, xp, q, bp, ab, ep, h, t, u , v;
    double			dble;
    int				ni, i, i2, j, k, jk, ijk, ij, l, ix, m, m1, isnou;
    
    vecalloc(&s, n0);
    a = 0.000000001;
    ni = 100;
    if (n0 == 1) {
	d[1] = w[1][1];
	w[1][1] = 1.0;
	*rang = 1;
	freevec (s);
	return;
    }
    
    for (i2=2;i2<=n0;i2++) {
	
	b=0.0;
	c=0.0;
	i=n0-i2+2;
	k=i-1;
	if (k < 2) goto Et1;
	for (l=1;l<=k;l++) {
	    c = c + fabs((double) w[i][l]);
	}
	if (c != 0.0) goto Et2;
	
    Et1:	s[i] = w[i][k];
	goto Etc;
	
    Et2:	for (l=1;l<=k;l++) {
	x = w[i][l] / c;
	w[i][l] = x;
	b = b + x * x;
    }
	xp = w[i][k];
	ix = 1;
	if (xp < 0.0) ix = -1;
		
/*		q = -sqrt(b) * ix; */
	dble = b;
	dble = -sqrt(dble);
	q = dble * ix;
	
	s[i] = c * q;
	b = b - xp * q;
	w[i][k] = xp - q;
	xp = 0;
	for (m=1;m<=k;m++) {
	    w[m][i] = w[i][m] / b / c;
	    q = 0;
	    for (l=1;l<=m;l++) {
		q = q + w[m][l] * w[i][l];
	    }
	    m1 = m + 1;
	    if (k < m1) goto Et3;
	    for (l=m1;l<=k;l++) {
		q = q + w[l][m] * w[i][l];
	    }
	    
	Et3:		s[m] = q / b;
	    xp = xp + s[m] * w[i][m];
	}
	bp = xp * 0.5 / b;
	for (m=1;m<=k;m++) {
	    xp = w[i][m];
	    q = s[m] - bp * xp;
	    s[m] = q;
	    for (l=1;l<=m;l++) {
		w[m][l] = w[m][l] - xp * s[l] - q * w[i][l];
	    }
	}
	for (l=1;l<=k;l++) {
	    w[i][l] = c * w[i][l];
	}
	
    Etc:	d[i] = b;
    } /* for (i2=2;i2<n0;i2++) */
    
    s[1] = 0.0;
    d[1] = 0.0;
    
    for (i=1;i<=n0;i++) {
	
	k = i - 1;
	if (d[i] == 0.0) goto Et4;
	for (m=1;m<=k;m++) {
	    q = 0.0;
	    for (l=1;l<=k;l++) {
		q = q + w[i][l] * w[l][m];
	    }
	    for (l=1;l<=k;l++) {
		w[l][m] = w[l][m] - q * w[l][i];
	    }
	}
	
    Et4:	d[i] = w[i][i];
	w[i][i] = 1.0;
	if (k < 1) goto Et5;
	for (m=1;m<=k;m++) {
	    w[i][m] = 0.0;
	    w[m][i] = 0.0;
	}
	
    Et5:;
    }
    
    for (i=2;i<=n0;i++) {
	s[i-1] = s[i];
    }
    s[n0] = 0.0;
    
    for (k=1;k<=n0;k++) {
	
	m = 0;
	
    Et6: 	for (j=k;j<=n0;j++) {
	if (j == n0) goto Et7;
	ab = fabs((double) s[j]);
	ep = a * (fabs((double) d[j]) + fabs((double) d[j+1]));
	if (ab < ep) goto Et7;
    }
	
    Et7: 	isnou = 1;
	h = d[k];
	if (j == k) goto Eta;
	if (m < ni) goto Etd;
	
	/* err_message("Error: can't compute matrix eigenvalues"); */
	
    Etd:	m = m + 1;
	q = (d[k+1]-h) * 0.5 / s[k];
	
/*		t = sqrt(q * q + 1.0); */
	dble = q * q + 1.0;
	dble = sqrt(dble);
	t = dble;
	
	if (q < 0.0) isnou = -1;
	q = d[j] - h + s[k] / (q + t * isnou);
	u = 1.0;
	v = 1.0;
	h = 0.0;
	jk = j-k;
	for (ijk=1;ijk<=jk;ijk++) {
	    i = j - ijk;
	    xp = u * s[i];
	    b = v * s[i];
	    if (fabs((double) xp) < fabs((double) q)) goto Et8;
	    u = xp / q;
	    
/*			t = sqrt(u * u + 1); */
	    dble = u * u + 1.0;
	    dble = sqrt(dble);
	    t = dble;
	    
	    s[i+1] = q * t;
	    v = 1 / t;
	    u = u * v;
	    goto Et9;
	    
	Et8:		v = q / xp;
	    
/*			t = sqrt(1 + v * v); */
	    dble = 1.0 + v * v;
	    dble = sqrt(dble);
	    t = dble;
	    
	    s[i+1] = t * xp;
	    u = 1 / t;
	    v = v * u;
	    
	Et9:
	    q = d[i+1] - h;
	    t = (d[i] - q) * u + 2.0 * v * b;
	    h = u * t;
	    d[i+1] = q + h;
	    q = v * t - b;
	    for (l=1;l<=n0;l++) {
		xp = w[l][i+1];
		w[l][i+1] = u * w[l][i] + v * xp;
		w[l][i] = v * w[l][i] - u * xp;
	    }
	}
	d[k] = d[k] - h;
	s[k] = q;
	s[j] = 0.0;
	
	goto Et6;
	
    Eta:;
    } /* for (k=1;k<=n0;k++) */
    
    for (ij=2;ij<=n0;ij++) {
	
	i = ij - 1;
	l = i;
	h = d[i];
	for (m=ij;m<=n0;m++) {
	    if (d[m] >= h) {
		l = m;
		h = d[m];
	    }
	}
	if (l == i) {
	    goto Etb;
	} else {
	    d[l] = d[i];
	    d[i] = h;
	}
	for (m=1;m<=n0;m++) {
	    h = w[m][i];
	    w[m][i] = w[m][l];
	    w[m][l] = h;
	}
	
    Etb:;
    } /* for (ij=2;ij<=n0;ij++) */
    
    /* final:; */
    *rang = 0;
    for (i=1;i<=n0;i++) {
	/*
	  if (d[i] / d[1] < 0.00001) d[i] = 0.0;
	  if (d[i] != 0.0) *rang = *rang + 1;
	*/
	if (d[i] > 0.0) *rang = *rang + 1;
    }
    freevec(s);
} /* DiagoCompbg */







/***********************************************************************/
void prodmatABC (double **a, double **b, double **c)
/*--------------------------------------------------
* Matrix product AB
--------------------------------------------------*/
{
    int j, k, i, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];
    
    for (i=1;i<=lig;i++) {
	for (k=1;k<=col2;k++) {
	    s = 0;
	    for (j=1;j<=col;j++) {
		s = s + a[i][j] * b[j][k];
	    }
	    c[i][k] = s;
	}		
    }
}

/***********************************************************************/
void prodmatAtAB (double **a, double **b)
/*--------------------------------------------------
* Matrix product AtA
--------------------------------------------------*/
{
    int j, k, i, lig, col;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    for (j=1;j<=col;j++) {
	for (k=j;k<=col;k++) {
	    s = 0;
	    for (i=1;i<=lig;i++) {
		s = s + a[i][k] * a[i][j];
	    }
	    b[j][k] = s;
	    b[k][j] = s;
	}		
    }
}

/***********************************************************************/
void prodmatAtBC (double **a, double **b, double **c)
/*--------------------------------------------------
 * Matrix product AtB
 --------------------------------------------------*/
{
    int j, k, i, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];
    
    for (j=1;j<=col;j++) {
	for (k=1;k<=col2;k++) {
	    s = 0;
	    for (i=1;i<=lig;i++) {
		s = s + a[i][j] * b[i][k];
	    }
	    c[j][k] = s;
	}		
    }
}


/***********************************************************************/
void prodmatAAtB (double **a, double **b)
/*--------------------------------------------------
 * Matrix product B = AAt
 --------------------------------------------------*/
{
    int j, k, i, lig, col;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    for (j=1;j<=lig;j++) {
	for (k=j;k<=lig;k++) {
	    s = 0;
	    for (i=1;i<=col;i++) {
		s = s + a[j][i] * a[k][i];
	    }
	    b[j][k] = s;
	    b[k][j] = s;
	}		
    }
}

/*******************/
void prodmatAtBrandomC (double **a, double **b, double **c, int*permut)
/*--------------------------------------------------
 * Produit matriciel AtB
 * les lignes de B sont permutees par la permutation permut
 --------------------------------------------------*/
{
    int j, k, i, i0, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];
    
    for (j=1;j<=col;j++) {
	for (k=1;k<=col2;k++) {
	    s = 0;
	    for (i=1;i<=lig;i++) {
		i0 = permut[i];
		s = s + a[i][j] * b[i0][k];
	    }
	    c[j][k] = s;
	}		
    }
}

/***********************************************************************/
void taballoc (double ***tab, int l1, int c1)
/*--------------------------------------------------
 * Dynamic Memory Allocation for a table (l1, c1)
 --------------------------------------------------*/
{
    int i;
    
    *tab = R_Calloc(l1+1, double *);
    for (i=0;i<=l1;i++) {
	*(*tab+i)= R_Calloc(c1+1, double);
    }
    
    **(*tab) = l1;
    **(*tab+1) = c1;
}

/***********************************************************************/
void vecalloc (double **vec, int n)
/*--------------------------------------------------
 * Memory Allocation for a vector of length n
 --------------------------------------------------*/
{
    *vec = R_Calloc(n+1, double);
    **vec = n;
    return;
}

/*****************/
void vecintalloc (int **vec, int n)
/*--------------------------------------------------
 * Memory allocation for an integer vector of length  n
 --------------------------------------------------*/
{
    *vec = R_Calloc(n+1, int);
    **vec = n;
}

/***********************************************************************/
void freetab (double **tab)
/*--------------------------------------------------
 * Free memory for a table
 --------------------------------------------------*/
{
    int 	i, n;
    
    n = *(*(tab));
    for (i=0;i<=n;i++) {
	R_Free(*(tab+i));
    }
    R_Free(tab);
}

/***********************************************************************/
void freevec (double *vec)
/*--------------------------------------------------
 * Free memory for a vector
 --------------------------------------------------*/
{
    R_Free(vec);	
}

/***********************************************************************/
void freeintvec (int *vec)
/*--------------------------------------------------
* Free memory for an integer  vector
--------------------------------------------------*/
{
    
    R_Free(vec);
    
}














/*********************************************************************
 *********************************************************************
 *********                                                       *****
 *********               The sources of adehabitat               *****
 *********               -------------------------               *****
 *********************************************************************
 *********************************************************************
 */



void dedans(double *pts, double *xc, double *yc, double *na,
	    double cs, double **asc)
{
    int nl, nc, i, ligne, colo;
    double x, y;
  
    x = pts[1];
    y = pts[2];
    
    nl = xc[0];
    nc = yc[0];
    
    ligne = 0;
    colo = 0;
    
    for (i = 1; i <= nl; i++) {
	if (((xc[i] - cs/2) <= x) && ((xc[i] + cs/2) > x))
	    ligne = i;
    }
    
    for (i = 1; i <= nc; i++) {
	if (((yc[i] - cs/2) <= y) && ((yc[i] + cs/2) > y))
	    colo = i;
    }
    *na = asc[ligne][colo]; 
}





/* ****************************************************************
   *                                                              *
   *   Randomization of a traject, based on the distances between *
   *   successive relocations and the time lag, as well as the    *
   *   angles between successive steps                            *
   *                                                              *
   **************************************************************** */

void rpath(double **xp, double *rcx, double *rcy, double **asc, 
	   double **tabdist, double *dt, 
	   double *angles, double *xc, double *yc,
	   double *cs, int r)
{
    /* Declaration */
    int i, j, k, l, m, nsam, nltd, *index, *indangles, nlocs;
    double *pts, na, interv, *dobs, dech, anglech, ang;
    
    /* Memory allocation */
    vecalloc(&pts, 2);
    nltd = tabdist[0][0];
    nlocs = xp[0][0];
    vecintalloc(&indangles, (nlocs-2));
    
    /* fills the vector indangle */
    for (i = 1; i <= (nlocs - 2); i++) {
	indangles[i] = i;
    }
    
    /* 1. First loc of the traject */
    k = 0;
    ang = 0;
    
    while (k==0) {
	
	/* Random draw of relocations coordinates */
	xp[1][1] = (alea())*(rcx[2]-rcx[1]) + rcx[1];
	xp[1][2] = (alea())*(rcy[2]-rcy[1]) + rcy[1];
	
	pts[1] = xp[1][1];
	pts[2] = xp[1][2];
	
	/* Verifies that the loc is inside the study area */
	dedans(pts, xc, yc, &na, *cs, asc);
	if (fabs(na + 9999) > 0.000000001)
	    k = 1;
	
    }
    
    
    /* loop for the following relocations */
    for (i = 1; i <= (nlocs-1); i++) {
	interv = dt[i];
	
	/* How many distances for the observed time lag ? */
	nsam = 0;
	for (j = 1; j <= nltd; j++) {
	    if (fabs(tabdist[j][1] - interv) < 0.000000001)
		nsam++;
	}
	
	/* Table of distances */
	vecalloc(&dobs, nsam);
	
	/* the vector index will be used to draw a random relocation */
	vecintalloc(&index, nsam);
	for (l = 1; l <= nsam; l++) {
	    index[l] = l;
	}
	
	/* In a first time, gets the corresponding distances */
	m = 1;
	for (j = 1; j <= nltd; j++) {
	    if (fabs(tabdist[j][1] - interv) < 0.000000001) {
		dobs[m] = tabdist[j][2];
		m++;
	    }
	}
	
	k = 0;
	while (k == 0) {
	    /* Sampled Distance */
	    r = (int) (alea() * 100);
	    getpermutation(index, j * r);
	    dech = dobs[index[1]];
	    
	    /* Sampled Angles */
	    getpermutation(indangles, j * r);
	    anglech = angles[indangles[1]];
	    
	    /* update the angles */
	    ang = (ang + anglech);
	    
	    /* The new coordinates */
	    xp[i+1][1] = xp[i][1] + dech * cos(ang);
	    xp[i+1][2] = xp[i][2] + dech * sin(ang);
	    
	    pts[1] = xp[i+1][1];
	    pts[2] = xp[i+1][2];
	    
	    dedans(pts, xc, yc, &na, *cs, asc);
	    if (fabs(na + 9999) > 0.000000001)
		k = 1;
	}
	freevec(dobs);
	freeintvec(index);
    }
    
    /* Free memory */
    freeintvec(indangles);
    freevec(pts);
}







/* Test of rpath with R */

void randpath(double *xpr, double *rcrx, double *rcry, double *ascr, 
	      double *xcr, double *ycr, double *csr,
	      double *tabdistr, double *dtr, double *anglesr, 
	      int *nlasc, int *ncasc, int *nltdr, int *nlocsr)
{
    /* declaration of the variables */
    int i, j, k, r, nlocs, nltd;
    double **xp, *rcx, *rcy, **asc, **tabdist, *dt, *angles;
    double *xc,*yc, cs;
    
    /* Memory allocation */
    nlocs = *nlocsr;
    nltd = *nltdr;
    cs = *csr;
    
    taballoc(&xp, nlocs, 2);
    vecalloc(&rcx, 2);
    vecalloc(&rcy, 2);
    taballoc(&asc, *nlasc, *ncasc);
    vecalloc(&xc, *nlasc);
    vecalloc(&yc, *ncasc);
    taballoc(&tabdist, nltd, 2);
    vecalloc(&dt, nlocs-1);
    vecalloc(&angles, nlocs-2);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= nlocs; i++) {
	for (j = 1; j <= 2; j++) {
	    xp[i][j] = xpr[k];
	    k++;
	}
    }
    rcx[1] = rcrx[0];
    rcx[2] = rcrx[1];
    rcy[1] = rcry[0];
    rcy[2] = rcry[1];
    
    k = 0;
    for (i = 1; i <= *nlasc; i++) {
	for (j = 1; j <= *ncasc; j++) {
	    asc[i][j] = ascr[k];
	    k++;
	}
    }
    
    k = 0;
    for (i = 1; i <= nltd; i++) {
	for (j = 1; j <= 2; j++) {
	    tabdist[i][j] = tabdistr[k];
	    k++;
	}
    }
    
    for (i = 1; i <= (nlocs - 1); i++) {
	dt[i] = dtr[i-1];
    }
    for (i = 1; i <= *nlasc; i++) {
	xc[i] = xcr[i-1];
    }
    for (i = 1; i <= *ncasc; i++) {
	yc[i] = ycr[i-1];
    }
    for (i = 1; i <= (nlocs - 2); i++) {
	angles[i] = anglesr[i-1];
    }
    
    /* test of randpath */
    r = 1;
    rpath(xp, rcx, rcy, asc, tabdist, dt, angles,
	  xc, yc, &cs, r);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= nlocs; i++) {
	for (j = 1; j <= 2; j++) {
	    xpr[k] = xp[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(xp);
    freevec(rcx);
    freevec(rcy);
    freetab(asc);
    freetab(tabdist);
    freevec(dt);
    freevec(angles);
    freevec(xc);
    freevec(yc);
}




/* *********************************************************************
 *                                                                     *
 *                   First passage time                                *
 *                                                                     *
 ***********************************************************************/

/* compute the distance cetween two points */
void dtmp(double x1, double x2, double y1, double y2, 
	  double *di)
{
  *di = sqrt(((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)));
}



/* compute the FPT for ONE relocation */
void fptt(double *x, double *y, double *t, int pos, double radius, double *fptto, int nlo)
{
    /* Declaration */
    int ok, pos2, naar, naav, na;
    double di, dt, dt2, di2, fptar, fptav;
    
    ok = 0;
    di = 0;
    di2 = 0;
    dt = 0;
    dt2 = 0;
    naar = 1;
    naav = 1;
    fptar = 0;
    fptav = 0;
  
  
    /* Search of the first loc outside the circle (before) */
    pos2 = pos;
    while (ok == 0) {
	pos2 = pos2 - 1;
	if (pos2 > 0) {
	    dtmp(x[pos2], x[pos], y[pos2], y[pos], &di);
	    if (di >= radius)
		ok = 1;
	} else {
	    ok = 1;
	    naar = 0;
	}
    }
    
    /* computes the linear approximation */
    if (naar > 0) {
	dt = fabs(t[pos] - t[pos2]);
	dt2 = fabs(t[pos] - t[(pos2+1)]);
	dtmp(x[(pos2+1)], x[pos], y[(pos2+1)], y[pos], &di2);
	fptar = dt2 + ( (dt - dt2) * (radius - di2) / (di - di2) );
    }
    
    
    /* Search of the first loc outside the circle (after) */
    pos2 = pos;
    ok = 0;
    while (ok == 0) {
	pos2 = pos2 + 1;
	if (pos2 <= nlo) {
	    dtmp(x[pos2], x[pos], y[pos2], y[pos], &di);
	    if (di >= radius)
		ok = 1;
	} else {
	    ok = 1;
	    naav = 0;
	}
    }
    
    /* Computes linear approximation */
    if (naav > 0) {
	dt = fabs(t[pos2] - t[pos]);
	dt2 = fabs(t[(pos2-1)] - t[pos]);
	dtmp(x[(pos2-1)], x[pos], y[(pos2-1)], y[pos], &di2);
	fptav = dt2 + ( (dt - dt2) * (radius - di2) / (di - di2) );
    }
    
    na = naar * naav;
    if (na > 0) {
	*fptto = fptar + fptav;
    } else {
	*fptto = -1;
    }
    
}



/* Computes the FPT for all relocations */
void fipati(double *x, double *y, double *t, 
	    int nlo, int nra, double *rad, 
	    double **fpt)
{
    /* Declaration */
    int i, j;
    double val;
    
    /* Computes the FPT */
    for (i=1; i<=nra; i++) {
	for (j=1; j<=nlo; j++) {
	    fptt(x, y, t, j, rad[i], &val, nlo);
	    fpt[j][i] = val;
	}
    }
}

/* for external call from within R */
void fipatir(double *xr, double *yr, double *tr, 
	     int *nlocs, double *radius, int *nrad, 
             double *fptr)
{
    /* Declaration */
    int i, j, k, nlo, nra;
    double *x, *y, *t, *rad, **fpt;
    
    /* Memory allocation */
    nlo = *nlocs;
    nra = *nrad;
    
    vecalloc(&x, nlo);
    vecalloc(&y, nlo);
    vecalloc(&t, nlo);
    vecalloc(&rad, nra);
    taballoc(&fpt, nlo, nra);
  
    /* R to C */
    for (i = 1; i <= nlo; i++) {
	x[i] = xr[i-1];
	y[i] = yr[i-1];
	t[i] = tr[i-1];
    }
    
    for (i = 1; i <= nra; i++) {
	rad[i] = radius[i-1];
    }
    
    /* main function */
    fipati(x,y,t, nlo, nra, rad, fpt);
    
    /* C to R */
    k = 0;
    for (i=1; i<= nlo; i++) {
	for (j = 1; j<=nra; j++) {
	    fptr[k]=fpt[i][j];
	    k++;
	}
    }
    
    /* free memory */
    freetab(fpt);
    freevec(x);
    freevec(y);
    freevec(t);
    freevec(rad);
}






/* *********************************************************************
 *                                                                     *
 *                   Percolation cluster                               *
 *                                                                     *
 ***********************************************************************/


void perclu(double **map, int nr, int nc, double *x, double *y,
	    int nmax, int *nreel, double *pm)
{
    /* Declaration */
    int i,j, encore, k, l, dir, *vois, *rvois, *cvois, choix, xt, yt, cons, *reord, len;
    double **rr, **cc, *cs, ptir;
    
    /* Memory allocation */
    xt = (int) x[1];
    yt = (int) y[1];
    len = nr;
    if (nc < nr)
	len = nc;
    encore = 1;
    i = 1;
    j = 1;
    l = 0;
    k = 2;
    ptir = 0;
    dir = 1;
    choix = 1;
    cons = 0;
    
    taballoc(&rr, nr, nc);
    taballoc(&cc, nr, nc);
    vecintalloc(&vois, 4);
    vecintalloc(&reord, 4);
    vecintalloc(&rvois, 4);
    vecintalloc(&cvois, 4);
    vecalloc(&cs, 4);
    
    
    /* Rows and columns matrices */
    for (i = 1; i <= nr; i++) {
	for (j = 1; j <= nc; j++) {
	    rr[i][j] = (double) i;
	    cc[i][j] = (double) j;
	}
    }
    
    
    cs[1] = pm[1];
    for (i = 2; i <= 4; i++) {
	cs[i] = cs[i-1] + pm[i];
    }
    
    /* Beginning of the loop */
    while (encore == 1) {
    
	/* Storage of the neighbouring */
	vois[1] = (int) map[xt][yt+1];
	vois[2] = (int) map[xt+1][yt];
	vois[3] = (int) map[xt][yt-1];
	vois[4] = (int) map[xt-1][yt];
	
	rvois[1] = (int) rr[xt][yt+1];
	rvois[2] = (int) rr[xt+1][yt];
	rvois[3] = (int) rr[xt][yt-1];
	rvois[4] = (int) rr[xt-1][yt];
	
	cvois[1] = (int) cc[xt][yt+1];
	cvois[2] = (int) cc[xt+1][yt];
	cvois[3] = (int) cc[xt][yt-1];
	cvois[4] = (int) cc[xt-1][yt];
	
	/* Re-order of the neighbouring according to the direction */
	l = 1;
	for (i = dir; i <= 4; i++) {
	    reord[l] = i;
	    l++;
	}
	i = 1;
	while (l != 4) {
	    reord[l] = i;
	    l++;
	    i++;
	}
	
	/* random draw of a direction */
	ptir = alea();
	choix = 4;
	if (ptir <= pm[1]) {
	    choix = 1;
	}
	if ((ptir > pm[1])&&(ptir <= pm[2])) {
	    choix = 2;
	}
	if ((ptir > pm[2])&&(ptir <= pm[3])) {
	    choix = 3;
	}
	
	/* And again, until the direction lead us into an accessible area */
	cons = reord[choix];
	
	while (vois[cons] == 1) {
	    ptir = alea();
	    choix = 4;
	    if (ptir <= pm[1]) {
		choix = 1;
	    }
	    if ((ptir > pm[1])&&(ptir <= pm[2])) {
		choix = 2;
	    }
	    if ((ptir > pm[2])&&(ptir <= pm[3])) {
		choix = 3;
	    }
	    cons = reord[choix];
	}
	
	/* Stores all information */
	xt = (int) (rvois[choix]);
	yt = (int) (cvois[choix]);
	
	x[k] = (double) xt;
	y[k] = (double) yt;
	
	if ((xt==1)|(xt==len)|(yt==1)|(yt==len))
	    encore = 0;
	if (k==nmax)
	    encore = 0;
	
	for (i = 1; i <= 4; i++) {
	    if ( ((int) rvois[i]) == xt) {
		if ( ((int) cvois[i]) == yt) {
		    dir = i;
		}
	    }
	}
	
	*nreel = k;
	k++;
    }
    
    /* free memory */
    freeintvec(vois);
    freeintvec(rvois);
    freeintvec(cvois);
    freeintvec(reord);
    freevec(cs);
    freetab(rr);
    freetab(cc);
}





/* For external call from within R */
void perclur(double *mapr, int *nrm, int *ncm, double *probamr,
	     double *xr, double *yr, int *nmaxr, int *nreel)
{
    /* Declaration */
    double **map, *pm, *x, *y;
    int i, j, k, nr, nc, nmax;
    
    /* Memory allocation */
    nr = *nrm;
    nc = *ncm;
    nmax = *nmaxr;
    
    taballoc(&map, nr, nc);
    vecalloc(&x, nmax);
    vecalloc(&y, nmax);
    vecalloc(&pm, 4);
    
    /* R to C */
    x[1] = xr[0];
    y[1] = yr[0];
    
    k = 0;
    for (i = 1; i <= nr; i++) {
	for (j = 1; j <= nc; j++) {
	    map[i][j] = mapr[k];
	    k++;
	}
    }
    
    for (i = 1; i <= 4; i++) {
	pm[i] = probamr[i-1];
    }
    
    /* Main function */
    perclu(map, nr, nc, x, y, nmax, nreel, pm);
    
    /* C to R */
    for (i = 1; i <= *nreel; i++) {
	xr[i-1] = x[i];
	yr[i-1] = y[i];
    }
    
    /* free memory */
    freevec(x);
    freevec(y);
    freevec(pm);
    freetab(map);
}



/* *********************************************************************
 *                                                                     *
 *         Rediscretization algorithm for a traject                    *
 *                                                                     *
 ***********************************************************************/



/* Resolves a quadratic equation of the type
 */

void resolpol(double a, double b, double c, double *x1, double *x2, int *warn)
{
    double delta;
    delta = (b * b) - 4 * a * c;
    *warn = 0;
    if (delta > 0) {
	*x1 = (-b - sqrt(delta)) / (2 * a);
	*x2 = (-b + sqrt(delta)) / (2 * a);
    } else {
	*warn = 1;
    }
}




void discretraj(double *x, double *y, double *dat, double *xn, 
		double *yn, int n, int nn, double *datn, 
		double u, int *neff)
{
    /* Declaration */
    double R, xt, yt, a, b, c, pente, ori, x1, x2, di1, di2;
    int k, m, p, fini, ok, warn, *dedans, lo, new, pp;
    
    /* memory allocation */
    fini = 0;
    k = 1;
    p = 2;
    m = 1;
    ok = 0;
    a = 0;
    b = 0;
    c = 0;
    pente = 0;
    ori = 0;
    x1 = 0;
    x2 = 0;
    lo = 0;
    di1 = 0;
    di2 = 0;
    *neff = 0;
    new = 0;
    pp = 1;
    
    
    vecintalloc(&dedans,2);
    
    /* Main algorithm */
    while (fini == 0) {
	
	dedans[1] = 0;
	dedans[2] = 0;
	ok = 0;
	xt = xn[k];
	yt = yn[k];
	k++;
	new = 0;
	
	/* Determines the "upper" point */
	while (ok == 0) {
	    if (new == 1)
		p++;
	    R = sqrt((((x[p] - xt) * (x[p] - xt)) + ((y[p] - yt) * (y[p] - yt))));
	    if (R > u) {
		ok = 1;
	    } else {
		if (p == n) {
		    fini = 1;
		    ok = 1;
		}
	    }
	    new = 1;
	}
	m = p-1;
    
    if (fini == 0) {
	/* Does the difference between x[p] and x[m] = 0? */
	if ((fabs(x[p] - x[m]) > 0.000000000001)) {
	    /* Computes the slope between m and p */
	    pente = (y[p] - y[m]) / (x[p] - x[m]); /* when diff(x) == 0 ? */
	    /* The intercept */
	    ori = y[p] - (pente * x[p]);
	    /* The parameters of the polynomial equation */
	    a = 1 + (pente * pente);
	    b = (-2 * xt) + (2 * pente * ori) - (2 * pente * yt);
	    c = (xt * xt) + (yt * yt) + (ori * ori) - (2 * ori * yt) - (u * u);
	    resolpol(a, b, c, &x1, &x2, &warn);
	    /* 
	       A line cuts a circle with radius u at two points. One has 
	       (i) to identify the point the closest from m,n and 
	       (ii) to keep the one on the segment m-p
	    */
	    
	    
	    /* Which one are in the interval ? */
	    if (x1 >= x[m]) {
		if (x1 < x[p]) {
		    dedans[1] = 1;
		    lo = 1;
		}
	    }
	    if (x1 >= x[p]) {
		if (x1 < x[m]) {
		    dedans[1] = 1;
		    lo = 1;
		}
	    }
	    if (x2 >= x[m]) {
		if (x2 < x[p]) {
		    dedans[2] = 1;
		    lo = 2;
		}
	    }
	    if (x2 >= x[p]) {
		if (x2 < x[m]) {
		    dedans[2] = 1;
		    lo = 2;
		}
	    }
	    
	    /* What is the minimum distance to m ? */
	    if ((dedans[1] + dedans[2]) > 1) {
		di1 = fabs((double) (x[p] - x1));
		di2 = fabs((double) (x[p] - x2));
		
		/* verify that xk-1 is not in the same interval. Otherwise one increase of 1 */
		if (di1 < di2) {
		    lo = 2;
		} else {
		    lo = 1;
		}
		if (pp == p) {
		    if (di1 < di2) {
			lo = 1;
		    } 
		    if (di2 < di1) {
			lo = 2;
		    } 
		}
	    }
	    
	    /* storage of the coordinates */
	    if (lo == 1) {
		xn[k] = x1;
		yn[k] = (pente * x1) + ori;
	    }
	    if (lo == 2) {
		xn[k] = x2;
		yn[k] = (pente * x2) + ori;
	    }

	} else { /* We change x and y coordinates */
	    
	    /* Computes the slope between m and p */
	    pente =  (x[p] - x[m]) / (y[p] - y[m]);
	    /* The intercept */
	    ori = x[p] - (pente * y[p]);
	    /* The parameters of the polynomial equation */
	    a = 1 + (pente * pente);
	    b = (-2 * yt) + (2 * pente * ori) - (2 * pente * xt);
	    c = (xt * xt) + (yt * yt) + (ori * ori) - (2 * ori * xt) - (u * u);
	    resolpol(a, b, c, &x1, &x2, &warn);
	    /* 
	       A line cuts a circle with radius u at two points. One has 
	       (i) to identify the point the closest from m,n and 
	       (ii) to keep the one on the segment m-p
	    */
	    
	    
	    /* Which one are in the interval ? */
	    if (x1 >= y[m]) {
		if (x1 < y[p]) {
		    dedans[1] = 1;
		    lo = 1;
		}
	    }
	    if (x1 >= y[p]) {
		if (x1 < y[m]) {
		    dedans[1] = 1;
		    lo = 1;
		}
	    }
	    if (x2 >= y[m]) {
		if (x2 < y[p]) {
		    dedans[2] = 1;
		    lo = 2;
		}
	    }
	    if (x2 >= y[p]) {
		if (x2 < y[m]) {
		    dedans[2] = 1;
		    lo = 2;
		}
	    }
	    
	    /* What is the minimum distance to m ? */
	    if ((dedans[1] + dedans[2]) > 1) {
		di1 = fabs((double) (y[p] - x1));
		di2 = fabs((double) (y[p] - x2));
		
		/* verify that yk-1 is not in the same interval. Otherwise one increase of 1 */
		if (di1 < di2) {
		    lo = 2;
		} else {
		    lo = 1;
		}
		if (pp == p) {
		    if (di1 < di2) {
			lo = 1;
		    } 
		    if (di2 < di1) {
			lo = 2;
		    } 
		}
	    }
	    
	    /* storage of the coordinates */
	    if (lo == 1) {
		yn[k] = x1;
		xn[k] = (pente * x1) + ori;
	    }
	    if (lo == 2) {
		yn[k] = x2;
		xn[k] = (pente * x2) + ori;
	    }
	}
	
	/* Computes the nnew date (linear approximation) */
	di1 = sqrt((((xn[k] - x[m]) * (xn[k] - x[m])) + ((yn[k] - y[m]) * (yn[k] - y[m]))));
	R = sqrt((((x[p] - x[m]) * (x[p] - x[m])) + ((y[p] - y[m]) * (y[p] - y[m]))));
	di2 = dat[p] - dat[m];
	datn[k] = dat[m] + (di1 * di2 / R);
    }
    if (k == nn) {
	fini = 1;
    }
    pp = p;
    }
    
    /* Free memory */
    *neff = k;
    freeintvec(dedans);
}



/* For external Call from within R */

void discretrajr(double *xr, double *yr, double *datr, double *xnr, 
		 double *ynr, int *nr, int *nnr, double *datnr, 
		 double *xdeb, double *ydeb, double *ur, double *dat0, int *neff)
{
    /* Declaration */
    int i, n, nn;
    double *x, *y, *xn, *yn, *dat, *datn, u;
    
    /* Memory allocation */
    n = *nr;
    nn = *nnr;
    u = *ur;
    
    vecalloc(&x, n);
    vecalloc(&y, n);
    vecalloc(&xn, nn);
    vecalloc(&yn, nn);
    vecalloc(&dat, n);
    vecalloc(&datn, nn);
    
    /* R to C */
    for (i = 1; i <= n; i++) {
	x[i] = xr[i-1];
	y[i] = yr[i-1];
	dat[i] = datr[i-1];
    }
    
    xn[1] = *xdeb;
    yn[1] = *ydeb;
    datn[1] = *dat0;
    
    /* Main function  */
    discretraj(x, y, dat, xn, yn, n, nn, datn, u, neff);
    
    /* C to R */
    for (i = 1; i <= nn; i++) {
	xnr[i-1] = xn[i];
	ynr[i-1] = yn[i];
	datnr[i-1] = datn[i];
    }
    
    /* Free memory */
    freevec(x);
    freevec(y);
    freevec(xn);
    freevec(yn);
    freevec(dat);
    freevec(datn);
}








void permutR2n(double *xyr, int *nro, int *nrepr, 
	       double *R2nr, double *dtr, double *dtsimr)
{
  double **xy, **R2n, *xp, *dt, **dtsim, tt;
  int n, i, j, k, *index, nr;
  
  n = *nro;
  nr = *nrepr;
  tt = 0;
  vecalloc(&xp, 2);
  vecalloc(&dt, n);
  taballoc(&R2n, n, nr);
  taballoc(&dtsim, n, nr);
  vecintalloc(&index, n);
  taballoc(&xy, n, 2);
  
  k = 0;
  for (i = 1; i <= n; i++) {
    dt[i] = dtr[i-1];
    for (j = 1; j <= 2; j++) {
      xy[i][j] = xyr[k];
      k++;
    }
  }
  
  for (k = 1; k <= nr; k++) {
    
    getpermutation(index, k);
    j = index[1];
    xp[1] = 0;
    xp[2] = 0;
    tt = 0;

    for (i = 1; i <= n; i++) {
      j = index[i];
      xp[1] = xp[1] + xy[j][1];
      xp[2] = xp[2] + xy[j][2];
      R2n[i][k] = (xp[2] * xp[2]) + (xp[1] * xp[1]);
      tt = tt + dt[j];
      dtsim[i][k] = tt;
    }
  }
  
  
  k = 0;
  for (i = 1; i <= n; i++) {
    for (j = 1; j<= nr; j++) {
      R2nr[k] = R2n[i][j];
      dtsimr[k] = dtsim[i][j];
      k++;
    }
  }
  
  
  freevec(xp);
  freevec(dt);
  freeintvec(index);
  freetab(xy);
  freetab(dtsim);
  freetab(R2n);
}




void runsltr(int *xr, int *nr, double *res, int *nrepr)
{
  int i, j, n, *x, *xb, nrep, nbsui, nz, nu, *numero;
  double m, s, n1, n2;
  
  n = *nr;
  nrep = *nrepr;
  
  vecintalloc(&x, n);
  vecintalloc(&xb, n);
  vecintalloc(&numero, n);

  
  for (i = 1; i <= n; i++) {
    x[i] = xr[i-1];
    numero[i] = i;
  }
  
  nbsui = 1;
  nz = 0;
  nu = 1;
  if (x[1] == 0)
    nz = 1;
  if (x[1] == 1)
    nu = 1;
  
  for (i = 2; i <= n; i++) {
    if (x[i] == 0)
      nz++;
    if (x[i] == 1)
      nu++;
    if (x[i-1] != x[i])
      nbsui++;
  }
  
  n1 = ((double) nz);
  n2 = ((double) nu);
  
  m = 1 + 2 * n1 * n2 / (n1 + n2);
  s = sqrt(2 * n1 * n2 * (2 * n1 * n2 - n1 - n2)/((n1 + n2) * (n1 + n2) * (n1 + n2 - 1)));
  
  res[0] = (((double) nbsui) - m) / s;
  
  
  for (j = 1; j <= nrep; j++) {
    nbsui = 1;
    getpermutation(numero, j);
    trirapideint(numero , x, 1, n);
    for (i = 2; i <= n; i++) {
      if (x[i-1] != x[i])
	nbsui++;
    }
    res[j] = (((double) nbsui) - m) / s;
  }
  
  freeintvec(x);
  freeintvec(xb);
  freeintvec(numero);

}




void testindepangl (double *sim, double *ang, int *nang, int *debut, int *fin, int *ndeb, int *ni){
  
  int i,j,k;
  double *angle;
  vecalloc(&angle, *nang);
  for (i=1; i<=*nang; i++) {
    angle[i] = ang[i-1];
  }
  for (k=0;k<=(*ndeb-1);k++){
    for (i=debut[k]; i<=(fin[k]-1); i++) {
      sim[0]=sim[0]+(1.0-cos(angle[i+1] - angle[i]));
    }
  }
  for (j=1; j<=*ni; j++) {
    aleapermutvec(angle);
    for (k=0;k<=(*ndeb-1);k++){
      for (i=debut[k]; i<=(fin[k]-1); i++) {
	sim[j]=sim[j]+(1.0-cos(angle[i+1] - angle[i]));
      }
    }
  }
  for (j=0; j<=*ni; j++) {
    sim[j]=2*sim[j];

  }
  freevec(angle);
}

void testindepdist (double *sim, double *di, int *ndi, int *debut, int *fin, int *ndeb, int *ni){
  
  int i,j,k;
  double *dist;
  vecalloc(&dist, *ndi);
  for (i=1; i<=*ndi; i++) {
    dist[i] = di[i-1];
  }
  for (k=0;k<=(*ndeb-1);k++){
    for (i=(debut[k]); i<=(fin[k]-1); i++) {
      sim[0]=sim[0]+pow(dist[i+1] - dist[i],2);
    }
  }
  for (j=1; j<=*ni; j++) {
    aleapermutvec(dist);
    for (k=0;k<=(*ndeb-1);k++){
      for (i=(debut[k]); i<=(fin[k]-1); i++) {
	sim[j]=sim[j]+pow(dist[i+1] - dist[i],2);
      }
    }

  }
  freevec(dist);
}


void prepquart (double *dtur, int *nur, double *dtrr, double *resr, int *nrr, 
		int *ncr, double *rur)
{
  int i,j,k, nu, nr, nc, *rt;
  double *dtu, **dtr, **res, **ru, tmp1;
  
  nu = *nur;
  nr = *nrr;
  nc = *ncr;
  
  vecalloc(&dtu, nu);
  vecintalloc(&rt, nc);
  taballoc(&dtr, nr, nc);
  taballoc(&res, nr, nc);
  taballoc(&ru, nu, nc);
  
  for (i = 1; i <= nu; i++) {
    dtu[i] = dtur[i-1];
  }
  
  k = 0;
  for (i = 1; i <= nr; i++) {
    for (j = 1; j<= nc; j++) {
      dtr[i][j] = dtrr[k];
      res[i][j] = resr[k];
      k++;
    }
  }
  
  for (i = 1; i <= nu; i++) {
    for (k = 1; k <= nc; k++) {
      rt[k] = 0;
    }
    tmp1 = dtu[i];
    for (j = 1; j <= nr; j++) {
      for (k = 1; k <= nc; k++) {
	if ((fabs((double) (dtr[j][k] - tmp1)))< 0.0000000001)
	  rt[k] = j;
      }
    }
    for (k = 1; k <= nc; k++) {
      if ((fabs((double) rt[k] )) < 0.00000000001) {
	ru[i][k] = -1;
      } else {
	j = rt[k];
	ru[i][k] = res[j][k];
      }
    }
  }
  
  
  k = 0;
  for (i = 1; i <= nu; i++) {
    for (j = 1; j<= nc; j++) {
      rur[k] = ru[i][j];
      k++;
    }
  }
  
  freevec(dtu);
  freeintvec(rt);
  freetab(dtr);
  freetab(res);
  freetab(ru);
}





/* *********************************************************************
 *                                                                     *
 *                   partition d'un trajet                             *
 *                                                                     *
 ***********************************************************************/

void optcut (double **Pid, double **mk, int *maxk) 
{
    /* Declaration of variables */
    int l, p, i, ii, Km, k, j;
    double **mkd, **mkdn, tmp, mi1, mi2, mi3;
    
    
    
    /* memory allocation */
    l = Pid[0][0];   /* The length of the sequence */
    p = Pid[1][0];   /* The number of models */
    Km = *maxk;     /* the max number of partition to compute */
    i = 0;     /* the sequence index used in the paper (from 0 to l-1) */
    ii = 0;    /* the position index used in the program (from 1 to l) */
    j = 0;    /* the position index for the model (from 1 to p) */
    k = 0;    /* the index for the partition */
    taballoc(&mkd, l, p); /* The table for the log-probabilities of a 
			     partition given a model 
			     for a k partition
			     (to be re-used for each k) */
    taballoc(&mkdn, l, p); /* The table for the log-probabilities 
			      of a partition given a model
			     for a k+1 partition
			     (to be re-used for each k) */
    mi1 = 0;
    mi2 = 0;
    mi3 = 0;
    
    
    /* 1. First computes the probability of one partition*/
    
    /* 1.1 Fills the table mkd   */
    
    for (j = 1; j <= p; j++) {
	mkd[1][j] = log(Pid[1][j]);
    }
    for (j = 1; j <= p; j++) {
	for (ii = 2; ii <= l; ii++) {
	    mkd[ii][j] = mkd[ii-1][j] + log(Pid[ii][j]);
	}
    }
    
    /* 1.2 Fills the table mk */
    for (ii = 1; ii <= l; ii++) {
	
        /* computes the max */
	tmp = mkd[ii][1];
	for (j = 1; j <= p; j++) {
	    if (mkd[ii][j] - tmp > 0.0000000000001) {
		tmp = mkd[ii][j];
	    }
	}
	
	/* removes the max */
	for (j = 1; j <= p; j++) {
	    mkd[ii][j] = mkd[ii][j] - tmp;
	}
	
	/* computes the mean */
	mk[ii][1] = exp(mkd[ii][1]) / ( (double) p ) ;
	for (j = 2; j <= p; j++) {
	    mk[ii][1] = mk[ii][1] + (exp(mkd[ii][j]) / ( (double) p ) );
	}
	
        /* takes the log and adds the max again */
	mk[ii][1] = log(mk[ii][1]) + tmp;
	
	/* adds the max again for mkd */
	for (j = 1; j <= p; j++) {
	    mkd[ii][j] = mkd[ii][j] + tmp;
	}
    }
    
    
    
    
    /* 2. computes mk for each possible r-partitions: a loop */
    for (k = 2; k <= Km; k++) {
	
	
        /* First start at the first step */
	i = k-1;
	ii = k;
	
	/* Computes mkdn for each d and ii */
	for (j = 1; j <= p; j++) {
	    
	    /* Removes the max */
	    tmp = mk[ii-1][k-1];
	    if (mk[ii-1][k-1] - mkd[ii-1][j] < -0.00000000001)
		tmp = mkd[ii-1][j];
	    mi1 = mk[ii-1][k-1] - tmp;
	    mi2 = mkd[ii-1][j] - tmp;
	    
	    /* Computes for i = k-1 */
	    mkdn[ii][j] = log(Pid[ii][j]) + 
		log((((double) k) - 1) / (((double) i) * (((double) p) - 1))) +
		log(((double) p) * exp(mi1) - exp(mi2));
	    
	    /* adds the max again */
	    mkdn[ii][j] = mkdn[ii][j] + tmp;
	}
	


	/* fills mk */
	/* computes the max */
	tmp = mkdn[ii][1];
	for (j = 1; j <= p; j++) {
	    if (mkdn[ii][j] - tmp > 0.00000000000000001) {
		tmp = mkdn[ii][j];
	    }
	}
	
	/* removes the max */
	for (j = 1; j <= p; j++) {
	    mkdn[ii][j] = mkdn[ii][j] - tmp;
	}
	
	/* computes the mean */
	mk[ii][k] = exp(mkdn[ii][1]) / ( (double) p ) ;
	for (j = 2; j <= p; j++) {
	    mk[ii][k] = mk[ii][k] + (exp(mkdn[ii][j]) / ( (double) p ) );
	}
	
	/* takes the log and adds the max again */
	mk[ii][k] = log(mk[ii][k]) + tmp;
	
	/* adds the max again for mkd */
	for (j = 1; j <= p; j++) {
	    mkdn[ii][j] = mkdn[ii][j] + tmp;
	}
	
	
	/* Then increase the sequence */
	for (ii = (k + 1) ; ii <= l; ii++) {
	    
	    i = ii - 1;
	    
	    /* again computes mkdn for each d */
	    for (j = 1; j <= p; j++) {
		
		tmp = mkdn[ii-1][j];
		if (tmp - mk[ii-1][k-1] < -0.00000000000001) {
		    tmp = mk[ii-1][k-1];
		}
		if (tmp - mkd[ii-1][j] < -0.000000000000000001) {
		    tmp = mkd[ii-1][j];
		} 
		mi1 = mkdn[ii-1][j] - tmp;
		mi2 = mk[ii-1][k-1] - tmp;
		mi3 = mkd[ii-1][j] - tmp;
		
		mkdn[ii][j] = log(Pid[ii][j]) +
		    log((( ((double) (i - k + 1)) / ((double) i)) * 
			exp(mi1)) + (((((double) k) - 1) / 
				    (((double) i) * 
				     (((double) p) - 1))) * 
			((((double) p) * exp(mi2)) - 
			 exp(mi3)))) + tmp;
	    }


	    
	    /* ... and again fills mk */
	    /* computes the max */
	    tmp = mkdn[ii][1];
	    for (j = 1; j <= p; j++) {
		if (mkdn[ii][j] - tmp > 0.000000000000001) {
		    tmp = mkdn[ii][j];
		}
	    }
	    
	    /* removes the max */
	    for (j = 1; j <= p; j++) {
		mkdn[ii][j] = mkdn[ii][j] - tmp;
	    }
	    
	    /* computes the mean */
	    mk[ii][k] = exp(mkdn[ii][1]) / ( (double) p ) ;
	    for (j = 2; j <= p; j++) {
		mk[ii][k] = mk[ii][k] + (exp(mkdn[ii][j]) / ( (double) p ) );
	    }
	    
	    /* takes the log and adds the max again */
	    mk[ii][k] = log(mk[ii][k]) + tmp;
	    
	    /* adds the max again for mkd */
	    for (j = 1; j <= p; j++) {
		mkdn[ii][j] = mkdn[ii][j] + tmp;
	    }
	    
	}
	
	/* and finally, replace mkd with mkdn for next k */    
	for (ii = 1; ii <= l; ii++) {
	    for (j = 1; j <= p; j++) {
		mkd[ii][j] = mkdn[ii][j];
	    }
	}
    }
    
    
    /* free memory */
    freetab(mkd);
    freetab(mkdn);

}



/* Now, the R version */

void optcutr (double *Pidr, double *mkr, int *maxk, int *lr, int *pr, 
	      double *baye) 
{
    /* Declaration */
    int i, j, k, l, p, Km;
    double **Pid, **mk, tmp, *mik, *msk, *kk;
    
    /* Memory allocation */
    l = *lr;
    Km = *maxk;
    p = *pr;
    taballoc(&Pid, l, p);
    taballoc(&mk, l, Km);
    vecalloc(&mik, Km);
    vecalloc(&msk, Km);
    vecalloc(&kk, Km);

    
    /* Fills local variables */
    k = 0;
    for (i = 1; i <= l; i++) {
	for (j = 1; j <= p; j++) {
	    Pid[i][j] = Pidr[k];
	    k++;
	}
    }
    
    optcut(Pid, mk, maxk);
    
    /* keeps the last row */
    for (k = 1; k <= Km; k++) {
	kk[k] = mk[l][k];
    }
    /* computes the max */
    tmp = kk[1];
    for (k = 1; k <= Km; k++) {
	if (tmp - kk[k] < -0.0000000001)
	    tmp = kk[k];
    }
    
    /* removes the max */
    for (k = 1; k <= Km; k++) {
	kk[k] = kk[k] - tmp;
    }
    
    /* Computes mik */
    mik[1] = exp(kk[1]);
    for (k = 2; k <= Km; k++) {
	mik[k] = mik[k-1] + (exp(kk[k]) / ((double) k ));
    }
    
    /* adds the max again */
    for (k = 1; k <= Km; k++) {
	mik[k] = log(mik[k]) + tmp;
    }
    
    /* Computes msk */
    msk[Km] = exp(kk[Km]) / ((double) (l - Km + 1)) ;
    for (k = Km-1; k >= 1; k--) {
	msk[k] = msk[k+1] + exp(kk[k]) / ((double) (l - k + 1)  );
    }

    /* adds the max again */
    for (k = 1; k <= Km; k++) {
	msk[k] = log(msk[k]) + tmp;
    }
    
    
    /* The bayes factor */
    for (k = 2; k <= (Km - 1); k++) {
	baye[k-2] = msk[k] - mik[k-1];
    }

    /* ... and back to R */
    k = 0;
    for (j = 1; j<=Km; j++) {
	mkr[k] = mk[l][j];
	k++;
    }
    


    /* Free memory */
    freevec(msk);
    freevec(mik);
    freevec(kk);
    freetab(Pid);
    freetab(mk);
}



/* *********************************************************************
 *                                                                     *
 *                   partitioning of a trajectory                      *
 *                                                                     *
 ***********************************************************************/

void partraj(double **Pid, int *maxk, double **Mk, double **Mkd, 
	     double **res)
{
    /* declaration of variables */
    int i, j, k, l, D, Km;
    double **Mkk, **cumPid, tmp;
    
    /* Memory allocation */
    l = Pid[0][0]; /* length of the sequence */
    D = Pid[1][0]; /* Number of models */
    Km = *maxk;    /* Partition size */
    tmp = 0;
    
    taballoc(&Mkk, Km, D); /* Contains mkd for i = k, for all models */
    taballoc(&cumPid, l, D); /* For the probability of the sequences */
    
    
    /* First Compute the prediction of the sequences */
    for (j = 1; j <= D; j++) {
	cumPid[1][j] = Pid[1][j];
    }
    for (i = 2; i <= l; i++) {
	for (j = 1; j <= D; j++) {
	    cumPid[i][j] = cumPid[i-1][j] + Pid[i][j];
	    /* fills res */
	    res[i][j] = cumPid[i][j];
	}
    }
    
    
    /* Computes M1 */    
    for (i = 1; i <= l; i++) {
	Mk[i][1] = cumPid[i][1];
	for (j = 2; j <= D; j++) {
	    if (Mk[i][1] < cumPid[i][j]) {
		Mk[i][1] = cumPid[i][j];
	    }
	}
    }
    
    /* Then, computes Mkk for all k <= Km */
    for (j = 1; j <= D; j++) {
	Mkk[1][j] = Pid[1][j];
    }
    
    for (k = 2; k <= Km; k++) {
	/* Computes Mkk */
	for (j = 1; j <= D; j++) {
	    Mkk[k][j] = Pid[k][j] + Mk[k-1][k-1];
	}
	/* Update Mk */
	Mk[k][k] = Mkk[k][1];
	for (j = 2; j <= D; j++) {
	    if (Mkk[k][j] - Mk[k][k] < 0.000000000001)
		Mk[k][k] = Mkk[k][j];
	}
    }
    
    
    /* Now, compute Mkd */
    for (k = 2; k <= Km; k++) {
	
	/* Deletes for i < k */
	for (i = 1; i < k; i++) {
	    for (j = 1; j <= D; j++) {
		Mkd[i][j] = 0;
	    }
	}
		
	/* first fill the first line of Mkd */
	for (j = 1; j <= D; j++) {
	    Mkd[k][j] = Mkk[k][j];
	}

	/* Computes Mkd */
	for (i = (k+1); i <= l; i++) {
	    for (j = 1; j <= D; j++) {
		tmp = Mk[i-1][k-1];
		if (Mkd[i-1][j] - tmp > 0.000000000000000001) {
		    tmp = Mkd[i-1][j];
		}
		Mkd[i][j] = Pid[i][j] + tmp;
	    }
	    
	    /* Update Mk */
	    Mk[i][k] = Mkd[i][1];
	    for (j = 1; j <= D; j++) {
		if (Mkd[i][j] - Mk[i][k] > 0.000000000000000001) {
		    Mk[i][k] = Mkd[i][j];
		}
	    }
	    
	}
	
	/* Fills res */
	for (i = 1; i <= l; i++) {
	    for (j = 1; j <= D; j++) {
		res[((k-1) * l)+i][j] = Mkd[i][j];
	    }
	}
    }

    /* free memory */
    freetab(Mkk);
    freetab(cumPid);
}



void partrajr(double *Pidr, double *curmar, int *curmodr, int *curlocr, 
	      int *lr, int *Dr, int *Kmr)
{
    /* Variable declaration */
    int l, D, Km, i, j, k, m, n, *curloc, *curmod;
    double **Mk, **Mkd, **Pid, **res, *curma, **grap;
    
    /* Memory allocation */
    l = *lr;
    D = *Dr;
    Km = *Kmr;
    m = 0;
    n = 0;
    taballoc(&Mk, l, Km);
    taballoc(&Mkd, l, D);
    taballoc(&Pid, l, D);
    taballoc(&res, (l * Km), D);
    taballoc(&grap, l, D);
    vecalloc(&curma, Km);
    vecintalloc(&curmod, Km);
    vecintalloc(&curloc, (Km+1));
    
    
    /* R to C */
    k = 0;
    for (i = 1; i <= l; i++) {
	for (j = 1; j <= D; j++) {
	    Pid[i][j] = log(Pidr[k]);
	    k++;
	}
    }
    
    /* The main algorithm */
    partraj(Pid, Kmr, Mk, Mkd, res);
    
    
    /* Backtracking */
    curloc[1] = l;
    for (k = Km; k>=1; k--) {
	if (k > 1) {
	    
	    m = Km - k + 2;
	    
            /* Stores the graph */
	    for (i = 1; i <= l; i++) {
		for (j = 1; j <= D; j++) {
		    if (res[(l * (k-1)) + i][j] > Mk[i][k-1]) {
			grap[i][j] = 1;
		    } else {
			grap[i][j] = 0;
		    }
		}
	    }
	    
	    /* The best model for the last step */
	    n = (l * (k-1)) + curloc[m-1];
	    curma[m-1] = res[n][1];
	    curmod[m-1] = 1;
	    for (j = 2; j <= D; j++) {
		if (res[n][j] > curma[m-1]) {
		    curma[m-1] = res[n][j];
		    curmod[m-1] = j;
		}
	    }
	    
	    /* Keep this model until ? */
	    n = curloc[m-1];
	    j = curmod[m-1];
	    while (grap[n][j] > 0.0000000001)
		n--;
	    curloc[m] = n;
	} else {
	    m = Km+1;
	    curloc[m] = 1;
	    curma[m-1] = res[curloc[m-1]][1];
	    curmod[m-1] = 1;
	    for (j = 1; j <= D; j++) {
		if (res[curloc[m-1]][j] > curma[m-1]) {
		    curma[m-1] = res[curloc[m-1]][j];
		    curmod[m-1] = j;
		}
	    }
	}
	
	
    }
    
    
    
    /* C to R */
    for (i = 1; i <= Km; i++) {
	curmodr[i-1] = curmod[i];
	curlocr[i-1] = curloc[i];
	curmar[i-1] = curma[i];
    }
    curlocr[Km] = curloc[Km+1];
    
    
    /* free memory */
    freetab(Mk);
    freetab(Mkd);
    freetab(res);
    freetab(grap);
    freevec(curma);
    freeintvec(curmod);
    freeintvec(curloc);
}






/* acfdist and acfang */


void acfdist (double *sim, double *di, int *ndi, double *diper, int *nbobs, int *indxnona, int *nsim, int *lag){
  
  int i,j, *nona, *numero;
  double nobslag;
  vecintalloc(&nona, *nbobs);
  vecintalloc(&numero, *nbobs);
  for (i=1; i<=*nbobs; i++) {
    nona[i] = indxnona[i-1]-1;
  }
  nobslag=0;
  /* compute observed statistic */
  for (i=1;i<=(*nbobs);i++){
    if((nona[i]+*lag)<=(*ndi-1)){
      if(!ISNAN(di[nona[i] + *lag])) {
	sim[0]=sim[0]+pow(di[nona[i]+*lag] - di[nona[i]],2);
	nobslag=nobslag+1; /* number of observed differences */

      }
    }
  }
 
  
  /* permute only observed values */
  for (j=1; j<=*nsim; j++) {
    getpermutation(numero,j);
    for (i=1;i<=(*nbobs);i++){
      diper[nona[i]]=di[nona[numero[i]]];
    }
    for (i=1;i<=(*nbobs);i++){
      if((nona[i]+*lag)<=(*ndi-1)){
	if(!ISNAN(di[nona[i] + *lag])) {
	  sim[j]=sim[j]+pow(diper[nona[i]+*lag] - diper[nona[i]],2);
	}
      }
    }
  }

  for (j=0; j<=*nsim; j++) {
    sim[j]=sim[j]/nobslag; /* rescale the stat by dividing by the number of observed differences */
  }
  freeintvec(nona);
  freeintvec(numero);
}

void acfangl (double *sim, double *ang, int *nang, double *angper, int *nbobs, int *indxnona, int *nsim, int *lag){
  
  int i,j, *nona, *numero;
  double nobslag;  
  vecintalloc(&nona, *nbobs);
  vecintalloc(&numero, *nbobs);
  for (i=1; i<=*nbobs; i++) {
    nona[i] = indxnona[i-1]-1;
  }
  nobslag=0;
  /* compute observed statistic */
  for (i=1;i<=(*nbobs);i++){
    if((nona[i]+*lag)<=(*nang-1)){
      if(!ISNAN(ang[nona[i] + *lag])) {
	sim[0]=sim[0]+2.0*(1.0-cos(ang[nona[i]+*lag] - ang[nona[i]]));
	nobslag=nobslag+1; /* number of observed differences */
      }
    }
  }
  /* permute only observed values */
  for (j=1; j<=*nsim; j++) {
    getpermutation(numero,j);
    for (i=1;i<=(*nbobs);i++){
      angper[nona[i]]=ang[nona[numero[i]]];
    }
    for (i=1;i<=(*nbobs);i++){
      if((nona[i]+*lag)<=(*nang-1)){
	if(!ISNAN(ang[nona[i] + *lag])) {
	  sim[j]=sim[j]+2.0*(1.0-cos(angper[nona[i]+*lag] - angper[nona[i]]));
	}
      }
    }
  }
  for (j=0; j<=*nsim; j++) {
    sim[j]=sim[j]/nobslag; /* rescale the stat by dividing by the number of observed differences */
  }
  freeintvec(nona);
  freeintvec(numero);

}

/* ***********************************************************************
 *                                                                       *
 *                          MKDE                                         *
 *                                                                       *
 * ********************************************************************* */



int HBT(double xt, double yt, SEXP hab, SEXP nrow, SEXP cs, double xll2, 
	double yll2)
{
    int hh, nl, nc;
    nl = (int) ftrunc(((xt - xll2)/REAL(cs)[0]) + REAL(cs)[0]*0.000001); 
    nc = (int) ftrunc(((yt - yll2)/REAL(cs)[0]) + REAL(cs)[0]*0.000001);
    hh = INTEGER(hab)[ nl + (nc * (INTEGER(nrow)[0])) ];
    return(hh);
}


int HBTl(SEXP xl, SEXP yl, SEXP PAtmp, SEXP hab, SEXP nrow, SEXP cs, double xll2, 
	double yll2, int k, int i)
{
    double t1, xt, yt;
    int n, j, hh, th;
    SEXP habp;
    
    PROTECT(habp = allocVector(INTSXP, k+1));
    

    t1 = REAL(PAtmp)[i+1] - REAL(PAtmp)[i];
    n = (int) round(t1);
    if (n <1)
	n = 1;
    for (j = 0; j < k+1; j++) {
	INTEGER(habp)[j] = 0;
    }
    
    /* identify the habitat at each step */
    for (j = 0; j <= n; j++) {
	xt = REAL(xl)[i] + (((double) j)/((double) n)) * (REAL(xl)[i+1] - REAL(xl)[i]);
	yt = REAL(yl)[i] + (((double) j)/((double) n)) * (REAL(yl)[i+1] - REAL(yl)[i]);
	hh = HBT(xt, yt, hab, nrow, cs, xll2, yll2);
	if (hh != NA_INTEGER) {
	    INTEGER(habp)[hh]++;
	} else {
	    INTEGER(habp)[k]++;
	}
    }
    hh=0;
    for (j = 0; j < k+1; j ++) {
	if (INTEGER(habp)[j] == (n+1)) {
	    th = j;
	    hh++;
	}
    }
    
    if (hh > 0) {
	UNPROTECT(1);
	return(th);
    } else {
	UNPROTECT(1);
	return(NA_INTEGER);
    }
}


/* df contient x,y,date en posix */
SEXP fillsegments(SEXP df, SEXP Tmaxr, SEXP taur, SEXP hminr, SEXP D, SEXP Lminr, 
		  SEXP b, SEXP hab, SEXP xll, SEXP yll, SEXP cs, SEXP nrowc, SEXP PA)
{
    int nrow, nnr, ni, i, m, k, nh, h, lp;
    double dt, dta, Tmax, tau, h2min, Lmin, dist, hmin, xll2, yll2, hm;
    SEXP x, y, date, resux, resuy, resuh, dfso, hmax, h2max, PAtmp, PA2;
    
    /* Le nombre de lignes de ce data.frame */
    nrow = length(VECTOR_ELT(df,0));
    nnr = 0;
    nh = length(D);
    if (nh > 1) {
	xll2 = REAL(xll)[0] - REAL(cs)[0]/2.0;
	yll2 = REAL(yll)[0] - REAL(cs)[0]/2.0;
    }
    
    Tmax = REAL(Tmaxr)[0];
    hmin = REAL(hminr)[0];
    h2min = R_pow(hmin, 2.0);
    
    PROTECT(hmax = allocVector(REALSXP, nh+1));
    PROTECT(h2max = allocVector(REALSXP, nh+1));

    REAL(hmax)[0] = sqrt(((1.0 - REAL(b)[0]) * h2min) + (REAL(D)[0] * Tmax / 2.0));
    hm = REAL(hmax)[0];
    if (nh > 1) {
	for (i = 1; i < nh; i++) {
	    REAL(hmax)[i] = sqrt(((1.0 - REAL(b)[0]) * h2min) + (REAL(D)[i] * Tmax / 2.0));
	    if (hm < REAL(hmax)[i])
		hm = REAL(hmax)[i];
	}
	REAL(hmax)[nh] = hm;
    }
    Lmin = REAL(Lminr)[0];
    tau = REAL(taur)[0];
    REAL(h2max)[0] = R_pow(REAL(hmax)[0], 2.0);
    if (nh > 1) {
	for (i = 1; i <= nh; i++) {
	    REAL(h2max)[i] = R_pow(REAL(hmax)[i], 2.0);
	}
    }

    /* Get the coordinates and date */
    PROTECT(x = coerceVector(VECTOR_ELT(df,0), REALSXP));
    PROTECT(y = coerceVector(VECTOR_ELT(df,1), REALSXP));
    PROTECT(date = coerceVector(VECTOR_ELT(df,2), REALSXP));
    lp = length(PA);
    PROTECT(PAtmp = allocVector(REALSXP, nrow));
    PROTECT(PA2 = coerceVector(PA, REALSXP));
    if (lp > 1) {
	REAL(PAtmp)[0] = 0.0;
	for (i = 1; i < nrow; i++) {
	    REAL(PAtmp)[i] = REAL(PAtmp)[i-1] + (REAL(PA2)[i-1] * (REAL(date)[i] - REAL(date)[i-1]));
	}	
    } else {
	for (i = 0; i < nrow; i++) {
	    REAL(PAtmp)[i] = REAL(date)[i];
	}
    }
    
    
    /* for each segment, calculates the number of points to add */
    for (i = 0; i < (nrow-1); i++) {
	dt = REAL(date)[i+1] - REAL(date)[i];
	dta = REAL(PAtmp)[i+1] - REAL(PAtmp)[i];
	dist = hypot(REAL(x)[i+1] - REAL(x)[i], REAL(y)[i+1] - REAL(y)[i]);
	if ((dt < Tmax)&&(dist > Lmin)) {
	    nnr = nnr + (int) round(dta/tau);
	}
    }
    nnr++;
    
    /* prepares the vector of output */
    PROTECT(resux = allocVector(REALSXP, nnr));
    PROTECT(resuy = allocVector(REALSXP, nnr));
    PROTECT(resuh = allocVector(REALSXP, nnr));
    
    /* and finds the coordinates of the points */
    k=0;
    for (i = 0; i < (nrow-1); i++) {

	dt = REAL(date)[i+1] - REAL(date)[i];
	dta = REAL(PAtmp)[i+1] - REAL(PAtmp)[i];
	dist = hypot(REAL(x)[i+1] - REAL(x)[i], REAL(y)[i+1] - REAL(y)[i]);

	if ((dt < Tmax)&&(dist>Lmin)&&(dta>0.0000001)) {
	    ni = (int) round(dta/tau);
	    for (m = 0; m < ni; m++) {
		REAL(resux)[k] = REAL(x)[i]+ 
			((double) m) * (REAL(x)[i+1] - REAL(x)[i])/((double) ni);
		REAL(resuy)[k] = REAL(y)[i]+ 
		    ((double) m) * (REAL(y)[i+1] - REAL(y)[i])/((double) ni);
		if (nh < 2) {
		    REAL(resuh)[k] = sqrt(h2min + 4.0*(((double) m)/((double) ni))*
					  (1 - ((double) m)/((double) ni))*(REAL(h2max)[0] - h2min)*dta/Tmax);
		} else {
		    h = HBT(REAL(resux)[k], REAL(resuy)[k], hab, nrowc, cs, xll2, 
			    yll2);
		    if (h == NA_INTEGER) {
			REAL(resuh)[k] = sqrt(h2min + 4.0*(((double) m)/((double) ni))*
					      (1 - ((double) m)/((double) ni))*
					      (REAL(h2max)[nh] - h2min)*dta/Tmax);
			
		    } else {
			REAL(resuh)[k] = sqrt(h2min + 4.0*(((double) m)/((double) ni))*
					      (1 - ((double) m)/((double) ni))*
					      (REAL(h2max)[h] - h2min)*dta/Tmax);
			
		    }
		}
		k++;
	    }
	}
    }
    
    REAL(resux)[k] = REAL(x)[nrow-1];
    REAL(resuy)[k] = REAL(y)[nrow-1];
    REAL(resuh)[k] = sqrt(h2min);
    
    PROTECT(dfso = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(dfso, 0, resux);
    SET_VECTOR_ELT(dfso, 1, resuy);
    SET_VECTOR_ELT(dfso, 2, resuh);

    UNPROTECT(11);

    return(dfso);
}


/* On calcule maintenant, sur la base d'une grille passee, l'estimation kernel */
SEXP mkde(SEXP xyh, SEXP grid)
{
    
    int n, nl, i, j;
    SEXP x, y, h, dens, xg, yg, gridso;
    double xmin, ymin, xmax, ymax, hmax, dist;
    
    /* on ajuste alors le noyau */
    n = length(VECTOR_ELT(grid,0));
    nl = length(VECTOR_ELT(xyh,0));


    PROTECT(x = coerceVector(VECTOR_ELT(xyh,0), REALSXP));
    PROTECT(y = coerceVector(VECTOR_ELT(xyh,1), REALSXP));
    PROTECT(h = coerceVector(VECTOR_ELT(xyh,2), REALSXP));
    PROTECT(xg = coerceVector(VECTOR_ELT(grid,0), REALSXP));
    PROTECT(yg = coerceVector(VECTOR_ELT(grid,1), REALSXP));
    PROTECT(dens = allocVector(REALSXP, n));
    
    
    xmin = REAL(x)[0];
    ymin = REAL(y)[0];
    xmax = REAL(x)[0];
    ymax = REAL(y)[0];
    hmax = REAL(h)[0];
    for (j = 1; j < nl; j++) {
	if (REAL(x)[j] < xmin)
	    xmin = REAL(x)[j];
	if (REAL(x)[j] > xmax)
	    xmax = REAL(x)[j];
	if (REAL(y)[j] < ymin)
	    ymin = REAL(y)[j];
	if (REAL(y)[j] > ymax)
	    ymax = REAL(y)[j];
	if (REAL(h)[j] > hmax)
	    hmax = REAL(h)[j];
    }
    hmax = hmax * 3.0;
    
    for (i = 0; i < n; i++) {
	R_CheckUserInterrupt();
	REAL(dens)[i] = 0.0;
	if ((xmin - REAL(xg)[i] < hmax)&&
	    (ymin - REAL(yg)[i] < hmax)&&
	    (REAL(xg)[i] - xmax < hmax)&&
	    (REAL(yg)[i] - ymax < hmax)) {

	    for (j = 0; j < nl; j++) {
		dist= hypot(REAL(x)[j] -REAL(xg)[i], REAL(y)[j] -REAL(yg)[i]);
		if (dist < 3.0*REAL(h)[j]) {
		    REAL(dens)[i] = REAL(dens)[i] + exp(-(R_pow(dist,2.0))/
							(2.0 * R_pow(REAL(h)[j], 2.0))) / 
			R_pow(REAL(h)[j], 2.0);
		}
	    }
	    REAL(dens)[i] = (1.0/(2.0 * M_PI * ((double) nl)))*(REAL(dens)[i]);
	}
    }

    PROTECT(gridso = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(gridso, 0, xg);
    SET_VECTOR_ELT(gridso, 1, yg);
    SET_VECTOR_ELT(gridso, 2, dens);

    UNPROTECT(7);
    return(gridso);
}




/* bis */
SEXP mkdeb(SEXP xyh, SEXP xll, SEXP yll, SEXP cs, SEXP nrow, SEXP ncol)
{
    
    int nl, nr, nc, nro, nco, i, j, l, c, hmaxdis;
    SEXP x, y, h, dens, xg, yg, gridso;
    double xlo, ylo, hmax, dist, xll2, yll2;
    
    /* on ajuste alors le noyau */
    nl = length(VECTOR_ELT(xyh,0));


    PROTECT(x = coerceVector(VECTOR_ELT(xyh,0), REALSXP));
    PROTECT(y = coerceVector(VECTOR_ELT(xyh,1), REALSXP));
    PROTECT(h = coerceVector(VECTOR_ELT(xyh,2), REALSXP));
    PROTECT(xg = allocVector(REALSXP, INTEGER(nrow)[0]*INTEGER(ncol)[0]));
    PROTECT(yg = allocVector(REALSXP, INTEGER(nrow)[0]*INTEGER(ncol)[0]));
    PROTECT(dens = allocVector(REALSXP, INTEGER(nrow)[0]*INTEGER(ncol)[0]));
    nro = INTEGER(nrow)[0];
    nco = INTEGER(ncol)[0];

    for (j = 0; j < nco; j++) {
	for (i = 0; i < nro; i++) {
	    REAL(xg)[ i + (j * (INTEGER(nrow)[0])) ] = REAL(xll)[0] + ((double) i)*REAL(cs)[0];
	    REAL(yg)[ i + (j * (INTEGER(nrow)[0])) ] = REAL(yll)[0] + ((double) j)*REAL(cs)[0];
	}
    }
	    
    for (i = 0; i < INTEGER(nrow)[0]*INTEGER(ncol)[0]; i++) {
	REAL(dens)[i] = 0.0;
    }
    
    hmax = REAL(h)[0];
    for (j = 1; j < nl; j++) {
	if (REAL(h)[j] > hmax)
	    hmax = REAL(h)[j];
    }
    hmax = hmax * 3.0;
    xll2 = REAL(xll)[0] - REAL(cs)[0]/2.0;
    yll2 = REAL(yll)[0] - REAL(cs)[0]/2.0;
    hmaxdis = (int) round(hmax / REAL(cs)[0]);
    
    for (j = 0; j < nl; j++) {
	R_CheckUserInterrupt();
	xlo = REAL(x)[j];
	ylo = REAL(y)[j];
	nr = (int) ftrunc(((xlo - xll2)/REAL(cs)[0]) + REAL(cs)[0]*0.000001); 
	nc = (int) ftrunc(((ylo - yll2)/REAL(cs)[0]) + REAL(cs)[0]*0.000001);
	for (l = (nr-hmaxdis-1); l <(nr+hmaxdis+1); l++) {
	    for (c = (nc-hmaxdis-1); c <(nc+hmaxdis+1); c++) {
		if ((l<nro)&&(l>0)) {
		    if ((c<nco)&&(c>0)) {
			dist= hypot(xlo -REAL(xg)[l + (c * (INTEGER(nrow)[0]))], 
				     ylo -REAL(yg)[l + (c * (INTEGER(nrow)[0]))]);
			REAL(dens)[ l + (c * (INTEGER(nrow)[0])) ] =
			    REAL(dens)[ l + (c * (INTEGER(nrow)[0])) ] +
			    exp(-(R_pow(dist,2.0))/
				(2.0 * R_pow(REAL(h)[j], 2.0))) / 
			    R_pow(REAL(h)[j], 2.0)/(2.0 * M_PI * ((double) nl));
		    }
		}
	    }
	}
    }
	

    PROTECT(gridso = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(gridso, 0, xg);
    SET_VECTOR_ELT(gridso, 1, yg);
    SET_VECTOR_ELT(gridso, 2, dens);

    UNPROTECT(7);
    return(gridso);
}



/* */
SEXP CalculD(SEXP tra, SEXP Tmaxr, SEXP Lmin, SEXP PA)
{
    double Tmax, t1, t2, xt, yt, delta2, D, l1, l2;
    int n, i, Nc, lp;
    SEXP x, y, date, Ds, PA2, PAtmp;
    
    Tmax = REAL(Tmaxr)[0];
    n = length(VECTOR_ELT(tra,0));
    xt = 0.0;
    yt = 0.0;

    /* Get the coordinates and date */
    PROTECT(x = coerceVector(VECTOR_ELT(tra,0), REALSXP));
    PROTECT(y = coerceVector(VECTOR_ELT(tra,1), REALSXP));
    PROTECT(date = coerceVector(VECTOR_ELT(tra,2), REALSXP));
    lp = length(PA);
    PROTECT(PAtmp = allocVector(REALSXP, n));
    PROTECT(PA2 = coerceVector(PA, REALSXP));
    if (lp > 1) {
	REAL(PAtmp)[0] = 0.0;
	for (i = 1; i < n; i++) {
	    REAL(PAtmp)[i] = REAL(PAtmp)[i-1] + (REAL(PA2)[i-1] * 
						 (REAL(date)[i] - REAL(date)[i-1]));
	}	
    } else {
	for (i = 0; i < n; i++) {
	    REAL(PAtmp)[i] = REAL(date)[i];
	}
    }


    D = 0.0;
    Nc = 0;
    for (i = 0; i < (n-2); i++) {	
	t1 = REAL(PAtmp)[i+1] - REAL(PAtmp)[i];
	t2 = REAL(PAtmp)[i+2] - REAL(PAtmp)[i+1];
	l1 = hypot(REAL(x)[i+1] - REAL(x)[i], 
		    REAL(y)[i+1] - REAL(y)[i]);
	l2 = hypot(REAL(x)[i+2] - REAL(x)[i+1],
		    REAL(y)[i+2] - REAL(y)[i+1]);
	if ((REAL(date)[i+2]-REAL(date)[i]) < Tmax) {
	    if (t1 > 0.0000000001) {
		if (t2 > 0.0000000001) {
		    if (t1 < 2.0 * t2) {
			if (t1 > t2/2.0) {
			    if (l1 <  2.0 * l2) {
				if (l1 >  l2/2.0) {
				    if (l1 > REAL(Lmin)[0]) {
					if (l2 > REAL(Lmin)[0]) {
					    xt = REAL(x)[i] + (REAL(x)[i+2] - REAL(x)[i])*(t1/(t1+t2));
					    yt = REAL(y)[i] + (REAL(y)[i+2] - REAL(y)[i])*(t1/(t1+t2));
					    delta2 = R_pow((xt - REAL(x)[i+1]), 2.0) + 
						R_pow((yt - REAL(y)[i+1]), 2.0);
					    D = D + (delta2*((1.0/t1) + (1.0/t2)));
					    Nc++;
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
    D = D/(4.0 * ((double) Nc));
    PROTECT(Ds = allocVector(REALSXP, 1));
    REAL(Ds)[0] = D;
    UNPROTECT(6);

    return(Ds);    
}





/* Composantes d'un trajet. Rasterisation d'un pas */
SEXP RasterPas(SEXP df, SEXP xllr, SEXP yllr, SEXP cs, SEXP type1)
{
    int npas, i, j, nso, k;
    SEXP xl, yl, resu, xso, yso, so, dfso;
    double x1, y1, x2, y2, 
	xc, yc, dist, xt, yt, csi, xll, yll;
    
    npas = length(VECTOR_ELT(df,0)) - 1;
    csi = REAL(cs)[0];
    nso = INTEGER(type1)[0];
    xll = REAL(xllr)[0];
    yll = REAL(yllr)[0];
    
    PROTECT(xl = coerceVector(VECTOR_ELT(df,0), REALSXP));
    PROTECT(yl = coerceVector(VECTOR_ELT(df,1), REALSXP));
    
    
    if (nso > 0) {
	PROTECT(xso = allocVector(REALSXP, nso));
	PROTECT(yso = allocVector(REALSXP, nso));
	PROTECT(so = allocVector(REALSXP, nso));
    }
    
    nso = 0;
    
    for (i = 0; i < npas; i++) {
	
	/* premier filtre */
	x1 = REAL(xl)[i];
	x2 = REAL(xl)[i+1];
	y1 = REAL(yl)[i];
	y2 = REAL(yl)[i+1];
	
	/* identification des coordonnees initiales et finales */	
	dist = hypot(x2-x1, y2-y1);
	k = (int) round(50.0*dist/csi);
	if (k == 0)
	    k++;
	xc = -230876.0;
	yc = -230876.0;
	
	for (j = 0; j<k; j++) {
	    yt = y1 + ((double) j)*(y2-y1)/((double) k);
	    xt = x1 + ((double) j)*(x2-x1)/((double) k);
	    xt = xll + round((xt - xll)/csi)*csi;
	    yt = yll + round((yt - yll)/csi)*csi;
	    if (hypot(yt - yc, xt - xc) > 0.00000000001) {
		xc = xt;
		yc = yt;
		if (INTEGER(type1)[0] != 0) {
		    REAL(xso)[nso] = xc;
		    REAL(yso)[nso] = yc;
		    REAL(so)[nso] = ((double) (i+1));
		} 
		nso++;
	    }
	}
    }
    
    if (INTEGER(type1)[0] != 0) {
	PROTECT(dfso = allocVector(VECSXP, 3));
    	SET_VECTOR_ELT(dfso, 0, xso);
	SET_VECTOR_ELT(dfso, 1, yso);
	SET_VECTOR_ELT(dfso, 2, so);
	UNPROTECT(6);
	
	return(dfso);
    } else {
	PROTECT(resu = allocVector(INTSXP, 1));
	INTEGER(resu)[0] = nso;
	PROTECT(dfso = RasterPas(df, xllr, yllr, cs, resu));
	UNPROTECT(4);
	return(dfso);
    }

}






/* Differences with the program of Simon Benhamou: trunc on a integer i stored as a double
   can return i or i-1
 */
SEXP calculDparhab(SEXP df, SEXP hab, SEXP xll, SEXP yll, SEXP cs, SEXP nrow,
		   SEXP Lmin, SEXP nombrehab, SEXP PA, SEXP Tmax)
{
    SEXP xl, yl, tem, typpas, habp, Nc, Dh, dfso, PAtmp, PA2;
    int i, k, nlocs, lp;
    double t1, t2, l1, l2, xt, yt, delta2, xll2, yll2;

    k = INTEGER(nombrehab)[0];
    nlocs = length(VECTOR_ELT(df,0));


    PROTECT(xl = coerceVector(VECTOR_ELT(df,0), REALSXP));
    PROTECT(yl = coerceVector(VECTOR_ELT(df,1), REALSXP));
    PROTECT(tem = coerceVector(VECTOR_ELT(df,2), REALSXP));
    PROTECT(typpas = allocVector(INTSXP, nlocs-1));
    PROTECT(habp = allocVector(INTSXP, k+1));
    lp = length(PA);
    PROTECT(PAtmp = allocVector(REALSXP, nlocs));
    PROTECT(PA2 = coerceVector(PA, REALSXP));

    /* From center to the corner of lower left pixel */
    xll2 = REAL(xll)[0] - (REAL(cs)[0]/2.0);
    yll2 = REAL(yll)[0] - (REAL(cs)[0]/2.0);

    /*  Take into account the proportion of activity time */
    if (lp > 1) {
	REAL(PAtmp)[0] = 0.0;
	for (i = 1; i < nlocs; i++) {
	    REAL(PAtmp)[i] = REAL(PAtmp)[i-1] + (REAL(PA2)[i-1] * (REAL(tem)[i] - REAL(tem)[i-1]));
	}	
    } else {
	for (i = 0; i < nlocs; i++) {
	    REAL(PAtmp)[i] = REAL(tem)[i];
	}
    }
    
    
    /* for each step */
    for (i = 0; i < nlocs-1; i++) {
	INTEGER(typpas)[i] = HBTl(xl, yl, PAtmp, hab, nrow, cs, xll2, yll2, k, i);
    }
    
    
    /* calculates the D coefficient */
    PROTECT(Nc = allocVector(INTSXP, k));
    PROTECT(Dh = allocVector(REALSXP, k));

    for (i = 0; i < k; i++) {
	REAL(Dh)[i] = 0.0;
	INTEGER(Nc)[i] = 0;
    }
    
    for (i = 0; i < (nlocs-2); i++) { 
	if ((INTEGER(typpas)[i+1] != NA_INTEGER)&&
	    (INTEGER(typpas)[i+1] == INTEGER(typpas)[i])) {
	    l2 = hypot(REAL(xl)[i+2] - REAL(xl)[i+1], REAL(yl)[i+2] - REAL(yl)[i+1]);
	    l1 = hypot(REAL(xl)[i+1] - REAL(xl)[i], REAL(yl)[i+1] - REAL(yl)[i]);
	    t2 = REAL(PAtmp)[i+2] - REAL(PAtmp)[i+1];
	    t1 = REAL(PAtmp)[i+1] - REAL(PAtmp)[i];
	    if (t1 > 0.0000000001) {
		if (t2 > 0.0000000001) {
		    if ((REAL(tem)[i+2]-REAL(tem)[i]) < REAL(Tmax)[0]) {
			if (t1 < 2.0 * t2) {
			    if (t1 > t2/2.0) {
				if (l1 <  2.0 * l2) {
				    if (l1 >  l2/2.0) {
					if (l1 > REAL(Lmin)[0]) {
					    if (l2 > REAL(Lmin)[0]) {
						xt = REAL(xl)[i] + (REAL(xl)[i+2] - 
								    REAL(xl)[i])*(t1/(t1+t2));
						yt = REAL(yl)[i] + (REAL(yl)[i+2] - 
								    REAL(yl)[i])*(t1/(t1+t2));
						delta2 = R_pow((xt - REAL(xl)[i+1]), 2.0) + 
						    R_pow((yt - REAL(yl)[i+1]), 2.0);
						REAL(Dh)[INTEGER(typpas)[i]] = 
						    REAL(Dh)[INTEGER(typpas)[i]] +
						    (delta2*((1.0/t1) + (1.0/t2)));
						INTEGER(Nc)[INTEGER(typpas)[i]]++;
					    }
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
    
    for (i = 0; i < k; i++) {
	REAL(Dh)[i] = 
	    REAL(Dh)[i] / 
	    (4.0 * ((double) INTEGER(Nc)[i]));	    
    }

    PROTECT(dfso = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(dfso, 0, Nc);
    SET_VECTOR_ELT(dfso, 1, Dh);


    UNPROTECT(10);
    return(dfso);    
}



/* **********************************************************


********************************************************** */

/* Calcul D par maximum de vraisemblance */

double calcv(SEXP xl, SEXP yl, SEXP da, double D, SEXP pc, int k)
{
    int n, i;
    double vrais, d, T, t;
    
    n = length(xl);
    vrais = 0.0;
    for (i = 1; i < n-1; i++) {
	if (k == 0) {
	    if (INTEGER(pc)[i] == 1) {
		T = REAL(da)[i+1]-REAL(da)[i-1];
		t = REAL(da)[i]-REAL(da)[i-1];
		d = hypot(REAL(xl)[i] - REAL(xl)[i-1] - (t/T)*(REAL(xl)[i+1] - REAL(xl)[i-1]),
			  REAL(yl)[i] - REAL(yl)[i-1] - (t/T)*(REAL(yl)[i+1] - REAL(yl)[i-1]));
		vrais = vrais + log(T/(4.0*M_PI*D*t*(T-t))) - R_pow(d,2.0)/(4.0*D*t*(T-t)/T);
		k++;
	    } 
	} else {
	    k=0;
	}
    }
    return(vrais);
}


double compteN(SEXP xl, SEXP pc, int k)
{
    int n, i,cons;
    
    n = length(xl);
    cons=0;
    for (i = 1; i < n-1; i++) {
	if (k == 0) {
	    if (INTEGER(pc)[i] == 1) {
		cons++;
		k++;
	    } 
	} else {
	    k=0;
	}
    }
    return((double) cons);
}


SEXP Dmv(SEXP df, SEXP Dr, SEXP pcr, SEXP kr)
{
    SEXP xl, yl, da, D, sor, pc;
    double fx2, fx4, x1, x2, x3, x4, phi;
    int conv, kk;
    

    PROTECT(D = coerceVector(Dr, REALSXP));
    PROTECT(pc = coerceVector(pcr, INTSXP));
    PROTECT(xl = coerceVector(VECTOR_ELT(df,0), REALSXP));
    PROTECT(yl = coerceVector(VECTOR_ELT(df,1), REALSXP));
    PROTECT(da = coerceVector(VECTOR_ELT(df,2), REALSXP));
    PROTECT(sor = allocVector(REALSXP, 2));
    kk = INTEGER(kr)[0];
    
    /* Golden section search */
    x1 = REAL(D)[0];
    x3 = REAL(D)[1];
    phi = (-1.0 + sqrt(5.0))/2.0;


    
    conv = 0;
    while (!conv) {
	x2 = x3 - phi*(x3-x1);
	x4 = x1 + phi*(x3-x1);
	fx2 = calcv(xl, yl, da, x2, pc, kk);
	fx4 = calcv(xl, yl, da, x4, pc, kk);

	if (fx2 < fx4) {
	    x1 = x2;
	} else {
	    x3 = x4;
	}
	if (fabs(x3-x1)<0.00000001) {
	    conv = 1;
	    x4 = (x3+x1)/2.0;
	}
    }
    REAL(sor)[0] = compteN(xl, pc, kk);
    REAL(sor)[1] = x4;


    UNPROTECT(6);
    return(sor);
}





double calcvb(SEXP xl, SEXP yl, SEXP da, double D, SEXP pc, SEXP nb, int k)
{
    int n, i;
    double vrais, d, T, t;
    
    n = length(xl);
    vrais = 0.0;
    for (i = 1; i < n-1; i++) {
	if (k == 0) {
	    if (INTEGER(pc)[i] == 1) {
		if (REAL(nb)[i] > 0.5) {
		    T = REAL(da)[i+1]-REAL(da)[i-1];
		    t = REAL(da)[i]-REAL(da)[i-1];
		    d = hypot(REAL(xl)[i] - REAL(xl)[i-1] - (t/T)*(REAL(xl)[i+1] - REAL(xl)[i-1]),
			      REAL(yl)[i] - REAL(yl)[i-1] - (t/T)*(REAL(yl)[i+1] - REAL(yl)[i-1]));
		    vrais = vrais + REAL(nb)[i] * (log(T/(4.0*M_PI*D*t*(T-t))) - 
						   R_pow(d,2.0)/(4.0*D*t*(T-t)/T));
		    k++;
		}
	    } 
	} else {
	    k=0;
	}
    }
    return(vrais);
}


SEXP Dmvb(SEXP df, SEXP Dr, SEXP pcr, SEXP nbr, SEXP kk)
{
    SEXP xl, yl, da, D, sor, pc, nb;
    double fx2, fx4, x1, x2, x3, x4, phi;
    int conv, k;
    

    PROTECT(D = coerceVector(Dr, REALSXP));
    PROTECT(pc = coerceVector(pcr, INTSXP));
    PROTECT(xl = coerceVector(VECTOR_ELT(df,0), REALSXP));
    PROTECT(yl = coerceVector(VECTOR_ELT(df,1), REALSXP));
    PROTECT(da = coerceVector(VECTOR_ELT(df,2), REALSXP));
    PROTECT(sor = allocVector(REALSXP, 2));
    PROTECT(nb = coerceVector(nbr, REALSXP));
    k = INTEGER(kk)[0];
    
    /* Golden section search */
    x1 = REAL(D)[0];
    x3 = REAL(D)[1];
    phi = (-1.0 + sqrt(5.0))/2.0;


    
    conv = 0;
    while (!conv) {
	x2 = x3 - phi*(x3-x1);
	x4 = x1 + phi*(x3-x1);
	fx2 = calcvb(xl, yl, da, x2, pc, nb, k);
	fx4 = calcvb(xl, yl, da, x4, pc, nb, k);

	if (fx2 < fx4) {
	    x1 = x2;
	} else {
	    x3 = x4;
	}
	if (fabs(x3-x1)<0.00000001) {
	    conv = 1;
	    x4 = (x3+x1)/2.0;
	}
    }
    REAL(sor)[0] = compteN(xl, pc, k);
    REAL(sor)[1] = x4;


    UNPROTECT(7);
    return(sor);
}


/* *************************************************** */

SEXP contrastM(SEXP serie, SEXP Lminr, SEXP type, SEXP ldr)
{
    int n, nd, i, j, Lmin, typei, ld, lmind;
    double *matr, moy, *serier, *x, *x2, *xi, *x2i;
    SEXP mat, seriec, Lminc, typec, ldc, xr, x2r, xir, x2ir;

    /* parameters */
    PROTECT(seriec = coerceVector(serie, REALSXP));
    PROTECT(Lminc = coerceVector(Lminr, INTSXP));
    PROTECT(typec = coerceVector(type, INTSXP));
    PROTECT(ldc = coerceVector(ldr, INTSXP));
    ld = INTEGER(ldc)[0];
    Lmin = INTEGER(Lminc)[0];
    typei = INTEGER(typec)[0];
    
    /* calculation of the grid */
    n = length(serie);
    nd = n;
    n = (int)  ((double) nd)/((double) (ld)); /* similar to floor for positive integer. Truncates the fractional part */
    if (Lmin < ld)
	error("'Lmin' should be greater than ld");
    if ((Lmin % ld) != 0)
	error("'Lmin' should be a multiple of ld");
    lmind = Lmin / ld;
    
    /* contrast matrix calculated for each point of the grid */
    PROTECT(mat = allocMatrix(REALSXP, n, n));
    serier = REAL(seriec);
    matr = REAL(mat);
    
    /* fill the contrast matrix with large values */
    for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++) {
	    matr[i+j*n] = 1000000000000000;
	}
    }    
    
    switch (typei) {
	
    case 1: /* Change in the mean */
	/* calculation of the new x and x2 */
	PROTECT(xr = allocVector(REALSXP, n));
	PROTECT(x2r = allocVector(REALSXP, n));
	PROTECT(xir = allocVector(REALSXP, n));
	PROTECT(x2ir = allocVector(REALSXP, n));
	x = REAL(xr);
	x2 = REAL(x2r);
	xi = REAL(xir);
	x2i = REAL(x2ir);
	
	/* sum of the series and squared series */
	for (i = 0; i < n; i++) {
	    x[i] = 0.0;
	    x2[i] = 0.0;
	    for (j = 0; j < ld; j++) {
		x[i] = x[i] + serier[(ld*i) + j];
		x2[i] = x2[i] + R_pow(serier[(ld*i) + j], 2.0);
	    }
	}
	
	/* cumsum of x2 and x */
	xi[0] = x[0];
	x2i[0] = x2[0];
	for (i = 1; i < n; i++) {
	    xi[i] = xi[i-1] + x[i];
	    x2i[i] = x2i[i-1] + x2[i];
	}
	
	/* fill the first row of the contrast matrix */
	for (i = (lmind-1); i < n; i++) {
	    matr[0 + i*n] = x2i[i] - (R_pow(xi[i], 2.0)/((double) ((i+1)*ld)));
	}
	
	/* fill the rest of the matrix */
	for (i = 1; i < (n - lmind + 1); i++) {
	    
	    for (j = 0; j < n; j++) {
		x2i[j] = x2i[j]-x2[i-1];
		xi[j] = xi[j]-x[i-1];
	    }
	    
	    for (j = (i + lmind - 1); j < n; j++) {
		matr[i + j*n] = x2i[j] - (R_pow(xi[j], 2.0)/((double) ( (j-i+1)*ld)));
	    }
	}
	UNPROTECT(4);
	
	break;


    case 2: /* Change in the variance */
	moy = 0.0;
	for (i = 0; i < nd; i++) {
	    moy = moy + serier[i];
	}
	moy = moy /((double) nd);

	/* calculation of the new x and x2 */
	PROTECT(x2r = allocVector(REALSXP, n));
	PROTECT(x2ir = allocVector(REALSXP, n));
	x2 = REAL(x2r);
	x2i = REAL(x2ir);
	
	/* sum of the squared series */
	for (i = 0; i < n; i++) {
	    x2[i] = 0;
	    for (j = 0; j < ld; j++) {
		x2[i] = x2[i] + R_pow(serier[(ld*i) + j] - moy, 2.0);
	    }
	}
	
	/* cumsum of x2 and x */
	x2i[0] = x2[0];
	for (i = 1; i < n; i++) {
	    x2i[i] = x2i[i-1] + x2[i];
	}
	
	/* fill the first row of the contrast matrix */
	for (i = (lmind-1); i < n; i++) {
	    matr[0 + i * n] = ((double) ((i+1)*ld)) * log(x2i[i]/((double) ((i+1)*ld)));
	}
	
	/* fill the rest of the matrix */
	for (i = 1; i < (n - lmind + 1); i++) {
	    
	    for (j = 0; j < n; j++) {
		x2i[j] = x2i[j]-x2[i-1];
	    }
	    
	    for (j = (i + lmind - 1); j < n; j++) {
		matr[i + j*n] = ((double) ( (j-i+1)*ld)) * log(x2i[j]/((double) ( (j-i+1)*ld)));
	    }
	}
	UNPROTECT(2);
	break;
	
	
	
    case 3: /* Change in the mean and variance */	
	/* calculation of the new x and x2 */
	PROTECT(xr = allocVector(REALSXP, n));
	PROTECT(x2r = allocVector(REALSXP, n));
	PROTECT(xir = allocVector(REALSXP, n));
	PROTECT(x2ir = allocVector(REALSXP, n));
	x = REAL(xr);
	x2 = REAL(x2r);
	xi = REAL(xir);
	x2i = REAL(x2ir);
	
	/* sum of the series and squared series */
	for (i = 0; i < n; i++) {
	    x[i] = 0.0;
	    x2[i] = 0.0;
	    for (j = 0; j < ld; j++) {
		x[i] = x[i] + serier[(ld*i) + j];
		x2[i] = x2[i] + R_pow(serier[(ld*i) + j], 2.0);
	    }
	}
	
	/* cumsum of x2 and x */
	xi[0] = x[0];
	x2i[0] = x2[0];
	for (i = 1; i < n; i++) {
	    xi[i] = xi[i-1] + x[i];
	    x2i[i] = x2i[i-1] + x2[i];
	}

	/* fill the first row of the contrast matrix */
	for (i = (lmind-1); i < n; i++) {
	    matr[0 + i*n] = ((double) ((i+1)*ld)) * 
		log((x2i[i] - (R_pow(xi[i], 2.0)/((double) ((i+1)*ld))))/((double) ((i+1)*ld)));
	}

	/* fill the rest of the matrix */
	for (i = 1; i < (n - lmind + 1); i++) {
	    
	    for (j = 0; j < n; j++) {
		x2i[j] = x2i[j]-x2[i-1];
		xi[j] = xi[j]-x[i-1];
	    }
	    
	    for (j = (i + lmind - 1); j < n; j++) {
		matr[i + j*n] = ((double) ((j-i+1)*ld)) * 
		    log((x2i[i] - (R_pow(xi[i], 2.0)/((double) ((j-i+1)*ld))))/((double) ((j-i+1)*ld)));
	    }
	}
	UNPROTECT(4);
	
	break;

    default:
	Rprintf("No value passed for type\n");
	break;

    }
    UNPROTECT(5);
    return(mat);
}


SEXP dynprog(SEXP mat, SEXP Kmaxr)
{
    int n2, n, i, j, k, L, Kmax, wmi, *ti;
    double mi, tmp, *matr, *Ir;
    SEXP I, t, Kmaxi, so;
    
    n2 = length(mat);
    n = sqrt(n2);
    PROTECT(Kmaxi = coerceVector(Kmaxr, INTSXP));
    Kmax = INTEGER(Kmaxi)[0];
    PROTECT(I = allocMatrix(REALSXP, Kmax, n));
    PROTECT(t = allocMatrix(INTSXP, Kmax, n));
    n2 = n * Kmax;
    Ir = REAL(I);
    matr = REAL(mat);
    ti = INTEGER(t);
    mi = 0.0;
    wmi = 0.0;
    
    for (i = 0; i < n2; i++) {
	Ir[i] = 1000000000000000;
	ti[i] = 0;
    }
    for (i = 0; i < n; i++) {
	Ir[Kmax*i] = matr[n*i];
    }
    tmp = 0.0;
    
    if (Kmax > 2) {
	for (k = 2; k <= Kmax-1; k++) {
	    for (L = k; L <= n; L++) {
		for (j = 1; j <= L-1; j++) {
		    tmp = Ir[k-2 + Kmax * (j-1)]+matr[j+n*(L-1)];
		    if (j==1) {
			mi = tmp;
			wmi = j;
		    } else {
			if (tmp < mi) {
			    mi = tmp;
			    wmi = j;
			}
		    }
		}
		Ir[k-1 + Kmax*(L-1)] = mi;
		ti[k-1 + Kmax*(L-1)] = wmi;
	    }
	}
    }
    
    PROTECT(so = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(so, 0, I);
    SET_VECTOR_ELT(so, 1, t);
    UNPROTECT(4);
    return(so);
}



SEXP findpathc(SEXP matr, SEXP Kr, SEXP Kmax)
{
    SEXP so, Kc, matc;
    int K, i, *sor, n2, n, *mat, Km;
    
    /* size of the matrix */
    n2 = length(matr);
    Km = INTEGER(Kmax)[0]+1;
    n = n2/Km;

    /* Coercion of the arguments -> integer */
    PROTECT(Kc = coerceVector(Kr, INTSXP));
    PROTECT(matc = coerceVector(matr, INTSXP));

    /* Usefull for further analysis */
    mat = INTEGER(matc);
    K = INTEGER(Kc)[0];
    
    /* Output vector */
    PROTECT(so = allocVector(INTSXP, K));
    sor = INTEGER(so);

    /* On rajoute 1 pour la sortie sous R */
    sor[0] = mat[K-1+Km*(n-1)]+1;
    for (i = 1; i < K; i++) {
	sor[i] = mat[K-1-i + Km*(sor[i-1]-1)]+1;
    }
    
    UNPROTECT(3);
    return(so);
}


/* Calculation of the residence time */


SEXP residtime(SEXP xyt, SEXP distr, SEXP maxt)
{
    /* declaring the variables */
    int n,i, j, *deds, sortie;
    double *xr, *yr, *tr, dist, maxtr, *resur, bti, fti, limitr, refti, a, u, v, p, lrb, lrf;
    SEXP x, y, t, dedsr, resu;
    
    /* coercing the arguments */
    PROTECT(x = coerceVector(VECTOR_ELT(xyt,0), REALSXP));
    PROTECT(y = coerceVector(VECTOR_ELT(xyt,1), REALSXP));
    PROTECT(t = coerceVector(VECTOR_ELT(xyt,2), REALSXP));

    n = length(x); /* number of relocations */
    PROTECT(dedsr = allocVector(INTSXP, n)); /* will be used to check whether the relocation j
						is within the distance distr of the relocation i */
    PROTECT(resu = allocVector(REALSXP, n)); /* the output vector */


    /* for the ease of manipulation: gets the pointers */
    resur=REAL(resu);
    xr = REAL(x);
    yr = REAL(y);
    tr = REAL(t);
    deds = INTEGER(dedsr);

    /* get the two constants passed as arguments */
    maxtr = REAL(maxt)[0];
    dist = REAL(distr)[0];

    /* Now, calculate the residence time for each relocation */
    for (i = 0; i < n; i++) {
	
	/* checks which relocations are within the distance dist from reloc j */
	for (j = 0; j < n; j++) {
	    if (hypot(xr[i]-xr[j], yr[i]-yr[j])<=dist) {
		deds[j] = 1;
	    } else {
		deds[j] = 0;
	    }
	}
	
	/* calculates the backward time */
	sortie = 0; /* = 0 when the animal is still inside the circle; =1 otherwise */
	limitr = -5.0; /* used to store the time point when the animal goes out of the circle.
			  Negative if the animal never goes out
			*/
	refti = tr[i]; /* reference time */
	bti = 0.0; /* The backward time */
	
	
	/* if this is not the first relocation (if it is, bti = 0) */
	if (i != 0) {
	    
	    /* for all previous relocations (because backward) */
	    for (j = i-1; j>=0; j--) {
		
		/* if the relocation is outside the circle */
		if (deds[j]==0) {
		    
		    /* if this is the first relocation outside */
		    if (sortie==0) {
			
			/* interpolating the time when the animal came in
			   of the circle */
			a = atan2(yr[j]-yr[j+1], xr[j]-xr[j+1]);
			u = ((xr[i]-xr[j+1])*cos(a))+((yr[i]-yr[j+1])*sin(a));
			v = ((yr[i]-yr[j+1])*cos(a))-((xr[i]-xr[j+1])*sin(a));
			p = (sqrt(R_pow(dist, 2.0) - R_pow(v, 2.0))-fabs(u))/
			    hypot(xr[j]-xr[j+1], yr[j]-yr[j+1]);
			limitr = tr[j+1] - p*(tr[j+1]-tr[j]);
			bti = bti + fabs(refti- limitr);
			sortie = 1;
		    } else {
			/* checks whether the
			   time is too long outside the circle. In this case,
			   break the loop. */
			if (fabs(limitr - tr[j]) > maxtr) {
			    break;
			}
		    }
		} else {
		    if (sortie>0) {
			/* interpolating the time when the animal came out of
			   the circle */
			a = atan2(yr[j+1]-yr[j], xr[j+1]-xr[j]);
			u = ((xr[i]-xr[j])*cos(a))+((yr[i]-yr[j])*sin(a));
			v = ((yr[i]-yr[j])*cos(a))-((xr[i]-xr[j])*sin(a));
			p = (sqrt(R_pow(dist, 2.0) - R_pow(v, 2.0))-fabs(u))/
			    hypot(xr[j]-xr[j+1], yr[j]-yr[j+1]);
			refti = tr[j] + p*(tr[j+1]-tr[j]);
			if (fabs(refti - limitr)> maxtr)
			    break;
			bti = bti + fabs(tr[j] - refti);
			refti = tr[j];
			sortie = 0;
		    } else {
			bti = bti + fabs(refti - tr[j]);
			refti = tr[j];
		    }
		}
	    }
	}

	/* stores the value of limitr (whether the animal came out of the circle */
	lrb = limitr;


	/* forward time */
	sortie = 0;
	limitr = -5.0;
	refti = tr[i];
	fti = 0.0; /* forward time */
	
	/* if this is not the last relocation (otherwise, fti = 0.0) */
	if (i < (n-1)) {
	    
	    /* for all next relocations */
	    for (j = i+1; j<n; j++) {
		
		/* if the relocation is outside the circle */
		if (deds[j]==0) {
		    
		    /* if this is the first relocation outside */
		    if (sortie==0) {
			/* interpolating the time when the animal come out
			   of the circle
			 */
			a = atan2(yr[j]-yr[j-1], xr[j]-xr[j-1]);
			u = ((xr[i]-xr[j-1])*cos(a))+((yr[i]-yr[j-1])*sin(a));
			v = ((yr[i]-yr[j-1])*cos(a))-((xr[i]-xr[j-1])*sin(a));
			p = (sqrt(R_pow(dist, 2.0) - R_pow(v, 2.0))-fabs(u))/
			    hypot(xr[j]-xr[j-1], yr[j]-yr[j-1]);
			limitr = tr[j-1] + p*(tr[j]-tr[j-1]);
			fti = fti + fabs(limitr - refti);
			sortie = 1;
		    } else {
			/* if it is not the first relocation: checks whether the
			   time is too long outside the circle. In this case,
			   break the loop. */
			if (fabs(tr[j] - limitr) > maxtr) {
			    break;
			}
		    }
		} else {
		    if (sortie>0) {
			/* interpolating the time when the animal came in
			   the circle */
			a = atan2(yr[j-1]-yr[j], xr[j-1]-xr[j]);
			u = ((xr[i]-xr[j])*cos(a))+((yr[i]-yr[j])*sin(a));
			v = ((yr[i]-yr[j])*cos(a))-((xr[i]-xr[j])*sin(a));
			p = (sqrt(R_pow(dist, 2.0) - R_pow(v, 2.0))-fabs(u))/
			    hypot(xr[j]-xr[j-1], yr[j]-yr[j-1]);
			refti = tr[j] - p*(tr[j]-tr[j-1]);
			if (fabs(refti - limitr)> maxtr)
			    break;
			fti = fti + fabs(tr[j] - refti);
			refti = tr[j];
			sortie = 0;
		    } else {
			fti = fti + fabs(tr[j] - refti);
			refti = tr[j];
		    }
		}
	    }
	}
	lrf = limitr;
	resur[i] = bti+fti;
	if ((lrb < 0)||(lrf < 0))
	    resur[i] = NA_REAL;
    }
    
    UNPROTECT(5);
    
    /* output */
    return(resu);
    
}


/* ******************************* */
/* Discretization of a trajectory with constant time lag */

SEXP redistime(SEXP xyt, SEXP tlr, SEXP sam0r)
{
    int n, i, nn, cur;
    SEXP x, y, t, xn, yn, tn, dfso;
    double *xr, *yr, *tr, *xnr, *ynr, *tnr, tp, tl, u, sam0, xp, yp, tlc;
    
    /* coercing the arguments */
    PROTECT(x = coerceVector(VECTOR_ELT(xyt,0), REALSXP));
    PROTECT(y = coerceVector(VECTOR_ELT(xyt,1), REALSXP));
    PROTECT(t = coerceVector(VECTOR_ELT(xyt,2), REALSXP));
    xr = REAL(x);
    yr = REAL(y);
    tr = REAL(t);
    tl = REAL(tlr)[0];
    sam0 = REAL(sam0r)[0];
    n = length(x);
    
    /* New number of relocations */
    nn = (int) (round((tr[n-1]-tr[0])/tl)+2);
    PROTECT(xn = allocVector(REALSXP, nn));
    PROTECT(yn = allocVector(REALSXP, nn));
    PROTECT(tn = allocVector(REALSXP, nn));
    xnr = REAL(xn);
    ynr = REAL(yn);
    tnr = REAL(tn);

    for (i = 0; i < nn; i++) {
	tnr[i] = -10.0;
	xnr[i] = -10.0;
	ynr[i] = -10.0;
    }
    
    /* Sample the first relocation ? */
    if (sam0>0.5) {
	GetRNGstate();
	u = unif_rand();
	PutRNGstate();
	xnr[0] = xr[0]+u*(xr[1]-xr[0]);
	ynr[0] = yr[0]+u*(yr[1]-yr[0]);
	tnr[0] = tr[0]+u*(tr[1]-tr[0]);
    } else {
	xnr[0] = xr[0];
	ynr[0] = yr[0];
	tnr[0] = tr[0];
    }
    cur = 0;
    tp = tnr[0];
    xp = xnr[0];
    yp = ynr[0];
    tlc = tl;
    
    /* */
    for (i = 1; i < n; i++) {
	while ((tr[i] - tp) > tlc) {
	    tnr[cur+1] = tp + tlc;
	    xnr[cur+1] = xp + (tlc/(tr[i] - tp))*(xr[i]-xp);
	    ynr[cur+1] = yp + (tlc/(tr[i] - tp))*(yr[i]-yp);
	    cur++;
	    tp = tnr[cur];
	    xp = xnr[cur];
	    yp = ynr[cur];
	    tlc = tl;
	}
	tlc = tlc - (tr[i]-tp);
	tp = tr[i];
	xp = xr[i];
	yp = yr[i];
    }
    if (cur < nn-1) {
	for (i=(cur+1); i<nn; i++) {
	    tnr[i] = -10;
	    xnr[i] = -10;
	    ynr[i] = -10;
	}
    }
    
    PROTECT(dfso = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(dfso, 0, xn);
    SET_VECTOR_ELT(dfso, 1, yn);
    SET_VECTOR_ELT(dfso, 2, tn);

    UNPROTECT(7);
    return(dfso);
}


/* ************************************************** */

/* Random direction herd */
SEXP tr_randomDirection(SEXP xyt, SEXP par1, SEXP par2, SEXP parcon, 
			SEXP traitement, SEXP constraint)
{
    int n, i, ok;
    SEXP x, y, alpha, xyso, xn, yn, t, resu, env, resucont;
    double *xr, *yr, *alphar, *ynr, *xnr, d2;
    
    /* checks the arguments */
    PROTECT(x = coerceVector(VECTOR_ELT(xyt,0), REALSXP));
    PROTECT(y = coerceVector(VECTOR_ELT(xyt,1), REALSXP));
    PROTECT(t = coerceVector(VECTOR_ELT(xyt,2), REALSXP));
    n = length(x);
    PROTECT(alpha = allocVector(REALSXP, n-1));
    PROTECT(xn = allocVector(REALSXP, n));
    PROTECT(yn = allocVector(REALSXP, n));
    PROTECT(env = VECTOR_ELT(par1,0));
    if(!isEnvironment(env)) error("'env' should be an environment");
    
    /* get the pointers */
    xnr = REAL(xn);
    ynr = REAL(yn);
    xr = REAL(x);
    yr = REAL(y);
    alphar = REAL(alpha);
    ok = 0;
    
    /* while the trajectory does not satisfy the constraints */
    while (ok==0) {

	/* sample random angles */
	GetRNGstate();
	for (i=0; i < (n-1); i++) {
	    alphar[i] = unif_rand()*2*M_PI;
	}
	PutRNGstate();
	
	/* define the first relocation */
	xnr[0] = xr[0];
	ynr[0] = yr[0];
	
	/* and the following */
	for (i=1; i < n; i++) {
	    d2 = hypot(xr[i]-xr[i-1], yr[i]-yr[i-1]);
	    xnr[i] = xnr[i-1] + d2*cos(alphar[i]);
	    ynr[i] = ynr[i-1] + d2*sin(alphar[i]);
	}
	
	/* prepares the data.frame storing the results */
	PROTECT(xyso = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(xyso, 0, xn);
	SET_VECTOR_ELT(xyso, 1, yn);
	SET_VECTOR_ELT(xyso, 2, t);
	
	/* checks whether the constraints are satisfied */
	defineVar(install("x"), xyso, env);
	defineVar(install("par"), parcon, env);
	PROTECT(resucont = coerceVector(eval(constraint, env), INTSXP));
	ok = INTEGER(resucont)[0];
	if (ok!=1) {
	    UNPROTECT(2);
	}
    }
    
    defineVar(install("x"), xyso, env);
    defineVar(install("par"), par2, env);
    PROTECT(resu = eval(traitement, env));
    UNPROTECT(10);
    return(resu);
    
}





SEXP tr_randomRotation(SEXP xyt, SEXP par1, SEXP par2, SEXP parcon,
		       SEXP traitement, SEXP constraint)
{
    int n, i, ok;
    SEXP x, y, xyso, xn, yn, t, resu, env, resucont;
    double *xr, *yr, alpha, alph2, *ynr, *xnr, d2, xm, ym;
        
    PROTECT(x = coerceVector(VECTOR_ELT(xyt,0), REALSXP));
    PROTECT(y = coerceVector(VECTOR_ELT(xyt,1), REALSXP));
    PROTECT(t = coerceVector(VECTOR_ELT(xyt,2), REALSXP));
    PROTECT(env = VECTOR_ELT(par1,0));
    n = length(x);
    ok = 0;
    
    PROTECT(xn = allocVector(REALSXP, n));
    PROTECT(yn = allocVector(REALSXP, n));
    if(!isEnvironment(env)) error("'env' should be an environment");
    
    xnr = REAL(xn);
    ynr = REAL(yn);
    xr = REAL(x);
    yr = REAL(y);
    
    /* random rotation */
    ok = 0;
    while (ok == 0) {
	GetRNGstate();
	alpha = unif_rand()*2*M_PI;
	PutRNGstate();
	
	/* the centroid */
	xm = 0;
	ym = 0;
	for (i = 0; i < n; i++) {
	    xm = xm + xr[i];
	    ym = ym + yr[i];	
	}
	xm = xm / ((double) n);
	ym = ym / ((double) n);
	
	/* rotations */
	for (i=0; i < n; i++) {
	    d2 = hypot(xr[i]-xm, yr[i]-ym);
	    alph2 = atan2(yr[i]-ym, xr[i]-xm);
	    alph2= alph2+alpha;
	    xnr[i] = xm + d2*cos(alph2);
	    ynr[i] = ym + d2*sin(alph2);
	}
	
	PROTECT(xyso = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(xyso, 0, xn);
	SET_VECTOR_ELT(xyso, 1, yn);
	SET_VECTOR_ELT(xyso, 2, t);
	
	/* checks whether the constraints are satisfied */
	defineVar(install("x"), xyso, env);
	defineVar(install("par"), parcon, env);
	
	PROTECT(resucont = coerceVector(eval(constraint, env), INTSXP));
	ok = INTEGER(resucont)[0];
	if (ok!=1) {
	    UNPROTECT(2);
	}
    }
    
    defineVar(install("x"), xyso, env);
    defineVar(install("par"), par2, env);
    PROTECT(resu = eval(traitement, env));
    
    UNPROTECT(9);
    return(resu);
    
}




SEXP tr_randomShiftRotation(SEXP xyt, SEXP par1, SEXP par2, SEXP parcon,
			    SEXP traitement, SEXP constraint)
{
    int n, i, ok;
    SEXP x, y, xyso, xn, yn, t, resu, env, rax, ray, resucont, rshift, rrot;
    SEXP namecol, namerow, classdf;
    double *xr, *yr, alpha, alph2, *ynr, *xnr, d2, xm, ym, xmn, ymn;
    
    /* coercion */
    PROTECT(x = coerceVector(VECTOR_ELT(xyt,0), REALSXP));
    PROTECT(y = coerceVector(VECTOR_ELT(xyt,1), REALSXP));
    PROTECT(t = coerceVector(VECTOR_ELT(xyt,2), REALSXP));
    PROTECT(env = VECTOR_ELT(par1,0));
    if(!isEnvironment(env)) error("'env' should be an environment");
    
    PROTECT(rshift = coerceVector(VECTOR_ELT(par1,1), INTSXP));
    PROTECT(rrot = coerceVector(VECTOR_ELT(par1,2), INTSXP));
    PROTECT(rax = coerceVector(VECTOR_ELT(par1,3), REALSXP));
    PROTECT(ray = coerceVector(VECTOR_ELT(par1,4), REALSXP));
    n = length(x);
    xmn = 0.0;
    ymn = 0.0;
    
    /* new vectors */
    PROTECT(xn = allocVector(REALSXP, n));
    PROTECT(yn = allocVector(REALSXP, n));
    
    /* prepares the attributes of the output data.frame */
    /* row.names */
    PROTECT(namerow = getAttrib(xyt, R_RowNamesSymbol)); 
    /* names */
    PROTECT(namecol = allocVector(STRSXP, 3)); 
    SET_STRING_ELT(namecol, 0, mkChar("x"));
    SET_STRING_ELT(namecol, 1, mkChar("y"));
    SET_STRING_ELT(namecol, 2, mkChar("date"));
    /* class */
    PROTECT(classdf = allocVector(STRSXP, 1));
    SET_STRING_ELT(classdf, 0, mkChar("data.frame"));
    
    
    xnr = REAL(xn);
    ynr = REAL(yn);
    xr = REAL(x);
    yr = REAL(y);
    
    /* random rotation */
    ok = 0;
    while (ok == 0) {
	R_CheckUserInterrupt();
	GetRNGstate();
	if (INTEGER(rrot)[0]>0) {
	    alpha = unif_rand()*2*M_PI;
	} else {
	    alpha = 0;
	}
	if (INTEGER(rshift)[0]>0) {
	    xmn = REAL(rax)[0] + unif_rand()*(REAL(rax)[1]-REAL(rax)[0]);
	    ymn = REAL(ray)[0] + unif_rand()*(REAL(ray)[1]-REAL(ray)[0]);
	}
	PutRNGstate();
	
	/* the centroid */
	xm = 0;
	ym = 0;
	for (i = 0; i < n; i++) {
	    xm = xm + xr[i];
	    ym = ym + yr[i];	
	}
	xm = xm / ((double) n);
	ym = ym / ((double) n);
	if (INTEGER(rshift)[0]==0) {
	    xmn = xm;
	    ymn = ym;
	}
	
	/* rotations */
	for (i=0; i < n; i++) {
	    d2 = hypot(xr[i]-xm, yr[i]-ym);
	    alph2 = atan2(yr[i]-ym, xr[i]-xm);
	    alph2= alph2+alpha;
	    xnr[i] = xmn + d2*cos(alph2);
	    ynr[i] = ymn + d2*sin(alph2);
	}
	
	PROTECT(xyso = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(xyso, 0, xn);
	SET_VECTOR_ELT(xyso, 1, yn);
	SET_VECTOR_ELT(xyso, 2, t);
	
	/* define class and attributes */
	classgets(xyso, classdf);
	setAttrib(xyso, R_NamesSymbol, namecol);
	setAttrib(xyso, R_RowNamesSymbol, namerow);
	
	/* checks whether the constraints are satisfied */
	defineVar(install("x"), xyso, env);
	defineVar(install("par"), parcon, env);

	PROTECT(resucont = coerceVector(eval(constraint, env), INTSXP));
	ok = INTEGER(resucont)[0];
	if (ok!=1) {
	    UNPROTECT(2);
	}
	
    }
    
    defineVar(install("x"), xyso, env);
    defineVar(install("par"), par2, env);
    PROTECT(resu = eval(traitement, env));
    
    UNPROTECT(16);
    return(resu);
    
}





SEXP tr_randomRotationCs(SEXP xyt, SEXP par1, SEXP par2, SEXP parcon,
			 SEXP traitement, SEXP constraint)
{
    int n, i, ok, nani, ncs;
    SEXP x, y, xyso, xn, yn, t, resu, env, resucont, rDistCs;
    SEXP rAngleCs, rCentroidAngle, cs, distances;
    SEXP namecol, namerow, classdf, samCs, rCs;
    double *xr, *yr, alpha, alpha2, alph2, *ynr, *xnr, d2, xm, ym, xmn, ymn, dist, *csr;
    
    /* coercion */
    PROTECT(x = coerceVector(VECTOR_ELT(xyt,0), REALSXP));
    PROTECT(y = coerceVector(VECTOR_ELT(xyt,1), REALSXP));
    PROTECT(t = coerceVector(VECTOR_ELT(xyt,2), REALSXP));
    PROTECT(env = VECTOR_ELT(par1,0));
    if(!isEnvironment(env)) error("'env' should be an environment");

    PROTECT(rDistCs = coerceVector(VECTOR_ELT(par1,1), INTSXP));
    PROTECT(rAngleCs = coerceVector(VECTOR_ELT(par1,2), INTSXP));
    PROTECT(rCentroidAngle = coerceVector(VECTOR_ELT(par1,3), INTSXP));
    PROTECT(cs = coerceVector(VECTOR_ELT(par1,4), REALSXP));

    PROTECT(distances = coerceVector(VECTOR_ELT(par1,5), REALSXP));
    PROTECT(rCs = coerceVector(VECTOR_ELT(par1, 6), INTSXP));
    samCs = VECTOR_ELT(par1, 7);

    ncs = length(samCs);
    n = length(x);
    nani = length(distances);
    xmn = 0.0;
    ymn = 0.0;
    
    /* new vectors */
    PROTECT(xn = allocVector(REALSXP, n));
    PROTECT(yn = allocVector(REALSXP, n));
    
    /* some pointers */
    xnr = REAL(xn);
    ynr = REAL(yn);
    xr = REAL(x);
    yr = REAL(y);
    xmn = 0.0;
    ymn = 0.0;
    
    /* the centroid */
    xm = 0;
    ym = 0;
    for (i = 0; i < n; i++) {
	xm = xm + xr[i];
	ym = ym + yr[i];	
    }
    xm = xm / ((double) n);
    ym = ym / ((double) n);
    

    /* prepares the attributes of the output data.frame */
    /* row.names */
    PROTECT(namerow = getAttrib(xyt, R_RowNamesSymbol)); 
    /* names */
    PROTECT(namecol = allocVector(STRSXP, 3)); 
    SET_STRING_ELT(namecol, 0, mkChar("x"));
    SET_STRING_ELT(namecol, 1, mkChar("y"));
    SET_STRING_ELT(namecol, 2, mkChar("date"));
    /* class */
    PROTECT(classdf = allocVector(STRSXP, 1));
    SET_STRING_ELT(classdf, 0, mkChar("data.frame"));
    
    
    
    ok = 0;
    while (ok == 0) {

	R_CheckUserInterrupt();
	
	/* random choice of the capture site */
	if (INTEGER(rCs)[0]>0) {
	    GetRNGstate();
	    i = (int) floor(unif_rand()*((double) ncs));
	    PutRNGstate();
	    csr = REAL(VECTOR_ELT(samCs, i));
	} else {
	    csr = REAL(cs);
	}

	
	/* random rotation around the centroid */   
	if (INTEGER(rCentroidAngle)[0]>0) {
	    GetRNGstate();
	    alpha = unif_rand()*2*M_PI;
	    PutRNGstate();
	} else {
	    alpha = 0;
	}
	
	/* random rotation around the capture site */
	if (INTEGER(rAngleCs)[0]>0) {
	    GetRNGstate();
	    alpha2 = unif_rand()*2*M_PI;
	    PutRNGstate();
	} else {
	    alpha2 = atan2(ym-csr[1], xm-csr[0]);
	}
	
	/* random choice of a distance to the capture site */
	if (INTEGER(rDistCs)[0]>0) {
	    GetRNGstate();
	    i = (int) floor(unif_rand()*((double) nani));
	    PutRNGstate();
	    dist = REAL(distances)[i];
	} else {
	    dist = hypot(xm-csr[0], ym-csr[1]);
	}
		
	/* update the coordinates of the centroid after randomization */
	xmn = csr[0] + dist*cos(alpha2);
	ymn = csr[1] + dist*sin(alpha2);
	
	/* update the coordinates of the relocations */
	for (i=0; i < n; i++) {
	    d2 = hypot(xr[i]-xm, yr[i]-ym);
	    alph2 = atan2(yr[i]-ym, xr[i]-xm);
	    alph2= alph2+alpha;
	    xnr[i] = xmn + d2*cos(alph2);
	    ynr[i] = ymn + d2*sin(alph2);
	}	
	
	/* build the object */
	PROTECT(xyso = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(xyso, 0, xn);
	SET_VECTOR_ELT(xyso, 1, yn);
	SET_VECTOR_ELT(xyso, 2, t);
	
	/* define class and attributes */
	classgets(xyso, classdf);
	setAttrib(xyso, R_NamesSymbol, namecol);
	setAttrib(xyso, R_RowNamesSymbol, namerow);
	
	/* checks whether the constraints are satisfied */
	defineVar(install("x"), xyso, env);
	defineVar(install("par"), parcon, env);

	PROTECT(resucont = coerceVector(eval(constraint, env), INTSXP));
	ok = INTEGER(resucont)[0];
	if (ok!=1) {
	    UNPROTECT(2);
	}
	
    }
    
    defineVar(install("x"), xyso, env);
    defineVar(install("par"), par2, env);
    PROTECT(resu = eval(traitement, env));
    
    UNPROTECT(18);
    return(resu);
    
}







SEXP tr_randomCRW(SEXP xyt, SEXP par1, SEXP par2, SEXP parcon,
		  SEXP traitement, SEXP constraint)
{
    int n, i, ok, *permr, *permdr, fs;
    SEXP x, y, xyso, xn, yn, t, resu, env, alphar, alphaa, dist, resucont, perm, alea, permd;
    SEXP pa, pd, alead, fixedStart, x0, rx, ry;
    SEXP namecol, namerow, classdf;
    double *xr, *yr, *ynr, *xnr, *alpharr, *alphaar, *distr, *alear, *aleadr;
    double alp, *x0r, *rxr, *ryr;
    
    /* the trajectory */
    PROTECT(x = coerceVector(VECTOR_ELT(xyt,0), REALSXP));
    PROTECT(y = coerceVector(VECTOR_ELT(xyt,1), REALSXP));
    PROTECT(t = coerceVector(VECTOR_ELT(xyt,2), REALSXP));
        
    /* The variables used in this prog */
    PROTECT(env = VECTOR_ELT(par1,0));
    if(!isEnvironment(env)) error("'env' should be an environment");
    
    /* The parameters of the randomization */
    PROTECT(pa = coerceVector(VECTOR_ELT(par1,1), INTSXP));
    PROTECT(pd = coerceVector(VECTOR_ELT(par1,2), INTSXP));
    PROTECT(fixedStart = coerceVector(VECTOR_ELT(par1,3), INTSXP));

    PROTECT(x0 = coerceVector(VECTOR_ELT(par1,4), REALSXP));
    PROTECT(rx = coerceVector(VECTOR_ELT(par1,5), REALSXP));
    PROTECT(ry = coerceVector(VECTOR_ELT(par1,6), REALSXP));

    n = length(x);
    PROTECT(alphar = allocVector(REALSXP,n-2));    
    PROTECT(perm = allocVector(INTSXP,n-2));    
    PROTECT(permd = allocVector(INTSXP,n-1));    

    PROTECT(alea = allocVector(REALSXP,n-2));    
    PROTECT(alead = allocVector(REALSXP,n-1));    
    PROTECT(alphaa = allocVector(REALSXP,n-1));

    PROTECT(dist = allocVector(REALSXP,n-1));    
    PROTECT(xn = allocVector(REALSXP, n));
    PROTECT(yn = allocVector(REALSXP, n));

    /* prepares the attributes of the output data.frame */
    /* row.names */
    PROTECT(namerow = getAttrib(xyt, R_RowNamesSymbol)); 
    /* names */
    PROTECT(namecol = allocVector(STRSXP, 3)); 
    SET_STRING_ELT(namecol, 0, mkChar("x"));
    SET_STRING_ELT(namecol, 1, mkChar("y"));
    SET_STRING_ELT(namecol, 2, mkChar("date"));
    /* class */
    PROTECT(classdf = allocVector(STRSXP, 1));
    SET_STRING_ELT(classdf, 0, mkChar("data.frame"));

    
    /* get the pointers */
    xnr = REAL(xn);
    ynr = REAL(yn);
    xr = REAL(x);
    yr = REAL(y);
    alphaar = REAL(alphaa);
    alpharr = REAL(alphar);
    distr = REAL(dist);
    permr = INTEGER(perm);
    alear = REAL(alea);
    permdr = INTEGER(permd);
    aleadr = REAL(alead);
    x0r=REAL(x0);
    rxr=REAL(rx);
    ryr=REAL(ry);
    fs = INTEGER(fixedStart)[0];

    /* calculate relative angles and distance */
    for (i=1; i<n; i++) {
	alphaar[i-1] = atan2(yr[i]-yr[i-1], xr[i]-xr[i-1]);
	distr[i-1] = hypot(yr[i]-yr[i-1], xr[i]-xr[i-1]);
	if (i>1) {
	    alpharr[i-2] = alphaar[i-1]-alphaar[i-2];
	}
    }
    
    /* random permutation of the relative angles and distance */
    ok = 0;
    while (ok == 0) {
	R_CheckUserInterrupt();
	
	GetRNGstate();
	for (i=0; i<n-1; i++) {
	    if (i<n-2) {
		if (INTEGER(pa)[0]>0) {
		    alear[i] = unif_rand();
		}
		permr[i] = i;
	    }
	    if (INTEGER(pd)[0]>0) {
		aleadr[i] = unif_rand();
	    }
	    permdr[i] = i;
	}
	PutRNGstate();
	if (INTEGER(pa)[0]>0) {
	    rsort_with_index(alear, permr, n-2);
	}
	if (INTEGER(pd)[0]>0) {
	    rsort_with_index(aleadr, permdr, n-1);
	}
	
	/* calculates the new trajectory */
	/* if the starting point should be drawn at random */
	if (fs<1) {
	    GetRNGstate();
	    x0r[0] = rxr[0] + unif_rand()*(rxr[1]-rxr[0]);
	    x0r[1] = ryr[0] + unif_rand()*(ryr[1]-ryr[0]);
	    PutRNGstate();
	}
	xnr[0] = x0r[0];
	ynr[0] = x0r[1];
	xnr[1] = xnr[0] + (xr[1]-xr[0]);
	ynr[1] = ynr[0] + (yr[1]-yr[0]);
	
	for (i = 2; i<n; i++) {
	    alp = atan2(ynr[i-1]-ynr[i-2], xnr[i-1]-xnr[i-2]);
	    xnr[i] = xnr[i-1] + distr[permdr[i-2]]*cos(alp+alpharr[permr[i-2]]);
	    ynr[i] = ynr[i-1] + distr[permdr[i-2]]*sin(alp+alpharr[permr[i-2]]);
	}
	
	/* stores the results */
	PROTECT(xyso = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(xyso, 0, xn);
	SET_VECTOR_ELT(xyso, 1, yn);
	SET_VECTOR_ELT(xyso, 2, t);
	
	/* define class and attributes */
	classgets(xyso, classdf);
	setAttrib(xyso, R_NamesSymbol, namecol);
	setAttrib(xyso, R_RowNamesSymbol, namerow);

	/* checks whether the constraints are satisfied */
	defineVar(install("x"), xyso, env);
	defineVar(install("par"), parcon, env);
	PROTECT(resucont = coerceVector(eval(constraint, env), INTSXP));
	ok = INTEGER(resucont)[0];
	if (ok!=1) {
	    UNPROTECT(2);
	}
    }
    
    /* evaluate the treatment */
    defineVar(install("x"), xyso, env);
    defineVar(install("par"), par2, env);
    PROTECT(resu = eval(traitement, env));
    
    /* Return the result */
    UNPROTECT(25);
    return(resu);
}



SEXP rwrpnorm(int n, double mu, double rho)
{
    SEXP so;
    int i;
    double *sor, sd;
    
    PROTECT(so = allocVector(REALSXP, n));
    sor = REAL(so);
    
    GetRNGstate();
    if (rho< 1e-12) {
	for (i=0; i<n; i++) {
	    sor[i] = (unif_rand()*2*M_PI);
	}
    } else {
	sd = sqrt(-2.0 * log(rho));
	for (i=0; i<n; i++) {
	    sor[i] = (mu + norm_rand()*sd);
	}
    }
    PutRNGstate();
    UNPROTECT(1);
    return(so);
}


SEXP rchi(int n, double h)
{
    int i;
    SEXP so;
    double *sor;
    
    PROTECT(so = allocVector(REALSXP, n));
    sor = REAL(so);
    
    GetRNGstate();
    for (i=0; i<n; i++) {
	sor[i] = sqrt(rchisq(2))*h;
    }
    PutRNGstate();
    
    UNPROTECT(1);
    return(so);
}

SEXP tr_CRW(SEXP xyt, SEXP par1, SEXP par2, SEXP parcon,
	    SEXP traitement, SEXP constraint)
{
    SEXP nn, env, rhor, hr, angl, dista, x0, xn, yn, xyso, t, resucont, resu;
    SEXP namecol, namerow, classdf, classPOSIX;
    double *angles, *distances, ang1, *xnr, *ynr, *tr;
    int n, i, ok;

    PROTECT(env = VECTOR_ELT(par1,0));
    if(!isEnvironment(env)) error("'env' should be an environment");
    PROTECT(nn = coerceVector(VECTOR_ELT(par1, 1), INTSXP));
    PROTECT(rhor = coerceVector(VECTOR_ELT(par1, 2), REALSXP));
    PROTECT(hr = coerceVector(VECTOR_ELT(par1, 3), REALSXP));
    PROTECT(x0 = coerceVector(VECTOR_ELT(par1, 4), REALSXP));
    n = INTEGER(nn)[0];

    PROTECT(t = allocVector(REALSXP, n));
    PROTECT(xn = allocVector(REALSXP, n));
    PROTECT(yn = allocVector(REALSXP, n));
    xnr = REAL(xn);
    ynr = REAL(yn);
    tr = REAL(t);
    ok = 0;
    
    for (i = 0; i < n; i++) {
	tr[i] = (double) i+1;
    }

    /* prepares the attributes of the output data.frame */
    /* row.names */
    PROTECT(namerow = allocVector(INTSXP, n)); 
    for (i = 0; i < n; i++) {
	INTEGER(namerow)[i] = (int) i+1;
    }	    
    /* names */
    PROTECT(namecol = allocVector(STRSXP, 3)); 
    SET_STRING_ELT(namecol, 0, mkChar("x"));
    SET_STRING_ELT(namecol, 1, mkChar("y"));
    SET_STRING_ELT(namecol, 2, mkChar("date"));
    /* class data.frame */
    PROTECT(classdf = allocVector(STRSXP, 1));
    SET_STRING_ELT(classdf, 0, mkChar("data.frame"));
    /* class POSIXct for the date */
    PROTECT(classPOSIX = allocVector(STRSXP, 2));
    SET_STRING_ELT(classPOSIX, 0, mkChar("POSIXct"));
    SET_STRING_ELT(classPOSIX, 1, mkChar("POSIXt"));
    classgets(t, classPOSIX);

    while (ok != 1) {

	/* random sampling of n-2 relative angles and n-1 distances */
	PROTECT(angl = rwrpnorm(n-2, 0.0, REAL(rhor)[0]));
	PROTECT(dista = rchi(n-1, REAL(hr)[0]));
	angles = REAL(angl);
	distances = REAL(dista);
	
    
	/* first location */
	xnr[0] = REAL(x0)[0];
	ynr[0] = REAL(x0)[1];

	/* sampling of the first angle */
	GetRNGstate();
	ang1 = unif_rand() * M_PI * 2.0;
	PutRNGstate();
	
	/* second relocation */
	xnr[1] = xnr[0] + cos(ang1)*distances[0];
	ynr[1] = ynr[0] + sin(ang1)*distances[0];
	
	/* following relocations */
	for (i = 2; i < n; i++) {
	    ang1 = atan2(ynr[i-1]-ynr[i-2], xnr[i-1]-xnr[i-2]);
	    xnr[i] = xnr[i-1] + distances[i-1]*cos(ang1+angles[i-2]);
	    ynr[i] = ynr[i-1] + distances[i-1]*sin(ang1+angles[i-2]);
	}
	
	/* stores the results */
	PROTECT(xyso = allocVector(VECSXP, 3));
	SET_VECTOR_ELT(xyso, 0, xn);
	SET_VECTOR_ELT(xyso, 1, yn);
	SET_VECTOR_ELT(xyso, 2, t);

	/* define class and attributes */
	classgets(xyso, classdf);
	setAttrib(xyso, R_NamesSymbol, namecol);
	setAttrib(xyso, R_RowNamesSymbol, namerow);

	/* checks whether the constraints are satisfied */
	defineVar(install("x"), xyso, env);
	defineVar(install("par"), parcon, env);

	PROTECT(resucont = coerceVector(eval(constraint, env), INTSXP));
	

	ok = INTEGER(resucont)[0];
	if (ok!=1) {
	    UNPROTECT(4);
	}
	ok = 1;
    }
    
    /* evaluate the treatment */
    defineVar(install("x"), xyso, env);
    defineVar(install("par"), par2, env);
    PROTECT(resu = eval(traitement, env));
    
    /* Return the result */
    UNPROTECT(17);
    return(resu);
}





void compteur(int i)
{
    if (i < 10)
	Rprintf("\b");
    if ((i > 9)&&(i<100))
	Rprintf("\b\b");
    if ((i > 99)&&(i<1000))
	Rprintf("\b\b\b");
    if ((i > 999)&&(i<10000))
	Rprintf("\b\b\b\b");
    if ((i > 9999)&&(i<100000))
	Rprintf("\b\b\b\b\b");
    if ((i > 99999)&&(i<1000000))
	Rprintf("\b\b\b\b\b\b");
    if ((i > 999999)&&(i<10000000))
	Rprintf("\b\b\b\b\b\b\b");
    if ((i > 9999999)&&(i<100000000))
	Rprintf("\b\b\b\b\b\b\b\b");
    Rprintf("%i", i+1);
}




SEXP simulmod(SEXP xyt, SEXP nrepr, SEXP type, SEXP par, SEXP countr)
{
    int typi, nrep, i, count;
    SEXP res, par1, par2, parcon, traitement, constraint;
    SEXP (*fct)(SEXP, SEXP, SEXP, SEXP,
		SEXP, SEXP) = NULL;
    
    typi = INTEGER(type)[0];
    nrep = INTEGER(nrepr)[0];
    count = INTEGER(countr)[0];
    PROTECT(res = allocVector(VECSXP, nrep));
    PROTECT(par1 = VECTOR_ELT(par,0));

    PROTECT(par2 = VECTOR_ELT(par,1));
    PROTECT(parcon = VECTOR_ELT(par,2));
    PROTECT(traitement = VECTOR_ELT(par,3));
    PROTECT(constraint = VECTOR_ELT(par,4));
    
    switch (typi) {
    case 0:
	fct = tr_randomShiftRotation;
	break;
    case 1:
	fct = tr_randomCRW;
	break;
    case 2:
	fct = tr_CRW;
	break;
    case 3:
	fct = tr_randomRotationCs;
	break;
    default:
	error("type of null model not specified");
    }
    
    if (count) {
	Rprintf("Iteration:             ");
    }
    for (i=0; i<nrep; i++) {
	if (count)
	    compteur(i);
	SET_VECTOR_ELT(res, i, fct(xyt, par1, par2, parcon, 
				   traitement, constraint));
    }
    if (count) 
	Rprintf("\n");
    
    UNPROTECT(6);
    return(res);
}



SEXP simulmodmv(SEXP lixyt, SEXP nrepr, SEXP litype, SEXP lipar, SEXP countr, 
		SEXP env2, SEXP convlt, SEXP na, SEXP nlo, SEXP traitement,
		SEXP treat_par, SEXP constraint, SEXP cons_par)
{
    int nli, i, count, nrep, r, ok;
    SEXP xyt, type, par, nreppr, reso, resook, res, countint, liso;
    SEXP resucont, resu;

    nli = length(lixyt);
    PROTECT(nreppr = allocVector(INTSXP, 1));
    PROTECT(countint = allocVector(INTSXP, 1));
    INTEGER(nreppr)[0] = 1;
    INTEGER(countint)[0] = 0;
    count = INTEGER(countr)[0];
    nrep = INTEGER(nrepr)[0];
    PROTECT(liso = allocVector(VECSXP, nrep));
    

    if(!isEnvironment(env2)) error("'env2' should be an environment");
    
    if (count) {
	Rprintf("Iteration:             ");
    }
    
    for (r = 0; r < nrep; r++) {
	
	ok = 0;
	
	while (ok!=1) {
	    
	    if (count)
		compteur(r);
	    
	    PROTECT(reso = allocVector(VECSXP, nli));
	    
	    for (i=0; i<nli; i++) {
		PROTECT(par = VECTOR_ELT(lipar,i));
		PROTECT(type = VECTOR_ELT(litype,i));
		PROTECT(xyt = VECTOR_ELT(lixyt,i));
		PROTECT(res = simulmod(xyt, nreppr, type, par, countint));
		SET_VECTOR_ELT(reso, i, VECTOR_ELT(res, 0));
		UNPROTECT(4);
	    }
	    defineVar(install("lixyt"), reso, env2);
	    defineVar(install("na"), na, env2);
	    defineVar(install("nlo"), nlo, env2);
	    PROTECT(resook = eval(convlt, env2));
	    
	    
	    /* checks whether the constraints are satisfied */
	    defineVar(install("x"), resook, env2);
	    defineVar(install("par"), cons_par, env2);
	    PROTECT(resucont = coerceVector(eval(constraint, env2), INTSXP));
	    
	    ok = INTEGER(resucont)[0];
	    if (ok!=1) {
		UNPROTECT(3);
	    }
	    ok = 1;
	}
	
	/* evaluate the treatment */
	defineVar(install("x"), resook, env2);
	defineVar(install("par"), treat_par, env2);
	PROTECT(resu = eval(traitement, env2));

	SET_VECTOR_ELT(liso, r, resu);
	UNPROTECT(4);
    }


    if (count) 
	Rprintf("\n");

    UNPROTECT(3);
    return(liso);
}



/* The number of visits */
SEXP nvisits(SEXP xyt, SEXP distr, SEXP maxt)
{
    /* declaring the variables */
    int n,i, j, *deds, sortie, *nvisi, tjdeds;
    double *xr, *yr, *tr, dist, maxtr;
    double refti, refti2, a, u, v, p;
    SEXP x, y, t, dedsr, nvisit;
    
    /* coercing the arguments */
    PROTECT(x = coerceVector(VECTOR_ELT(xyt,0), REALSXP));
    PROTECT(y = coerceVector(VECTOR_ELT(xyt,1), REALSXP));
    PROTECT(t = coerceVector(VECTOR_ELT(xyt,2), REALSXP));

    n = length(x); /* number of relocations */
    PROTECT(dedsr = allocVector(INTSXP, n)); /* will be used to check whether 
						the relocation j
						is within the distance distr 
						of the relocation i */
    PROTECT(nvisit = allocVector(INTSXP, n)); /* the output vector */

    /* for the ease of manipulation: gets the pointers */
    xr = REAL(x);
    yr = REAL(y);
    tr = REAL(t);
    deds = INTEGER(dedsr);
    nvisi = INTEGER(nvisit);

    /* get the three constants passed as arguments */
    maxtr = REAL(maxt)[0];
    dist = REAL(distr)[0];
    
    /* Now, calculate the residence time for each relocation */
    for (i = 0; i < n; i++) {
	/* At least one visit in the relocation */
	nvisi[i] = 1;
	
	/* checks which relocations are within the distance dist from reloc j */
	for (j = 0; j < n; j++) {
	    if (hypot(xr[i]-xr[j], yr[i]-yr[j])<=dist) {
		deds[j] = 1;
	    } else {
		deds[j] = 0;
	    }
	}
	
	/* calculates the backward number of visits */
	sortie = 0; /* = 0 when the animal is still 
		       inside the circle; =1 otherwise */
	refti = tr[i]; /* last backward time */
	refti2 = tr[i]; /* current time of circle crossing */
	
	
	/* if this is not the first relocation (otherwise, the number of backward visits is 0) */
	if (i != 0) {
	    
	    /* not outside the circle since relocation i */
	    tjdeds = 1;
	    
	    /* for all previous relocations (because backward) */
	    for (j = i-1; j>=0; j--) {
		
		/* if the relocation is outside the circle */
		if (deds[j]==0) {

		    /* if this is the end of a movement coming from inside */
		    if (sortie==0) {
			
			/* we are outside */
			sortie = 1;

			/* interpolation of the current time of crossing.
			   use the parallelogram with a basis equal to the distance dist,
			   the diagonal corresponding to the vector connecting the
			   the relocations i and j+1 and the "other side" corresponding
			   to the segment having a length to be estimated.

			   This parallelogram is included in a rectangle having the
			   same diagonal as the parallelogram, one side defined by the 
			   projection of i on the segment (j, j+1) and by the reloc 
			   j+1 (its length is denoted u below); and the perpendicular 
			   one defined by the reloc j+1 and the projection of i on 
			   the vector orthogonal to the segment (j, j+1).
			   
			   The space in this rectangle that is not part of the parallelogram
			   forms two rectangle triangles allowing to solve this problem, 
			   estimating the length of the "other side" of the parallelogram,
			   actually equal to p*u. p is the required proportion
			*/
			a = atan2(yr[j]-yr[j+1], xr[j]-xr[j+1]);
			u = ((xr[i]-xr[j+1])*cos(a))+((yr[i]-yr[j+1])*sin(a));
			v = ((yr[i]-yr[j+1])*cos(a))-((xr[i]-xr[j+1])*sin(a));
			p = (sqrt(R_pow(dist, 2.0) - R_pow(v, 2.0))-fabs(u))/
			    hypot(xr[j]-xr[j+1], yr[j]-yr[j+1]);
			refti2 = tr[j+1] - p*(tr[j+1]-tr[j]);
			
			/* if this is the first time that the animal 
			   cross the circle since relocation i */
			if (tjdeds) {
			    /* then defines the "last backward time" =
			       current time of crossing. This means that for the 
			       first time of crossing, there is no chance that
			       this crossing results into an increment of 
			       the number of visits
			    */
			    tjdeds = 0;
			    refti = refti2;
			}
			
			/* If the difference between the current time of
			   crossing and the last backward time is larger
			   than maxtr, increase the number of visits
			 */
			if (fabs(refti2 - refti) > maxtr) {
			    nvisi[i]++;
			}
			/* and defines the new backward time */
			refti = refti2;
		    }
		} else {
		    /* We note that the animal is now out */
		    if (sortie>0) {
			sortie = 0;
		    }
		}
	    }
	}


	/* forward time */
	sortie = 0;
	refti = tr[i];
	refti2 = tr[i];
	
	/* if this is not the last relocation (otherwise, n forward visit = 0.0) */
	if (i < (n-1)) {
	    
	    /* not outside the circle since relocation i */
	    tjdeds = 1;

	    /* for all next relocations */
	    for (j = i+1; j<n; j++) {
		
		/* if the relocation is outside the circle */
		if (deds[j]==0) {
		    
		    /* if this is the first relocation outside */
		    if (sortie==0) {
			sortie = 1;
			/* interpolating the time when the animal come out
			   of the circle. Same procedure as above
			 */
			a = atan2(yr[j]-yr[j-1], xr[j]-xr[j-1]);
			u = ((xr[i]-xr[j-1])*cos(a))+((yr[i]-yr[j-1])*sin(a));
			v = ((yr[i]-yr[j-1])*cos(a))-((xr[i]-xr[j-1])*sin(a));
			p = (sqrt(R_pow(dist, 2.0) - R_pow(v, 2.0))-fabs(u))/
			    hypot(xr[j]-xr[j-1], yr[j]-yr[j-1]);
			refti2 = tr[j-1] + p*(tr[j]-tr[j-1]);
			
			if (tjdeds) {
			    tjdeds = 0;
			    refti = refti2;
			}
			if (fabs(refti2 - refti) > maxtr) {
			    nvisi[i]++;
			}
			refti = refti2;
		    }
		} else {
		    if (sortie>0) {
			sortie = 0;
		    }
		}
	    }
	}
    }
    
    UNPROTECT(5);
    
    /* output */
    return(nvisit);
    
}


