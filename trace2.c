/*
trace2.c  :  Perform a counting of components (segments) of a 
            PGM - file in ASCII (non-raw!) format
             (Note: use pnmnoraw file1.PGM file2.PGM to transform into
                       ASCII format.)

	     This version of the program (the original one is called `trace.c')
	     has an accelerated algorithm.
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define YES 1
#define NO 0

char *strcpy(char*, char*);
char *strcat(char*, char*);

output_f(float *w, int B, int H, char *filename)
{
  int i, j, temp;
  FILE *fpo;

  if((fpo=fopen(filename,"w"))==NULL)
    {
      printf("Couldn't open %s for writing !\n", filename);
      exit(-1);
    }

  fprintf(fpo,"P2\n%d %d\n%d\n", B, H, 255);

  for( j=0; j<H; j++ )
    {
      for(i = 0; i<B; i++)
	{
	  temp = (int)w[j*B+i];
	  fprintf(fpo, "%d ", temp);
	}
      fprintf(fpo, "\n");
    }

  fclose(fpo);
}


main(int argc, char **argv)
{
  int B, H, N, NN, temp1, temp2, n, r, t, s, m, i, j, k, l;
  int ngreyvalues, length_runs;
  FILE *fp;
  char line[100], outfile[120], outfile1[120];
  int *f, *range, *newrange, *deg,  *kl;
  float *u, *v, *g, *w, W, T;
  int *nr, *nu, *nl, *nd, *runlen;
  int darkfields=YES;
  
  if(argc!=4 && argc!=6)
    {
      printf("Usage: trace2 <file> n s [-b <r>] \n\n`file' is plain PGM (Magic number P2)\
\n\nWrites smoothed file in PGM ASCII format to `file.out' or `outfile', if specified.\
  \n\n n is an integer greater or equal to 1 = number of iterates \n\
     in the random walk (e.g. 3)\n\n\
 s is the minimal difference in greyvalue between two neighboring greyvalues \n\
   (range: 0 - 255) to be sufficiently different so that the edge connecting \n\
   the two vertices is removed from the lattice and no smoothing across it \n\
   takes place.\n\n -b <r> :  only count components that are lighter or equal to`r' in {0...254}.\n\n");
      exit(-1);
    }

  t = atoi(argv[2]);   /* Number of smoothing steps */
  s = atoi(argv[3]);   /* threshold in greyvalues, when to cut an edge */
  if(argc==6)          /* then last argument is output filename.  */
    {          
      if(strcmp(argv[4],"-b")==0)
	{
	  r = atoi(argv[5]);
	  darkfields = NO;
	}
      else
	{
	  printf("Usage: trace2 <file> n s [-b <r>] \n\n");
	  exit(0);
	}
    }

  strcpy(outfile, argv[1]);  
  strcat(outfile,".out");
    
  if((fp=fopen(argv[1],"r"))==NULL)
    {
      printf("Couldn't open %s for reading !\n", argv[1]);
      exit(-1);
    }

  fscanf(fp,"%s", line);
  if(strcmp(line, "P2")!=0)
    {
      printf("'%s' is not a proper pgm-file.\n\n", argv[1]);
      exit(-1);
      }  

  fscanf(fp,"%s", line);
  if(line[0]=='#')
    {
      printf("Remove comment line (2nd line) from '%s'.\n\n", argv[1]);
      exit(-1);
      }
  fclose(fp);          /* Initial checking of inputfile finished. */

  if((fp=fopen(argv[1],"r"))==NULL)
    {
      printf("Couldn't open %s for reading !\n", argv[1]);
      exit(-1);
    }                  /* file is reopened, so that fp* jumps to beginning of it. */

  fscanf(fp,"%s", line);
  temp1 = fscanf(fp,"%d %d", &B, &H);
  temp2 = fscanf(fp,"%d", &ngreyvalues);
  if(ngreyvalues!=255)
    { 
      printf("File doesn't have 0-255 greyvalues! (B:%d, H:%d)\n\n", B, H);
      exit(-1);
    }
  if(temp1!=2 || temp2!=1)  /* in this case fscanf didn't read in 3 values successfully, s.a. */
    {
      printf("'%s' is not a proper file.\n\n", argv[1]);
      exit(-1);
    }         /* Second checking of inputfile finished. If no `exit(-1)' so far, file ok. */
  

  N = H*B;    /* Number of Pixels */

  f = (int *)malloc((sizeof(int)*N+1));       /* The input image - stays fixed */
  nr = (int *)malloc((sizeof(int)*N+1));      /* edge to right present? */
  nu = (int *)malloc((sizeof(int)*N+1));      /* edge up present? */
  nl = (int *)malloc((sizeof(int)*N+1));      /* edge to left present? */
  nd = (int *)malloc((sizeof(int)*N+1));      /* edge downward present? */
  u = (float *)malloc((sizeof(float)*N*N+N)); /* The image which is changed from step to step */
  v = (float *)malloc((sizeof(float)*N+1));   /* The image one time step back  -  the last image */
  g = (float *)malloc((sizeof(float)*N*N+N)); /* The image one time step back  -  the last image */  
  w = (float *)malloc((sizeof(float)*N+1));   /* Eventually the value of the normalized trace */
  deg = (int *)malloc((sizeof(int)*N+1));     /* number of neighbors of each of the vertices */
  runlen = (int *)malloc((sizeof(int)*N+1));  /* list of runs */
  kl = (int *)malloc((sizeof(int)*N+1));      /* list ofinitial positions of runs */

  for(j = 0; j<H; j++)
    {
      for(i = 0; i<B; i++)
	{
	  fscanf(fp, "%d", &(f[j*B+i]));
	}
    }
  fclose(fp);   
  /* Input finished. f[0..N] should have values between 0 and 255, now */

  
  for(i=0; i<N; i++)  /* Initializing deg[], and nr[],nu[],nl[],nd[] */
    {                 /* This is the abbreviation for neighborright, neighborup, 			                      neighborleft, neighbordown */
      deg[i] = 0;
      nr[i] = 0;
      nu[i] = 0;
      nl[i] = 0;
      nd[i] = 0;
      if((i+1)%B!=0)  /* not on right edge */
	if(f[i]>=f[i+1]-s && f[i]<=f[i+1]+s)                      
	  {         
	         /* a large value of s here will make it */
	         /* more likely that smoothing occurs between f[i] and f[i+1] etc. */
	    deg[i]++;
	    nr[i] = 1;
	  }
      if(i%B !=0)    /* not on left edge */
	if(f[i]>=f[i-1]-s && f[i]<=f[i-1]+s)
	  {
	    deg[i]++;
	    nl[i]=1;
	  }
      if(i<N-B)     /* not on upper edge */
	if(f[i]>=f[i+B]-s && f[i]<=f[i+B]+s)
	  {
	    deg[i]++;
	    nu[i]=1;
	  }
      if(i>B-1)     /* not on lower edge */
	if(f[i]>=f[i-B]-s && f[i]<=f[i-B]+s)
	  {
	    deg[i]++;
	    nd[i]=1;
	  }
    } 
      /* Degree-vector deg[] set up  -  Note: maximum is 4; */
      /* Also nr[], nu[], nl[], nd[] also set up, now*/

      /*  vertex around which smoothing occurs: */
      /*     m = j*B + i  - goes through whole picture: */

  for(k=0; k<N; k++)
    {
      v[k] = f[k]; 	  /*  v[i]=f[i];  /* Picture in the beginning - no smoothing yet */
      for(i=0; i<N; i++)
	{
	  g[k*N+i] = 0.0; /* Initialization */
	  u[k*N+i] = 0.0;
	}
      g[k*N+k] = 1.0;     /* The k'th g[] is g[k*N + i ] = 0 unless k=i, then it is equal to 1.0 */
    }


  /* Setting up runlen[.]: runlen[i] = inf{ k>0 : | f[kl[i] + k] - f[kl[i] + k + 1 | > s } */

  k=0;  /* in {0, ..., N-1 } : goes through all pixels */
  i=0;  /* in {0, ..., length_runs-1} : goes only through the start-pixels of runs */
  while(k<N)
    {
      if(darkfields==NO) /* If option `-b <r>' was set, then k is increased until f[k] >= r. */	
	while(f[k]<r)    /* Count only the `bright spots' -- r  defines what bright means. */
	  k++;
      if(k<N)
	{
	  runlen[i]=1;        /* each run has at least one pixel */
	  kl[i] = k;          /* coordinate in {0, ..., N-1} of the i'th start-pixel */
	  l=k+1;              /* coordinate of the pixel to the right of k */
	  while(nr[l-1]==1)   /* if there is an edge between l-1 and l */
	    {                 /* Note: using nr[l-1] circumvents segm. faults! */
	      runlen[i]++;    /* In this case the run has at least length 2. */
	      l++;            /* go to the next pixel on the right */
	    }	
	  k = k+runlen[i];    /* Jump to the start-pixel of the next run */
	  i++;                /* Increase the index of the list of runs */
	}
    } 


  length_runs = i;            /* This is the largest index + 1 (in the list of runs) */

  m=0;                        /* This is a checksum, to see if the sum of the runs */
                              /* is equal to the number of pixels. */   
  for(j=0;j<length_runs;j++)  /* adding up the runs' lengths */
    m=m+runlen[j];
  
  if(darkfields==YES && m!=N) /* Error check if runlen[] has been determined correctly */
    {
      printf("Warning: Sum %d of runs is not equal to number of pixels %d !\n\n", m ,N);
      exit(-1);
    }
      
         /* Strip runlen[]  of it's runs across dark < r (in {1, ...., 255})  fields */
  
         /* Main Loop: */

  system("rm /home/florian/Pictures/experiment/out/*");

  T=0.0;  /* This is the value of the trace */

  for(l=0; l<length_runs; l++)     /* Run through all starting-pixels of runs. */
  {
    for(n = 0; n<t; n++)           /* Time loop */
      {	  
	for(i=0; i<N; i++)         /* Go through ALL the pixels of the pictures. */
	  {
	    u[kl[l]*N+i] = 0.0;    /* Initialize i'th pixel of l'th picture. */

	    if(nr[i]==1)           /* If there is an edge between i and i+1 ...*/
	      u[kl[l]*N+i] += g[kl[l]*N+i+1];   /* ... add the value of g at i+1 */
	    if(nl[i]==1)
	      u[kl[l]*N+i] += g[kl[l]*N+i-1];   /* ... same for other neighbors ...*/
	    if(nu[i]==1)
	      u[kl[l]*N+i] += g[kl[l]*N+i+B];
	    if(nd[i]==1)
	      u[kl[l]*N+i] += g[kl[l]*N+i-B];
	    
	    u[kl[l]*N+i] += g[kl[l]*N+i]*(4-deg[i]);  /* Finally add (4-deg(i)) times g[i] */
	                                              /* Corresponds to `loops' */
	                                              /* of delayed random walk  (DRW) */
	    u[kl[l]*N+i] = 0.25*u[kl[l]*N+i];         /* Normalize (Stochastic matrix!) */
	  }
	
	for(i=0; i<N; i++)  /* Put new picture into old */
	  g[kl[l]*N+i] = u[kl[l]*N+i];      /*  after every step */
      }    
    
    T += runlen[l]*u[kl[l]*N+kl[l]];   /* Add (length of l'th run) times P^n_xx where x=kl[l] */
    
    for(i=0;i<N;i++)
      v[i] = 10000*g[N*kl[l]+i];         /* Store final l'th annealed picture into v[] */
    
    sprintf(outfile1, "out/%s_%d", outfile, l);   /* Give it name `outfile1' */
    
    output_f(v, B, H, outfile1);   /* Output l'th picture into out/ */    
  }

  printf(" %9.2f \n", T );         /* Report approximation of Trace[P^n] */
}
