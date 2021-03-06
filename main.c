#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "kd.h"


void usage(void)
{
	fprintf(stderr,"USAGE:\n");
	fprintf(stderr,"fof -e <Linking Length>\n");
	fprintf(stderr,"   [-m <nMinMembers>] [-dgs] [-v]\n");
	fprintf(stderr,"   [-o <Output Name>] [-p <xyzPeriod>] [-c <xyzCenter>]\n");
	fprintf(stderr,"   [-px <xPeriod>] [-py <yPeriod>] [-pz <zPeriod>]\n");
	fprintf(stderr,"   [-cx <xCenter>] [-cy <yCenter>] [-cz <zCenter>]\n");
	fprintf(stderr,"   [-std]\n");
	fprintf(stderr,"Input taken from stdin in tipsy binary format.\n");
    fprintf(stderr,"Default is to link all types of particles.\n");
	fprintf(stderr,"SEE MAN PAGE: fof(1) for more information.\n");
	exit(1);
	}

int main(int argc,char **argv)
{
	KD kd;
	int nBucket,i,j;
	char ach[80],tmp[80];
	float fPeriod[3],fCenter[3],fEps;
	int bDark,bGas,bStar;
	int nMembers,nGroup,bVerbose,bStandard;
	int sec,usec;
	char *p;
    long lStart,lStartAll;
	
   	nBucket = 32;
	nMembers = 8;
	bDark = 1;
	bGas = 1;
	bStar = 1;
	bVerbose = 0;
	bStandard = 0;
	strcpy(ach,"fof");
	i = 1;
	for (j=0;j<3;++j) {
		fPeriod[j] = FLT_MAX;
		fCenter[j] = 0.0;
		}
	while (i < argc) {
		if (!strcmp(argv[i],"-e")) {
			++i;
			fEps = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-m")) {
			++i;
			nMembers = atoi(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-o")) {
			++i;
			strcpy(ach,argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-p")) {
			++i;
			fPeriod[0] = atof(argv[i]);
			fPeriod[1] = atof(argv[i]);
			fPeriod[2] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-px")) {
			++i;
			fPeriod[0] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-py")) {
			++i;
			fPeriod[1] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-pz")) {
			++i;
		    fPeriod[2] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-c")) {
			++i;
			if (i >= argc) usage();
			fCenter[0] = atof(argv[i]);
			fCenter[1] = atof(argv[i]);
			fCenter[2] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-cx")) {
			++i;
			if (i >= argc) usage();
			fCenter[0] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-cy")) {
			++i;
			if (i >= argc) usage();
			fCenter[1] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-cz")) {
			++i;
			if (i >= argc) usage();
		    fCenter[2] = atof(argv[i]);
			++i;
			}
		else if (!strcmp(argv[i],"-v")) {
			bVerbose = 1;
			++i;
			}
		else if (!strcmp(argv[i],"-std")) {
		        bStandard = 1;
			++i;
		        }
		else if (*argv[i] == '-') {
			p = argv[i];
			++p;
			if (*p == 'd' || *p == 'g' || *p == 's') {
				bDark = 0;
				bGas = 0;
				bStar = 0;
				}
			else usage();
			while (isalpha(*p)) {
				switch (*p) {
				case 'd':
					bDark = 1;
					break;
				case 'g':
					bGas = 1;
					break;
				case 's':
					bStar = 1;
					break;
				default:
					usage();
					}
				++p;
				}
			++i;
			}
		else usage();
		}

#ifdef _OPENMP 
#pragma omp parallel
        if (omp_get_thread_num() == 0)
            printf("FOF (OpenMP) running on %d threads.\n",omp_get_num_threads());
#else
        printf("FOF (serial) running on 1 thread.\n");
#endif
	kdInit(&kd,nBucket,fPeriod,fCenter);
    /* Time the input phase */
    lStart = time(0);
    lStartAll = lStart;
	kdTime(kd,&sec,&usec);
	kdReadTipsy(kd,stdin,bDark,bGas,bStar,bStandard);
	kdTime(kd,&sec,&usec);
    if (bVerbose) {
        printf("Read input file of %d total particles.\n",kd->nParticles);
        printf("READ CPU TIME:  %d.%06d secs\n",sec,usec);
        printf("READ WALL TIME: %d secs\n",(int)(time(0)-lStart));
        }
	kdBuildTree(kd);
    if (bVerbose) printf("Tree built.  Beginning FOF...\n");
    lStart = time(0);
	kdTime(kd,&sec,&usec);
	nGroup = kdFoF(kd,fEps);
	kdTime(kd,&sec,&usec);
	if (bVerbose) {
        printf("Number of initial groups: %d\n",nGroup);
#ifdef _OPENMP
        printf("Size of particle lock array: %d\n",kd->nHash);
#endif
        }
	nGroup = kdTooSmall(kd,nMembers);
	if (bVerbose) {
		printf("Number of groups: %d\n",nGroup);
		printf("FOF CPU TIME:  %d.%06d secs\n",sec,usec);
        printf("FOF WALL TIME: %d\n",(int)(time(0)-lStart));
		}
	strcpy(tmp,ach);
	strcat(tmp,".gtp");
	kdOutGTP(kd,tmp,bStandard);
	/*
	 ** Now re-order the particles.
	 */
	kdOrder(kd);
	strcpy(tmp,ach);
	strcat(tmp,".grp");
	kdOutGroup(kd,tmp);
	kdFinish(kd);
    printf("FOF completed in %d wallclock seconds\n.",(int)(time(0)-lStartAll));
	}

