#ifndef KD_HINCLUDED
#define KD_HINCLUDED

#define ROOT		1
#define LOWER(i)	(i<<1)
#define UPPER(i)	((i<<1)+1)
#define PARENT(i)	(i>>1)
#define SIBLING(i) 	((i&1)?i-1:i+1)
#define SETNEXT(i)\
{\
	while (i&1) i=i>>1;\
	++i;\
	}

#define DARK	1
#define GAS		2
#define STAR	4

#define KD_ORDERTEMP	256

typedef struct Particle {
	float fMass;
	float r[3];
	float v[3];
	int iGroup;
	int iOrder;
	} PARTICLE;

typedef struct bndBound {
	float fMin[3];
	float fMax[3];
	} BND;

typedef struct kdNode {
	float fSplit;
	BND bnd;
	int iDim;
	int pLower;
	int pUpper;
	} KDN;

typedef struct kdContext {
	int nBucket;
	int nParticles;
	int nDark;
	int nGas;
	int nStar;
	int bDark;
	int bGas;
	int bStar;
	int nActive;
	float fTime;
	float fPeriod[3];
	float fCenter[3];
	int nLevels;
	int nNodes;
	int nSplit;
	PARTICLE *p;
	KDN *kdNodes;
	int nGroup;
	int uSecond;
	int uMicro;
	} * KD;


#define INTERSECT(c,cp,fBall2,lx,ly,lz,x,y,z,sx,sy,sz)\
{\
	float INTRSCT_dx,INTRSCT_dy,INTRSCT_dz;\
	float INTRSCT_dx1,INTRSCT_dy1,INTRSCT_dz1;\
	float INTRSCT_fDist2,INTRSCT_fMax2;\
	INTRSCT_dx = c[cp].bnd.fMin[0]-x;\
	INTRSCT_dx1 = x-c[cp].bnd.fMax[0];\
	INTRSCT_dy = c[cp].bnd.fMin[1]-y;\
	INTRSCT_dy1 = y-c[cp].bnd.fMax[1];\
	INTRSCT_dz = c[cp].bnd.fMin[2]-z;\
	INTRSCT_dz1 = z-c[cp].bnd.fMax[2];\
	if (INTRSCT_dx > 0.0) {\
		if (INTRSCT_dx1+lx < INTRSCT_dx) {\
			INTRSCT_dx1 += lx;\
			INTRSCT_dx -= lx;\
			sx = x+lx;\
			INTRSCT_fDist2 = INTRSCT_dx1*INTRSCT_dx1;\
			INTRSCT_fMax2 = INTRSCT_dx*INTRSCT_dx;\
			}\
		else {\
			sx = x;\
			INTRSCT_fDist2 = INTRSCT_dx*INTRSCT_dx;\
			INTRSCT_fMax2 = INTRSCT_dx1*INTRSCT_dx1;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (INTRSCT_dx1 > 0.0) {\
		if (INTRSCT_dx+lx < INTRSCT_dx1) {\
		    INTRSCT_dx += lx;\
			INTRSCT_dx1 -= lx;\
			sx = x-lx;\
			INTRSCT_fDist2 = INTRSCT_dx*INTRSCT_dx;\
			INTRSCT_fMax2 = INTRSCT_dx1*INTRSCT_dx1;\
			}\
		else {\
			sx = x;\
			INTRSCT_fDist2 = INTRSCT_dx1*INTRSCT_dx1;\
			INTRSCT_fMax2 = INTRSCT_dx*INTRSCT_dx;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sx = x;\
		INTRSCT_fDist2 = 0.0;\
		if (INTRSCT_dx < INTRSCT_dx1) INTRSCT_fMax2 = INTRSCT_dx*INTRSCT_dx;\
		else INTRSCT_fMax2 = INTRSCT_dx1*INTRSCT_dx1;\
		}\
	if (INTRSCT_dy > 0.0) {\
		if (INTRSCT_dy1+ly < INTRSCT_dy) {\
		    INTRSCT_dy1 += ly;\
			INTRSCT_dy -= ly;\
			sy = y+ly;\
			INTRSCT_fDist2 += INTRSCT_dy1*INTRSCT_dy1;\
			INTRSCT_fMax2 += INTRSCT_dy*INTRSCT_dy;\
			}\
		else {\
			sy = y;\
			INTRSCT_fDist2 += INTRSCT_dy*INTRSCT_dy;\
			INTRSCT_fMax2 += INTRSCT_dy1*INTRSCT_dy1;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (INTRSCT_dy1 > 0.0) {\
		if (INTRSCT_dy+ly < INTRSCT_dy1) {\
		    INTRSCT_dy += ly;\
			INTRSCT_dy1 -= ly;\
			sy = y-ly;\
			INTRSCT_fDist2 += INTRSCT_dy*INTRSCT_dy;\
			INTRSCT_fMax2 += INTRSCT_dy1*INTRSCT_dy1;\
			}\
		else {\
			sy = y;\
			INTRSCT_fDist2 += INTRSCT_dy1*INTRSCT_dy1;\
			INTRSCT_fMax2 += INTRSCT_dy*INTRSCT_dy;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sy = y;\
		if (INTRSCT_dy < INTRSCT_dy1) INTRSCT_fMax2 += INTRSCT_dy*INTRSCT_dy;\
		else INTRSCT_fMax2 += INTRSCT_dy1*INTRSCT_dy1;\
		}\
	if (INTRSCT_dz > 0.0) {\
		if (INTRSCT_dz1+lz < INTRSCT_dz) {\
		    INTRSCT_dz1 += lz;\
            INTRSCT_dz -= lz;\
			sz = z+lz;\
			INTRSCT_fDist2 += INTRSCT_dz1*INTRSCT_dz1;\
			INTRSCT_fMax2 += INTRSCT_dz*INTRSCT_dz;\
			}\
		else {\
			sz = z;\
			INTRSCT_fDist2 += INTRSCT_dz*INTRSCT_dz;\
			INTRSCT_fMax2 += INTRSCT_dz1*INTRSCT_dz1;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else if (INTRSCT_dz1 > 0.0) {\
		if (INTRSCT_dz+lz < INTRSCT_dz1) {\
			INTRSCT_dz += lz;\
		    INTRSCT_dz1 -= lz;\
			sz = z-lz;\
			INTRSCT_fDist2 += INTRSCT_dz*INTRSCT_dz;\
			INTRSCT_fMax2 += INTRSCT_dz1*INTRSCT_dz1;\
			}\
		else {\
			sz = z;\
			INTRSCT_fDist2 += INTRSCT_dz1*INTRSCT_dz1;\
			INTRSCT_fMax2 += INTRSCT_dz*INTRSCT_dz;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto GetNextCell;\
		}\
	else {\
		sz = z;\
		if (INTRSCT_dz < INTRSCT_dz1) INTRSCT_fMax2 += INTRSCT_dz*INTRSCT_dz;\
		else INTRSCT_fMax2 += INTRSCT_dz1*INTRSCT_dz1;\
		}\
	if (INTRSCT_fMax2 < fBall2) goto ContainedCell;\
	}


void kdTime(KD,int *,int *);
int kdInit(KD *,int,float *,float *);
void kdReadTipsy(KD,FILE *,int,int,int,int);
void kdBuildTree(KD);
int kdFoF(KD,float);
int kdTooSmall(KD,int);
void kdOrder(KD);
void kdOutGroup(KD,char *);
void kdOutGTP(KD,char *,int);
void kdFinish(KD);

#endif











