#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>
#include <rpc/types.h>
#include <rpc/xdr.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "kd.h"
#include "tipsydefs.h"

int xdrHeader(XDR *pxdrs,struct dump *ph)
{
	int pad = 0;
	
	if (!xdr_double(pxdrs,&ph->time)) return 0;
	if (!xdr_int(pxdrs,&ph->nbodies)) return 0;
	if (!xdr_int(pxdrs,&ph->ndim)) return 0;
	if (!xdr_int(pxdrs,&ph->nsph)) return 0;
	if (!xdr_int(pxdrs,&ph->ndark)) return 0;
	if (!xdr_int(pxdrs,&ph->nstar)) return 0;
	if (!xdr_int(pxdrs,&pad)) return 0;
	return 1;
	}

void kdTime(KD kd,int *puSecond,int *puMicro)
{
	struct rusage ru;

	getrusage(0,&ru);
	*puMicro = ru.ru_utime.tv_usec - kd->uMicro;
	*puSecond = ru.ru_utime.tv_sec - kd->uSecond;
	if (*puMicro < 0) {
		*puMicro += 1000000;
		*puSecond -= 1;
		}
	kd->uSecond = ru.ru_utime.tv_sec;
	kd->uMicro = ru.ru_utime.tv_usec;
	}


int kdInit(KD *pkd,int nBucket,float *fPeriod,float *fCenter)
{
	KD kd;
	int j;

	kd = (KD)malloc(sizeof(struct kdContext));
	assert(kd != NULL);
	kd->nBucket = nBucket;
	for (j=0;j<3;++j) {
		kd->fPeriod[j] = fPeriod[j];
		kd->fCenter[j] = fCenter[j];
		}
	kd->p = NULL;
	kd->kdNodes = NULL;
	*pkd = kd;
	return(1);
	}


void kdReadTipsy(KD kd,FILE *fp,int bDark,int bGas,int bStar,int bStandard)
{
	int i,j,nCnt;
	struct dump h;
	struct gas_particle gp;
	struct dark_particle dp;
	struct star_particle sp;
	XDR xdrs;

	if (bStandard) {
	    assert(sizeof(Real)==sizeof(float)); /* Otherwise, this XDR stuff
						    ain't gonna work */
	    xdrstdio_create(&xdrs, fp, XDR_DECODE);
	    xdrHeader(&xdrs,&h);
	} else {
	    fread(&h,sizeof(struct dump),1,fp);
	}
	kd->nParticles = h.nbodies;
	kd->nDark = h.ndark;
	kd->nGas = h.nsph;
	kd->nStar = h.nstar;
	kd->fTime = h.time;
	kd->nActive = 0;
	if (bDark) kd->nActive += kd->nDark;
	if (bGas) kd->nActive += kd->nGas;
	if (bStar) kd->nActive += kd->nStar;
	kd->bDark = bDark;
	kd->bGas = bGas;
	kd->bStar = bStar;
	/*
	 ** Allocate particles.
	 */
	kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));
	assert(kd->p != NULL);
	/*
	 ** Read Stuff!
	 */
	nCnt = 0;
	for (i=0;i<h.nsph;++i) {
		if (bStandard) {
			xdr_vector(&xdrs, (char *) &gp,
				   sizeof(struct gas_particle)/sizeof(Real),
				   sizeof(Real), (xdrproc_t)xdr_float);
		} else {
			fread(&gp,sizeof(struct gas_particle),1,fp);
		}
		if (bGas) {
			kd->p[nCnt].iOrder = nCnt;
			kd->p[nCnt].fMass = gp.mass;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = gp.pos[j];
			for (j=0;j<3;++j) kd->p[nCnt].v[j] = gp.vel[j];
			++nCnt;
		}
	}
	for (i=0;i<h.ndark;++i) {
		if (bStandard) {
			xdr_vector(&xdrs, (char *) &dp,
				   sizeof(struct dark_particle)/sizeof(Real),
                       sizeof(Real), (xdrproc_t)xdr_float);
		} else {
		    fread(&dp,sizeof(struct dark_particle),1,fp);
		}
		if (bDark) {
			kd->p[nCnt].iOrder = nCnt;
			kd->p[nCnt].fMass = dp.mass;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = dp.pos[j];
			for (j=0;j<3;++j) kd->p[nCnt].v[j] = dp.vel[j];
			++nCnt;
		}
	}
	for (i=0;i<h.nstar;++i) {
		if (bStandard) {
			xdr_vector(&xdrs, (char *) &sp,
				   sizeof(struct star_particle)/sizeof(Real),
				   sizeof(Real), (xdrproc_t)xdr_float);
		} else {
			fread(&sp,sizeof(struct star_particle),1,fp);
		}
		if (bStar) {
			kd->p[nCnt].iOrder = nCnt;
			kd->p[nCnt].fMass = sp.mass;
			for (j=0;j<3;++j) kd->p[nCnt].r[j] = sp.pos[j];
			for (j=0;j<3;++j) kd->p[nCnt].v[j] = sp.vel[j];
			++nCnt;
		}
	}
	if (bStandard) xdr_destroy(&xdrs);
}


void kdSelect(KD kd,int d,int k,int l,int r)
{
	PARTICLE *p,t;
	double v;
	int i,j;

	p = kd->p;
	while (r > l) {
		v = p[k].r[d];
		t = p[r];
		p[r] = p[k];
		p[k] = t;
		i = l - 1;
		j = r;
		while (1) {
			while (i < j) if (p[++i].r[d] >= v) break;
			while (i < j) if (p[--j].r[d] <= v) break;
			t = p[i];
			p[i] = p[j];
			p[j] = t;
			if (j <= i) break;
			}
		p[j] = p[i];
		p[i] = p[r];
		p[r] = t;
		if (i >= k) r = i - 1;
		if (i <= k) l = i + 1;
		}
	}


void kdCombine(KDN *p1,KDN *p2,KDN *pOut)
{
	int j;

	/*
	 ** Combine the bounds.
	 */
	for (j=0;j<3;++j) {
		if (p2->bnd.fMin[j] < p1->bnd.fMin[j])
			pOut->bnd.fMin[j] = p2->bnd.fMin[j];
		else
			pOut->bnd.fMin[j] = p1->bnd.fMin[j];
		if (p2->bnd.fMax[j] > p1->bnd.fMax[j])
			pOut->bnd.fMax[j] = p2->bnd.fMax[j];
		else
			pOut->bnd.fMax[j] = p1->bnd.fMax[j];
		}
	}


void kdUpPass(KD kd,int iCell)
{
	KDN *c;
	int l,u,pj,j;

	c = kd->kdNodes;
	if (c[iCell].iDim != -1) {
		l = LOWER(iCell);
		u = UPPER(iCell);
		kdUpPass(kd,l);
		kdUpPass(kd,u);
		kdCombine(&c[l],&c[u],&c[iCell]);
		}
	else {
		l = c[iCell].pLower;
		u = c[iCell].pUpper;
		for (j=0;j<3;++j) {
			c[iCell].bnd.fMin[j] = kd->p[u].r[j];
			c[iCell].bnd.fMax[j] = kd->p[u].r[j];
			}
		for (pj=l;pj<u;++pj) {
			for (j=0;j<3;++j) {
				if (kd->p[pj].r[j] < c[iCell].bnd.fMin[j])
					c[iCell].bnd.fMin[j] = kd->p[pj].r[j];
				if (kd->p[pj].r[j] > c[iCell].bnd.fMax[j])
					c[iCell].bnd.fMax[j] = kd->p[pj].r[j];
				}
			}
		}
	}

void kdBuildTree(KD kd)
{
	int l,n,i,d,m,j,diff;
	KDN *c;
	BND bnd;

	n = kd->nActive;
	kd->nLevels = 1;
	l = 1;
	while (n > kd->nBucket) {
		n = n>>1;
		l = l<<1;
		++kd->nLevels;
		}
	kd->nSplit = l;
	kd->nNodes = l<<1;
	if (kd->kdNodes != NULL) free(kd->kdNodes);
	kd->kdNodes = (KDN *)malloc(kd->nNodes*sizeof(KDN));
	assert(kd->kdNodes != NULL);
	/*
	 ** Calculate Bounds.
	 */
	for (j=0;j<3;++j) {
		bnd.fMin[j] = kd->p[0].r[j];
		bnd.fMax[j] = kd->p[0].r[j];
		}
	for (i=1;i<kd->nActive;++i) {
		for (j=0;j<3;++j) {
			if (bnd.fMin[j] > kd->p[i].r[j]) 
				bnd.fMin[j] = kd->p[i].r[j];
			else if (bnd.fMax[j] < kd->p[i].r[j])
				bnd.fMax[j] = kd->p[i].r[j];
			}
		}
	/*
	 ** Set up ROOT node
	 */
	c = kd->kdNodes;
	c[ROOT].pLower = 0;
	c[ROOT].pUpper = kd->nActive-1;
	c[ROOT].bnd = bnd;
	i = ROOT;
	while (1) {
		assert(c[i].pUpper - c[i].pLower + 1 > 0);
		if (i < kd->nSplit && (c[i].pUpper - c[i].pLower) > 0) {
			d = 0;
			for (j=1;j<3;++j) {
				if (c[i].bnd.fMax[j]-c[i].bnd.fMin[j] > 
					c[i].bnd.fMax[d]-c[i].bnd.fMin[d]) d = j;
				}
			c[i].iDim = d;

			m = (c[i].pLower + c[i].pUpper)/2;
			kdSelect(kd,d,m,c[i].pLower,c[i].pUpper);

			c[i].fSplit = kd->p[m].r[d];
			c[LOWER(i)].bnd = c[i].bnd;
			c[LOWER(i)].bnd.fMax[d] = c[i].fSplit;
			c[LOWER(i)].pLower = c[i].pLower;
			c[LOWER(i)].pUpper = m;
			c[UPPER(i)].bnd = c[i].bnd;
			c[UPPER(i)].bnd.fMin[d] = c[i].fSplit;
			c[UPPER(i)].pLower = m+1;
			c[UPPER(i)].pUpper = c[i].pUpper;
			diff = (m-c[i].pLower+1)-(c[i].pUpper-m);
			assert(diff == 0 || diff == 1);
			i = LOWER(i);
			}
		else {
			c[i].iDim = -1;
			SETNEXT(i);
			if (i == ROOT) break;
			}
		}
	kdUpPass(kd,ROOT);
	}

#ifdef _OPENMP
/* Returns the lock ID associated with particle pid. */
int _hashLock(KD kd,int pid)
{
    return pid % (kd->nHash);
}
#endif

int kdFoF(KD kd,float fEps)
{
	PARTICLE *p;
	KDN *c;
	int pi,pj,pn,cp;

	int iGroup;

	int *Fifo,iHead,iTail,nFifo;
	float fEps2;
	float dx,dy,dz,x,y,z,lx,ly,lz,sx,sy,sz,fDist2;
#ifdef _OPENMP
    int idSelf;
    omp_lock_t *locks;

	for (pn=0;pn<kd->nActive;++pn) kd->p[pn].iTouched = -1;
    /* We really want to make an independent lock for each particle.  However, each lock
     * seems to use a buttload of memory (something like 312 bytes per lock).  Therefore,
     * to ensure that we don't use too much memory, only use 1 lock per 100 particles.
     * This should still provide very low lock contention while not using oodles of 
     * memory at the same time, since it is extremely rare that two threads will be looking
     * two particles that map to the same lock at the same time.*/
    kd->nHash = (int)(kd->nActive/100);
    locks = (omp_lock_t *)malloc(kd->nHash*sizeof(omp_lock_t));
    assert(locks != NULL);
    for (pn=0;pn<kd->nHash;++pn) omp_init_lock(&locks[pn]);
#endif

	p = kd->p;
	c = kd->kdNodes;
	lx = kd->fPeriod[0];
	ly = kd->fPeriod[1];
	lz = kd->fPeriod[2];
	fEps2 = fEps*fEps;
	for (pn=0;pn<kd->nActive;++pn) p[pn].iGroup = 0;
#pragma omp parallel default(none) shared(kd,locks,p,c,lx,ly,lz,fEps2) \
    private(pi,pj,pn,cp,iGroup,Fifo,iHead,iTail,dx,dy,dz,x,y,z,sx,sy,sz,fDist2,idSelf,nFifo)
  {
#ifdef _OPENMP
    nFifo = kd->nActive/omp_get_num_threads();
    idSelf = omp_get_thread_num();
#else
    nFifo = kd->nActive;
#endif
	Fifo = (int *)malloc(nFifo*sizeof(int));
	assert(Fifo != NULL);
	iHead = 0;
	iTail = 0;
	iGroup = 0;
#pragma omp for schedule(runtime)
	for (pn=0;pn<kd->nActive;++pn) {
		if (p[pn].iGroup) continue;
		/*
		 ** Mark it and add to the do-fifo.
		 */
#ifdef _OPENMP
        omp_set_lock(&locks[_hashLock(kd,pn)]);
        if (p[pn].iTouched >= 0 && p[pn].iTouched < idSelf ) {
            assert(p[pn].iGroup > 0);
            omp_unset_lock(&locks[_hashLock(kd,pn)]);
            continue;
        }
        p[pn].iTouched = idSelf;
        iGroup = pn+1;
		p[pn].iGroup = iGroup;
        omp_unset_lock(&locks[_hashLock(kd,pn)]);
#else
		++iGroup;
		p[pn].iGroup = iGroup;
#endif
		Fifo[iTail++] = pn;
		if (iTail == nFifo) iTail = 0;
		while (iHead != iTail) {
			pi = Fifo[iHead++];
			if (iHead == nFifo) iHead = 0;
			/*
			 ** Now do an fEps-Ball Gather!
			 */
			x = p[pi].r[0];
			y = p[pi].r[1];
			z = p[pi].r[2];
			cp = ROOT;
			while (1) {
				INTERSECT(c,cp,fEps2,lx,ly,lz,x,y,z,sx,sy,sz);
				/*
				 ** We have an intersection to test.
				 */
				if (c[cp].iDim >= 0) {
					cp = LOWER(cp);
					continue;
					}
				else {
					for (pj=c[cp].pLower;pj<=c[cp].pUpper;++pj) {
#ifdef _OPENMP
                        if (p[pj].iGroup == iGroup) {
                            /* We have already looked at this particle */
                            //assert(p[pj].iTouched == idSelf);  particle is not locked.
                            continue;
                        }
                        if (p[pj].iTouched >= 0 && p[pj].iTouched < idSelf) {
                            /* Somebody more important than us is already looking at this
                             * particle.  However, we do not yet know if this particle belongs
                             * in our group, so just skip it to save time but don't restart the
                             * entire group. */
                            // assert(p[pj].iGroup > 0); particle is not locked
                            continue;
                        }
#else
						if (p[pj].iGroup) continue;
#endif
						dx = sx - p[pj].r[0];
						dy = sy - p[pj].r[1];
						dz = sz - p[pj].r[2];
						fDist2 = dx*dx + dy*dy + dz*dz;
						if (fDist2 < fEps2) {
							/*
							 ** Mark it and add to the do-fifo.
							 */
#ifdef _OPENMP
                            omp_set_lock(&locks[_hashLock(kd,pj)]);
                            if (p[pj].iTouched >= 0 && p[pj].iTouched < idSelf) {
                                /* Now we know this particle should be in our group.  If somebody more
                                 * important than us touched it, about the entire group. */
                                assert(p[pj].iGroup > 0);
                                omp_unset_lock(&locks[_hashLock(kd,pj)]);
                                iHead = iTail;
                                /*printf("Thread %d: Aborting group %d. p[%d].iOrder  p.iGroup=%d  p.iTouched=%d (Per-Particle2)\n",
                                  idSelf, iGroup, pj, p[pj].iOrder, p[pj].iGroup, p[pj].iTouched);*/
                                goto RestartSnake;
                            }
                            p[pj].iTouched = idSelf;
							p[pj].iGroup = iGroup;
                            omp_unset_lock(&locks[_hashLock(kd,pj)]);
#else
							p[pj].iGroup = iGroup;
#endif
							Fifo[iTail++] = pj;
							if (iTail == nFifo) iTail = 0;
							}
						}
					SETNEXT(cp);
					if (cp == ROOT) break;
					continue;
					}
			ContainedCell:
				for (pj=c[cp].pLower;pj<=c[cp].pUpper;++pj) {
#ifdef _OPENMP
                    if (p[pj].iGroup == iGroup) continue;
                    if (p[pj].iTouched >= 0 && p[pj].iTouched < idSelf) {
                        /* Somebody more important that us is already looking at this
                         * group.  Abort this entire group! */
                        //assert(p[pj].iGroup > 0); particle is not locked
                        iHead = iTail;
                        /*printf("Thread %d: Aborting group %d. p[%d].iOrder=%d  p.iGroup=%d  p.iTouched=%d (Per-Cell1)\n",
                          idSelf, iGroup, pj, p[pj].iOrder, p[pj].iGroup, p[pj].iTouched);*/
                        goto RestartSnake;
                    }
#else
					if (p[pj].iGroup) continue;
#endif                    
					/*
					 ** Mark it and add to the do-fifo.
					 */
#ifdef _OPENMP
                    omp_set_lock(&locks[_hashLock(kd,pj)]);
                    if (p[pj].iTouched >= 0 && p[pj].iTouched < idSelf) {
                        /* Check again in case somebody touched it before the lock. */
                        assert(p[pj].iGroup > 0);
                        omp_unset_lock(&locks[_hashLock(kd,pj)]);
                        iHead = iTail;
                        /*printf("Thread %d: Aborting group %d.  p[%d].iGroup=%d  p[%d].iTouched=%d (Per-Cell2)\n",
                          idSelf, iGroup, pj, p[pj].iGroup, pj, p[pj].iTouched);*/
                        goto RestartSnake;
                    }
                    p[pj].iTouched = idSelf;
                    p[pj].iGroup = iGroup;
                    omp_unset_lock(&locks[_hashLock(kd,pj)]);
#else
					p[pj].iGroup = iGroup;
#endif
					Fifo[iTail++] = pj;
					if (iTail == nFifo) iTail = 0;
					}
			GetNextCell:
				SETNEXT(cp);
				if (cp == ROOT) break;
            }
        } /* End while(iHead != iTail) */
#ifdef _OPENMP
    RestartSnake:
#endif
        assert(iHead == iTail);
    }
	free(Fifo);
  }  /* End of the OpenMP PARALLEL section */

#ifdef _OPENMP
    /* Now we have count how many groups there are.  This is straightforward,
     * since the number of groups is the number of particles whose groupID equals
     * their particleID+1. */
    pj = 0;
	for (pn=0;pn<kd->nActive;++pn)
        if (p[pn].iGroup == pn+1) ++pj;
    kd->nGroup = (kd->nActive)+1;
    free(locks);
#else
	kd->nGroup = iGroup+1;
#endif
	return(kd->nGroup-1);
	}


int kdTooSmall(KD kd,int nMembers)
{
	int *pnMembers,*pMap;
	int i,pi,nGroup;

	pnMembers = (int *)malloc(kd->nGroup*sizeof(int));
	assert(pnMembers != NULL);
	pMap = (int *)malloc(kd->nGroup*sizeof(int));
	assert(pMap != NULL);
	for (i=0;i<kd->nGroup;++i) pnMembers[i] = 0;
	for (pi=0;pi<kd->nActive;++pi) {
		++pnMembers[kd->p[pi].iGroup];
		}
	for (i=1;i<kd->nGroup;++i) {
		if (pnMembers[i] < nMembers) {
			pnMembers[i] = 0;
			}
		}
	/*
	 ** Create a remapping!
	 */
	pMap[0] = 0;
	nGroup = 1;
	for (i=1;i<kd->nGroup;++i) {
		pMap[i] = nGroup;
		if (pnMembers[i] == 0) {
			pMap[i] = 0;
			}
		else {
			++nGroup;
			}
		}
	/*
	 ** Remap the groups.
	 */
	for (pi=0;pi<kd->nActive;++pi) {
		kd->p[pi].iGroup = pMap[kd->p[pi].iGroup];
		}
	free(pMap);
	free(pnMembers);
	kd->nGroup = nGroup;
	return(nGroup-1);
	}


int CmpParticles(const void *v1,const void *v2)
{
	PARTICLE *p1 = (PARTICLE *)v1;
	PARTICLE *p2 = (PARTICLE *)v2;
	return(p1->iOrder - p2->iOrder);
	}

void kdOrder(KD kd)
{
	qsort(kd->p,kd->nActive,sizeof(PARTICLE),CmpParticles);
	}


void kdOutGroup(KD kd,char *pszFile)
{
	FILE *fp;
	int i,iCnt;

	fp = fopen(pszFile,"w");
	assert(fp != NULL);
	fprintf(fp,"%d\n",kd->nParticles);
	iCnt = 0;
	for (i=0;i<kd->nGas;++i) {
		if (kd->bGas) fprintf(fp,"%d\n",kd->p[iCnt++].iGroup);
		else fprintf(fp,"0\n");
		}
	for (i=0;i<kd->nDark;++i) {
		if (kd->bDark) fprintf(fp,"%d\n",kd->p[iCnt++].iGroup);
		else fprintf(fp,"0\n");
		}
	for (i=0;i<kd->nStar;++i) {
		if (kd->bStar) fprintf(fp,"%d\n",kd->p[iCnt++].iGroup);
		else fprintf(fp,"0\n");
		}
	fclose(fp);
	}


typedef struct GroupStats {
	double m;
	double r[3];
	double v[3];
	double rel[3];
	double rm;
	} GROUP_STAT;


void kdOutGTP(KD kd,char *pszFile,int bStandard)
{
	FILE *fp;
	GROUP_STAT *grp;
	int pi,i,j;
	struct dump h;
	struct star_particle sp;
	double d,d2;
	XDR xdrs;

	fp = fopen(pszFile,"w");
	assert(fp != NULL);
	grp = malloc(kd->nGroup*sizeof(GROUP_STAT));
	assert(grp != NULL);
	for (i=1;i<kd->nGroup;++i) {
		for (j=0;j<3;++j) grp[i].r[j] = 0.0;
		for (j=0;j<3;++j) grp[i].v[j] = 0.0;
		grp[i].m = 0.0;
		grp[i].rm = 0.0;
		}
	for (pi=0;pi<kd->nActive;++pi) {
		i = kd->p[pi].iGroup;
		if (!i) continue;
		for (j=0;j<3;++j) {
		    grp[i].rel[j] = kd->p[pi].r[j];
			}
		}
	for (pi=0;pi<kd->nActive;++pi) {
		i = kd->p[pi].iGroup;
		if (!i) continue;
		grp[i].m += kd->p[pi].fMass;
		for (j=0;j<3;++j) {
			d = kd->p[pi].r[j] - grp[i].rel[j];
			if (d > 0.5*kd->fPeriod[j]) d -= kd->fPeriod[j];
			if (d <= -0.5*kd->fPeriod[j]) d += kd->fPeriod[j];
			grp[i].r[j] += kd->p[pi].fMass*d;
			}
		for (j=0;j<3;++j) grp[i].v[j] += kd->p[pi].fMass*kd->p[pi].v[j];
		}
	for (i=1;i<kd->nGroup;++i) {
		for (j=0;j<3;++j) {
			grp[i].r[j] /= grp[i].m;
			grp[i].r[j] += grp[i].rel[j];
			if (grp[i].r[j] > kd->fCenter[j]+0.5*kd->fPeriod[j])
				grp[i].r[j] -= kd->fPeriod[j];
			if (grp[i].r[j] <= kd->fCenter[j]-0.5*kd->fPeriod[j])
				grp[i].r[j] += kd->fPeriod[j];
			}
		for (j=0;j<3;++j) grp[i].v[j] /= grp[i].m;
		}
	for (pi=0;pi<kd->nActive;++pi) {
		i = kd->p[pi].iGroup;
		if (!i) continue;
		d2 = 0.0;
		for (j=0;j<3;++j) {
			 d = kd->p[pi].r[j] - grp[i].r[j]; 
			 if (d > 0.5*kd->fPeriod[j]) d -= kd->fPeriod[j];
			 if (d <= -0.5*kd->fPeriod[j]) d += kd->fPeriod[j];
			 d2 += d*d;
			 }
		if (d2 > grp[i].rm) grp[i].rm = d2;
		}
	h.time = kd->fTime;
	h.nbodies = kd->nGroup-1;
	h.nsph = 0;
	h.ndark = 0;
	h.nstar = h.nbodies;
	h.ndim = 3;
	
	if (bStandard) {
		assert(sizeof(Real)==sizeof(float)); /* Else this XDR stuff
							ain't gonna work */
		xdrstdio_create(&xdrs, fp, XDR_ENCODE);
		xdrHeader(&xdrs,&h);
	        }
	else {
		fwrite(&h,sizeof(struct dump),1,fp);
	        }
	for (i=1;i<kd->nGroup;++i) {
		sp.mass = grp[i].m;
		for (j=0;j<3;++j) sp.pos[j] = grp[i].r[j];
		for (j=0;j<3;++j) sp.vel[j] = grp[i].v[j];
		sp.eps = sqrt(grp[i].rm);
		sp.tform = kd->fTime;
		sp.metals = 0.0;
		sp.phi = 0.0;
		if (bStandard) {
			xdr_vector(&xdrs, (char *) &sp,
				   sizeof(struct star_particle)/sizeof(Real),
				   sizeof(Real), (xdrproc_t)xdr_float);
		        }
		else {
			fwrite(&sp,sizeof(struct star_particle),1,fp);
		        }
		}
	if (bStandard) xdr_destroy(&xdrs);
	free(grp);
	fclose(fp);
	}


void kdFinish(KD kd)
{
	free(kd->p);
	free(kd->kdNodes);
	free(kd);
	}

