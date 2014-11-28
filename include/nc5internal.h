#ifndef _NC5INTERNAL_
#define _NC5INTERNAL_

/*
 *	netcdf library 'private' data structures, objects and interfaces
 */
#include <config.h>
#include <stddef.h>	/* size_t */
#ifndef HAVE_STDINT_H
#  include "pstdint.h"	/* attempts to define uint32_t etc portably */
#else
#  include <stdint.h>
#endif /* HAVE_STDINT_H */
#include <sys/types.h>	/* off_t */
#ifdef USE_PARALLEL
#include <netcdf_par.h>
#else
#include <netcdf.h>
#endif /* USE_PARALLEL */

/* Always needed */
#include "nc.h"

#include <nc3internal.h>

/*#ifndef HAVE_SSIZE_T
#define ssize_t int
#endif*/

#ifndef NC_ARRAY_GROWBY
#define NC_ARRAY_GROWBY 4
#endif

/*
 * The extern size of an empty
 * netcdf version 1 file.
 * The initial value of ncp->xsz.
 */
#ifndef MIN_NC_XSZ
#define MIN_NC_XSZ 32
#endif
#define MIN_NC5_XSZ 48

/* Forward */
typedef struct NC5_INFO NC5_INFO;

/* Begin defined in nc5dim.c */

extern void
nc5i_free_NC_dim(NC_dim *dimp);

extern NC_dim *
nc5i_new_x_NC_dim(NC_string *name);

extern int
nc5i_find_NC_Udim(const NC_dimarray *ncap, NC_dim **dimpp);

/* dimarray */

extern void
nc5i_free_NC_dimarrayV(NC_dimarray *ncap);

extern int
nc5i_dup_NC_dimarrayV(NC_dimarray *ncap, const NC_dimarray *ref);

extern NC_dim *
nc5i_elem_NC_dimarray(const NC_dimarray *ncap, size_t elem);

/* End defined in nc5dim.c */

/* Begin defined in nc5attr.c */

extern void
nc5i_free_NC_attr(NC_attr *attrp);

extern NC_attr *
nc5i_new_x_NC_attr(NC_string *strp, nc_type type, size_t nelems);

extern NC_attr **
NC5_findattr(const NC_attrarray *ncap, const char *name);

/* attrarray */

extern void
nc5i_free_NC_attrarrayV0(NC_attrarray *ncap);

extern void
nc5i_free_NC_attrarrayV(NC_attrarray *ncap);

extern int
nc5i_dup_NC_attrarrayV(NC_attrarray *ncap, const NC_attrarray *ref);

/* End defined in nc5attr.c */


/*
 * NC variable: description and data
 */

/* Begin defined in nc5var.c */

extern void
nc5i_free_NC_var(NC_var *varp);

extern NC_var *
nc5i_new_x_NC_var(NC_string *strp, size_t ndims);

/* vararray */

extern void
nc5i_free_NC_vararrayV(NC_vararray *ncap);

extern int
nc5i_dup_NC_vararrayV(NC_vararray *ncap, const NC_vararray *ref);

extern int
NC5_var_shape(NC_var *varp, const NC_dimarray *dims);

extern int
NC5_check_vlen(NC_var *varp, size_t vlen_max);

extern NC_var *
NC5_lookupvar(NC5_INFO* ncp, int varid);

/* End defined in nc5var.c */

struct NC5_INFO {
	/* contains the previous NC during redef. */
	NC5_INFO *old;
	/* flags */
#define NC_CREAT 2	/* in create phase, cleared by ncendef */
#define NC_INDEF 8	/* in define mode, cleared by ncendef */
#define NC_NSYNC 0x10	/* synchronise numrecs on change */
#define NC_HSYNC 0x20	/* synchronise whole header on change */
#define NC_NDIRTY 0x40	/* numrecs has changed */
#define NC_HDIRTY 0x80  /* header info has changed */
/*	NC_NOFILL in netcdf.h, historical interface */
	int flags;
	struct ncio* nciop;
	size_t chunk;	/* largest extent this layer will request from ncio->get() */
	size_t xsz;	/* external size of this header, == var[0].begin */
	off_t begin_var; /* position of the first (non-record) var */
	off_t begin_rec; /* position of the first 'record' */
        /* Don't constrain maximum size of record unnecessarily */
#if SIZEOF_OFF_T > SIZEOF_SIZE_T
        off_t recsize;   /* length of 'record' */
#else
	size_t recsize;  /* length of 'record' */
#endif
	/* below gets xdr'd */
	size_t numrecs; /* number of 'records' allocated */
	NC_dimarray dims;
	NC_attrarray attrs;
	NC_vararray vars;
#ifdef LOCKNUMREC
/* size and named indexes for the lock array protecting NC.numrecs */
#  define LOCKNUMREC_DIM	4
#  define LOCKNUMREC_VALUE	0
#  define LOCKNUMREC_LOCK	1
#  define LOCKNUMREC_SERVING	2
#  define LOCKNUMREC_BASEPE	3
	/* Used on Cray T3E MPP to maintain the
	 * integrity of numrecs for an unlimited dimension
	 */
	ushmem_t lock[LOCKNUMREC_DIM];
#endif
   int pnetcdf_access_mode;
   int use_parallel;
};

#ifndef LOCKNUMREC
#  define NC5_get_numrecs(nc5i) \
	((nc5i)->numrecs)

#  define NC5_set_numrecs(nc5i, nrecs) \
	{(nc5i)->numrecs = (nrecs);}

#  define NC5_increase_numrecs(nc5i, nrecs) \
	{if((nrecs) > (nc5i)->numrecs) ((nc5i)->numrecs = (nrecs));}
#else
	size_t NC5_get_numrecs(const NC5_INFO *nc5i);
	void   NC5_set_numrecs(NC5_INFO *nc5i, size_t nrecs);
	void   NC5_increase_numrecs(NC5_INFO *nc5i, size_t nrecs);
#endif

/* Begin defined in nc5internal.c */

extern int
nc5_cktype(int mode, nc_type datatype);

extern size_t
nc5i_howmany(nc_type type, size_t xbufsize);

extern int
nc5i_read_numrecs(NC5_INFO* ncp);

extern int
nc5i_write_numrecs(NC5_INFO* ncp);

extern int
nc5i_sync(NC5_INFO* ncp);

/* End defined in nc5internal.c */

/* Begin defined in nc5v1hpg.c */

extern size_t
nc5i_len_NC(const NC5_INFO* ncp, size_t sizeof_off_t);

extern int
nc5i_put_NC(const NC5_INFO* ncp, void **xpp, off_t offset, size_t extent);

extern int
nc5_get_NC(NC5_INFO* ncp);

/* End defined in nc5v1hpg.c */
/* Begin defined in nc5putget.c */

extern int
nc5i_fill_NC_var(NC5_INFO* ncp, const NC_var *varp, size_t varsize, size_t recno);

/* End defined in nc5putget.c */

/* Define accessors for the dispatchdata */
#define NC5_DATA(nc) ((NC5_INFO*)(nc)->dispatchdata)
#define NC5_DATA_SET(nc,data) ((nc)->dispatchdata = (void*)(data))

#endif /* _NC5INTERNAL_ */
