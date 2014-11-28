/* WARNING: Order of mpi.h, nc.h, and pnetcdf.h is important */
#include "config.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#if defined(LOCKNUMREC) /* && _CRAYMPP */
#  include <mpp/shmem.h>
#  include <intrinsics.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef USE_PARALLEL
#include <mpi.h>
#endif

#include "nc5internal.h"
#include "ncdispatch.h"
#include "nc5dispatch.h"
#include "rnd.h"
#include "nc5ncx.h"

#ifdef USE_PNETCDF
/* Must follow netcdf.h */
#include <pnetcdf.h>
#endif

/* Define accessors for the dispatchdata */
#define NC5_DATA(nc) ((NC5_INFO*)(nc)->dispatchdata)
#define NC5_DATA_SET(nc,data) ((nc)->dispatchdata = (void*)(data))

/* Cannot have NC_MPIPOSIX flag, ignore NC_MPIIO as PnetCDF use MPIIO */
#define LEGAL_CREATE_FLAGS (NC_NOCLOBBER | NC_64BIT_OFFSET | NC_CLASSIC_MODEL | NC_SHARE | NC_MPIIO | NC_LOCK | NC_PNETCDF | NC_64BIT_DATA)

#define LEGAL_OPEN_FLAGS (NC_WRITE | NC_NOCLOBBER | NC_SHARE | NC_MPIIO | NC_LOCK | NC_PNETCDF | NC_CLASSIC_MODEL | NC_64BIT_OFFSET | NC_64BIT_DATA)


/* These have to do with version numbers. */
#define MAGIC_NUM_LEN 4
#define VER_CLASSIC 1
#define VER_64BIT_OFFSET 2
#define VER_HDF5 3

static void
free_NC5INFO(NC5_INFO *nc5)
{
	if(nc5 == NULL)
		return;
	nc5i_free_NC_dimarrayV(&nc5->dims);
	nc5i_free_NC_attrarrayV(&nc5->attrs);
	nc5i_free_NC_vararrayV(&nc5->vars);
#if _CRAYMPP && defined(LOCKNUMREC)
	shfree(nc5);
#else
	free(nc5);
#endif /* _CRAYMPP && LOCKNUMREC */
}

static NC5_INFO *
new_NC5INFO(const size_t *chunkp)
{
	NC5_INFO *ncp;
	ncp = (NC5_INFO*)calloc(1,sizeof(NC5_INFO));
	if(ncp == NULL) return ncp;
/* delay this setting
	ncp->xsz = MIN_NC_XSZ;
	assert(ncp->xsz == nc5i_len_NC(ncp,0));
*/
        ncp->chunk = chunkp != NULL ? *chunkp : NC_SIZEHINT_DEFAULT;
	return ncp;
}

static NC5_INFO *
dup_NC5INFO(const NC5_INFO *ref)
{
	NC5_INFO *ncp;
	ncp = (NC5_INFO*)calloc(1,sizeof(NC5_INFO));
	if(ncp == NULL) return ncp;

	if(nc5i_dup_NC_dimarrayV(&ncp->dims, &ref->dims) != NC_NOERR)
		goto err;
	if(nc5i_dup_NC_attrarrayV(&ncp->attrs, &ref->attrs) != NC_NOERR)
		goto err;
	if(nc5i_dup_NC_vararrayV(&ncp->vars, &ref->vars) != NC_NOERR)
		goto err;

	ncp->xsz = ref->xsz;
	ncp->begin_var = ref->begin_var;
	ncp->begin_rec = ref->begin_rec;
	ncp->recsize = ref->recsize;
	NC5_set_numrecs(ncp, NC5_get_numrecs(ref));
	return ncp;
err:
	free_NC5INFO(ncp);
	return NULL;
}


/*
 *  Verify that this is a user nc_type
 * Formerly
NCcktype()
 * Sense of the return is changed.
 */
int
nc5_cktype(int mode, nc_type type)
{
    if (mode & NC_64BIT_OFFSET) { /* CDF-1 and CDF-2 format */
        if (type >= NC_BYTE && type <= NC_DOUBLE) return NC_NOERR;
	return NC_EBADTYPE;
    }
    else {  /* CDF-5 format */
        if (type >= NC_BYTE && type <= NC_STRING) return NC_NOERR;
	return NC_EBADTYPE;
    }
}


/*
 * How many objects of 'type'
 * will fit into xbufsize?
 */
size_t
nc5i_howmany(nc_type type, size_t xbufsize)
{
	switch(type){
	case NC_BYTE:
	case NC_CHAR:
	case NC_UBYTE:
		return xbufsize;
	case NC_SHORT:
		return xbufsize/X_SIZEOF_SHORT;
	case NC_INT:
		return xbufsize/X_SIZEOF_INT;
	case NC_FLOAT:
		return xbufsize/X_SIZEOF_FLOAT;
	case NC_DOUBLE:
		return xbufsize/X_SIZEOF_DOUBLE;
	case NC_USHORT:
		return xbufsize/X_SIZEOF_USHORT;
	case NC_UINT:
		return xbufsize/X_SIZEOF_UINT;
	case NC_INT64:
		return xbufsize/X_SIZEOF_LONGLONG;
	case NC_UINT64:
		return xbufsize/X_SIZEOF_ULONGLONG;
	default:
	        assert("nc5i_howmany: Bad type" == 0);
		return(0);
	}
}

#define	D_RNDUP(x, align) _RNDUP(x, (off_t)(align))

/*
 * Compute each variable's 'begin' offset,
 * update 'begin_rec' as well.
 */
static int
NC_begins(NC5_INFO* ncp,
	size_t h_minfree, size_t v_align,
	size_t v_minfree, size_t r_align)
{
	size_t ii, j;
	int sizeof_off_t;
	off_t index = 0;
	NC_var **vpp;
	NC_var *last = NULL;
	NC_var *first_var = NULL;       /* first "non-record" var */


	if(v_align == NC_ALIGN_CHUNK)
		v_align = ncp->chunk;
	if(r_align == NC_ALIGN_CHUNK)
		r_align = ncp->chunk;

	if (fIsSet(ncp->flags, NC_64BIT_OFFSET) || fIsSet(ncp->flags, NC_64BIT_DATA)) {
	  sizeof_off_t = 8;
	} else {
	  sizeof_off_t = 4;
	}

	ncp->xsz = nc5i_len_NC(ncp,sizeof_off_t);

	if(ncp->vars.nelems == 0)
		return NC_NOERR;

	/* only (re)calculate begin_var if there is not sufficient space in header
	   or start of non-record variables is not aligned as requested by valign */
	if (ncp->begin_var < ncp->xsz + h_minfree ||
	    ncp->begin_var != D_RNDUP(ncp->begin_var, v_align) )
	{
	  index = (off_t) ncp->xsz;
	  ncp->begin_var = D_RNDUP(index, v_align);
	  if(ncp->begin_var < index + h_minfree)
	  {
	    ncp->begin_var = D_RNDUP(index + (off_t)h_minfree, v_align);
	  }
	}

	if (ncp->old != NULL) {
            /* check whether the new begin_var is smaller */
            if (ncp->begin_var < ncp->old->begin_var)
                ncp->begin_var = ncp->old->begin_var;
	}

	index = ncp->begin_var;

	/* loop thru vars, first pass is for the 'non-record' vars */
	j = 0;
	vpp = ncp->vars.value;
	for(ii = 0; ii < ncp->vars.nelems ; ii++, vpp++)
	{
		if( IS_RECVAR(*vpp) )
		{
			/* skip record variables on this pass */
			continue;
		}
		if (first_var == NULL) first_var = *vpp;

#if 0
fprintf(stderr, "    VAR %d %s: %ld\n", ii, (*vpp)->name->cp, (long)index);
#endif
                if( sizeof_off_t == 4 && (index > X_OFF_MAX || index < 0) )
		{
		    return NC_EVARSIZE;
                }
		(*vpp)->begin = index;

		if (ncp->old != NULL) {
		    /* move to the next fixed variable */
		    for (; j<ncp->old->vars.nelems; j++)
		        if (!IS_RECVAR(ncp->old->vars.value[j]))
		            break;
		    if (j < ncp->old->vars.nelems) {
		        if ((*vpp)->begin < ncp->old->vars.value[j]->begin)
		            /* the first ncp->vars.nelems fixed variables
                               should be the same. If the new begin is smaller,
                               reuse the old begin */
                            (*vpp)->begin = ncp->old->vars.value[j]->begin;
                        j++;
		    }
		}

		index += (*vpp)->len;
	}

	if (ncp->old != NULL) {
	    /* check whether the new begin_rec is smaller */
	    if (ncp->begin_rec < ncp->old->begin_rec)
	        ncp->begin_rec = ncp->old->begin_rec;
	}

	/* only (re)calculate begin_rec if there is not sufficient
	   space at end of non-record variables or if start of record
	   variables is not aligned as requested by r_align */
	if (ncp->begin_rec < index + v_minfree ||
	    ncp->begin_rec != D_RNDUP(ncp->begin_rec, r_align) )
	{
	  ncp->begin_rec = D_RNDUP(index, r_align);
	  if(ncp->begin_rec < index + v_minfree)
	  {
	    ncp->begin_rec = D_RNDUP(index + (off_t)v_minfree, r_align);
	  }
	}

	if (first_var != NULL)
	    ncp->begin_var = first_var->begin;
	else
	    ncp->begin_var = ncp->begin_rec;

	index = ncp->begin_rec;

	ncp->recsize = 0;

	/* loop thru vars, second pass is for the 'record' vars */
	j = 0;
	vpp = (NC_var **)ncp->vars.value;
	for(ii = 0; ii < ncp->vars.nelems; ii++, vpp++)
	{
		if( !IS_RECVAR(*vpp) )
		{
			/* skip non-record variables on this pass */
			continue;
		}

#if 0
fprintf(stderr, "    REC %d %s: %ld\n", ii, (*vpp)->name->cp, (long)index);
#endif
                if( sizeof_off_t == 4 && (index > X_OFF_MAX || index < 0) )
		{
		    return NC_EVARSIZE;
                }
		(*vpp)->begin = index;

                if (ncp->old != NULL) {
                    /* move to the next record variable */
                    for (; j<ncp->old->vars.nelems; j++)
                        if (IS_RECVAR(ncp->old->vars.value[j]))
                            break;
                    if (j < ncp->old->vars.nelems) {
                        if ((*vpp)->begin < ncp->old->vars.value[j]->begin)
                            /* if the new begin is smaller, use the old begin */
                            (*vpp)->begin = ncp->old->vars.value[j]->begin;
                        j++;
                    }
                }

		index += (*vpp)->len;
		/* check if record size must fit in 32-bits */
#if SIZEOF_OFF_T == SIZEOF_SIZE_T && SIZEOF_SIZE_T == 4
		if( ncp->recsize > X_UINT_MAX - (*vpp)->len )
		{
		    return NC_EVARSIZE;
		}
#endif
		if((*vpp)->len != UINT32_MAX) /* flag for vars >= 2**32 bytes */
		    ncp->recsize += (*vpp)->len;
		last = (*vpp);
	}

	/*
	 * for special case of
	 */
	if(last != NULL) {
	    if(ncp->recsize == last->len) { /* exactly one record variable, pack value */
		ncp->recsize = *last->dsizes * last->xsz;
	    } else if(last->len == UINT32_MAX) { /* huge last record variable */
		ncp->recsize += *last->dsizes * last->xsz;
	    }
	}
	if(NC_IsNew(ncp))
		NC5_set_numrecs(ncp, 0);
	return NC_NOERR;
}


/*
 * Read just the numrecs member.
 * (A relatively expensive way to do things.)
 */
int
nc5i_read_numrecs(NC5_INFO *ncp)
{
	int status = NC_NOERR;
	const void *xp = NULL;
	size_t new_nrecs, old_nrecs = NC5_get_numrecs(ncp);
	size_t nc_numrecs_extent=4; /* CDF-1 and CDF-2 */

	assert(!NC_indef(ncp));

	if (fIsSet(ncp->flags, NC_64BIT_DATA))
		nc_numrecs_extent = 8; /* CDF-5 */

#define NC_NUMRECS_OFFSET 4
	status = ncio_get(ncp->nciop,
		 NC_NUMRECS_OFFSET, nc_numrecs_extent, 0, (void **)&xp);
					/* cast away const */
	if(status != NC_NOERR)
		return status;

	if (fIsSet(ncp->flags, NC_64BIT_DATA)) {
		long long tmp=0;
		status = nc5x_get_int64(&xp, &tmp);
		new_nrecs = tmp;
        }
	else
		status = nc5x_get_size_t(&xp, &new_nrecs);

	(void) ncio_rel(ncp->nciop, NC_NUMRECS_OFFSET, 0);

	if(status == NC_NOERR && old_nrecs != new_nrecs)
	{
		NC5_set_numrecs(ncp, new_nrecs);
		fClr(ncp->flags, NC_NDIRTY);
	}

	return status;
}


/*
 * Write out just the numrecs member.
 * (A relatively expensive way to do things.)
 */
int
nc5i_write_numrecs(NC5_INFO *ncp)
{
	int status = NC_NOERR;
	void *xp = NULL;
	size_t nc_numrecs_extent=4; /* CDF-1 and CDF-2 */

	assert(!NC_readonly(ncp));
	assert(!NC_indef(ncp));

	if (fIsSet(ncp->flags, NC_64BIT_DATA))
		nc_numrecs_extent = 8; /* CDF-5 */

	status = ncio_get(ncp->nciop,
		 NC_NUMRECS_OFFSET, nc_numrecs_extent, RGN_WRITE, &xp);
	if(status != NC_NOERR)
		return status;

	{
		const size_t nrecs = NC5_get_numrecs(ncp);
		if (fIsSet(ncp->flags, NC_64BIT_DATA))
			status = nc5x_put_int64(&xp, nrecs);
		else
			status = nc5x_put_size_t(&xp, &nrecs);
	}

	(void) ncio_rel(ncp->nciop, NC_NUMRECS_OFFSET, RGN_MODIFIED);

	if(status == NC_NOERR)
		fClr(ncp->flags, NC_NDIRTY);

	return status;
}


/*
 * Read in the header
 * It is expensive.
 */
static int
read_NC(NC5_INFO *ncp)
{
	int status = NC_NOERR;

	nc5i_free_NC_dimarrayV(&ncp->dims);
	nc5i_free_NC_attrarrayV(&ncp->attrs);
	nc5i_free_NC_vararrayV(&ncp->vars);

	status = nc5_get_NC(ncp);

	if(status == NC_NOERR)
		fClr(ncp->flags, NC_NDIRTY | NC_HDIRTY);

	return status;
}


/*
 * Write out the header
 */
static int
write_NC(NC5_INFO *ncp)
{
	int status = NC_NOERR;

	assert(!NC_readonly(ncp));

	status = nc5i_put_NC(ncp, NULL, 0, 0);

	if(status == NC_NOERR)
		fClr(ncp->flags, NC_NDIRTY | NC_HDIRTY);

	return status;
}


/*
 * Write the header or the numrecs if necessary.
 */
int
nc5i_sync(NC5_INFO *ncp)
{
	assert(!NC_readonly(ncp));

	if(NC_hdirty(ncp))
	{
		return write_NC(ncp);
	}
	/* else */

	if(NC_ndirty(ncp))
	{
		return nc5i_write_numrecs(ncp);
	}
	/* else */

	return NC_NOERR;
}


/*
 * Initialize the 'non-record' variables.
 */
static int
fillerup(NC5_INFO *ncp)
{
	int status = NC_NOERR;
	size_t ii;
	NC_var **varpp;

	assert(!NC_readonly(ncp));
	assert(NC_dofill(ncp));

	/* loop thru vars */
	varpp = ncp->vars.value;
	for(ii = 0; ii < ncp->vars.nelems; ii++, varpp++)
	{
		if(IS_RECVAR(*varpp))
		{
			/* skip record variables */
			continue;
		}

		status = nc5i_fill_NC_var(ncp, *varpp, (*varpp)->len, 0);
		if(status != NC_NOERR)
			break;
	}
	return status;
}

/* Begin endef */

/*
 */
static int
fill_added_recs(NC5_INFO *gnu, NC5_INFO *old)
{
	NC_var ** const gnu_varpp = (NC_var **)gnu->vars.value;

	const int old_nrecs = (int) NC5_get_numrecs(old);
	int recno = 0;
	NC_var **vpp = gnu_varpp;
	NC_var *const *const end = &vpp[gnu->vars.nelems];
	int numrecvars = 0;

	/* Determine if there is only one record variable.  If so, we
	   must treat as a special case because there's no record padding */
	for(; vpp < end; vpp++) {
	    if(IS_RECVAR(*vpp)) {
		numrecvars++;
	    }
	}

	for(; recno < old_nrecs; recno++)
	    {
		int varid = (int)old->vars.nelems;
		for(; varid < (int)gnu->vars.nelems; varid++)
		    {
			const NC_var *const gnu_varp = *(gnu_varpp + varid);
			if(!IS_RECVAR(gnu_varp))
			    {
				/* skip non-record variables */
				continue;
			    }
			/* else */
			{
			    size_t varsize = numrecvars == 1 ? gnu->recsize :  gnu_varp->len;
			    const int status = nc5i_fill_NC_var(gnu, gnu_varp, varsize, recno);
			    if(status != NC_NOERR)
				return status;
			}
		    }
	    }
	return NC_NOERR;
}

/*
 */
static int
fill_added(NC5_INFO *gnu, NC5_INFO *old)
{
	NC_var ** const gnu_varpp = (NC_var **)gnu->vars.value;
	int varid = (int)old->vars.nelems;

	for(; varid < (int)gnu->vars.nelems; varid++)
	{
		const NC_var *const gnu_varp = *(gnu_varpp + varid);
		if(IS_RECVAR(gnu_varp))
		{
			/* skip record variables */
			continue;
		}
		/* else */
		{
		const int status = nc5i_fill_NC_var(gnu, gnu_varp, gnu_varp->len, 0);
		if(status != NC_NOERR)
			return status;
		}
	}

	return NC_NOERR;
}


/*
 * Move the records "out".
 * Fill as needed.
 */
static int
move_recs_r(NC5_INFO *gnu, NC5_INFO *old)
{
	int status;
	int recno;
	int varid;
	NC_var **gnu_varpp = (NC_var **)gnu->vars.value;
	NC_var **old_varpp = (NC_var **)old->vars.value;
	NC_var *gnu_varp;
	NC_var *old_varp;
	off_t gnu_off;
	off_t old_off;
	const size_t old_nrecs = NC5_get_numrecs(old);

	/* Don't parallelize this loop */
	for(recno = (int)old_nrecs -1; recno >= 0; recno--)
	{
	/* Don't parallelize this loop */
	for(varid = (int)old->vars.nelems -1; varid >= 0; varid--)
	{
		gnu_varp = *(gnu_varpp + varid);
		if(!IS_RECVAR(gnu_varp))
		{
			/* skip non-record variables on this pass */
			continue;
		}
		/* else */

		/* else, a pre-existing variable */
		old_varp = *(old_varpp + varid);
		gnu_off = gnu_varp->begin + (off_t)(gnu->recsize * recno);
		old_off = old_varp->begin + (off_t)(old->recsize * recno);

		if(gnu_off == old_off)
			continue; 	/* nothing to do */

		assert(gnu_off > old_off);

		status = ncio_move(gnu->nciop, gnu_off, old_off,
			 old_varp->len, 0);

		if(status != NC_NOERR)
			return status;

	}
	}

	NC5_set_numrecs(gnu, old_nrecs);

	return NC_NOERR;
}


/*
 * Move the "non record" variables "out".
 * Fill as needed.
 */
static int
move_vars_r(NC5_INFO *gnu, NC5_INFO *old)
{
	int err, status=NC_NOERR;
	int varid;
	NC_var **gnu_varpp = (NC_var **)gnu->vars.value;
	NC_var **old_varpp = (NC_var **)old->vars.value;
	NC_var *gnu_varp;
	NC_var *old_varp;
	off_t gnu_off;
	off_t old_off;

	/* Don't parallelize this loop */
	for(varid = (int)old->vars.nelems -1;
		 varid >= 0; varid--)
	{
		gnu_varp = *(gnu_varpp + varid);
		if(IS_RECVAR(gnu_varp))
		{
			/* skip record variables on this pass */
			continue;
		}
		/* else */

		old_varp = *(old_varpp + varid);
		gnu_off = gnu_varp->begin;
		old_off = old_varp->begin;

		if (gnu_off > old_off) {
		    err = ncio_move(gnu->nciop, gnu_off, old_off,
			               old_varp->len, 0);
		    if (status == NC_NOERR) status = err;
		}
	}
	return status;
}


/*
 * Given a valid ncp, return NC_EVARSIZE if any variable has a bad len
 * (product of non-rec dim sizes too large), else return NC_NOERR.
 */
static int
NC_check_vlens(NC5_INFO *ncp)
{
    NC_var **vpp;
    /* maximum permitted variable size (or size of one record's worth
       of a record variable) in bytes.  This is different for format 1
       and format 2. */
    size_t vlen_max;
    size_t ii;
    size_t large_vars_count;
    size_t rec_vars_count;
    int last = 0;

    if(ncp->vars.nelems == 0)
	return NC_NOERR;

    if (ncp->flags & NC_64BIT_DATA) {
	/* CDF5 format allows many large vars */
        return NC_NOERR;
    }
    else if ((ncp->flags & NC_64BIT_OFFSET) && sizeof(off_t) > 4) {
	/* CDF2 format and LFS */
	vlen_max = X_UINT_MAX - 3; /* "- 3" handles rounded-up size */
    } else {
	/* CDF1 format */
	vlen_max = X_INT_MAX - 3;
    }
    /* Loop through vars, first pass is for non-record variables.   */
    large_vars_count = 0;
    rec_vars_count = 0;
    vpp = ncp->vars.value;
    for (ii = 0; ii < ncp->vars.nelems; ii++, vpp++) {
	if( !IS_RECVAR(*vpp) ) {
	    last = 0;
	    if( NC5_check_vlen(*vpp, vlen_max) == 0 ) {
		large_vars_count++;
		last = 1;
	    }
	} else {
	  rec_vars_count++;
	}
    }
    /* OK if last non-record variable size too large, since not used to
       compute an offset */
    if( large_vars_count > 1) { /* only one "too-large" variable allowed */
      return NC_EVARSIZE;
    }
    /* and it has to be the last one */
    if( large_vars_count == 1 && last == 0) {
      return NC_EVARSIZE;
    }
    if( rec_vars_count > 0 ) {
	/* and if it's the last one, there can't be any record variables */
	if( large_vars_count == 1 && last == 1) {
	    return NC_EVARSIZE;
	}
	/* Loop through vars, second pass is for record variables.   */
	large_vars_count = 0;
	vpp = ncp->vars.value;
	for (ii = 0; ii < ncp->vars.nelems; ii++, vpp++) {
	    if( IS_RECVAR(*vpp) ) {
		last = 0;
		if( NC5_check_vlen(*vpp, vlen_max) == 0 ) {
		    large_vars_count++;
		    last = 1;
		}
	    }
	}
	/* OK if last record variable size too large, since not used to
	   compute an offset */
	if( large_vars_count > 1) { /* only one "too-large" variable allowed */
	    return NC_EVARSIZE;
	}
	/* and it has to be the last one */
	if( large_vars_count == 1 && last == 0) {
	    return NC_EVARSIZE;
	}
    }
    return NC_NOERR;
}


/*
 *  End define mode.
 *  Common code for ncendef, ncclose(endef)
 *  Flushes I/O buffers.
 */
static int
NC_endef(NC5_INFO *ncp,
	size_t h_minfree, size_t v_align,
	size_t v_minfree, size_t r_align)
{
	int status = NC_NOERR;

	assert(!NC_readonly(ncp));
	assert(NC_indef(ncp));

	status = NC_check_vlens(ncp);
	if(status != NC_NOERR)
	    return status;
	status = NC_begins(ncp, h_minfree, v_align, v_minfree, r_align);
	if(status != NC_NOERR)
	    return status;

	if(ncp->old != NULL)
	{
		/* a plain redef, not a create */
		assert(!NC_IsNew(ncp));
		assert(fIsSet(ncp->flags, NC_INDEF));
		assert(ncp->begin_rec >= ncp->old->begin_rec);
		assert(ncp->begin_var >= ncp->old->begin_var);

		if(ncp->vars.nelems != 0)
		{
		if(ncp->begin_rec > ncp->old->begin_rec)
		{
			status = move_recs_r(ncp, ncp->old);
			if(status != NC_NOERR)
				return status;
			if(ncp->begin_var > ncp->old->begin_var)
			{
				status = move_vars_r(ncp, ncp->old);
				if(status != NC_NOERR)
					return status;
			}
			/* else if (ncp->begin_var == ncp->old->begin_var) { NOOP } */
		}
		else
                {
			/* due to fixed variable alignment, it is possible that header
                           grows but begin_rec did not change */
			if(ncp->begin_var > ncp->old->begin_var)
			{
				status = move_vars_r(ncp, ncp->old);
				if(status != NC_NOERR)
					return status;
			}
		 	/* Even if (ncp->begin_rec == ncp->old->begin_rec)
			   and     (ncp->begin_var == ncp->old->begin_var)
			   might still have added a new record variable */
		        if(ncp->recsize > ncp->old->recsize)
			{
			        status = move_recs_r(ncp, ncp->old);
				if(status != NC_NOERR)
				      return status;
			}
		}
		}
	}

	status = write_NC(ncp);
	if(status != NC_NOERR)
		return status;

	if(NC_dofill(ncp))
	{
		if(NC_IsNew(ncp))
		{
			status = fillerup(ncp);
			if(status != NC_NOERR)
				return status;

		}
		else if(ncp->old == (void*)0 ?
                0 : (ncp->vars.nelems > ncp->old->vars.nelems))
          {
            status = fill_added(ncp, ncp->old);
            if(status != NC_NOERR)
              return status;
            status = fill_added_recs(ncp, ncp->old);
            if(status != NC_NOERR)
              return status;
          }
	}

	if(ncp->old != NULL)
	{
		free_NC5INFO(ncp->old);
		ncp->old = NULL;
	}

	fClr(ncp->flags, NC_CREAT | NC_INDEF);

	return ncio_sync(ncp->nciop);
}

#ifdef LOCKNUMREC
static int
NC_init_pe(NC *ncp, int basepe) {
	if (basepe < 0 || basepe >= _num_pes()) {
		return NC_EINVAL; /* invalid base pe */
	}
	/* initialize common values */
	ncp->lock[LOCKNUMREC_VALUE] = 0;
	ncp->lock[LOCKNUMREC_LOCK] = 0;
	ncp->lock[LOCKNUMREC_SERVING] = 0;
	ncp->lock[LOCKNUMREC_BASEPE] =  basepe;
	return NC_NOERR;
}
#endif

/*
 * Compute the expected size of the file.
 */
static int
nc5i_NC_calcsize(const NC5_INFO *ncp, off_t *calcsizep)
{
	NC_var **vpp = (NC_var **)ncp->vars.value;
	NC_var *const *const end = &vpp[ncp->vars.nelems];
	NC_var *last_fix = NULL;	/* last "non-record" var */
	int numrecvars = 0;	/* number of record variables */

	if(ncp->vars.nelems == 0) { /* no non-record variables and
				       no record variables */
	    *calcsizep = ncp->xsz; /* size of header */
	    return NC_NOERR;
	}

	for( /*NADA*/; vpp < end; vpp++) {
	    if(IS_RECVAR(*vpp)) {
		numrecvars++;
	    } else {
		last_fix = *vpp;
	    }
	}

	if(numrecvars == 0) {
	    off_t varsize;
	    assert(last_fix != NULL);
	    varsize = last_fix->len;
	    if(last_fix->len == X_UINT_MAX) { /* huge last fixed var */
          int i;
          varsize = 1;
          for(i = 0; i < last_fix->ndims && last_fix->shape != NULL; i++ ) {
            varsize *= last_fix->shape[i];
          }
	    }
	    *calcsizep = last_fix->begin + varsize;
	    /*last_var = last_fix;*/
	} else {       /* we have at least one record variable */
      *calcsizep = ncp->begin_rec + ncp->numrecs * ncp->recsize;
	}

	return NC_NOERR;
}

/**************************************************/

static int
nc5_create_file(const char *path, int ioflags,
		size_t initialsz, int basepe,
		size_t *chunksizehintp,
		int use_parallel, void* parameters,
                NC_Dispatch* dispatch, NC* nc)
{
	int status, default_format;
	void *xp = NULL;
	int sizeof_off_t = 0;
	NC5_INFO* nc5;

	/* Create our specific NC5_INFO instance */
	nc5 = new_NC5INFO(chunksizehintp);

#if ALWAYS_NC_SHARE /* DEBUG */
	fSet(ioflags, NC_SHARE);
#endif

#if defined(LOCKNUMREC) /* && _CRAYMPP */
	if (status = NC_init_pe(nc5, basepe)) {
		return status;
	}
#else
	/*
	 * !_CRAYMPP, only pe 0 is valid
	 */
	if(basepe != 0) {
		free(nc5);
		return NC_EINVAL;
	}
#endif

	assert(nc5->flags == 0);

	if (!fIsSet(ioflags, NC_64BIT_DATA) && !fIsSet(ioflags, NC_64BIT_OFFSET)) {
            /* if no format set in ioflag, apply default create format. */
            default_format = nc_get_default_format();
            if (default_format == NC_FORMAT_CDF2)
              ioflags |= NC_64BIT_OFFSET;
            else if (default_format == NC_FORMAT_CDF5)
              ioflags |= NC_64BIT_DATA;
        }

	nc5->xsz = MIN_NC_XSZ;
	if (fIsSet(ioflags, NC_64BIT_DATA)) {
		fSet(nc5->flags, NC_64BIT_DATA);
		sizeof_off_t = 8;
		nc5->xsz = MIN_NC5_XSZ; /* CDF-5 has minimum 16 extra bytes */
	} else if (fIsSet(ioflags, NC_64BIT_OFFSET)) {
		fSet(nc5->flags, NC_64BIT_OFFSET);
		sizeof_off_t = 8;
	} else {
		sizeof_off_t = 4;
	}

	assert(nc5->xsz == nc5i_len_NC(nc5,sizeof_off_t));

        status = ncio_create(path, ioflags, initialsz,
			     0, nc5->xsz, &nc5->chunk,
			     &nc5->nciop, &xp);
	if(status != NC_NOERR)
	{
		/* translate error status */
		if(status == EEXIST)
			status = NC_EEXIST;
		goto unwind_alloc;
	}

	fSet(nc5->flags, NC_CREAT);

	if(fIsSet(nc5->nciop->ioflags, NC_SHARE))
	{
		/*
		 * NC_SHARE implies sync up the number of records as well.
		 * (File format version one.)
		 * Note that other header changes are not shared
		 * automatically.  Some sort of IPC (external to this package)
		 * would be used to trigger a call to nc_sync().
		 */
		fSet(nc5->flags, NC_NSYNC);
	}

	status = nc5i_put_NC(nc5, &xp, sizeof_off_t, nc5->xsz);
	if(status != NC_NOERR)
		goto unwind_ioc;

	if(chunksizehintp != NULL)
		*chunksizehintp = nc5->chunk;

	/* Link nc5 and nc */
        NC5_DATA_SET(nc,nc5);
	nc->int_ncid = nc5->nciop->fd;

	return NC_NOERR;

unwind_ioc:
	if(nc5 != NULL) {
	    (void) ncio_close(nc5->nciop, 1); /* N.B.: unlink */
	    nc5->nciop = NULL;
	}
	/*FALLTHRU*/
unwind_alloc:
	free_NC5INFO(nc5);
	if(nc)
            NC5_DATA_SET(nc,NULL);
	return status;
}

int
NC5_create(const char *path, int cmode,
	  size_t initialsz, int basepe, size_t *chunksizehintp,
	  int use_parallel, void* mpidata,
	  struct NC_Dispatch* table, NC* nc)
{
#ifdef USE_PNETCDF
    if (use_parallel) {
        int res, default_format;
        NC5_INFO* nc5;

        /* Check the cmode for only valid flags*/
        if(cmode & ~LEGAL_CREATE_FLAGS)
	    return NC_EINVAL;

        /* Cannot have both NC_64BIT_OFFSET & NC_64BIT_DATA */
        if((cmode & (NC_64BIT_OFFSET|NC_64BIT_DATA)) == (NC_64BIT_OFFSET|NC_64BIT_DATA))
	    return NC_EINVAL;

	if (!fIsSet(cmode, NC_64BIT_DATA) && !fIsSet(cmode, NC_64BIT_OFFSET)) {
            /* if no format set in ioflag, apply default create format. */
            default_format = nc_get_default_format();
            if (default_format == NC_FORMAT_CDF2)
                cmode |= NC_64BIT_OFFSET;
            else if (default_format == NC_FORMAT_CDF5)
                cmode |= NC_64BIT_DATA;
        }

        /* Create our specific NC5_INFO instance */
	nc5 = new_NC5INFO(chunksizehintp);
        if(nc5 == NULL) return NC_ENOMEM;

        nc5->use_parallel = use_parallel;

        /* Link nc5 and nc */
        NC5_DATA_SET(nc,nc5);

        /* PnetCDF recognizes the flags below for create and ignores NC_LOCK and
         * NC_SHARE */
        cmode &= (NC_WRITE | NC_NOCLOBBER | NC_SHARE | NC_64BIT_OFFSET | NC_64BIT_DATA);

        /* No MPI environment initialized */
        if (mpidata == NULL) return NC_ENOPAR;

        res = ncmpi_create(((NC_MPI_INFO *)mpidata)->comm, path, cmode,
                           ((NC_MPI_INFO *)mpidata)->info, &(nc->int_ncid));

        if(res != NC_NOERR && nc5 != NULL) free_NC5INFO(nc5); /* reclaim allocated space */
        return res;
    }
#endif
    /* sequential program to create a file */
    return nc5_create_file(path, cmode, initialsz, basepe, chunksizehintp,
			   use_parallel, mpidata, table, nc);

}

static int
nc5_open_file(const char * path, int ioflags,
               int basepe, size_t *chunksizehintp,
	       int use_parallel,void* parameters,
               NC_Dispatch* dispatch, NC* nc)
{
	int status;
	NC5_INFO* nc5;

	/* Create our specific NC5_INFO instance */
	nc5 = new_NC5INFO(chunksizehintp);

#if ALWAYS_NC_SHARE /* DEBUG */
	fSet(ioflags, NC_SHARE);
#endif

#if defined(LOCKNUMREC) /* && _CRAYMPP */
	if (status = NC_init_pe(nc5, basepe)) {
		return status;
	}
#else
	/*
	 * !_CRAYMPP, only pe 0 is valid
	 */
	if(basepe != 0) {
      free(nc5);
      return NC_EINVAL;
    }
#endif

	status = ncio_open(path, ioflags, 0, 0, &nc5->chunk, &nc5->nciop, 0);
	if(status)
		goto unwind_alloc;

	assert(nc5->flags == 0);

	if(fIsSet(nc5->nciop->ioflags, NC_SHARE))
	{
		/*
		 * NC_SHARE implies sync up the number of records as well.
		 * (File format version one.)
		 * Note that other header changes are not shared
		 * automatically.  Some sort of IPC (external to this package)
		 * would be used to trigger a call to nc_sync().
		 */
		fSet(nc5->flags, NC_NSYNC);
	}

	status = nc5_get_NC(nc5);
	if(status != NC_NOERR)
		goto unwind_ioc;

	if(chunksizehintp != NULL)
		*chunksizehintp = nc5->chunk;

	/* Link nc5 and nc */
        NC5_DATA_SET(nc,nc5);
	nc->int_ncid = nc5->nciop->fd;

	return NC_NOERR;

unwind_ioc:
	if(nc5) {
    	    (void) ncio_close(nc5->nciop, 0);
	    nc5->nciop = NULL;
	}
	/*FALLTHRU*/
unwind_alloc:
	free_NC5INFO(nc5);
	if(nc)
            NC5_DATA_SET(nc,NULL);
	return status;
}

int
NC5_open(const char *path, int cmode,
	    int basepe, size_t *chunksizehintp,
	    int use_parallel, void* mpidata,
	    struct NC_Dispatch* table, NC* nc)
{
#ifdef USE_PNETCDF
    if (use_parallel) {
        int res;
        NC5_INFO* nc5;

        /* Check the cmode for only valid flags*/
        if(cmode & ~LEGAL_OPEN_FLAGS)
	    return NC_EINVAL;

        /* Create our specific NC5_INFO instance */
	nc5 = new_NC5INFO(chunksizehintp);
        if(nc5 == NULL) return NC_ENOMEM;

        nc5->use_parallel = use_parallel;

        /* Link nc5 and nc */
        NC5_DATA_SET(nc,nc5);

        /* PnetCDF recognizes the flags NC_WRITE and NC_NOCLOBBER for file open
         * and ignores NC_LOCK, NC_SHARE, NC_64BIT_OFFSET, and NC_64BIT_DATA.
         * Ignoring the NC_64BIT_OFFSET and NC_64BIT_DATA flags is because the
         * file is already in one of the CDF-formats, and setting these 2 flags
         * will not change the format of that file.
         */
        cmode &= (NC_WRITE | NC_NOCLOBBER);

        /* No MPI environment initialized */
        if (mpidata == NULL) return NC_ENOPAR;

        res = ncmpi_open(((NC_MPI_INFO *)mpidata)->comm, path, cmode,
                         ((NC_MPI_INFO *)mpidata)->info, &(nc->int_ncid));

        /* Default to independent access, like netCDF-4/HDF5 files. */
        if (res == NC_NOERR) {
	    res = ncmpi_begin_indep_data(nc->int_ncid);
	    nc5->pnetcdf_access_mode = NC_INDEPENDENT;
        }
        if(res != NC_NOERR && nc5 != NULL) free_NC5INFO(nc5); /* reclaim allocated space */
        return res;
    }
#endif
    /* sequential program to open a file */
    return nc5_open_file(path, cmode, basepe, chunksizehintp, use_parallel,
		         mpidata, table, nc);
}

int
NC5_redef(int ncid)
{
	int status;
	NC *nc;
	NC5_INFO* nc5;

	status = NC_check_id(ncid, &nc);
	if(status != NC_NOERR)
		return status;
	nc5 = NC5_DATA(nc);

#ifdef USE_PNETCDF
        if (nc5->use_parallel)
            return ncmpi_redef(nc->int_ncid);
#endif
	if(NC_readonly(nc5))
		return NC_EPERM;

	if(NC_indef(nc5))
		return NC_EINDEFINE;


	if(fIsSet(nc5->nciop->ioflags, NC_SHARE))
	{
		/* read in from disk */
		status = read_NC(nc5);
		if(status != NC_NOERR)
			return status;
	}

	nc5->old = dup_NC5INFO(nc5);
	if(nc5->old == NULL)
		return NC_ENOMEM;

	fSet(nc5->flags, NC_INDEF);

	return NC_NOERR;
}

int
NC5__enddef(int ncid,
            size_t h_minfree, size_t v_align,
            size_t v_minfree, size_t r_align)
{
    int status;
    NC* nc;
    NC5_INFO* nc5;

    status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR)
	return status;
    nc5 = NC5_DATA(nc);

#ifdef USE_PNETCDF
    if (nc5->use_parallel) {
        status = ncmpi__enddef(nc->int_ncid, h_minfree, v_align, v_minfree, r_align);
        if(status == NC_NOERR) {
	    if (nc5->pnetcdf_access_mode == NC_INDEPENDENT)
	        status = ncmpi_begin_indep_data(nc->int_ncid);
        }
        return status;
    }
#endif
    if(!NC_indef(nc5))
        return(NC_ENOTINDEFINE);

    return NC_endef(nc5, h_minfree, v_align, v_minfree, r_align);
}

int
NC5_sync(int ncid)
{
	int status;
	NC *nc;
	NC5_INFO* nc5;

	status = NC_check_id(ncid, &nc);
	if(status != NC_NOERR)
		return status;
	nc5 = NC5_DATA(nc);

#ifdef USE_PNETCDF
	if (nc5->use_parallel)
		return ncmpi_sync(nc->int_ncid);
#endif

	if(NC_indef(nc5))
		return NC_EINDEFINE;

	if(NC_readonly(nc5))
	{
		return read_NC(nc5);
	}
	/* else, read/write */

	status = nc5i_sync(nc5);
	if(status != NC_NOERR)
		return status;

	status = ncio_sync(nc5->nciop);
	if(status != NC_NOERR)
		return status;

#ifdef USE_FSYNC
	/* may improve concurrent access, but slows performance if
	 * called frequently */
#ifndef WIN32
	status = fsync(nc5->nciop->fd);
#else
	status = _commit(nc5->nciop->fd);
#endif	/* WIN32 */
#endif	/* USE_FSYNC */

	return status;
}

/*
 * In data mode, same as ncclose.
 * In define mode, restore previous definition.
 * In create, remove the file.
 */
int
NC5_abort(int ncid)
{
	int status;
	NC *nc;
	NC5_INFO* nc5;
	int doUnlink = 0;

	status = NC_check_id(ncid, &nc);
	if(status != NC_NOERR)
	    return status;
	nc5 = NC5_DATA(nc);

#ifdef USE_PNETCDF
	if (nc5->use_parallel) {
		status = ncmpi_abort(nc->int_ncid);
		if(nc5 != NULL) free_NC5INFO(nc5); /* reclaim allocated space */
		return status;
	}
#endif
	doUnlink = NC_IsNew(nc5);

	if(nc5->old != NULL)
	{
		/* a plain redef, not a create */
		assert(!NC_IsNew(nc5));
		assert(fIsSet(nc5->flags, NC_INDEF));
		free_NC5INFO(nc5->old);
		nc5->old = NULL;
		fClr(nc5->flags, NC_INDEF);
	}
	else if(!NC_readonly(nc5))
	{
		status = nc5i_sync(nc5);
		if(status != NC_NOERR)
			return status;
	}


	(void) ncio_close(nc5->nciop, doUnlink);
	nc5->nciop = NULL;

	free_NC5INFO(nc5);
	if(nc)
            NC5_DATA_SET(nc,NULL);

	return NC_NOERR;
}

int
NC5_close(int ncid)
{
	int status = NC_NOERR;
	NC *nc;
	NC5_INFO* nc5;

	status = NC_check_id(ncid, &nc);
	if(status != NC_NOERR)
	    return status;
	nc5 = NC5_DATA(nc);

#ifdef USE_PNETCDF
	if (nc5->use_parallel) {
		status = ncmpi_close(nc->int_ncid);
		if(nc5 != NULL) free_NC5INFO(nc5); /* reclaim allocated space */
		return status;
	}
#endif
	if(NC_indef(nc5))
	{
		status = NC_endef(nc5, 0, 1, 0, 1); /* TODO: defaults */
		if(status != NC_NOERR )
		{
			(void) NC5_abort(ncid);
			return status;
		}
	}
	else if(!NC_readonly(nc5))
	{
		status = nc5i_sync(nc5);
		/* flush buffers before any filesize comparisons */
		(void) ncio_sync(nc5->nciop);
	}

	/*
	 * If file opened for writing and filesize is less than
	 * what it should be (due to previous use of NOFILL mode),
	 * pad it to correct size, as reported by nc5i_NC_calcsize().
	 */
	if (status == ENOERR) {
	    off_t filesize; 	/* current size of open file */
	    off_t calcsize;	/* calculated file size, from header */
	    status = ncio_filesize(nc5->nciop, &filesize);
	    if(status != ENOERR)
		return status;
	    status = nc5i_NC_calcsize(nc5, &calcsize);
	    if(status != NC_NOERR)
		return status;
	    if(filesize < calcsize && !NC_readonly(nc5)) {
		status = ncio_pad_length(nc5->nciop, calcsize);
		if(status != ENOERR)
		    return status;
	    }
	}

	(void) ncio_close(nc5->nciop, 0);
	nc5->nciop = NULL;

	free_NC5INFO(nc5);
        NC5_DATA_SET(nc,NULL);

	return status;
}


int
NC5_set_fill(int ncid,
	int fillmode, int *old_mode_ptr)
{
	int status;
	NC *nc;
	NC5_INFO* nc5;
	int oldmode;

	status = NC_check_id(ncid, &nc);
	if(status != NC_NOERR)
		return status;
	nc5 = NC5_DATA(nc);

#ifdef USE_PNETCDF
	if (nc5->use_parallel)
		return ncmpi_set_fill(nc->int_ncid,fillmode,old_mode_ptr);
#endif

	if(NC_readonly(nc5))
		return NC_EPERM;

	oldmode = fIsSet(nc5->flags, NC_NOFILL) ? NC_NOFILL : NC_FILL;

	if(fillmode == NC_NOFILL)
	{
		fSet(nc5->flags, NC_NOFILL);
	}
	else if(fillmode == NC_FILL)
	{
		if(fIsSet(nc5->flags, NC_NOFILL))
		{
			/*
			 * We are changing back to fill mode
			 * so do a sync
			 */
			status = nc5i_sync(nc5);
			if(status != NC_NOERR)
				return status;
		}
		fClr(nc5->flags, NC_NOFILL);
	}
	else
	{
		return NC_EINVAL; /* Invalid fillmode */
	}

	if(old_mode_ptr != NULL)
		*old_mode_ptr = oldmode;

	return NC_NOERR;
}

#ifdef LOCKNUMREC

/* create function versions of the NC_*_numrecs macros */
size_t
NC5_get_numrecs(const NC *nc5) {
	shmem_t numrec;
	shmem_short_get(&numrec, (shmem_t *) nc5->lock + LOCKNUMREC_VALUE, 1,
		nc5->lock[LOCKNUMREC_BASEPE]);
	return (size_t) numrec;
}

void
NC5_set_numrecs(NC *nc5, size_t nrecs)
{
    shmem_t numrec = (shmem_t) nrecs;
    /* update local value too */
    nc5->lock[LOCKNUMREC_VALUE] = (ushmem_t) numrec;
    shmem_short_put((shmem_t *) nc5->lock + LOCKNUMREC_VALUE, &numrec, 1,
    nc5->lock[LOCKNUMREC_BASEPE]);
}

void NC5_increase_numrecs(NC *nc5, size_t nrecs)
{
    /* this is only called in one place that's already protected
     * by a lock ... so don't worry about it */
    if (nrecs > NC5_get_numrecs(nc5))
	NC5_set_numrecs(nc5, nrecs);
}

#endif /* LOCKNUMREC */

/* everyone in communicator group will be executing this */
/*ARGSUSED*/
int
NC5_set_base_pe(int ncid, int pe)
{
#ifdef USE_PNETCDF
    NC* nc;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;

    if (((NC5_INFO*)nc->dispatchdata)->use_parallel)
        return NC_NOERR;
#endif

#if _CRAYMPP && defined(LOCKNUMREC)
	int status;
	NC *nc;
	NC5_INFO* nc5;
	shmem_t numrecs;

	if ((status = NC_check_id(ncid, &nc) != NC_NOERR) {
		return status;
	}
	if (pe < 0 || pe >= _num_pes()) {
		return NC_EINVAL; /* invalid base pe */
	}
	nc5 = NC5_DATA(nc);

	numrecs = (shmem_t) NC5_get_numrecs(nc5);

	nc5->lock[LOCKNUMREC_VALUE] = (ushmem_t) numrecs;

	/* update serving & lock values for a "smooth" transition */
	/* note that the "real" server will being doing this as well */
	/* as all the rest in the group */
	/* must have syncronization before & after this step */
	shmem_short_get(
		(shmem_t *) nc5->lock + LOCKNUMREC_SERVING,
		(shmem_t *) nc5->lock + LOCKNUMREC_SERVING,
		1, nc5->lock[LOCKNUMREC_BASEPE]);

	shmem_short_get(
		(shmem_t *) nc5->lock + LOCKNUMREC_LOCK,
		(shmem_t *) nc5->lock + LOCKNUMREC_LOCK,
		1, nc5->lock[LOCKNUMREC_BASEPE]);

	/* complete transition */
	nc5->lock[LOCKNUMREC_BASEPE] = (ushmem_t) pe;

#endif /* _CRAYMPP && LOCKNUMREC */
	return NC_NOERR;
}

/*ARGSUSED*/
int
NC5_inq_base_pe(int ncid, int *pe)
{
#ifdef USE_PNETCDF
    NC* nc;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;

    if (((NC5_INFO*)nc->dispatchdata)->use_parallel) {
        if(pe) *pe = 0;
        return NC_NOERR;
    }
#endif

#if _CRAYMPP && defined(LOCKNUMREC)
	int status;
	NC *nc;
	NC5_INFO* nc5;

	if ((status = NC_check_id(ncid, &nc)) != NC_NOERR) {
		return status;
	}

	*pe = (int) nc5->lock[LOCKNUMREC_BASEPE];
	nc5 = NC5_DATA(nc);
#else
	/*
	 * !_CRAYMPP, only pe 0 is valid
	 */
	*pe = 0;
#endif /* _CRAYMPP && LOCKNUMREC */
	return NC_NOERR;
}

int
NC5_inq_format(int ncid, int* formatp)
{
        int status;
        NC *nc;
        NC5_INFO* nc5;

        status = NC_check_id(ncid, &nc);
        if(status != NC_NOERR)
                return status;

        nc5 = NC5_DATA(nc);

        /* only need to check for netCDF-3 variants, since this is never
           called for netCDF-4 files */
        if (fIsSet(nc5->flags, NC_64BIT_DATA))
                *formatp = NC_FORMAT_CDF5;
        else if (fIsSet(nc5->flags, NC_64BIT_OFFSET))
                *formatp = NC_FORMAT_CDF2;
        else
                *formatp = NC_FORMAT_CLASSIC;

        return NC_NOERR;
}

int
NC5_inq_format_extended(int ncid, int* formatp, int *modep)
{
    NC* nc;
    int status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;

    if(modep) *modep = nc->mode;

#ifdef USE_PNETCDF
    if (((NC5_INFO*)nc->dispatchdata)->use_parallel) {
        if(formatp) *formatp = NC_FORMAT_PNETCDF;
    }
    else
#endif
    {
        if(formatp) *formatp = NC_FORMAT_NC5;
    }
    return NC_NOERR;
}

int
NC5_inq(int ncid,
	int *ndimsp,
	int *nvarsp,
	int *nattsp,
	int *xtendimp)
{
	int status;
	NC *nc;
	NC5_INFO* nc5;

	status = NC_check_id(ncid, &nc);
	if(status != NC_NOERR)
		return status;
	nc5 = NC5_DATA(nc);

#ifdef USE_PNETCDF
	if (nc5->use_parallel)
		return ncmpi_inq(nc->int_ncid,ndimsp,nvarsp,nattsp,xtendimp);
#endif
	if(ndimsp != NULL)
		*ndimsp = (int) nc5->dims.nelems;
	if(nvarsp != NULL)
		*nvarsp = (int) nc5->vars.nelems;
	if(nattsp != NULL)
		*nattsp = (int) nc5->attrs.nelems;
	if(xtendimp != NULL)
		*xtendimp = nc5i_find_NC_Udim(&nc5->dims, NULL);

	return NC_NOERR;
}

/* The sizes of types may vary from platform to platform, but within
 * netCDF files, type sizes are fixed. */
#define NC_CHAR_LEN sizeof(char)
#define NC_STRING_LEN sizeof(char *)
#define NC_BYTE_LEN 1
#define NC_SHORT_LEN 2
#define NC_INT_LEN 4
#define NC_FLOAT_LEN 4
#define NC_DOUBLE_LEN 8
#define NC_INT64_LEN 8

static int atomic_size[12] = {NC_BYTE_LEN,  NC_CHAR_LEN,  NC_SHORT_LEN,
                              NC_INT_LEN,   NC_FLOAT_LEN, NC_DOUBLE_LEN,
                              NC_BYTE_LEN,  NC_SHORT_LEN, NC_INT_LEN,
                              NC_INT64_LEN, NC_INT64_LEN, NC_STRING_LEN};

static char
atomic_name[12][NC_MAX_NAME + 1] = {"byte", "char", "short",
                                    "int", "float", "double",
                                    "ubyte", "ushort","uint",
                                    "int64", "uint64", "string"};

int
NC5_inq_type(int ncid, nc_type typeid, char* name, size_t* size)
{
    NC* nc;
    int format, status;

    status = NC_check_id(ncid, &nc);
    if (status != NC_NOERR) return status;

#ifdef USE_PNETCDF
    if (((NC5_INFO*)nc->dispatchdata)->use_parallel)
        status = ncmpi_inq_format(nc->int_ncid, &format);
    else
#endif
        status = NC5_inq_format(ncid, &format);

    if (status != NC_NOERR) return status;

    if (format == NC_FORMAT_CDF5) {
        if (typeid < NC_BYTE || typeid >= NC_STRING)
           return NC_EBADTYPE;
    }
    else { /* CDF-1 and CDF-2 */
        if (typeid < NC_BYTE || typeid > NC_DOUBLE)
           return NC_EBADTYPE;
    }

    /* Give the user the values they want. Subtract one because types
     * are numbered starting at 1, not 0. */
    if (name)
        strcpy(name, atomic_name[typeid - 1]);
    if (size)
        *size = atomic_size[typeid - 1];

    return NC_NOERR;
}

int
NC5_inq_var_all(int ncid, int varid, char *name, nc_type *xtypep,
               int *ndimsp, int *dimidsp, int *nattsp,
               int *shufflep, int *deflatep, int *deflate_levelp,
               int *fletcher32p, int *contiguousp, size_t *chunksizesp,
               int *no_fill, void *fill_valuep, int *endiannessp,
	       int *options_maskp, int *pixels_per_blockp)
{
    int status;

#ifdef USE_PNETCDF
    NC* nc;
    status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;

    if (((NC5_INFO*)nc->dispatchdata)->use_parallel)
        status = ncmpi_inq_var(nc->int_ncid, varid, name, xtypep, ndimsp, dimidsp, nattsp);
    else
#endif
        status = NC5_inq_var(ncid, varid, name, xtypep, ndimsp, dimidsp, nattsp);

    if(status != NC_NOERR) return status;
    if(shufflep) *shufflep = 0;
    if(deflatep) *deflatep = 0;
    if(fletcher32p) *fletcher32p = 0;
    if(contiguousp) *contiguousp = NC_CONTIGUOUS;
    if(no_fill) *no_fill = 1;
    if(endiannessp) return NC_ENOTNC4;
    if(options_maskp) return NC_ENOTNC4;
    return NC_NOERR;
}

int
NC5_var_par_access(int ncid, int varid, int par_access)
{
    NC *nc;
    int status;

    status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;
#ifdef USE_PNETCDF
    if (par_access != NC_INDEPENDENT && par_access != NC_COLLECTIVE)
	return NC_EINVAL;

    if(par_access == ((NC5_INFO*)nc->dispatchdata)->pnetcdf_access_mode)
	return NC_NOERR;
    ((NC5_INFO*)nc->dispatchdata)->pnetcdf_access_mode = par_access;

    if (par_access == NC_INDEPENDENT)
	status = ncmpi_begin_indep_data(nc->int_ncid);
    else
	status = ncmpi_end_indep_data(nc->int_ncid);
#endif
    return status;
}

