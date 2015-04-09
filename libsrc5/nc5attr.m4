dnl This is m4 source.
dnl Process using m4 to produce 'C' language file.
dnl
dnl If you see this line, you can ignore the next one.
/* Do not edit this file. It is produced from the corresponding .m4 source */
dnl

#include "nc5internal.h"
#include "ncdispatch.h"
#include "nc5dispatch.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "nc5ncx.h"
#include "fbits.h"
#include "rnd.h"
#include "utf8proc.h"

#ifdef USE_PNETCDF
/* Must follow netcdf.h */
/* In PnetCDF version 1.5.0 and prior, NC_64BIT_DATA is defined in conflict
 * with 1.6.0 and netCDF 4.3.4.
 * In PnetCDF version 1.4.1 and prior, NC_FILL_STRING and NC4_LAST_ERROR are
 * either undefined or in conflict with netCDF.
 * Redefine these constants to avoid compilation warning.
 */
#undef NC_64BIT_DATA
#undef NC_FILL_STRING
#undef NC4_LAST_ERROR
#include <pnetcdf.h>
#ifdef NC_64BIT_DATA
#undef NC_64BIT_DATA
#endif
#ifdef NC_FILL_STRING
#undef NC_FILL_STRING
#endif
#ifdef NC4_LAST_ERROR
#undef NC4_LAST_ERROR
#endif
#define NC_64BIT_DATA 0x0020
#define NC_FILL_STRING  ((char *)"")
#define NC4_LAST_ERROR  (-131)
#endif

/*
 * Free attr
 * Formerly
NC_free_attr()
 */
void
nc5i_free_NC_attr(NC_attr *attrp)
{

	if(attrp == NULL)
		return;
	free_NC_string(attrp->name);
	free(attrp);
}


/*
 * How much space will 'nelems' of 'type' take in
 *  external representation (as the values of an attribute)?
 */
static size_t
nc5i_len_NC_attrV(nc_type type, size_t nelems)
{
	switch(type) {
	case NC_BYTE:
	case NC_CHAR:
		return nc5x_len_char(nelems);
	case NC_SHORT:
		return nc5x_len_short(nelems);
	case NC_INT:
		return nc5x_len_int(nelems);
	case NC_FLOAT:
		return nc5x_len_float(nelems);
	case NC_DOUBLE:
		return nc5x_len_double(nelems);
	case NC_UBYTE:
		return nc5x_len_ubyte(nelems);
	case NC_USHORT:
		return nc5x_len_ushort(nelems);
	case NC_UINT:
		return nc5x_len_uint(nelems);
	case NC_INT64:
		return nc5x_len_int64(nelems);
	case NC_UINT64:
		return nc5x_len_uint64(nelems);
	default:
	        assert("nc5i_len_NC_attr bad type" == 0);
	}
	return 0;
}


NC_attr *
nc5i_new_x_NC_attr(
	NC_string *strp,
	nc_type type,
	size_t nelems)
{
	NC_attr *attrp;
	const size_t xsz = nc5i_len_NC_attrV(type, nelems);
	size_t sz = M_RNDUP(sizeof(NC_attr));

	assert(!(xsz == 0 && nelems != 0));

	sz += xsz;

	attrp = (NC_attr *) malloc(sz);
	if(attrp == NULL )
		return NULL;

	attrp->xsz = xsz;

	attrp->name = strp;
	attrp->type = type;
	attrp->nelems = nelems;
	if(xsz != 0)
		attrp->xvalue = (char *)attrp + M_RNDUP(sizeof(NC_attr));
	else
		attrp->xvalue = NULL;

	return(attrp);
}


/*
 * Formerly
NC_new_attr(name,type,count,value)
 */
static NC_attr *
new_NC_attr(
	const char *uname,
	nc_type type,
	size_t nelems)
{
	NC_string *strp;
	NC_attr *attrp;

	char *name = (char *)utf8proc_NFC((const unsigned char *)uname);
	if(name == NULL)
	    return NULL;
	assert(name != NULL && *name != 0);

	strp = new_NC_string(strlen(name), name);
	free(name);
	if(strp == NULL)
		return NULL;
	
	attrp = nc5i_new_x_NC_attr(strp, type, nelems);
	if(attrp == NULL)
	{
		free_NC_string(strp);
		return NULL;
	}

	return(attrp);
}


static NC_attr *
dup_NC_attr(const NC_attr *rattrp)
{
	NC_attr *attrp = new_NC_attr(rattrp->name->cp,
		 rattrp->type, rattrp->nelems);
	if(attrp == NULL)
		return NULL;
	(void) memcpy(attrp->xvalue, rattrp->xvalue, rattrp->xsz);
	return attrp;
}

/* attrarray */

/*
 * Free the stuff "in" (referred to by) an NC_attrarray.
 * Leaves the array itself allocated.
 */
void
nc5i_free_NC_attrarrayV0(NC_attrarray *ncap)
{
	assert(ncap != NULL);

	if(ncap->nelems == 0)
		return;

	assert(ncap->value != NULL);

	{
		NC_attr **app = ncap->value;
		NC_attr *const *const end = &app[ncap->nelems];
		for( /*NADA*/; app < end; app++)
		{
			nc5i_free_NC_attr(*app);
			*app = NULL;
		}
	}
	ncap->nelems = 0;
}


/*
 * Free NC_attrarray values.
 * formerly
NC_free_array()
 */
void
nc5i_free_NC_attrarrayV(NC_attrarray *ncap)
{
	assert(ncap != NULL);
	
	if(ncap->nalloc == 0)
		return;

	assert(ncap->value != NULL);

	nc5i_free_NC_attrarrayV0(ncap);

	free(ncap->value);
	ncap->value = NULL;
	ncap->nalloc = 0;
}


int
nc5i_dup_NC_attrarrayV(NC_attrarray *ncap, const NC_attrarray *ref)
{
	int status = NC_NOERR;

	assert(ref != NULL);
	assert(ncap != NULL);

	if(ref->nelems != 0)
	{
		const size_t sz = ref->nelems * sizeof(NC_attr *);
		ncap->value = (NC_attr **) malloc(sz);
		if(ncap->value == NULL)
			return NC_ENOMEM;

		(void) memset(ncap->value, 0, sz);
		ncap->nalloc = ref->nelems;
	}

	ncap->nelems = 0;
	{
		NC_attr **app = ncap->value;
		const NC_attr **drpp = (const NC_attr **)ref->value;
		NC_attr *const *const end = &app[ref->nelems];
		for( /*NADA*/; app < end; drpp++, app++, ncap->nelems++)
		{
			*app = dup_NC_attr(*drpp);
			if(*app == NULL)
			{
				status = NC_ENOMEM;
				break;
			}
		}
	}

	if(status != NC_NOERR)
	{
		nc5i_free_NC_attrarrayV(ncap);
		return status;
	}

	assert(ncap->nelems == ref->nelems);

	return NC_NOERR;
}


/*
 * Add a new handle on the end of an array of handles
 * Formerly
NC_incr_array(array, tail)
 */
static int
incr_NC_attrarray(NC_attrarray *ncap, NC_attr *newelemp)
{
	NC_attr **vp;

	assert(ncap != NULL);

	if(ncap->nalloc == 0)
	{
		assert(ncap->nelems == 0);
		vp = (NC_attr **) malloc(NC_ARRAY_GROWBY * sizeof(NC_attr *));
		if(vp == NULL)
			return NC_ENOMEM;

		ncap->value = vp;
		ncap->nalloc = NC_ARRAY_GROWBY;
	}
	else if(ncap->nelems +1 > ncap->nalloc)
	{
		vp = (NC_attr **) realloc(ncap->value,
			(ncap->nalloc + NC_ARRAY_GROWBY) * sizeof(NC_attr *));
		if(vp == NULL)
			return NC_ENOMEM;
	
		ncap->value = vp;
		ncap->nalloc += NC_ARRAY_GROWBY;
	}

	if(newelemp != NULL)
	{
		ncap->value[ncap->nelems] = newelemp;
		ncap->nelems++;
	}
	return NC_NOERR;
}


static NC_attr *
nc5i_elem_NC_attrarray(const NC_attrarray *ncap, size_t elem)
{
	assert(ncap != NULL);
	/* cast needed for braindead systems with signed size_t */
	if(ncap->nelems == 0 || (unsigned long) elem >= ncap->nelems)
		return NULL;

	assert(ncap->value != NULL);

	return ncap->value[elem];
}

/* End attarray per se */

/*
 * Given ncp and varid, return ptr to array of attributes
 *  else NULL on error
 */
static NC_attrarray *
NC_attrarray0(NC5_INFO* ncp, int varid)
{
	NC_attrarray *ap;

	if(varid == NC_GLOBAL) /* Global attribute, attach to cdf */
	{
		ap = &ncp->attrs;
	}
	else if(varid >= 0 && (size_t) varid < ncp->vars.nelems)
	{
		NC_var **vpp;
		vpp = (NC_var **)ncp->vars.value;
		vpp += varid;
		ap = &(*vpp)->attrs;
	} else {
		ap = NULL;
	}
	return(ap);
}


/*
 * Step thru NC_ATTRIBUTE array, seeking match on name.
 *  return match or NULL if Not Found or out of memory.
 */
NC_attr **
NC5_findattr(const NC_attrarray *ncap, const char *uname)
{
	NC_attr **attrpp;
	size_t attrid;
	size_t slen;
	char *name;

	assert(ncap != NULL);

	if(ncap->nelems == 0)
		return NULL;

	attrpp = (NC_attr **) ncap->value;

	/* normalized version of uname */
	name = (char *)utf8proc_NFC((const unsigned char *)uname);
	if(name == NULL)
	    return NULL; /* TODO: need better way to indicate no memory */
	slen = strlen(name);

	for(attrid = 0; attrid < ncap->nelems; attrid++, attrpp++)
	{
		if(strlen((*attrpp)->name->cp) == slen &&
			strncmp((*attrpp)->name->cp, name, slen) == 0)
		{
		        free(name);
			return(attrpp); /* Normal return */
		}
	}
	free(name);
	return(NULL);
}


/*
 * Look up by ncid, varid and name, return NULL if not found
 */
static int 
NC_lookupattr(int ncid,
	int varid,
	const char *name, /* attribute name */
	NC_attr **attrpp) /* modified on return */
{
	int status;
	NC* nc;
	NC5_INFO *ncp;
	NC_attrarray *ncap;
	NC_attr **tmp;

	status = NC_check_id(ncid, &nc);
	if(status != NC_NOERR)
		return status;
	ncp = NC5_DATA(nc);

	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;

	tmp = NC5_findattr(ncap, name);
	if(tmp == NULL)
		return NC_ENOTATT;

	if(attrpp != NULL)
		*attrpp = *tmp;

	return ENOERR;
}

/* Public */

int
NC5_inq_attname(int ncid, int varid, int attnum, char *name)
{
	int status;
	NC* nc;
	NC5_INFO *ncp;
	NC_attrarray *ncap;
	NC_attr *attrp;

	status = NC_check_id(ncid, &nc);
	if(status != NC_NOERR)
		return status;
	ncp = NC5_DATA(nc);

#ifdef USE_PNETCDF
	if (ncp->use_parallel)
		return ncmpi_inq_attname(nc->int_ncid,varid,attnum,name);
#endif
	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;

	attrp = nc5i_elem_NC_attrarray(ncap, (size_t)attnum);
	if(attrp == NULL)
		return NC_ENOTATT;

	(void) strncpy(name, attrp->name->cp, attrp->name->nchars);
	name[attrp->name->nchars] = 0;

	return NC_NOERR;
}


int 
NC5_inq_attid(int ncid, int varid, const char *name, int *attnump)
{
	int status;
	NC *nc;
	NC5_INFO* ncp;
	NC_attrarray *ncap;
	NC_attr **attrpp;

	status = NC_check_id(ncid, &nc);
	if(status != NC_NOERR)
		return status;
	ncp = NC5_DATA(nc);

#ifdef USE_PNETCDF
	if (ncp->use_parallel)
		return ncmpi_inq_attid(nc->int_ncid,varid,name,attnump);
#endif
	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;
	

	attrpp = NC5_findattr(ncap, name);
	if(attrpp == NULL)
		return NC_ENOTATT;

	if(attnump != NULL)
		*attnump = (int)(attrpp - ncap->value);

	return NC_NOERR;
}

int
NC5_inq_att(int ncid,
	int varid,
	const char *name, /* input, attribute name */
	nc_type *datatypep,
	size_t *lenp)
{
	int status;
	NC_attr *attrp;

#ifdef USE_PNETCDF
	NC* nc;
	status = NC_check_id(ncid, &nc);
	if(status != NC_NOERR) return status;

	if (((NC5_INFO*)nc->dispatchdata)->use_parallel) {
		MPI_Offset mpilen;
		status = ncmpi_inq_att(nc->int_ncid,varid,name,datatypep,&mpilen);
		if(status == NC_NOERR && lenp) *lenp = mpilen;
		return status;
	}
#endif
	status = NC_lookupattr(ncid, varid, name, &attrp);
	if(status != NC_NOERR)
		return status;

	if(datatypep != NULL)
		*datatypep = attrp->type;
	if(lenp != NULL)
		*lenp = attrp->nelems;

	return NC_NOERR;
}


int
NC5_rename_att( int ncid, int varid, const char *name, const char *unewname)
{
	int status;
	NC *nc;
	NC5_INFO* ncp;
	NC_attrarray *ncap;
	NC_attr **tmp;
	NC_attr *attrp;
	NC_string *newStr, *old;
	char *newname;  /* normalized version */

	/* sortof inline clone of NC_lookupattr() */
	status = NC_check_id(ncid, &nc);
	if(status != NC_NOERR)
		return status;
	ncp = NC5_DATA(nc);

#ifdef USE_PNETCDF
	if (ncp->use_parallel)
		return ncmpi_rename_att(nc->int_ncid,varid,name,unewname);
#endif
	if(NC_readonly(ncp))
		return NC_EPERM;

	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;

	status = NC_check_name(unewname);
	if(status != NC_NOERR)
		return status;

	tmp = NC5_findattr(ncap, name);
	if(tmp == NULL)
		return NC_ENOTATT;
	attrp = *tmp;
			/* end inline clone NC_lookupattr() */

	if(NC5_findattr(ncap, unewname) != NULL)
	{
		/* name in use */
		return NC_ENAMEINUSE;
	}

	old = attrp->name;
	newname = (char *)utf8proc_NFC((const unsigned char *)unewname);
	if(newname == NULL)
	    return NC_EBADNAME;
	if(NC_indef(ncp))
	{
		newStr = new_NC_string(strlen(newname), newname);
		free(newname);
		if( newStr == NULL)
			return NC_ENOMEM;
		attrp->name = newStr;
		free_NC_string(old);
		return NC_NOERR;
	}
	/* else */
	status = set_NC_string(old, newname);
	free(newname);
	if( status != NC_NOERR)
		return status;

	set_NC_hdirty(ncp);

	if(NC_doHsync(ncp))
	{
		status = nc5i_sync(ncp);
		if(status != NC_NOERR)
			return status;
	}

	return NC_NOERR;
}

int
NC5_del_att(int ncid, int varid, const char *uname)
{
	int status;
	NC *nc;
	NC5_INFO* ncp;
	NC_attrarray *ncap;
	NC_attr **attrpp;
	NC_attr *old = NULL;
	int attrid;
	size_t slen;

	status = NC_check_id(ncid, &nc);
	if(status != NC_NOERR)
		return status;
	ncp = NC5_DATA(nc);

#ifdef USE_PNETCDF
	if (ncp->use_parallel)
		return ncmpi_del_att(nc->int_ncid,varid,uname);
#endif
	if(!NC_indef(ncp))
		return NC_ENOTINDEFINE;

	ncap = NC_attrarray0(ncp, varid);
	if(ncap == NULL)
		return NC_ENOTVAR;

	{
	char *name = (char *)utf8proc_NFC((const unsigned char *)uname);
	if(name == NULL)
	    return NC_ENOMEM;
	
			/* sortof inline NC5_findattr() */
	slen = strlen(name);

	attrpp = (NC_attr **) ncap->value;
	for(attrid = 0; (size_t) attrid < ncap->nelems; attrid++, attrpp++)
	    {
		if( slen == (*attrpp)->name->nchars &&
			strncmp(name, (*attrpp)->name->cp, slen) == 0)
		{
			old = *attrpp;
			break;
		}
	    }
	free(name);
	}
	if( (size_t) attrid == ncap->nelems )
		return NC_ENOTATT;
			/* end inline NC5_findattr() */

	/* shuffle down */
	for(attrid++; (size_t) attrid < ncap->nelems; attrid++)
	{
		*attrpp = *(attrpp + 1);
		attrpp++;
	}
	*attrpp = NULL;
	/* decrement count */
	ncap->nelems--;

	nc5i_free_NC_attr(old);

	return NC_NOERR;
}

dnl
dnl XNCX_PAD_PUTN(Type)
dnl
define(`XNCX_PAD_PUTN',dnl
`dnl
static int
nc5x_pad_putn_I$1(void **xpp, size_t nelems, const $1 *tp, nc_type type)
{
	switch(type) {
	case NC_CHAR:
		return NC_ECHAR;
	case NC_BYTE:
		return nc5x_pad_putn_schar_$1(xpp, nelems, tp);
	case NC_SHORT:
		return nc5x_pad_putn_short_$1(xpp, nelems, tp);
	case NC_INT:
		return nc5x_putn_int_$1(xpp, nelems, tp);
	case NC_FLOAT:
		return nc5x_putn_float_$1(xpp, nelems, tp);
	case NC_DOUBLE:
		return nc5x_putn_double_$1(xpp, nelems, tp);
	case NC_UBYTE:
		return nc5x_pad_putn_uchar_$1(xpp, nelems, tp);
	case NC_USHORT:
		return nc5x_putn_ushort_$1(xpp, nelems, tp);
	case NC_UINT:
		return nc5x_putn_uint_$1(xpp, nelems, tp);
	case NC_INT64:
		return nc5x_putn_int64_$1(xpp, nelems, tp);
	case NC_UINT64:
		return nc5x_putn_uint64_$1(xpp, nelems, tp);
	default:
                assert("nc5x_pad_putn_I$1 invalid type" == 0);
	}
	return NC_EBADTYPE;
}
')dnl
dnl
dnl XNCX_PAD_GETN(Type)
dnl
define(`XNCX_PAD_GETN',dnl
`dnl
static int
nc5x_pad_getn_I$1(const void **xpp, size_t nelems, $1 *tp, nc_type type)
{
	switch(type) {
	case NC_CHAR:
		return NC_ECHAR;
	case NC_UBYTE:
		return nc5x_pad_getn_uchar_$1(xpp, nelems, tp);
	case NC_BYTE:
		return nc5x_pad_getn_schar_$1(xpp, nelems, tp);
	case NC_SHORT:
		return nc5x_pad_getn_short_$1(xpp, nelems, tp);
	case NC_INT:
		return nc5x_getn_int_$1(xpp, nelems, tp);
	case NC_FLOAT:
		return nc5x_getn_float_$1(xpp, nelems, tp);
	case NC_DOUBLE:
		return nc5x_getn_double_$1(xpp, nelems, tp);
	case NC_USHORT:
		return nc5x_getn_ushort_$1(xpp, nelems, tp);
	case NC_UINT:
		return nc5x_getn_uint_$1(xpp, nelems, tp);
	case NC_INT64:
		return nc5x_getn_int64_$1(xpp, nelems, tp);
	case NC_UINT64:
		return nc5x_getn_uint64_$1(xpp, nelems, tp);
	default:
	        assert("nc5x_pad_getn_I$1 invalid type" == 0);
	}
	return NC_EBADTYPE;
}
')dnl
dnl Implement

XNCX_PAD_PUTN(uchar)
XNCX_PAD_GETN(uchar)

XNCX_PAD_PUTN(schar)
XNCX_PAD_GETN(schar)

XNCX_PAD_PUTN(short)
XNCX_PAD_GETN(short)

XNCX_PAD_PUTN(int)
XNCX_PAD_GETN(int)

XNCX_PAD_PUTN(float)
XNCX_PAD_GETN(float)

XNCX_PAD_PUTN(double)
XNCX_PAD_GETN(double)

#ifdef IGNORE
XNCX_PAD_PUTN(long)
XNCX_PAD_GETN(long)
#endif

XNCX_PAD_PUTN(longlong)
XNCX_PAD_GETN(longlong)

XNCX_PAD_PUTN(ushort)
XNCX_PAD_GETN(ushort)

XNCX_PAD_PUTN(uint)
XNCX_PAD_GETN(uint)

XNCX_PAD_PUTN(ulonglong)
XNCX_PAD_GETN(ulonglong)


/* Common dispatcher for put cases */
static int
dispatchput(void **xpp, size_t nelems, const void* tp,
	    nc_type atype, nc_type memtype)
{
    switch (memtype) {
    case NC_CHAR:
        return nc5x_pad_putn_text(xpp,nelems, (char *)tp);
    case NC_BYTE:
        return nc5x_pad_putn_Ischar(xpp, nelems, (schar*)tp, atype);
    case NC_SHORT:
        return nc5x_pad_putn_Ishort(xpp, nelems, (short*)tp, atype);
    case NC_INT:
          return nc5x_pad_putn_Iint(xpp, nelems, (int*)tp, atype);
    case NC_FLOAT:
        return nc5x_pad_putn_Ifloat(xpp, nelems, (float*)tp, atype);
    case NC_DOUBLE:
        return nc5x_pad_putn_Idouble(xpp, nelems, (double*)tp, atype);
    case NC_UBYTE: /*Synthetic*/
        return nc5x_pad_putn_Iuchar(xpp,nelems, (uchar *)tp, atype);
    case NC_USHORT:
          return nc5x_pad_putn_Iushort(xpp, nelems, (ushort*)tp, atype);
    case NC_UINT:
          return nc5x_pad_putn_Iuint(xpp, nelems, (uint*)tp, atype);
    case NC_INT64:
          return nc5x_pad_putn_Ilonglong(xpp, nelems, (longlong*)tp, atype);
    case NC_UINT64:
          return nc5x_pad_putn_Iulonglong(xpp, nelems, (ulonglong*)tp, atype);
    case NC_NAT:
        return NC_EBADTYPE;
    default:
        break;
    }
    return NC_EBADTYPE;
}

int
NC5_put_att(
	int ncid,
	int varid,
	const char *name,
	nc_type type,
	size_t nelems,
	const void *value,
	nc_type memtype)
{
    int status;
    NC *nc;
    NC5_INFO* ncp;
    NC_attrarray *ncap;
    NC_attr **attrpp;
    NC_attr *old = NULL;
    NC_attr *attrp;

    status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR)
	return status;
    ncp = NC5_DATA(nc);

#ifdef USE_PNETCDF
    if (ncp->use_parallel) {
        if (!name || (strlen(name) > NC_MAX_NAME))
	    return NC_EBADNAME;

        switch (memtype) {
        case NC_CHAR:
            return ncmpi_put_att_text(nc->int_ncid, varid, name, nelems, (char*)value);
        case NC_BYTE:
            return ncmpi_put_att_schar(nc->int_ncid, varid, name, type, nelems, (signed char*)value);
        case NC_SHORT:
            return ncmpi_put_att_short(nc->int_ncid, varid, name, type, nelems, (short*)value);
        case NC_INT:
            return ncmpi_put_att_int(nc->int_ncid, varid, name, type, nelems, (int*)value);
        case NC_FLOAT:
            return ncmpi_put_att_float(nc->int_ncid, varid, name, type, nelems, (float*)value);
        case NC_DOUBLE:
            return ncmpi_put_att_double(nc->int_ncid, varid, name, type, nelems, (double*)value);
        case NC_UBYTE:
            return ncmpi_put_att_uchar(nc->int_ncid, varid, name, type, nelems, (unsigned char*)value);
        case NC_USHORT:
            return ncmpi_put_att_ushort(nc->int_ncid, varid, name, type, nelems, (unsigned short*)value);
        case NC_UINT:
            return ncmpi_put_att_uint(nc->int_ncid, varid, name, type, nelems, (unsigned int*)value);
        case NC_INT64:
            return ncmpi_put_att_longlong(nc->int_ncid, varid, name, type, nelems, (long long*)value);
        case NC_UINT64:
            return ncmpi_put_att_ulonglong(nc->int_ncid, varid, name, type, nelems, (unsigned long long*)value);
        default:
            return NC_EBADTYPE;
        }
    }
#endif

    if(NC_readonly(ncp))
	return NC_EPERM;

    ncap = NC_attrarray0(ncp, varid);
    if(ncap == NULL)
	return NC_ENOTVAR;

    status = nc5_cktype(nc->mode, type);
    if(status != NC_NOERR)
	return status;

    if(memtype == NC_NAT) memtype = type;

    if(memtype != NC_CHAR && type == NC_CHAR)
	return NC_ECHAR;
    if(memtype == NC_CHAR && type != NC_CHAR)
	return NC_ECHAR;

    /* cast needed for braindead systems with signed size_t */
    if((unsigned long) nelems > X_INT_MAX) /* backward compat */
	return NC_EINVAL; /* Invalid nelems */

    if(nelems != 0 && value == NULL)
	return NC_EINVAL; /* Null arg */

    attrpp = NC5_findattr(ncap, name);

    /* 4 cases: exists X indef */

    if(attrpp != NULL) { /* name in use */
        if(!NC_indef(ncp)) {
	    const size_t xsz = nc5i_len_NC_attrV(type, nelems);
            attrp = *attrpp; /* convenience */
    
	    if(xsz > attrp->xsz) return NC_ENOTINDEFINE;
	    /* else, we can reuse existing without redef */
                    
	    attrp->xsz = xsz;
            attrp->type = type;
            attrp->nelems = nelems;

            if(nelems != 0) {
                void *xp = attrp->xvalue;
                status = dispatchput(&xp, nelems, (const void*)value, type, memtype);
            }
                       
            set_NC_hdirty(ncp);

            if(NC_doHsync(ncp)) {
	        const int lstatus = nc5i_sync(ncp);
                /*
                 * N.B.: potentially overrides NC_ERANGE
                 * set by nc5x_pad_putn_I$1
                 */
                if(lstatus != ENOERR) return lstatus;
            }

            return status;
        }
        /* else, redefine using existing array slot */
        old = *attrpp;
    } else {
        if(!NC_indef(ncp)) return NC_ENOTINDEFINE;

        if(ncap->nelems >= NC_MAX_ATTRS) return NC_EMAXATTS;
    }

    status = NC_check_name(name);
    if(status != NC_NOERR) return status;

    attrp = new_NC_attr(name, type, nelems);
    if(attrp == NULL) return NC_ENOMEM;

    if(nelems != 0) {
        void *xp = attrp->xvalue;
        status = dispatchput(&xp, nelems, (const void*)value, type, memtype);
    }

    if(attrpp != NULL) {
        assert(old != NULL);
        *attrpp = attrp;
        nc5i_free_NC_attr(old);
    } else {
        const int lstatus = incr_NC_attrarray(ncap, attrp);
        /*
         * N.B.: potentially overrides NC_ERANGE
         * set by nc5x_pad_putn_I$1
         */
        if(lstatus != NC_NOERR) {
           nc5i_free_NC_attr(attrp);
           return lstatus;
        }
    }
    return status;
}

int
NC5_get_att(
	int ncid,
	int varid,
	const char *name,
	void *value,
	nc_type memtype)
{
    int status;
    NC_attr *attrp;
    const void *xp;

#ifdef USE_PNETCDF
    NC* nc;
    nc_type xtype;

    status = NC_check_id(ncid, &nc);
    if(status != NC_NOERR) return status;

    if (((NC5_INFO*)nc->dispatchdata)->use_parallel) {

        status = NC5_inq_att(ncid,varid,name,&xtype,NULL);

        if(memtype == NC_NAT) memtype = xtype;

        switch (memtype) {
        case NC_CHAR:
            return ncmpi_get_att_text(nc->int_ncid, varid, name, (char*)value);
        case NC_BYTE:
            return ncmpi_get_att_schar(nc->int_ncid, varid, name, (signed char*)value);
        case NC_SHORT:
            return ncmpi_get_att_short(nc->int_ncid, varid, name, (short*)value);
        case NC_INT:
            return ncmpi_get_att_int(nc->int_ncid, varid, name, (int*)value);
        case NC_FLOAT:
            return ncmpi_get_att_float(nc->int_ncid, varid, name, (float*)value);
        case NC_DOUBLE:
            return ncmpi_get_att_double(nc->int_ncid, varid, name, (double*)value);
        case NC_UBYTE:
            return ncmpi_get_att_uchar(nc->int_ncid, varid, name, (unsigned char*)value);
        case NC_USHORT:
            return ncmpi_get_att_ushort(nc->int_ncid, varid, name, (unsigned short*)value);
        case NC_UINT:
            return ncmpi_get_att_uint(nc->int_ncid, varid, name, (unsigned int*)value);
        case NC_INT64:
            return ncmpi_get_att_longlong(nc->int_ncid, varid, name, (long long*)value);
        case NC_UINT64:
            return ncmpi_get_att_ulonglong(nc->int_ncid, varid, name, (unsigned long long*)value);
        default:
	    break;
        }
        return NC_EBADTYPE;
    }
#endif

    status = NC_lookupattr(ncid, varid, name, &attrp);
    if(status != NC_NOERR) return status;

    if(attrp->nelems == 0) return NC_NOERR;

    if(memtype == NC_NAT) memtype = attrp->type;

    if(memtype != NC_CHAR && attrp->type == NC_CHAR)
	return NC_ECHAR;
    if(memtype == NC_CHAR && attrp->type != NC_CHAR)
	return NC_ECHAR;

    xp = attrp->xvalue;
    switch (memtype) {
    case NC_CHAR:
        return nc5x_pad_getn_text(&xp, attrp->nelems , (char *)value);
    case NC_BYTE:
        return nc5x_pad_getn_Ischar(&xp,attrp->nelems,(schar*)value,attrp->type);
    case NC_SHORT:
        return nc5x_pad_getn_Ishort(&xp,attrp->nelems,(short*)value,attrp->type);
    case NC_INT:
          return nc5x_pad_getn_Iint(&xp,attrp->nelems,(int*)value,attrp->type);
    case NC_FLOAT:
        return nc5x_pad_getn_Ifloat(&xp,attrp->nelems,(float*)value,attrp->type);
    case NC_DOUBLE:
        return nc5x_pad_getn_Idouble(&xp,attrp->nelems,(double*)value,attrp->type);
    case NC_UBYTE: /* Synthetic */
        return nc5x_pad_getn_Iuchar(&xp, attrp->nelems , (uchar *)value, attrp->type);
    case NC_USHORT:
          return nc5x_pad_getn_Iushort(&xp,attrp->nelems,(ushort*)value,attrp->type);
    case NC_UINT:
          return nc5x_pad_getn_Iuint(&xp,attrp->nelems,(uint*)value,attrp->type);
    case NC_INT64:
          return nc5x_pad_getn_Ilonglong(&xp,attrp->nelems,(longlong*)value,attrp->type);
    case NC_UINT64:
          return nc5x_pad_getn_Iulonglong(&xp,attrp->nelems,(ulonglong*)value,attrp->type);
    case NC_NAT:
        return NC_EBADTYPE;
    default:
        break;
    }
    status =  NC_EBADTYPE;
    return status;
}

