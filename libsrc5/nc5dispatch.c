/*********************************************************************
 *   Copyright 1993, UCAR/Unidata
 *   See netcdf/COPYRIGHT file for copying and redistribution conditions.
 *   $Header: /upc/share/CVS/netcdf-3/libsrc4/NC5dispatch.c,v 1.5 2010/05/27 02:19:37 dmh Exp $
 *********************************************************************/

/* WARNING: Order of mpi.h, nc.h, and pnetcdf.h is important */
#include "config.h"
#include <stdlib.h>

#ifdef USE_PARALLEL
#include <mpi.h>
#endif

#include "netcdf.h"
#include "nc5internal.h"
#include "nc5dispatch.h"

#ifndef NC_CONTIGUOUS
#define NC_CONTIGUOUS 1
#endif

#ifndef NC_ENOTNC4
#define NC_ENOTNC4 (-111)
#endif

#ifndef NC_ENOGRP
#define NC_ENOGRP (-125)
#endif

#ifndef NC_STRING
#define NC_STRING (12)
#endif


#ifdef USE_NETCDF4
static int NC5_show_metadata(int);
static int NC5_inq_unlimdims(int,int*,int*);
static int NC5_inq_ncid(int,const char*,int*);
static int NC5_inq_grps(int,int*,int*);
static int NC5_inq_grpname(int,char*);
static int NC5_inq_grpname_full(int,size_t*,char*);
static int NC5_inq_grp_parent(int,int*);
static int NC5_inq_grp_full_ncid(int,const char*,int*);
static int NC5_inq_varids(int,int* nvars,int*);
static int NC5_inq_dimids(int,int* ndims,int*,int);
static int NC5_inq_typeids(int,int* ntypes,int*);
static int NC5_inq_type_equal(int,nc_type,int,nc_type,int*);
static int NC5_def_grp(int,const char*,int*);
static int NC5_rename_grp(int,const char*);
static int NC5_inq_user_type(int,nc_type,char*,size_t*,nc_type*,size_t*,int*);
static int NC5_inq_typeid(int,const char*,nc_type*);
static int NC5_def_compound(int,size_t,const char*,nc_type*);
static int NC5_insert_compound(int,nc_type,const char*,size_t,nc_type);
static int NC5_insert_array_compound(int,nc_type,const char*,size_t,nc_type,int,const int*);
static int NC5_inq_compound_field(int,nc_type,int,char*,size_t*,nc_type*,int*,int*);
static int NC5_inq_compound_fieldindex(int,nc_type,const char*,int*);
static int NC5_def_vlen(int,const char*,nc_type base_typeid,nc_type*);
static int NC5_put_vlen_element(int,int,void*,size_t,const void*);
static int NC5_get_vlen_element(int,int,const void*,size_t*,void*);
static int NC5_def_enum(int,nc_type,const char*,nc_type*);
static int NC5_insert_enum(int,nc_type,const char*,const void*);
static int NC5_inq_enum_member(int,nc_type,int,char*,void*);
static int NC5_inq_enum_ident(int,nc_type,long long,char*);
static int NC5_def_opaque(int,size_t,const char*,nc_type*);
static int NC5_def_var_deflate(int,int,int,int,int);
static int NC5_def_var_fletcher32(int,int,int);
static int NC5_def_var_chunking(int,int,int,const size_t*);
static int NC5_def_var_fill(int,int,int,const void*);
static int NC5_def_var_endian(int,int,int);
static int NC5_set_var_chunk_cache(int,int,size_t,size_t,float);
static int NC5_get_var_chunk_cache(int,int,size_t*,size_t*,float*);
#endif /*USE_NETCDF4*/


NC_Dispatch NC5_dispatcher = {

NC_DISPATCH_NC5,

NC5_create,
NC5_open,

NC5_redef,
NC5__enddef,
NC5_sync,
NC5_abort,
NC5_close,
NC5_set_fill,
NC5_inq_base_pe,
NC5_set_base_pe,
NC5_inq_format,
NC5_inq_format_extended,

NC5_inq,
NC5_inq_type,

NC5_def_dim,
NC5_inq_dimid,
NC5_inq_dim,
NC5_inq_unlimdim,
NC5_rename_dim,

NC5_inq_att,
NC5_inq_attid,
NC5_inq_attname,
NC5_rename_att,
NC5_del_att,
NC5_get_att,
NC5_put_att,

NC5_def_var,
NC5_inq_varid,
NC5_rename_var,
NC5_get_vara,
NC5_put_vara,
NC5_get_vars,
NC5_put_vars,
NC5_get_varm,
NC5_put_varm,

NC5_inq_var_all,

NC5_var_par_access,

#ifdef USE_NETCDF4
NC5_show_metadata,
NC5_inq_unlimdims,

NC5_inq_ncid,
NC5_inq_grps,
NC5_inq_grpname,
NC5_inq_grpname_full,
NC5_inq_grp_parent,
NC5_inq_grp_full_ncid,
NC5_inq_varids,
NC5_inq_dimids,
NC5_inq_typeids,
NC5_inq_type_equal,
NC5_def_grp,
NC5_rename_grp,
NC5_inq_user_type,
NC5_inq_typeid,

NC5_def_compound,
NC5_insert_compound,
NC5_insert_array_compound,
NC5_inq_compound_field,
NC5_inq_compound_fieldindex,
NC5_def_vlen,
NC5_put_vlen_element,
NC5_get_vlen_element,
NC5_def_enum,
NC5_insert_enum,
NC5_inq_enum_member,
NC5_inq_enum_ident,
NC5_def_opaque,
NC5_def_var_deflate,
NC5_def_var_fletcher32,
NC5_def_var_chunking,
NC5_def_var_fill,
NC5_def_var_endian,
NC5_set_var_chunk_cache,
NC5_get_var_chunk_cache,

#endif /*USE_NETCDF4*/

};

NC_Dispatch* NC5_dispatch_table = NULL; /* moved here from ddispatch.c */

int
NC5_initialize(void)
{
    NC5_dispatch_table = &NC5_dispatcher;
    return NC_NOERR;
}


#ifdef USE_NETCDF4

static int
NC5_show_metadata(int ncid)
{
    return NC_NOERR;
}

static int
NC5_inq_unlimdims(int ncid, int *ndimsp, int *unlimdimidsp)
{
    int retval;
    int unlimid;

    if((retval = NC5_inq_unlimdim(ncid, &unlimid)))
        return retval;
    if(unlimid != -1) {
        if(ndimsp) *ndimsp = 1;
        if(unlimdimidsp)
            unlimdimidsp[0] = unlimid;
    } else
        if(ndimsp) *ndimsp = 0;
    return NC_NOERR;
}

static int
NC5_inq_type_equal(int ncid1, nc_type typeid1, int ncid2, nc_type typeid2, int* equalp)
{
    /* Check input. */
    if(equalp == NULL) return NC_NOERR;

    if (typeid1 <= NC_NAT || typeid2 <= NC_NAT)
       return NC_EINVAL;

    *equalp = 0; /* assume */

    /* If one is atomic, and the other user-defined, the types are not equal */
    if ((typeid1 <= NC_STRING && typeid2 > NC_STRING) ||
        (typeid2 <= NC_STRING && typeid1 > NC_STRING)) {
        if (equalp) *equalp = 0;
        return NC_NOERR;
    }

    /* If both are atomic types, the answer is easy. */
    if (typeid1 <= ATOMICTYPEMAX) {
        if (equalp) {
            if (typeid1 == typeid2)
                *equalp = 1;
            else
                *equalp = 0;
        }
        return NC_NOERR;
    }
    return NC_NOERR;
}

static int
NC5_def_grp(int parent_ncid, const char *name, int *new_ncid)
{
    return NC_ENOTNC4;
}

static int
NC5_rename_grp(int ncid, const char *name)
{
    return NC_ENOTNC4;
}

static int
NC5_inq_ncid(int ncid, const char *name, int *grp_ncid)
{
    if(grp_ncid) *grp_ncid = ncid;
    return NC_NOERR;
}

static int
NC5_inq_grps(int ncid, int *numgrps, int *ncids)
{
    if(numgrps)
       *numgrps = 0;
    return NC_NOERR;
}

static int
NC5_inq_grpname(int ncid, char *name)
{
    if(name)
        strcpy(name, "/");
    return NC_NOERR;
}

static int
NC5_inq_grpname_full(int ncid, size_t *lenp, char *full_name)
{
    if(full_name)
        strcpy(full_name, "/");
    if(lenp) *lenp = 1;
    return NC_NOERR;
}

static int
NC5_inq_grp_parent(int ncid, int *parent_ncid)
{
    return NC_ENOGRP;
}

static int
NC5_inq_grp_full_ncid(int ncid, const char *full_name, int *grp_ncid)
{
    return NC_ENOGRP;
}

static int
NC5_inq_varids(int ncid, int *nvarsp, int *varids)
{
    int retval,v,nvars;
    /* This is, effectively, a netcdf-3 file, there is only one group, the root
        group, and its vars have ids 0 thru nvars - 1. */
    if((retval = NC5_inq(ncid, NULL, &nvars, NULL, NULL)))
        return retval;
    if(nvarsp) *nvarsp = nvars;
    if(varids)
        for (v = 0; v < nvars; v++)
            varids[v] = v;
    return NC_NOERR;
}

static int
NC5_inq_dimids(int ncid, int *ndimsp, int *dimids, int include_parents)
{
    int retval,d,ndims;
    /* If this is like a netcdf-3 file, then the dimids are going to be 0
       thru ndims-1, so just provide them. */
    if((retval = NC5_inq(ncid, &ndims,  NULL, NULL, NULL)))
        return retval;
    if(ndimsp) *ndimsp = ndims;
    if(dimids)
        for (d = 0; d < ndims; d++)
            dimids[d] = d;
    return NC_NOERR;
}

static int
NC5_inq_typeid(int ncid, const char *name, nc_type *typeidp)
{
    int i;
    for (i = 0; i <= ATOMICTYPEMAX; i++)
        if(!strcmp(name, NC_atomictypename(i))) {
            if(typeidp) *typeidp = i;
                return NC_NOERR;
        }
    return NC_ENOTNC4;
}

static int
NC5_inq_typeids(int ncid, int *ntypes, int *typeids)
{
    if(ntypes) *ntypes = 0;
    return NC_NOERR;
}

static int
NC5_inq_user_type(int ncid, nc_type typeid, char *name, size_t *size,
		 nc_type *base_nc_typep, size_t *nfieldsp, int *classp)
{
    return NC_ENOTNC4;
}

static int
NC5_def_compound(int ncid, size_t size, const char *name, nc_type *typeidp)
{
    return NC_ENOTNC4;
}

static int
NC5_insert_compound(int ncid, nc_type typeid, const char *name, size_t offset,
                    nc_type field_typeid)
{
    return NC_ENOTNC4;
}

static int
NC5_insert_array_compound(int ncid, nc_type typeid, const char *name,
			 size_t offset, nc_type field_typeid,
			 int ndims, const int *dim_sizes)
{
    return NC_ENOTNC4;
}


static int
NC5_inq_compound_field(int ncid, nc_type typeid, int fieldid, char *name,
		      size_t *offsetp, nc_type *field_typeidp, int *ndimsp,
		      int *dim_sizesp)
{
    return NC_ENOTNC4;
}

static int
NC5_inq_compound_fieldindex(int ncid, nc_type typeid, const char *name, int *fieldidp)
{
    return NC_ENOTNC4;
}

static int
NC5_def_opaque(int ncid, size_t datum_size, const char *name, nc_type* xtypep)
{
    return NC_ENOTNC4;
}

static int
NC5_def_vlen(int ncid, const char *name, nc_type base_typeid, nc_type* xtypep)
{
    return NC_ENOTNC4;
}

static int
NC5_def_enum(int ncid, nc_type base_typeid, const char *name,
	    nc_type *typeidp)
{
    return NC_ENOTNC4;
}

static int
NC5_inq_enum_ident(int ncid, nc_type xtype, long long value, char *identifier)
{
    return NC_ENOTNC4;
}

static int
NC5_inq_enum_member(int ncid, nc_type typeid, int idx, char *identifier,
		   void *value)
{
    return NC_ENOTNC4;
}

static int
NC5_insert_enum(int ncid, nc_type typeid, const char *identifier,
	       const void *value)
{
    return NC_ENOTNC4;
}

static int
NC5_put_vlen_element(int ncid, int typeid, void *vlen_element,
		    size_t len, const void *data)
{
    return NC_ENOTNC4;
}

static int
NC5_get_vlen_element(int ncid, int typeid, const void *vlen_element,
		    size_t *len, void *data)
{
    return NC_ENOTNC4;
}

static int
NC5_set_var_chunk_cache(int ncid, int varid, size_t size, size_t nelems, float preemption)
{
    return NC_ENOTNC4;
}

static int
NC5_get_var_chunk_cache(int ncid, int varid, size_t *sizep, size_t *nelemsp, float *preemptionp)
{
    return NC_ENOTNC4;
}

static int
NC5_def_var_deflate(int ncid, int varid, int shuffle, int deflate,
		   int deflate_level)
{
    return NC_ENOTNC4;
}

static int
NC5_def_var_fletcher32(int ncid, int varid, int fletcher32)
{
    return NC_ENOTNC4;
}

static int
NC5_def_var_chunking(int ncid, int varid, int contiguous, const size_t *chunksizesp)
{
    return NC_ENOTNC4;
}

static int
NC5_def_var_fill(int ncid, int varid, int no_fill, const void *fill_value)
{
    return NC_ENOTNC4;
}

static int
NC5_def_var_endian(int ncid, int varid, int endianness)
{
    return NC_ENOTNC4;
}

#endif /*USE_NETCDF4*/

