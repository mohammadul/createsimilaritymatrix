/* ---------------------------------------------------------------------------
** This software is furnished "as is", without technical support,
** and with no warranty, express or implied, as to its usefulness for
** any purpose.
**
** nnfield2similarity.cpp
** Converts a nearest neighbour field to its equivalent similarity matrix
**
** Author: Sk. Mohammadul Haque
** Copyright (c) 2013 Sk. Mohammadul Haque
** For more details and updates, visit http://mohammadulhaque.alotspace.com
** -------------------------------------------------------------------------*/

#include <mex.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
#include <climits>
#include <stdint.h>
#ifdef __PARALLEL__
#include <omp.h>
#endif

#define DATA_TYPE int64_t

#define IS_REAL_2D_FULL_DOUBLE(P) (!mxIsComplex(P) && \
mxGetNumberOfDimensions(P)==2 && !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_3D_OR_4D_FULL_INT32(P) (!mxIsComplex(P) && \
(mxGetNumberOfDimensions(P)==3 || mxGetNumberOfDimensions(P)==4) && !mxIsEmpty(P) && !mxIsSparse(P) && mxIsInt32(P))
#define IS_REAL_SCALAR(P) (IS_REAL_2D_FULL_DOUBLE(P) && mxGetNumberOfElements(P)==1)
#define SIMILARITY_OUT plhs[0]
#define NNF_IN prhs[0]
#define WTTYPE_IN prhs[1]
#define SIGMA_IN prhs[2]
#define NORMALIZED_IN prhs[3]
#define SYMMETRIC_IN prhs[4]
#define PATCHSIZE_IN prhs[5]
#define ALPHA_IN prhs[6]
#define THRESHOLD_IN prhs[7]
#define MASK_IN prhs[8]
#define REFIMG_SIZE_IN prhs[9]

using namespace std;
typedef vector<DATA_TYPE>::const_iterator vecintiter;
struct ordering
{
    bool operator ()(pair<DATA_TYPE, vecintiter> const& a, pair<DATA_TYPE, vecintiter> const& b)
    {
        return *(a.second) < *(b.second);
    }
};

typedef vector<double>::const_iterator vecdoubleiter;

template <typename T>vector<T> sort_from_ref( vector<T> const& in, vector<pair<DATA_TYPE, vecintiter> > const& reference)
{
    vector<T> ret(in.size());
    DATA_TYPE const size = in.size();
    for (DATA_TYPE i = 0; i < size; ++i)ret[i] = in[reference[i].first];
    return ret;
}

template <typename T>vector<T> sort_from_ref( vector<T> const& in, vector<pair<DATA_TYPE, vecdoubleiter> > const& reference)
{
    vector<T> ret(in.size());
    DATA_TYPE const size = in.size();
    for(DATA_TYPE i = 0; i < size; ++i) ret[i] = in[reference[i].first];
    return ret;
}

template <typename T> pair<bool,double*> get_element(const T& row, const T& col, const mxArray *spmat)
{
    bool flag = 0;
    double *valmat = mxGetPr(spmat), *val = NULL;
    mwIndex *is = mxGetIr(spmat), *js = mxGetJc(spmat);
    for(T i = js[col]; i<js[col+1]; i++)
    {
        if(is[i]==row)
        {
            flag = true;
            val = &(valmat[i]);
            break;
        }
    }
    pair<bool, double*> ret = make_pair(flag, val);
    return ret;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double sigmasquared = 65536.0, *similarity = NULL, *texturemat = NULL, *refimg_size = NULL, alpha = 1.0, threshold = static_cast<double>(DBL_MAX);
    double **mask = NULL, *rmask = NULL;
    int32_t *nnf = NULL;
    bool normalized = false, texture = false, symmetric = false, same_image = false;
    int ndims, NNN = 0, wttype = 0;
    DATA_TYPE NNF_DIMS[4] = {0,0,0,0}, NELEMS = 0, winoffset = 0, patchsize_minus_one = 0, patchsize = 1;
    DATA_TYPE REFIMG_SIZE[2] = {0,0}, REFIMG_SIZE_REDUCED[2] = {0,0};
    vector<vector<DATA_TYPE> > idx_x;
    vector<vector<DATA_TYPE> > idx_y;
    vector<vector<double> > idx_d;
    const mwSize *NNF_DIMS_ORI = NULL;
    mwIndex *is = NULL, *js = NULL;

    /* check number of arguments */

    if(nrhs<1 || nrhs>9) mexErrMsgTxt("Wrong number of input arguments.");
    else if(nlhs>1) mexErrMsgTxt("Too many output arguments.");

    /* check all arguments one by one if valid or set default if empty */

    if(!IS_REAL_3D_OR_4D_FULL_INT32(NNF_IN)||mxIsEmpty(NNF_IN)) mexErrMsgTxt("NNF must be a real 3D or 4D full int32 array.");

    if(nrhs<8) threshold = static_cast<double>(DBL_MAX);
    else
    {
        if(mxIsEmpty(THRESHOLD_IN)) threshold = static_cast<double>(DBL_MAX);
        else
        {
            if(!IS_REAL_SCALAR(THRESHOLD_IN) ||(mxGetScalar(THRESHOLD_IN)<=0.0)) mexErrMsgTxt("THRESHOLD must be real positive scalar or empty.");
            threshold = mxGetScalar(THRESHOLD_IN);
        }
    }

    if(nrhs<7) alpha = 1.0;
    else
    {
        if(mxIsEmpty(ALPHA_IN)) alpha = 1.0;
        else
        {
            if(!IS_REAL_SCALAR(ALPHA_IN)) mexErrMsgTxt("ALPHA must be real scalar or empty.");
            alpha = mxGetScalar(ALPHA_IN);
        }
    }

    if(nrhs<6) winoffset = 0;
    else
    {
        if(mxIsEmpty(PATCHSIZE_IN)) winoffset = 0;
        else
        {
            if(!IS_REAL_SCALAR(PATCHSIZE_IN)||(mxGetScalar(PATCHSIZE_IN)<1.0)) mexErrMsgTxt("PATCHSIZE must be real positive scalar or empty.");
            winoffset = static_cast<DATA_TYPE>(((mxGetScalar(PATCHSIZE_IN))-1)/2);
        }
    }
    patchsize_minus_one = 2*winoffset;
    patchsize = patchsize_minus_one+1;

    if(nrhs<5) symmetric = false;
    else
    {
        if(mxIsEmpty(SYMMETRIC_IN)) symmetric = false;
        else
        {
            if(!IS_REAL_SCALAR(SYMMETRIC_IN)) mexErrMsgTxt("SYMMETRIC must be 0, 1 or empty.");
            symmetric = static_cast<bool>(mxGetScalar(SYMMETRIC_IN));
        }
    }

    if(nrhs<4) normalized = false;
    else
    {
        if(mxIsEmpty(NORMALIZED_IN)) normalized = false;
        else
        {
            if(!IS_REAL_SCALAR(NORMALIZED_IN)) mexErrMsgTxt("NORMALIZED must be real scalar or empty.");
            normalized = static_cast<bool>(mxGetScalar(NORMALIZED_IN));
        }
    }

    if(nrhs<3) sigmasquared = 32768.0; /* 65536.0/2 */
    else
    {
        if((!IS_REAL_SCALAR(SIGMA_IN))&&(!IS_REAL_2D_FULL_DOUBLE(SIGMA_IN))) mexErrMsgTxt("SIGMA must be real scalar, matrix or empty.");
        if(IS_REAL_SCALAR(SIGMA_IN)) texture = false;
        else texture = true;

        if(mxIsEmpty(SIGMA_IN))
        {
            texture = false;
            sigmasquared = 32768.0; /* 65536.0/2 */
        }
        else if(texture==false)
        {
            sigmasquared = mxGetScalar(SIGMA_IN);
            sigmasquared *= sigmasquared;
            sigmasquared *= 2;
        }
        else texturemat = mxGetPr(SIGMA_IN);
    }

    if(nrhs<2) wttype = 0;
    else
    {
        if(mxIsEmpty(WTTYPE_IN)) wttype = 0;
        else
        {
            if(!IS_REAL_SCALAR(WTTYPE_IN)) mexErrMsgTxt("WTTYPE must be real scalar or empty.");
            wttype = mxGetScalar(WTTYPE_IN);
        }
    }

    if((wttype==1) && (texture ==1)) wttype = 2;

    /* get pointer to NN-Field and its dimensions */

    nnf = reinterpret_cast<int32_t*>(mxGetPr(NNF_IN));
    NNF_DIMS_ORI = mxGetDimensions(NNF_IN);
    ndims = mxGetNumberOfDimensions(NNF_IN);
    switch(ndims)
    {
    case 3:
        NNF_DIMS[0] = static_cast<DATA_TYPE>(NNF_DIMS_ORI[0]);
        NNF_DIMS[1] = static_cast<DATA_TYPE>(NNF_DIMS_ORI[1]);
        NNF_DIMS[2] = static_cast<DATA_TYPE>(NNF_DIMS_ORI[2]);
        NNF_DIMS[3] = 1;
        NNN = 1;
        break;
    case 4:
        NNF_DIMS[0] = static_cast<DATA_TYPE>(NNF_DIMS_ORI[0]);
        NNF_DIMS[1] = static_cast<DATA_TYPE>(NNF_DIMS_ORI[1]);
        NNF_DIMS[2] = static_cast<DATA_TYPE>(NNF_DIMS_ORI[2]);
        NNF_DIMS[3] = static_cast<DATA_TYPE>(NNF_DIMS_ORI[3]);
        NNN = NNF_DIMS[3];
        break;
    default:
        mexErrMsgTxt("NNF size error.");
        break;
    }
    size_t matelem = NNF_DIMS[0]*NNF_DIMS[1];
    size_t refmatelem;

    idx_x.resize(NNN);
    idx_y.resize(NNN);
    idx_d.resize(NNN);

    if(nrhs<10)
    {
        REFIMG_SIZE[0] = NNF_DIMS[0];
        REFIMG_SIZE[1] = NNF_DIMS[1];
    }
    else
    {
        if(mxIsEmpty(REFIMG_SIZE_IN))
        {
            REFIMG_SIZE[0] = NNF_DIMS[0];
            REFIMG_SIZE[1] = NNF_DIMS[1];
        }
        else
        {
            if(!IS_REAL_2D_FULL_DOUBLE(REFIMG_SIZE_IN) || (mxGetM(REFIMG_SIZE_IN)*mxGetN(REFIMG_SIZE_IN)!=2)) mexErrMsgTxt("REFIMG_SIZE must be real positive 2X1 matrix or empty.");
            refimg_size = mxGetPr(REFIMG_SIZE_IN);
            REFIMG_SIZE[0] = static_cast<DATA_TYPE>(refimg_size[0]);
            REFIMG_SIZE[1] = static_cast<DATA_TYPE>(refimg_size[1]);
        }
    }
    refmatelem = REFIMG_SIZE[0]*REFIMG_SIZE[1];
    if((REFIMG_SIZE[0] == NNF_DIMS[0]) && (REFIMG_SIZE[1] == NNF_DIMS[1])) same_image = true;
    else same_image = false;

    REFIMG_SIZE_REDUCED[0] = REFIMG_SIZE[0]-2*winoffset;
    REFIMG_SIZE_REDUCED[1] = REFIMG_SIZE[1]-2*winoffset;

    if(nrhs<9||mxIsEmpty(MASK_IN))
    {
        mask = new double*[patchsize];
        for(int i=0; i<patchsize; ++i)
        {
            mask[i] = new double[patchsize];
            for(int j=0; j<patchsize; ++j)
            {
                if((i==winoffset) && (j==winoffset)) mask[i][j] = 1.0;
                else mask[i][j] = 0.0;
            }
        }
    }
    else
    {
        if(!IS_REAL_2D_FULL_DOUBLE(MASK_IN)) mexErrMsgTxt("MASK must be real full array or empty.");
        rmask = mxGetPr(MASK_IN);
        mask = new double*[patchsize];
        for(int i=0; i<patchsize; ++i)
        {
            mask[i] = new double[patchsize];
            for(int j=0; j<patchsize; ++j) mask[i][j] = rmask[i*patchsize+j];
        }
    }

    /* based on the conditions */

    switch(wttype)
    {
    case 0: /* uniform weights without texture */
    {
#ifdef __PARALLEL__
#pragma omp parallel for shared(NNF_DIMS, nnf, idx_x, idx_y, idx_d)
#endif
        for(DATA_TYPE k=0; k<NNN; ++k)
        {
            DATA_TYPE tsz = NNF_DIMS[0]*NNF_DIMS[1];
            DATA_TYPE idx_start_y = tsz*NNF_DIMS[2]*k;
            DATA_TYPE idx_start_x = tsz + idx_start_y;
            DATA_TYPE idx_start_d = tsz + idx_start_x;

            for(DATA_TYPE j=0; j<static_cast<DATA_TYPE>(NNF_DIMS[1]-patchsize_minus_one); ++j)
            {
                DATA_TYPE idx_curr = (NNF_DIMS[0]*j);
                DATA_TYPE idx_start_xx = idx_start_x+idx_curr;
                DATA_TYPE idx_start_yy = idx_start_y+idx_curr;
                DATA_TYPE idx_start_dd = idx_start_d+idx_curr;
                for(DATA_TYPE i=0; i<static_cast<DATA_TYPE>(NNF_DIMS[0]-patchsize_minus_one); ++i)
                {
                    DATA_TYPE curr_x, curr_y, curr_d;
                    curr_d = nnf[idx_start_dd+i];
                    for (DATA_TYPE mj=0; mj<=patchsize_minus_one; ++mj)
                    {
                        curr_y = nnf[idx_start_yy+i]+mj;
                        for (DATA_TYPE mi=0; mi<=patchsize_minus_one; ++mi)
                        {
                            curr_x = nnf[idx_start_xx+i]+mi;
                            if(mask[mi][mj]>0.0 && (curr_x>=0 && curr_x<REFIMG_SIZE[0] && curr_y>=0 && curr_y<REFIMG_SIZE[1] )&& (!same_image ||((idx_curr+i+(mj*NNF_DIMS[0]+mi))!=(curr_x+curr_y*REFIMG_SIZE[0]))) && (static_cast<double>(curr_d)<threshold))
                            {
                                idx_x[k].push_back(idx_curr+i+(mj*NNF_DIMS[0]+mi));
                                idx_y[k].push_back(curr_x+curr_y*REFIMG_SIZE[0]);
                                idx_d[k].push_back(mask[mi][mj]);

                                if(symmetric && same_image)
                                {
                                    idx_y[k].push_back(idx_curr+i+(mj*NNF_DIMS[0]+mi));
                                    idx_x[k].push_back(curr_x+curr_y*REFIMG_SIZE[0]);
                                    idx_d[k].push_back(mask[mi][mj]);
                                }
                            }
                        }
                    }
                }
            }
            if(k==0 && same_image)
            {
                for(DATA_TYPE j=winoffset; j<static_cast<DATA_TYPE>(NNF_DIMS[1]-winoffset); ++j)
                {
                    DATA_TYPE idx_curr = (NNF_DIMS[0]*j);
                    for(DATA_TYPE i=winoffset; i<static_cast<DATA_TYPE>(NNF_DIMS[0]-winoffset); ++i)
                    {
                        idx_x[k].push_back(idx_curr+i);
                        idx_y[k].push_back(idx_curr+i);
                        idx_d[k].push_back(alpha);
                    }
                }
            }
        }
        break;
    }

    case 1: /* gaussian weights without texture */
    {
#ifdef __PARALLEL__
#pragma omp parallel for shared(NNF_DIMS, nnf, idx_x, idx_y, idx_d)
#endif
        for(DATA_TYPE k=0; k<NNN; ++k)
        {
            DATA_TYPE tsz = NNF_DIMS[0]*NNF_DIMS[1];
            DATA_TYPE idx_start_y = tsz*NNF_DIMS[2]*k;
            DATA_TYPE idx_start_x = tsz + idx_start_y;
            DATA_TYPE idx_start_d = tsz + idx_start_x;
            double expval;

            for(DATA_TYPE j=0; j<static_cast<DATA_TYPE>(NNF_DIMS[1]-patchsize_minus_one); ++j)
            {
                DATA_TYPE idx_curr = (NNF_DIMS[0]*j);
                DATA_TYPE idx_start_xx = idx_start_x+idx_curr;
                DATA_TYPE idx_start_yy = idx_start_y+idx_curr;
                DATA_TYPE idx_start_dd = idx_start_d+idx_curr;
                for(DATA_TYPE i=0; i<static_cast<DATA_TYPE>(NNF_DIMS[0]-patchsize_minus_one); ++i)
                {
                    DATA_TYPE curr_x, curr_y, curr_d;
                    curr_d = nnf[idx_start_dd+i];
                    for (DATA_TYPE mj=0; mj<=patchsize_minus_one; ++mj)
                    {
                        curr_y = nnf[idx_start_yy+i]+mj;
                        for (DATA_TYPE mi=0; mi<=patchsize_minus_one; ++mi)
                        {
                            curr_x = nnf[idx_start_xx+i]+mi;
                            if(mask[mi][mj]>0.0 && (curr_x>=0 && curr_x<REFIMG_SIZE[0] && curr_y>=0 && curr_y<REFIMG_SIZE[1] )&& (!same_image ||((idx_curr+i+(mj*NNF_DIMS[0]+mi))!=(curr_x+curr_y*REFIMG_SIZE[0]))) && (static_cast<double>(curr_d)<threshold))
                            {
                                expval = exp(-static_cast<double>(curr_d)/sigmasquared)*mask[mi][mj];
                                idx_x[k].push_back(idx_curr+i+(mj*NNF_DIMS[0]+mi));
                                idx_y[k].push_back(curr_x+curr_y*REFIMG_SIZE[0]);
                                idx_d[k].push_back(expval);

                                if(symmetric && same_image)
                                {
                                    idx_y[k].push_back(idx_curr+i+(mj*NNF_DIMS[0]+mi));
                                    idx_x[k].push_back(curr_x+curr_y*REFIMG_SIZE[0]);
                                    idx_d[k].push_back(expval);
                                }
                            }
                        }
                    }
                }
            }
            if(k==0 && same_image)
            {
                for(DATA_TYPE j=winoffset; j<static_cast<DATA_TYPE>(NNF_DIMS[1]-winoffset); ++j)
                {
                    DATA_TYPE idx_curr = (NNF_DIMS[0]*j);
                    for(DATA_TYPE i=winoffset; i<static_cast<DATA_TYPE>(NNF_DIMS[0]-winoffset); ++i)
                    {
                        idx_x[k].push_back(idx_curr+i);
                        idx_y[k].push_back(idx_curr+i);
                        idx_d[k].push_back(alpha);
                    }
                }
            }
        }
        break;
    }

    case 2: /* gaussian weight with texture */
    {
        /* check if texturemat is valid */

        if((NNF_DIMS[0]!=static_cast<DATA_TYPE>(mxGetM(SIGMA_IN))) || (NNF_DIMS[1]!=static_cast<DATA_TYPE>(mxGetN(SIGMA_IN)))) mexErrMsgTxt("SIGMA TEXTURE size error.");
#ifdef __PARALLEL__
#pragma omp parallel for shared(NNF_DIMS, nnf, idx_x, idx_y, idx_d)
#endif
        for(DATA_TYPE k=0; k<NNN; ++k)
        {
            DATA_TYPE curr_x, curr_y, curr_d;
            DATA_TYPE tsz = NNF_DIMS[0]*NNF_DIMS[1];
            DATA_TYPE idx_start_y = tsz*NNF_DIMS[2]*k;
            DATA_TYPE idx_start_x = tsz + idx_start_y;
            DATA_TYPE idx_start_d = tsz + idx_start_x;
            double expval;

            for(DATA_TYPE j=0; j<static_cast<DATA_TYPE>(NNF_DIMS[1]-patchsize_minus_one); ++j)
            {
                DATA_TYPE idx_curr = (NNF_DIMS[0]*j);
                DATA_TYPE idx_start_xx = idx_start_x+idx_curr;
                DATA_TYPE idx_start_yy = idx_start_y+idx_curr;
                DATA_TYPE idx_start_dd = idx_start_d+idx_curr;
                for(DATA_TYPE i=0; i<static_cast<DATA_TYPE>(NNF_DIMS[0]-patchsize_minus_one); ++i)
                {
                    DATA_TYPE curr_x, curr_y, curr_d;
                    curr_d = nnf[idx_start_dd+i];
                    for (DATA_TYPE mj=0; mj<=patchsize_minus_one; ++mj)
                    {
                        curr_y = nnf[idx_start_yy+i]+mj;
                        for (DATA_TYPE mi=0; mi<=patchsize_minus_one; ++mi)
                        {
                            curr_x = nnf[idx_start_xx+i]+mi;
                            if(mask[mi][mj]>0.0 && (curr_x>=0 && curr_x<REFIMG_SIZE[0] && curr_y>=0 && curr_y<REFIMG_SIZE[1] )&& (!same_image ||((idx_curr+i+(mj*NNF_DIMS[0]+mi))!=(curr_x+curr_y*REFIMG_SIZE[0]))) && (static_cast<double>(curr_d)<threshold))
                            {
                                expval = exp(-static_cast<double>(curr_d)*texturemat[idx_curr+i]/sigmasquared)*mask[mi][mj]; /* texture_mat */
                                idx_x[k].push_back(idx_curr+i+(mj*NNF_DIMS[0]+mi));
                                idx_y[k].push_back(curr_x+curr_y*REFIMG_SIZE[0]);
                                idx_d[k].push_back(expval);

                                if(symmetric && same_image)
                                {
                                    idx_y[k].push_back(idx_curr+i+(mj*NNF_DIMS[0]+mi));
                                    idx_x[k].push_back(curr_x+curr_y*REFIMG_SIZE[0]);
                                    idx_d[k].push_back(expval);
                                }
                            }
                        }
                    }
                }
            }
            if(k==0 && same_image)
            {
                for(DATA_TYPE j=winoffset; j<static_cast<DATA_TYPE>(NNF_DIMS[1]-winoffset); ++j)
                {
                    DATA_TYPE idx_curr = (NNF_DIMS[0]*j);
                    for(DATA_TYPE i=winoffset; i<static_cast<DATA_TYPE>(NNF_DIMS[0]-winoffset); ++i)
                    {
                        idx_x[k].push_back(idx_curr+i);
                        idx_y[k].push_back(idx_curr+i);
                        idx_d[k].push_back(alpha);
                    }
                }
            }
        }
        break;
    }

    default:
        mexErrMsgTxt("WTTYPE is invalid.");
        break;
    }

    /* collect all the entries in NNF in 3 single vectors */

    for (int k=1; k<NNN; ++k)
    {
        idx_x[0].insert(idx_x[0].end(), idx_x[k].begin(), idx_x[k].end());
        idx_y[0].insert(idx_y[0].end(), idx_y[k].begin(), idx_y[k].end());
        idx_d[0].insert(idx_d[0].end(), idx_d[k].begin(), idx_d[k].end());
        idx_x[k].clear();
        idx_y[k].clear();
        idx_d[k].clear();
    }

    /* sort the 3 vectors in order of location in matrix */

    vector<DATA_TYPE> tempvec(idx_y[0]);
    transform(tempvec.begin(), tempvec.end(), tempvec.begin(), bind1st(multiplies<DATA_TYPE>(), matelem));
    transform(tempvec.begin(), tempvec.end(), idx_x[0].begin(), tempvec.begin(), plus<DATA_TYPE>());

    vector<pair<DATA_TYPE, vecintiter> > order(tempvec.size());

    DATA_TYPE n = 0;
    for (vecintiter it = tempvec.begin(); it != tempvec.end(); ++it, ++n) order[n] = make_pair(n, it);
    sort(order.begin(), order.end(), ordering());

    vector<DATA_TYPE> sorted_idx_x = sort_from_ref(idx_x[0], order), trimmed_idx_x;
    vector<DATA_TYPE> sorted_idx_y = sort_from_ref(idx_y[0], order), trimmed_idx_y;
    vector<double> sorted_idx_d = sort_from_ref(idx_d[0], order), trimmed_idx_d;

    DATA_TYPE oldidx = -1, tmp_x, tmp_y;
    double maxval = 0;

    size_t m;
    for(m=0; m<order.size(); ++m)
    {
        if(oldidx==*(order[m].second))
        {
            /*             if(maxval<sorted_idx_d[m])
             *             {
             *                 tmp_x = sorted_idx_x[m];
             *                 tmp_y = sorted_idx_y[m];
             *             }
             */
            maxval += sorted_idx_d[m];

        }
        else
        {
            if(m>0)
            {
                trimmed_idx_x.push_back(tmp_x);
                trimmed_idx_y.push_back(tmp_y);
                trimmed_idx_d.push_back(maxval);
                maxval = 0;
            }
            maxval = sorted_idx_d[m];
            tmp_x = sorted_idx_x[m];
            tmp_y = sorted_idx_y[m];

        }
        oldidx = *(order[m].second);
    }
    if(m>0)
    {
        trimmed_idx_x.push_back(tmp_x);
        trimmed_idx_y.push_back(tmp_y);
        trimmed_idx_d.push_back(maxval);
    }

    NELEMS = trimmed_idx_d.size();

    /* first check if ref image size is large enough */

    if((NELEMS>0) && ((trimmed_idx_y[NELEMS-1])>=(REFIMG_SIZE[1]*REFIMG_SIZE[0]))) mexErrMsgTxt("REF IMG SIZE is smaller than required.");

    /* create the output similarity matrix */

    SIMILARITY_OUT = mxCreateSparse(matelem, refmatelem, NELEMS, mxREAL);
    if(SIMILARITY_OUT == NULL) mexErrMsgTxt("Cannot allocate enough space.");
    similarity = mxGetPr(SIMILARITY_OUT);
    is = mxGetIr(SIMILARITY_OUT);
    js = mxGetJc(SIMILARITY_OUT);
    DATA_TYPE tmp = -1, col_counter = 0;

    size_t k;

    for(k=0; k<trimmed_idx_d.size(); ++k)
    {
        is[k] = trimmed_idx_x[k];
        if((trimmed_idx_y[k])>tmp)
        {
            for(DATA_TYPE i=1+tmp; i<=trimmed_idx_y[k] ; ++i)
            {
                js[col_counter] = k;
                ++col_counter;
            }
            tmp = trimmed_idx_y[k];
        }
        similarity[k] = trimmed_idx_d[k];
    }
    for(size_t i=col_counter; i<=refmatelem; ++i) js[i] = k;

    /* normalize if any */

    if(normalized)
    {
        vector<double> rowsum(matelem, 0.0) ;
        for(DATA_TYPE i=0; i<NELEMS; ++i)
        {
            rowsum[is[i]] += trimmed_idx_d[i];
        }
        for(DATA_TYPE i=0; i<NELEMS; ++i)
        {
            if(rowsum[is[i]]>0.0) similarity[i] /= rowsum[is[i]];
        }
    }
}

