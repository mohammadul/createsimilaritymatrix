%NNFIELD2SIMILARITY - Converts NN-field to similarity matrix
%   sm = nnfield2similarity(nnf, wttype, sigma, normalized, make_symm, patchsize, alpha, threshold, mask, ref_img_size)
%       nnf - nearest neighbour field
%       wttype - 0 - one-type (default) , 1 - gaussian(exp)
%       sigma - exponential wt sigma (default - 256)
%       normalized - 0/1 (default - 0)
%       make_symm - 0/1 (default - 0)
%       patchsize - size of search patch (default - [1 1])
%       alpha - for diagonal elements (default - 1.0)
%       threshold - for limiting similar patches (default - DBL_MAX)
%       mask - patch weight mask (default - center dot)
%       ref_img_size - reference image size (2X1) (default - same as nnf)
%       sm - output similarity matrix
%
%   Author: Sk. Mohammadul Haque
%   Copyright (c) 2013 Sk. Mohammadul Haque
%   Website: http://mohammadulhaque.alotspace.com
%
