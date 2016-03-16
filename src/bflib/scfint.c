/* scfint.c (interface to Schur-complement-based factorization) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Copyright (C) 2013-2014 Andrew Makhorin, Department for Applied
*  Informatics, Moscow Aviation Institute, Moscow, Russia. All rights
*  reserved. E-mail: <mao@gnu.org>.
*
*  GLPK is free software: you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  GLPK is distributed in the hope that it will be useful, but WITHOUT
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
*  License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with GLPK. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#include "env.h"
#include "scfint.h"

SCFINT *scfint_create(int type)
{     /* create interface to SC-factorization */
      SCFINT *fi;
      fi = talloc(1, SCFINT);
      memset(fi, 0, sizeof(SCFINT));
      switch ((fi->scf.type = type))
      {  case 1:
            fi->u.lufi = lufint_create();
            break;
         case 2:
            fi->u.btfi = btfint_create();
            break;
         default:
            xassert(type != type);
      }
      return fi;
}

int scfint_factorize(SCFINT *fi, int n, int (*col)(void *info, int j,
      int ind[], double val[]), void *info)
{     /* compute SC-factorization of specified matrix A */
      int nn_max, old_n0_max, n0_max, k, ret;
      xassert(n > 0);
      fi->valid = 0;
      /* get required value of nn_max */
      nn_max = fi->nn_max;
      if (nn_max == 0)
         nn_max = 100;
      xassert(nn_max > 0);
      /* compute factorization of specified matrix A */
      switch (fi->scf.type)
      {  case 1:
            old_n0_max = fi->u.lufi->n_max;
            fi->u.lufi->sva_n_max = 4 * n + 2 * nn_max;
            ret = lufint_factorize(fi->u.lufi, n, col, info);
            n0_max = fi->u.lufi->n_max;
            fi->scf.sva = fi->u.lufi->sva;
            fi->scf.a0.luf = fi->u.lufi->luf;
            break;
         case 2:
            old_n0_max = fi->u.btfi->n_max;
            fi->u.btfi->sva_n_max = 6 * n + 2 * nn_max;
            ret = btfint_factorize(fi->u.btfi, n, col, info);
            n0_max = fi->u.btfi->n_max;
            fi->scf.sva = fi->u.btfi->sva;
            fi->scf.a0.btf = fi->u.btfi->btf;
            break;
         default:
            xassert(fi != fi);
      }
      /* allocate/reallocate arrays, if necessary */
      if (old_n0_max < n0_max)
      {  if (fi->w1 != NULL)
            tfree(fi->w1);
         if (fi->w2 != NULL)
            tfree(fi->w2);
         if (fi->w3 != NULL)
            tfree(fi->w3);
         fi->w1 = talloc(1+n0_max, double);
         fi->w2 = talloc(1+n0_max, double);
         fi->w3 = talloc(1+n0_max, double);
      }
      if (fi->scf.nn_max != nn_max)
      {  if (fi->scf.ifu.f != NULL)
            tfree(fi->scf.ifu.f);
         if (fi->scf.ifu.u != NULL)
            tfree(fi->scf.ifu.u);
         fi->scf.ifu.f = talloc(nn_max * nn_max, double);
         fi->scf.ifu.u = talloc(nn_max * nn_max, double);
      }
      if (old_n0_max < n0_max || fi->scf.nn_max != nn_max)
      {  if (fi->scf.pp_ind != NULL)
            tfree(fi->scf.pp_ind);
         if (fi->scf.pp_inv != NULL)
            tfree(fi->scf.pp_inv);
         if (fi->scf.qq_ind != NULL)
            tfree(fi->scf.qq_ind);
         if (fi->scf.qq_inv != NULL)
            tfree(fi->scf.qq_inv);
         if (fi->w4 != NULL)
            tfree(fi->w4);
         if (fi->w5 != NULL)
            tfree(fi->w5);
         fi->scf.pp_ind = talloc(1+n0_max+nn_max, int);
         fi->scf.pp_inv = talloc(1+n0_max+nn_max, int);
         fi->scf.qq_ind = talloc(1+n0_max+nn_max, int);
         fi->scf.qq_inv = talloc(1+n0_max+nn_max, int);
         fi->w4 = talloc(1+n0_max+nn_max, double);
         fi->w5 = talloc(1+n0_max+nn_max, double);
      }
      /* initialize SC-factorization */
      fi->scf.n = n;
      fi->scf.n0 = n;
      fi->scf.nn_max = nn_max;
      fi->scf.nn = 0;
      fi->scf.rr_ref = sva_alloc_vecs(fi->scf.sva, nn_max);
      fi->scf.ss_ref = sva_alloc_vecs(fi->scf.sva, nn_max);
      fi->scf.ifu.n_max = nn_max;
      fi->scf.ifu.n = 0;
      for (k = 1; k <= n; k++)
      {  fi->scf.pp_ind[k] = k;
         fi->scf.pp_inv[k] = k;
         fi->scf.qq_ind[k] = k;
         fi->scf.qq_inv[k] = k;
      }
      /* set validation flag */
      if (ret == 0)
         fi->valid = 1;
      return ret;
}

int scfint_update(SCFINT *fi, int upd, int j, int len, const int ind[],
      const double val[])
{     /* update SC-factorization after replacing j-th column of A */
      int n = fi->scf.n;
      int n0 = fi->scf.n0;
      int nn = fi->scf.nn;
      int *pp_ind = fi->scf.pp_ind;
      int *qq_ind = fi->scf.qq_ind;
      int *qq_inv = fi->scf.qq_inv;
      double *bf = fi->w4;
      double *dg = fi->w5;
      int k, t, ret;
      xassert(fi->valid);
      xassert(0 <= n && n <= n0+nn);
      /* (b, f) := inv(P) * (beta, 0) */
      for (k = 1; k <= n0+nn; k++)
         bf[k] = 0.0;
      for (t = 1; t <= len; t++)
      {  k = ind[t];
         xassert(1 <= k && k <= n);
#if 1 /* FIXME: currently P = I */
         xassert(pp_ind[k] == k);
#endif
         xassert(bf[k] == 0.0);
         xassert(val[t] != 0.0);
         bf[k] = val[t];
      }
      /* (d, g) := Q * (cj, 0) */
      for (k = 1; k <= n0+nn; k++)
         dg[k] = 0.0;
      xassert(1 <= j && j <= n);
      dg[fi->scf.qq_inv[j]] = 1;
      /* update factorization of augmented matrix */
      ret = scf_update_aug(&fi->scf, &bf[0], &dg[0], &bf[n0], &dg[n0],
         0.0, upd, fi->w1, fi->w2, fi->w3);
      if (ret == 0)
      {  /* swap j-th and last columns of new matrix Q */
         scf_swap_q_cols(j, n0+nn+1);
      }
      else
      {  /* updating failed */
         fi->valid = 0;
      }
      return ret;
}

void scfint_ftran(SCFINT *fi, double x[])
{     /* solve system A * x = b */
      xassert(fi->valid);
      scf_a_solve(&fi->scf, x, fi->w4, fi->w5, fi->w1, fi->w2);
      return;
}

void scfint_btran(SCFINT *fi, double x[])
{     /* solve system A'* x = b */
      xassert(fi->valid);
      scf_at_solve(&fi->scf, x, fi->w4, fi->w5, fi->w1, fi->w2);
      return;
}

double scfint_estimate(SCFINT *fi)
{     /* estimate 1-norm of inv(A) */
      double norm;
      xassert(fi->valid);
      xassert(fi->scf.n == fi->scf.n0);
      switch (fi->scf.type)
      {  case 1:
            norm = luf_estimate_norm(fi->scf.a0.luf, fi->w1, fi->w2);
            break;
         case 2:
            norm = btf_estimate_norm(fi->scf.a0.btf, fi->w1, fi->w2,
               fi->w3, fi->w4);
            break;
         default:
            xassert(fi != fi);
      }
      return norm;
}

void scfint_delete(SCFINT *fi)
{     /* delete interface to SC-factorization */
      switch (fi->scf.type)
      {  case 1:
            lufint_delete(fi->u.lufi);
            break;
         case 2:
            btfint_delete(fi->u.btfi);
            break;
         default:
            xassert(fi != fi);
      }
      if (fi->scf.ifu.f != NULL)
         tfree(fi->scf.ifu.f);
      if (fi->scf.ifu.u != NULL)
         tfree(fi->scf.ifu.u);
      if (fi->scf.pp_ind != NULL)
         tfree(fi->scf.pp_ind);
      if (fi->scf.pp_inv != NULL)
         tfree(fi->scf.pp_inv);
      if (fi->scf.qq_ind != NULL)
         tfree(fi->scf.qq_ind);
      if (fi->scf.qq_inv != NULL)
         tfree(fi->scf.qq_inv);
      if (fi->w1 != NULL)
         tfree(fi->w1);
      if (fi->w2 != NULL)
         tfree(fi->w2);
      if (fi->w3 != NULL)
         tfree(fi->w3);
      if (fi->w4 != NULL)
         tfree(fi->w4);
      if (fi->w5 != NULL)
         tfree(fi->w5);
      tfree(fi);
      return;
}

void scfint_copy(SCFINT *dst, SCFINT *src)
{     /* copy interface to SC-factorization */
      int n0_max, old_n0_max, nn_max, n0, nn, k, n;
      xassert(src->scf.type == dst->scf.type);
      switch (src->scf.type)
      {  case 1:
            old_n0_max = dst->u.lufi->n_max;
            lufint_copy(dst->u.lufi, src->u.lufi);
            n0_max = dst->u.lufi->n_max;
            n0 = dst->u.lufi->luf->n;
            dst->scf.sva = dst->u.lufi->sva;
            dst->scf.a0.luf = dst->u.lufi->luf;
            break;
         case 2:
            old_n0_max = dst->u.btfi->n_max;
            btfint_copy(dst->u.btfi, src->u.btfi);
            n0_max = dst->u.btfi->n_max;
            n0 = dst->u.btfi->btf->n;
            dst->scf.sva = dst->u.btfi->sva;
            dst->scf.a0.btf = dst->u.btfi->btf;
            break;
         default:
            xassert(dst != dst);
      }
      /* allocate/reallocate arrays, if necessary */
      if (old_n0_max < n0_max)
      {  if (dst->w1 != NULL)
            tfree(dst->w1);
         if (dst->w2 != NULL)
            tfree(dst->w2);
         if (dst->w3 != NULL)
            tfree(dst->w3);
         dst->w1 = talloc(1+n0_max, double);
         dst->w2 = talloc(1+n0_max, double);
         dst->w3 = talloc(1+n0_max, double);
      }
      nn_max = dst->scf.nn_max;
      nn = src->scf.nn;
      if (nn > nn_max)
      {  nn_max = src->scf.nn_max;
         if (dst->scf.ifu.f != NULL)
            tfree(dst->scf.ifu.f);
         if (dst->scf.ifu.u != NULL)
            tfree(dst->scf.ifu.u);
         dst->scf.ifu.f = talloc(nn_max * nn_max, double);
         dst->scf.ifu.u = talloc(nn_max * nn_max, double);
      }
      n = src->scf.ifu.n;
      if (src->scf.ifu.f != NULL)
         for (k = 0; k < n; k++)
            memcpy(dst->scf.ifu.f + k, src->scf.ifu.f + k,
               n * sizeof(double));
      if (src->scf.ifu.u != NULL)
         for (k = 0; k < n; k++)
            memcpy(dst->scf.ifu.u + k, src->scf.ifu.u + k,
               n * sizeof(double));
      if (old_n0_max < n0_max || dst->scf.nn_max != nn_max)
      {  if (dst->scf.pp_ind != NULL)
            tfree(dst->scf.pp_ind);
         if (dst->scf.pp_inv != NULL)
            tfree(dst->scf.pp_inv);
         if (dst->scf.qq_ind != NULL)
            tfree(dst->scf.qq_ind);
         if (dst->scf.qq_inv != NULL)
            tfree(dst->scf.qq_inv);
         if (dst->w4 != NULL)
            tfree(dst->w4);
         if (dst->w5 != NULL)
            tfree(dst->w5);
         dst->scf.pp_ind = talloc(1+n0_max+nn_max, int);
         dst->scf.pp_inv = talloc(1+n0_max+nn_max, int);
         dst->scf.qq_ind = talloc(1+n0_max+nn_max, int);
         dst->scf.qq_inv = talloc(1+n0_max+nn_max, int);
         dst->w4 = talloc(1+n0_max+nn_max, double);
         dst->w5 = talloc(1+n0_max+nn_max, double);
      }
      k = src->scf.n0 + src->scf.nn;
      if (src->scf.pp_ind != NULL)
         memcpy(dst->scf.pp_ind, src->scf.pp_ind, (1+k) * sizeof(int));
      if (src->scf.pp_inv != NULL)
         memcpy(dst->scf.pp_inv, src->scf.pp_inv, (1+k) * sizeof(int));
      if (src->scf.qq_ind != NULL)
         memcpy(dst->scf.qq_ind, src->scf.qq_ind, (1+k) * sizeof(int));
      if (src->scf.qq_inv != NULL)
         memcpy(dst->scf.qq_inv, src->scf.qq_inv, (1+k) * sizeof(int));
      /* copy SC-factorization */
      dst->scf.n = src->scf.n;
      dst->scf.n0 = src->scf.n0;
      dst->scf.nn_max = nn_max;
      dst->scf.nn = src->scf.nn;
      dst->scf.rr_ref = src->scf.rr_ref;
      dst->scf.ss_ref = src->scf.ss_ref;
      dst->scf.ifu.n_max = nn_max;
      dst->scf.ifu.n = src->scf.ifu.n;
      dst->nn_max = src->nn_max;
      dst->valid = src->valid;
      return;
}

/* eof */
