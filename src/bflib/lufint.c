/* lufint.c (interface to LU-factorization) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Copyright (C) 2012-2013 Andrew Makhorin, Department for Applied
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
#include "lufint.h"

LUFINT *lufint_create(void)
{     /* create interface to LU-factorization */
      LUFINT *fi;
      fi = talloc(1, LUFINT);
      fi->n_max = 0;
      fi->valid = 0;
      fi->sva = NULL;
      fi->luf = NULL;
      fi->sgf = NULL;
      fi->sva_n_max = fi->sva_size = 0;
      fi->delta_n0 = fi->delta_n = 0;
      fi->sgf_updat = 0;
      fi->sgf_piv_tol = 0.10;
      fi->sgf_piv_lim = 4;
      fi->sgf_suhl = 1;
      fi->sgf_eps_tol = DBL_EPSILON;
      return fi;
}

int lufint_factorize(LUFINT *fi, int n, int (*col)(void *info, int j,
      int ind[], double val[]), void *info)
{     /* compute LU-factorization of specified matrix A */
      SVA *sva;
      LUF *luf;
      SGF *sgf;
      int k;
      xassert(n > 0);
      fi->valid = 0;
      /* create sparse vector area (SVA), if necessary */
      sva = fi->sva;
      if (sva == NULL)
      {  int sva_n_max = fi->sva_n_max;
         int sva_size = fi->sva_size;
         if (sva_n_max == 0)
            sva_n_max = 4 * n;
         if (sva_size == 0)
            sva_size = 10 * n;
         sva = fi->sva = sva_create_area(sva_n_max, sva_size);
      }
      /* allocate/reallocate underlying objects, if necessary */
      if (fi->n_max < n)
      {  int n_max = fi->n_max;
         if (n_max == 0)
            n_max = fi->n_max = n + fi->delta_n0;
         else
            n_max = fi->n_max = n + fi->delta_n;
         xassert(n_max >= n);
         /* allocate/reallocate LU-factorization (LUF) */
         luf = fi->luf;
         if (luf == NULL)
         {  luf = fi->luf = talloc(1, LUF);
            memset(luf, 0, sizeof(LUF));
            luf->sva = sva;
         }
         else
         {  tfree(luf->vr_piv);
            tfree(luf->pp_ind);
            tfree(luf->pp_inv);
            tfree(luf->qq_ind);
            tfree(luf->qq_inv);
         }
         luf->vr_piv = talloc(1+n_max, double);
         luf->pp_ind = talloc(1+n_max, int);
         luf->pp_inv = talloc(1+n_max, int);
         luf->qq_ind = talloc(1+n_max, int);
         luf->qq_inv = talloc(1+n_max, int);
         /* allocate/reallocate factorizer workspace (SGF) */
         sgf = fi->sgf;
         if (sgf == NULL)
         {  sgf = fi->sgf = talloc(1, SGF);
            memset(sgf, 0, sizeof(SGF));
            sgf->luf = luf;
         }
         else
         {  tfree(sgf->rs_head);
            tfree(sgf->rs_prev);
            tfree(sgf->rs_next);
            tfree(sgf->cs_head);
            tfree(sgf->cs_prev);
            tfree(sgf->cs_next);
            tfree(sgf->vr_max);
            tfree(sgf->flag);
            tfree(sgf->work);
         }
         sgf->rs_head = talloc(1+n_max, int);
         sgf->rs_prev = talloc(1+n_max, int);
         sgf->rs_next = talloc(1+n_max, int);
         sgf->cs_head = talloc(1+n_max, int);
         sgf->cs_prev = talloc(1+n_max, int);
         sgf->cs_next = talloc(1+n_max, int);
         sgf->vr_max = talloc(1+n_max, double);
         sgf->flag = talloc(1+n_max, char);
         sgf->work = talloc(1+n_max, double);
      }
      luf = fi->luf;
      sgf = fi->sgf;
#if 1 /* FIXME */
      /* initialize SVA */
      sva->n = 0;
      sva->m_ptr = 1;
      sva->r_ptr = sva->size + 1;
      sva->head = sva->tail = 0;
#endif
      /* allocate sparse vectors in SVA */
      luf->n = n;
      luf->fr_ref = sva_alloc_vecs(sva, n);
      luf->fc_ref = sva_alloc_vecs(sva, n);
      luf->vr_ref = sva_alloc_vecs(sva, n);
      luf->vc_ref = sva_alloc_vecs(sva, n);
      /* store matrix V = A in column-wise format */
      luf_store_v_cols(luf, col, info, sgf->rs_prev, sgf->work);
      /* setup factorizer control parameters */
      sgf->updat = fi->sgf_updat;
      sgf->piv_tol = fi->sgf_piv_tol;
      sgf->piv_lim = fi->sgf_piv_lim;
      sgf->suhl = fi->sgf_suhl;
      sgf->eps_tol = fi->sgf_eps_tol;
      /* compute LU-factorization of specified matrix A */
      k = sgf_factorize(sgf, 1);
      if (k == 0)
         fi->valid = 1;
      return k;
}

void lufint_delete(LUFINT *fi)
{     /* delete interface to LU-factorization */
      SVA *sva = fi->sva;
      LUF *luf = fi->luf;
      SGF *sgf = fi->sgf;
      if (sva != NULL)
         sva_delete_area(sva);
      if (luf != NULL)
      {  tfree(luf->vr_piv);
         tfree(luf->pp_ind);
         tfree(luf->pp_inv);
         tfree(luf->qq_ind);
         tfree(luf->qq_inv);
         tfree(luf);
      }
      if (sgf != NULL)
      {  tfree(sgf->rs_head);
         tfree(sgf->rs_prev);
         tfree(sgf->rs_next);
         tfree(sgf->cs_head);
         tfree(sgf->cs_prev);
         tfree(sgf->cs_next);
         tfree(sgf->vr_max);
         tfree(sgf->flag);
         tfree(sgf->work);
         tfree(sgf);
      }
      tfree(fi);
      return;
}

void lufint_copy(LUFINT *dst, LUFINT *src)
{     /* copy interface to LU-factorization */
      int n, n_max;
      n = src->luf->n;
      n_max = dst->n_max;
      if (n > n_max)
      {  n_max = dst->n_max = src->n_max;
         if (dst->luf != NULL)
         {  tfree(dst->luf->vr_piv);
            tfree(dst->luf->pp_ind);
            tfree(dst->luf->pp_inv);
            tfree(dst->luf->qq_ind);
            tfree(dst->luf->qq_inv);
         }
         else
            dst->luf = talloc(1, LUF);
         dst->luf->vr_piv = talloc(1+n_max, double);
         dst->luf->pp_ind = talloc(1+n_max, int);
         dst->luf->pp_inv = talloc(1+n_max, int);
         dst->luf->qq_ind = talloc(1+n_max, int);
         dst->luf->qq_inv = talloc(1+n_max, int);
         if (dst->sgf != NULL)
         {  tfree(dst->sgf->rs_head);
            tfree(dst->sgf->rs_prev);
            tfree(dst->sgf->rs_next);
            tfree(dst->sgf->cs_head);
            tfree(dst->sgf->cs_prev);
            tfree(dst->sgf->cs_next);
            tfree(dst->sgf->vr_max);
            tfree(dst->sgf->flag);
            tfree(dst->sgf->work);
         }
         else
            dst->sgf = talloc(1, SGF);
         dst->sgf->rs_head = talloc(1+n_max, int);
         dst->sgf->rs_prev = talloc(1+n_max, int);
         dst->sgf->rs_next = talloc(1+n_max, int);
         dst->sgf->cs_head = talloc(1+n_max, int);
         dst->sgf->cs_prev = talloc(1+n_max, int);
         dst->sgf->cs_next = talloc(1+n_max, int);
         dst->sgf->vr_max = talloc(1+n_max, double);
         dst->sgf->flag = talloc(1+n_max, char);
         dst->sgf->work = talloc(1+n_max, double);
      }
      dst->valid = src->valid;
      dst->sva_n_max = src->sva_n_max;
      dst->sva_size = src->sva_size;
      dst->delta_n0 = src->delta_n0;
      dst->delta_n = src->delta_n;
      dst->sgf_updat = src->sgf_updat;
      dst->sgf_piv_tol = src->sgf_piv_tol;
      dst->sgf_piv_lim  = src->sgf_piv_lim;
      dst->sgf_suhl = src->sgf_suhl;
      dst->sgf_eps_tol = src->sgf_eps_tol;
      if (dst->sva == NULL)
         dst->sva = sva_create_area(src->sva->n_max, src->sva->size);
      sva_copy_area(dst->sva, src->sva);
      dst->luf->n = n;
      dst->luf->fr_ref = src->luf->fr_ref;
      dst->luf->fc_ref = src->luf->fc_ref;
      dst->luf->vr_ref = src->luf->vr_ref;
      dst->luf->vc_ref = src->luf->vc_ref;
      memcpy(dst->luf->vr_piv, src->luf->vr_piv,
         (1+n) * sizeof(double));
      memcpy(dst->luf->pp_ind, src->luf->pp_ind, (1+n) * sizeof(int));
      memcpy(dst->luf->pp_inv, src->luf->pp_inv, (1+n) * sizeof(int));
      memcpy(dst->luf->qq_ind, src->luf->qq_ind, (1+n) * sizeof(int));
      memcpy(dst->luf->qq_inv, src->luf->qq_inv, (1+n) * sizeof(int));
      dst->luf->sva = dst->sva;
      dst->sgf->luf = dst->luf;
      dst->sgf->updat = src->sgf->updat;
      dst->sgf->piv_tol = src->sgf->piv_tol;
      dst->sgf->piv_lim = src->sgf->piv_lim;
      dst->sgf->suhl = src->sgf->suhl;
      dst->sgf->eps_tol = src->sgf->eps_tol;
      return;
}

/* eof */
