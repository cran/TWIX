/*
 *  split-summary.h
 *  TWIX
 *
 *  Created by Sergej Potapov on 03.02.11.
 *  Copyright 2011 __IMBE__. All rights reserved.
 *
 */


SEXP split_summary( SEXP TBASE, SEXP tol );
SEXP split_summary_padj( SEXP SPLITS, SEXP tol );
SEXP split_summary_dev( SEXP SPLITS, SEXP tol );
SEXP padjust( SEXP PVALUE );


