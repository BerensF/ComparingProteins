@echo off

REM This batch file compiles the program with the GNU C/C++ compiler under Windows

set src=../lp_solve_5.5/shared/commonlib.c ../lp_solve_5.5/shared/mmio.c ../lp_solve_5.5/shared/myblas.c ../lp_solve_5.5/ini.c ../lp_solve_5.5/lp_rlp.c ../lp_solve_5.5/lp_crash.c ../lp_solve_5.5/bfp/bfp_LUSOL/lp_LUSOL.c ../lp_solve_5.5/bfp/bfp_LUSOL/LUSOL/lusol.c ../lp_solve_5.5/lp_Hash.c ../lp_solve_5.5/lp_lib.c ../lp_solve_5.5/lp_wlp.c ../lp_solve_5.5/lp_matrix.c ../lp_solve_5.5/lp_mipbb.c ../lp_solve_5.5/lp_MPS.c ../lp_solve_5.5/lp_params.c ../lp_solve_5.5/lp_presolve.c ../lp_solve_5.5/lp_price.c ../lp_solve_5.5/lp_pricePSE.c ../lp_solve_5.5/lp_report.c ../lp_solve_5.5/lp_scale.c ../lp_solve_5.5/lp_simplex.c ../lp_solve_5.5/lp_SOS.c ../lp_solve_5.5/lp_utils.c ../lp_solve_5.5/yacc_read.c ../lp_solve_5.5/lp_MDO.c ../lp_solve_5.5/colamd/colamd.c
set c=g++

%c% -I../lp_solve_5.5 -I../lp_solve_5.5/bfp -I../lp_solve_5.5/bfp/bfp_LUSOL -I../lp_solve_5.5/bfp/bfp_LUSOL/LUSOL -I../lp_solve_5.5/colamd -I../lp_solve_5.5/shared -DYY_NEVER_INTERACTIVE -DPARSER_LP -DINVERSE_ACTIVE=INVERSE_LUSOL -DRoleIsExternalInvEngine main.cpp %src% -w -o main -lpthread
