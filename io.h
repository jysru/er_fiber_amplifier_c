#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "generic.h"

//int GetNbLines(char* s) ;
//void wtodbm(double* wp, double** ws, double* dbmp, double** dbms) ;
void ecriture_fichier(double* tab1, double* tab2, char* fichier1,char* fichier2) ;

void ecriture_fichier_Z(double* tab1, double* tab2, double* tab3, double* tab4, double* tab5, char* nom_fichier) ;
void ecriture_fichier_lambda(double* tab1, double* tab2, double* tab3, double* tab4, double* tab5, double* tab6, double* tab7, char* nom_fichier) ;
//void ecriture_fichier_Matrix_Ps(double** tab, char* nom_fichier);
void lecture_fichier(double* tab1, double* tab2, double* tab3, char* nom_fichier) ;
//void trace_gnuplot() ;
void Ecriture2D(double **X, char* nom);
