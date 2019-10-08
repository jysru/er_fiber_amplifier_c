#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "generic.h"

/*
int GetNbLines(char* s)
{
    int n=0 ;
    char c ;
    FILE* fp=NULL ;
    fp = fopen(s,"r") ;
    if(fp==NULL)
    {
        printf("Erreur d'ouverture du fichier %s !\n",s);
        exit(-1) ;
    }
    while(c != EOF)
    {
        c=fgetc(fp) ;
        if(c=='\n')
        {
            n++ ;
        }
    }
    fclose(fp) ;
    return n ;
}
*/


/*void wtodbm(double* wp, double** ws, double* dbmp, double** dbms)
{
    int i=0, j=0 ;
    for(i=0 ; i<Lz ; i++)
    {
        dbmp[i] = 10*log10(wp[i]*1000) ;
        for(j=0 ; j<Ll ; j++)
        {
            dbms[i][j]= 10*log10(ws[i][j]*1000) ;
        }
    }
}*/


void lecture_fichier(double* tab1, double* tab2, double* tab3, char* nom_fichier)
{
    int i=0 ;
    FILE* fp=NULL ;
    fp = fopen(nom_fichier,"r") ;
    for(i=0 ; i<Ll ; i++)
    {
        fscanf(fp, "%lf\t%lf\t%lf\n", &(tab1[i]), &(tab2[i]), &(tab3[i])) ;
        tab1[i]=tab1[i]*1e-9 ;
    }
    fclose(fp) ;
}

void ecriture_fichier(double* tab1, double* tab2, char* fichier1,char* fichier2)
{
    int i,j ;
    FILE* fp=NULL ;
    FILE* fp2=NULL ;
    fp = fopen(fichier1, "w");
    fp2 = fopen(fichier2, "w");

    if(fp==NULL||fp2==NULL)
    {
        printf("Erreur : fichier inexistant \n") ;
        exit(-1) ;
    }

    for(i=0 ; i < Lz ; i=i+1)
    {
        fprintf(fp, "%7.6e\n", tab1[i]) ;
    }
    for(j=0 ; j < Ll ; j=j+1)
    {
        fprintf(fp2, "%7.6e\n", (tab2[j]*1e6)) ;
    }

    fclose(fp) ;
    fclose(fp2) ;
}

void ecriture_fichier_Z(double* tab1, double* tab2, double* tab3, double* tab4, double* tab5, char* nom_fichier)
{
    int i=0 ;
    FILE* fp=NULL ;
    fp = fopen(nom_fichier, "w");

    if(fp==NULL)
    {
        printf("Erreur : fichier inexistant \n") ;
        exit(-1) ;
    }

    fprintf(fp, "z,m\t\t Pp,W\t Pp,dBm\t N1,m-3\t N2,m-3\n") ;
    for(i=0 ; i < Lz ; i++)
    {
        fprintf(fp, "%3.2e\t %3.2e\t %3.2e\t %3.2e\t %3.2e\n", tab1[i], tab2[i], tab3[i], tab4[i], tab5[i]) ;
    }

    fclose(fp) ;
}


void ecriture_fichier_lambda(double* tab1, double* tab2, double* tab3, double* tab4, double* tab5, double* tab6, double* tab7, char* nom_fichier)
{
    int j=0 ;
    FILE* fp=NULL ;
    fp = fopen(nom_fichier, "w");

    if(fp==NULL)
    {
        printf("Erreur : fichier inexistant \n") ;
        exit(-1) ;
    }

    fprintf(fp, "l,nm\t\t s_A,m-2\t s_E,m-2\t V_S\t W0_S\t AeffS\t\t g_S\n") ;
    for(j=0 ; j < Ll ; j++)
    {
        tab1[j]=tab1[j]*1e9 ;
        fprintf(fp, "%3.2e\t %3.2e\t %3.2e\t %3.2e\t %3.2e\t %3.2e\t %3.2e\n", tab1[j], tab2[j], tab3[j], tab4[j], tab5[j], tab6[j], tab7[j]) ;
    }

    fclose(fp) ;
}

void Ecriture2D(double X[][Ll], char* nom)
{
    FILE* fp=fopen(nom, "w+") ;
    int i,j;
    if(fp==NULL)
    {
        printf("Erreur : fichier inexistant \n") ;
        exit(-1) ;
    }
    for(i=0 ; i<Lz ; i=i+1)
    {
        for(j=0 ; j<Ll ; j=j+1)
        {
            fprintf(fp, "%3.2e\t", X[i][j]) ;
        }
        fprintf(fp, "\n") ;
    }
    fclose(fp) ;
}

/*void ecriture_fichier_Matrix_Ps(double** tab, char* nom_fichier)
{
    int i=0, j=0 ;
    FILE* fp=NULL ;
    fp = fopen(nom_fichier, "w");

    if(fp==NULL)
    {
        printf("Erreur : fichier inexistant \n") ;
        exit(-1) ;
    }

    //fprintf(fp, "l,nm\t\t s_A,m-2\t s_E,m-2\t V_S\t W0_S\t AeffS\t\t g_S\n") ;
    for(i=0 ; i<Lz ; i++)
    {
        for(j=0 ; j<Ll ; j++)
        {
            fprintf(fp, "%3.2e\t", tab[i][j]) ;
        }
        fprintf(fp, "\n") ;
    }

    fclose(fp) ;
}*/



/*
void trace_gnuplot()
{
    FILE* gp=NULL;
    gp=popen("gnuplot -persist","w");
    if (gp==NULL)
    {
        printf("Erreur pour le trace sur Gnuplot ! \n");
        exit(-1);
    }

    fprintf(gp, "load 'plots.plt'\n");
    fclose(gp);
}
*/
