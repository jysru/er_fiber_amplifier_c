#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "generic.h"
#include "io.h"

void main(void)
{
    /************* Définition des variables *************/
    int i=0, j=0, k=0 ;
    double Z[Lz]={0}, lambda[Ll]={0}, sigma_E[Ll]={0}, sigma_A[Ll]={0};
    double ad=0, K1=0, K2=0, K3=0, K4=0 ;
    double Ps[Lz][Ll]={{0}}, Pp[Lz]={0}, Pp_dbm[Lz]={0};
    double Pase_p[Lz][Ll]={{1e-9}}, Pase_m[Lz][Ll]={{1e-9}};
    double N1[Lz]={0}, N2[Lz]={0}, R_Pa=0, W_Sa=0, W_Se=0 ;
    double V_P=0, V_S[Ll]={0}, W0_P=0, W0_S[Ll]={0}, AeffP=0, AeffS[Ll]={0}, gamma_P=0, gamma_S[Ll]={0} ;
    double Ps2[Lz][Nc]={{0}} ;

    /************* Conditions initiales *************/
    ad=a*1 ;
    // Calculs à lambda_P
    V_P=2.0*M_PI/lambda_P*a*ON ;
    W0_P=a*(0.65+1.619*pow(V_P,(-1.5))+2.879*pow(V_P,(-6))) ;
    AeffP=M_PI*pow(W0_P,2) ;
    gamma_P=(1-exp(-2*pow(ad,2)/pow(W0_P,2))) ;
    // Remplissage des tableaux lambda, sigma_A et sigma_E via le fichier fourni
    lecture_fichier(lambda, sigma_E, sigma_A, "sigma.txt") ;
    // Calculs à lambda_S
    for(j=0 ; j<Ll ; j++)
    {
        V_S[j]=2.0*M_PI/lambda[j]*a*ON ;
        W0_S[j]=a*(0.65+1.619*pow(V_S[j],(-1.5))+2.879*pow(V_S[j],(-6))) ;
        AeffS[j]=M_PI*pow(W0_S[j],2) ;
        gamma_S[j]=(1-exp(-2*pow(ad,2)/pow(W0_S[j],2))) ;
    }
    // Remplissage du tableau en Z
    for(i=0 ; i<Lz ; i++)
    {
        Z[i]=i*dz ;
    }
    // Initialisations
    N2[0] = N0*(19.0/20.0) ;
    N1[0] = N0 - N2[0] ;
    Pp[0] = 30e-3 ;

    for(i=0 ; i< Lz ; i++)
        for(j=0 ; j<Ll ; j++)
        {
            Pase_p[i][j]=1e-9 ;
            Pase_m[i][j]=1e-9 ;
        }

    for(j=280 ; j<481 ; j=j+4) // Ajout de puissance seulement dans les canaux
    {
        Ps[0][j]=1e-6 ;
    }


    /************* Méthode de Runge-Kutta d'ordre 4 *************/
    for(i=0 ; i<(Lz-1) ; i++)
    {
        K1 = -sigma_AP*N1[i]*gamma_P*Pp[i] ;
        K2 = -sigma_AP*N1[i]*gamma_P*(Pp[i]+dz*K1/2) ;
        K3 = -sigma_AP*N1[i]*gamma_P*(Pp[i]+dz*K2/2) ;
        K4 = -sigma_AP*N1[i]*gamma_P*(Pp[i]+dz*K3) ;
        Pp[i+1] = Pp[i] + (dz/6)*(K1+2*K2+2*K3+K4) ;
        R_Pa = lambda_P*sigma_AP*gamma_P*Pp[i]/(h*c*AeffP) ;

        for(j=0 ; j<Ll ; j++) // Calcul de Pase+ et Pase-
        {
            // Calcul de Pase+
            K1 = (2*h*c/lambda[j]*sigma_E[j]*N2[i] + sigma_E[j]*N2[i] - sigma_A[j]*N1[i])*gamma_S[j]*Pase_p[i][j] ;
            K2 = (2*h*c/lambda[j]*sigma_E[j]*N2[i] + sigma_E[j]*N2[i] - sigma_A[j]*N1[i])*gamma_S[j]*(Pase_p[i][j]+dz*K1/2) ;
            K3 = (2*h*c/lambda[j]*sigma_E[j]*N2[i] + sigma_E[j]*N2[i] - sigma_A[j]*N1[i])*gamma_S[j]*(Pase_p[i][j]+dz*K2/2) ;
            K4 = (2*h*c/lambda[j]*sigma_E[j]*N2[i] + sigma_E[j]*N2[i] - sigma_A[j]*N1[i])*gamma_S[j]*(Pase_p[i][j]+dz*K3) ;
            Pase_p[i+1][j] = Pase_p[i][j] + (dz/6)*(K1+2*K2+2*K3+K4) ;

            /*// Calcul de Pase-
            K1 = -(2*h*c/lambda[j]*sigma_E[j]*N2[i] + sigma_E[j]*N2[i] - sigma_A[j]*N1[i])*gamma_S[j]*Pase_m[i][j] ;
            K2 = -(2*h*c/lambda[j]*sigma_E[j]*N2[i] + sigma_E[j]*N2[i] - sigma_A[j]*N1[i])*gamma_S[j]*(Pase_m[i][j]+dz*K1/2) ;
            K3 = -(2*h*c/lambda[j]*sigma_E[j]*N2[i] + sigma_E[j]*N2[i] - sigma_A[j]*N1[i])*gamma_S[j]*(Pase_m[i][j]+dz*K2/2) ;
            K4 = -(2*h*c/lambda[j]*sigma_E[j]*N2[i] + sigma_E[j]*N2[i] - sigma_A[j]*N1[i])*gamma_S[j]*(Pase_m[i][j]+dz*K3) ;
            Pase_m[i+1][j] = Pase_m[i][j] + (dz/6)*(K1+2*K2+2*K3+K4) ;*/
        }

        for(j=0 ; j<Ll ; j++) // Méthode des rectangles pour WSa dépendant de Pase
        {
            W_Sa = W_Sa + lambda[j]*sigma_A[j]*gamma_S[j]*(Pase_p[i+1][j]+Pase_m[i][j])/(h*c*AeffS[j]) ;
            W_Se = W_Se + lambda[j]*sigma_E[j]*gamma_S[j]*(Pase_p[i+1][j]+Pase_m[i][j])/(h*c*AeffS[j]) ;
        }

        for(j=280 ; j<481 ; j=j+4) // Somme pour WSa au niveau des canaux, et calcul de Ps au niveau des canaux
        {
            K1 = (sigma_E[j]*N2[i] - sigma_A[j]*N1[i])*gamma_S[j]*Ps[i][j] ;
            K2 = (sigma_E[j]*N2[i] - sigma_A[j]*N1[i])*gamma_S[j]*(Ps[i][j]+dz*K1/2) ;
            K3 = (sigma_E[j]*N2[i] - sigma_A[j]*N1[i])*gamma_S[j]*(Ps[i][j]+dz*K2/2) ;
            K4 = (sigma_E[j]*N2[i] - sigma_A[j]*N1[i])*gamma_S[j]*(Ps[i][j]+dz*K3) ;
            Ps[i+1][j] = Ps[i][j] + (dz/6)*(K1+2*K2+2*K3+K4) ;

            W_Sa = W_Sa + lambda[j]*sigma_A[j]*gamma_S[j]*Ps[i+1][j]/(h*c*AeffS[j]) ;
            W_Se = W_Se + lambda[j]*sigma_E[j]*gamma_S[j]*Ps[i+1][j]/(h*c*AeffS[j]) ;
        }

        N2[i+1] = N0*(R_Pa+W_Sa)/(R_Pa+W_Sa+W_Se+(1/tau)) ;
        N1[i+1] = N0-N2[i+1] ;
        W_Sa=0 ; W_Se=0 ;
        R_Pa=0 ;
    }
    /************* Conversions, écriture et affichage *************/
    // Conversion des puissances en dbm
    for(i=0 ; i<Lz ; i++)
    {
        Pp_dbm[i] = 10*log10(Pp[i]*1000) ;
    }
    ecriture_fichier(Z,lambda,"z.txt","lambda.txt");
    ecriture_fichier_Z(Z, Pp, Pp_dbm, N1, N2, "dataZ.dat");
    ecriture_fichier_lambda(lambda, sigma_A, sigma_E, V_S, W0_S, AeffS, gamma_S, "dataL.dat");

    Ecriture2D(Ps, "Ps.txt");
    Ecriture2D(Pase_p, "Pase.txt");
    for(i=0 ; i<Lz ; i++)
    {
        for(j=0 ; j<Ll ; j++)
        {
            Ps[i][j]= 10*log10((Ps[i][j]+Pase_p[i][j])*1000) ;
            //Ps[i][j]= 10*log10((Pase_p[i][j])*1000) ;
        }
    }
    Ecriture2D(Ps, "Psdbm.txt");

}
