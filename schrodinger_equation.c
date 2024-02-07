#include <stdio.h>
#include <math.h>
#include "complex.h"
#define N 100
#define PI 3.141592

void calculo_potencial (double V[N+1], double k_0, double lambda)
{
    int j;

    for(j=0; j<=N+1; j++)
    {
        V[j]=0.;
        if(j>= (2.*N)/5 && j<=(3.*N)/5)
        {
            V[j]=lambda*k_0*k_0;
        }
    }

    return;
}

void inicio_psi (fcomplex psi[N+1], double k_0, double s)
{
    int j;
    double modulo, fase, norm;
    fcomplex aux[N+1];

    norm = 0.;

    aux[0].r = 0;
    aux[N].r = 0;
    aux[0].i = 0;
    aux[N].i = 0;
    for(j=1; j<=N; j++)
    {
        modulo = exp(-8.*(pow((4.*j-N), 2))/(N*N));
        fase = k_0*j;

        aux[j]= Cgauss(fase, modulo);

        norm += Cabs(aux[j])*Cabs(aux[j]);
    }


    for(j=1; j<=N; j++)
    {
        psi[j] = RCmul((1./sqrt(norm)), aux[j]);
        //psi[j] = aux[j];

    }

    return;
}

void calculo_auxb (double s, fcomplex psi[N+1], fcomplex b[N+1])
{
    int j;
    fcomplex im;

    im.r= 0.;
    im.i = 1.;
    for(j=0; j<N+1; j++)
    {
        b[j] = Cmul(RCmul(4./s, psi[j]), im);
    }

    return;
}

void calculo_alpha (double s, double V[N+1], fcomplex alpha[N+1])
{
    int j;
    fcomplex gamma[N+1], num1;
 
    num1.r=-1.;
    num1.i=0.;

    for(j=N-1; j>0; j--)
    {
        gamma[j].r = -2. -V[j] + alpha[j].r;
        gamma[j].i = 2./s + alpha[j].i;

        alpha[j-1] = Cdiv(num1, gamma[j]);      
    }

    return;
}

void calculo_beta (double s, double V[N+1], fcomplex b[N+1], fcomplex alpha[N+1], fcomplex beta[N+1])
{
    int j;
    fcomplex gamma[N+1], resta, num1, aux;

    
    num1.r = 1.;
    num1.i = 0.;

    for(j=N-1; j>0; j--)
    {
        gamma[j].r = -2. -V[j] +alpha[j].r;
        gamma[j].i = 2./s + alpha[j].i;

        resta = Csub(b[j], beta[j]);
        aux = Cdiv(num1, gamma[j]);

        beta[j-1] = Cmul(resta, aux); 
    }

    return;
}

void calculo_chi (fcomplex alpha[N+1], fcomplex beta[N+1], fcomplex chi[N+1])
{
    int j;
    fcomplex aux;
    double re, im, auxre, auxim;

    for(j=0; j<N; j++)
    {
        aux = Cmul(alpha[j], chi[j]);

        chi[j+1] = Cadd(aux, beta[j]);
    }

    return;
}

void calculo_psi (fcomplex chi[N+1], fcomplex psi[N+1])
{
    int j;

    for(j=0; j<=N; j++)
    {
        psi[j] = Csub(chi[j], psi[j]);
    }

    return;
}

double calculo_norma(fcomplex psi[N+1])
{
    int j;
    double norm = 0.;

    for(j=0; j<=N; j++)
    {
        norm += Cabs(psi[j])*Cabs(psi[j]);
    }

    return norm;
}

int main ()
{
    int n;
    int i;
    int nciclos = N/4;
    double lambda = 0.3;
    double k_0, s, V[N+1], mod[N+1], muestro;
    fcomplex alpha[N+2], beta[N+1], chi[N+1], psi[N+1], b[N+1];

    FILE *f1, *f2;

    f1 = fopen("DatosSch.txt", "w");
    f2 = fopen("Norma.txt", "w");

    //Obtengo los parámetros del pasos
    k_0=(2.*PI*nciclos)/N;
    s=1/(4*k_0*k_0);

    //Inicio los valores de los vectores
    psi[0].r = 0.;
    psi[N+1].r = 0.;
    psi[0].i = 0.;
    psi[N+1].i = 0.;

    alpha[N-1].r = 0.;
    alpha[N-1].i = 0.;

    beta[N-1].r=0.;
    beta[N-1].i=0.;

    chi[N+1].r=0.;
    chi[N+1].i=0.;    
    chi[0].r=0.;
    chi[0].i=0.;   

    calculo_potencial(V, k_0, lambda);
    inicio_psi(psi, k_0, s);
    calculo_alpha(s, V, alpha);

    for(n=0; n<500; n++)
    {

        muestro = calculo_norma(psi);
        fprintf(f2, "La norma del paso %i es: %lf\n", n, muestro);
        for(i=0; i<=N; i++)
        {
            mod[i] = Cabs(psi[i]);

            fprintf(f1, "%i, %e, %lf,\n", i, mod[i]*mod[i], V[i]);
        }

        fprintf(f1, "\n");

        calculo_auxb(s, psi, b);
        calculo_beta(s, V, b, alpha, beta);
        calculo_chi(alpha, beta, chi);
        calculo_psi(chi, psi);       
    }

    fclose(f1);
    fclose(f2);

    return 0;
}