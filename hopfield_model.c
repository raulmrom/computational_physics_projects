#include <stdio.h>
#include <math.h>
#include "gsl_rng.h"
#define N 20
#define n 3


gsl_rng *tau;

void calculo_amu (int Xi[N+2][N+2][n], double amu[n])
{
    int i, j, k;

    for(k=0; k<n; k++)
    {
        amu[k] = 0;
        for(i=1; i<N+1; i++)
        {
            for(j=1; j<N+1; j++)
            {
                amu[k] += Xi[i][j][k];
            }
        }
        amu[k] = amu[k]/(N*N);
    }
    
    return;  
}

double vdoble (int i, int j, int k, int l, double amu[n], int Xi[N+2][N+2][n])
{
    int m;
    double aux =0;
  
    if(i == k && j == l)
    {
        return 0.;
    }
    else
    {
        for(m=0; m<n; m++)
        {
            aux += 1.*(Xi[i][j][m]-amu[m])*(Xi[k][l][m]-amu[m]);
        }
        return aux/(1.*N*N);
    }
}

double func_theta (int i, int j, int Xi[N+2][N+2][n], double amu[n])
{
    int k, l;
    double sol, w;

    sol =0.;

    for(k=1; k<N+1; k++)
    {
        for(l=1; l<N+1; l++)
        {
            w = vdoble(i, j, k, l, amu, Xi);
            sol +=  w;
        }
    }

    return 0.5*sol;
}

double deltaE (int Xi[N+2][N+2][n], int s[N+2][N+2], double amu[n], int x, int y, double T)
{
    int i, j, k, l;
    double E, sol, sol2, w, rec;
    E=0.;
    sol=0.;
    sol2=0.;


    for(i=1; i<N+1; i++)
    {
        for(j=1; j<N+1; j++)
        {
            if(i!= x && j != y)
            {
                w = vdoble(i, j, x, y, amu, Xi);
                sol2 += w*s[i][j];
            }
        }
    }    

    E = pow(-1, s[x][y])*(-sol2 + func_theta(x, y, Xi, amu));

    return exp(-E/T);
}

void cambioSigno (int s[N+2][N+2], int x, int y)
{
    if(s[x][y] ==0)
    {
        s[x][y] = 1;
    }
    else
    {
        s[x][y] = 0;
    }

    return;
}

void func_solapamiento (double amu[n], int Xi[N+2][N+2][n], int s[N+2][N+2], double sol[n])
{
    int i, j, k;

    for(k=0; k<n; k++)
    {
        sol[k]=0;
        for(i=1; i<N+1; i++)
        {
            for(j=1; j<N+1; j++)
            {
                sol[k] += (Xi[i][j][k] - amu[k])*(s[i][j] - amu[k]);
            }
        }
        sol[k] = sol[k]/(N*N*amu[k]*(1-amu[k]));
    }
    
    return;
}

int main()
{
    int i, j, k, l, s[N+2][N+2], Xi[N+2][N+2][n], x, y, paso;
    double p, min, xi, variable_a, solapamiento, amu[n], sol[n], rec;
    extern gsl_rng *tau;
    int semilla = 39249354;

    FILE *f1, *f2, *f3, *f4, *f5;

    f1=fopen("cruz.txt", "r");
    f2=fopen("luna.txt", "r");
    f3=fopen("flecha.txt", "r");
    f4=fopen("Datos_Hopefield_patrones.txt", "w");
    f5=fopen("Solapamiento.txt", "w");
    

    tau=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau, semilla);

    //Fijo temperatura
    double T=1e-4;
    //printf("La temperatura generada es: %lf\n", T);

    //Leo las imágenes que se desean recordar
    for(i=1; i<N+1; i++)
    {
        for(j=1; j<N+1; j++)
        {
            fscanf(f1, "%i ", &Xi[i][j][0]);
        }
    }

    for(i=1; i<N+1; i++)
    {
        for(j=1; j<N+1; j++)
        {
            fscanf(f2, "%i ", &Xi[i][j][1]);
        }
    }

    for(i=1; i<N+1; i++)
    {
        for(j=1; j<N+1; j++)
        {
            fscanf(f3, "%i ", &Xi[i][j][2]);
        }
    }

    //Genero patrones aleatorios que se desean recordar
    /*for(i=1; i<N+1; i++)
    {
        for(j=1; j<N+1; j++)
        {
            for(k=0; k<n; k++)
            {
                Xi[i][j][k] = gsl_rng_uniform_int(tau, 2);
            }
        }
    }*/

    //Desordeno uno de los patrones un 50% (para N=20)
    for(i=1; i<201; i++)
    {
        x=gsl_rng_uniform_int(tau, N)+1;
        y=gsl_rng_uniform_int(tau, N)+1;

        Xi[x][y][2] = gsl_rng_uniform_int(tau, 2);
    }


    //Inicializo la matriz de trabajo con números aleatorios
    int cont=0;
    for(i=0; i<=N+1; i++)
    {
        for(j=0; j<=N+1; j++)
        {
            //Inicio el spin en 0 y 1 con igual probabilidad (todos)
            s[i][j]= gsl_rng_uniform_int(tau, 2);
        } 
    }

    //Muestro en un fichero la matriz de la que se parte (ya sea deformada o no)
    for(i=1; i<N+1; i++)
    {
        for(j=1; j<N+1; j++)
        {
            fprintf(f4, "%i,", Xi[i][j][1]);
        }
        fprintf(f4, "\n");
    }
    
    calculo_amu(Xi, amu);    
    paso=0;
    fprintf(f4, "\n");
    for(k=0; k<100*N*N; k++)
    {
        //Obtengo un punto aleatorio de spin
        x=gsl_rng_uniform_int(tau, N)+1;
        y=gsl_rng_uniform_int(tau, N)+1;

        //Obtengo la energía
        p=deltaE(Xi, s, amu, x, y, T);

        //Veo si la energía es menor que 1
        if(p >= 1)
        {
            min=1.;
        }
        else 
        {
            min =p;
        }

        //Genero número entero aleatorio entre 0 y 1
        xi=gsl_rng_uniform(tau);

        //Si el número generado es menor que la energía cambio el signo
        if(xi < min)
        {
            cambioSigno(s, x, y);
        }
        
        //Impongo que se muestran los valores en cada paso Montecarlo (N*N)
        if(paso%(N*N) == 0)
        {
            for(i=1; i<N+1; i++)
            {
                for(j=1; j<N+1; j++)
                    {
                        fprintf(f4, "%i,", s[i][j]);
                    }
                fprintf(f4, "\n");
            }   
            fprintf(f4, "\n");

            func_solapamiento(amu, Xi, s, sol);

            fprintf(f5, "%i\t%lf\t%lf\t%lf\n", paso/(N*N), sol[0], sol[1], sol[2]);

        }
        paso++;
    }

    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);
    fclose(f5);

    return 0;
}