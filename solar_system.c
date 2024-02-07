#include <stdio.h>
#include <math.h>
#define c 1.496e11
#define Ms 1.989e30
#define N 9
#define G 6.67e-11


double modulo (double v1, double v2)
{
    return sqrt(v1*v1+v2*v2);
}

void aceleracion (double m[N], double x[N], double y[N], double ax[N], double ay[N])
{
    int i, j;
    double w1[N], w2[N];

    for(i=0; i<N; i++)
    {
        //Inicializo las variables a 0
        ax[i]=0.;
        ay[i]=0.;
        w1[i]=0.;
        w2[i]=0.;
        for(j=0; j<N; j++)
        {
            //Obtengo las aceleraciones teneindo en cuenta las interacciones de cada planeta con el resto
            if(i != j)
            {
                w1[j]=x[i]-x[j];
                w2[j]=y[i]-y[j];
                ax[i] = ax[i] - m[j]*(w1[j])/(pow(modulo(w1[j],w2[j]) ,3));
                ay[i] = ay[i] - m[j]*(w2[j])/(pow(modulo(w1[j],w2[j]) ,3));
                //printf("ax: %lf ay: %lf\n", ax[j], ay[j]);
            }
        }
    }

    return;
}


void verlet (double m[N], double x[N], double y[N], double vx[N], double vy[N], double ax[N], double ay[N], double h)
{
    FILE *f4, *f7, *f9;
    double wx[N], wy[N], Ec, Ep, w1[N], w2[N], Em, aux[N], periodo[N];
    int T=1040;
    double t=0.;
    int cont =0;
    double vuelta[N];
    int i, j, l, k;

    //Inicializo los vectores auxiliares a 0
    for(k=0; k<N; k++)
    {
        wx[k]=0.;
        wy[k]=0.;
        w1[k]=0.;
        w2[k]=0.;
        aux[k]=0.;
        vuelta[k]=0.;
    }
    
    //Obtengo la aceleración inicial
    aceleracion(m, x, y, ax, ay);

    f4=fopen("Planetas1.txt", "w");
    f7=fopen("Energia.txt", "w");
    f9=fopen("Periodos.txt", "w");
    while (t<T)
    {

        for(i=0; i<N; i++)
        {
            aux[i]=y[i];
        }

        //Calculo las nuevas posiciones y las copio en un fichero
        for(i=0; i<N; i++)
        {            
            fprintf(f4, "%lf , \t %lf\n", x[i], y[i]);

            if(i>0)
            {
                x[i] += vx[i]*h + ((h*h)/2.)*ax[i];
                y[i] += vy[i]*h + ((h*h)/2.)*ay[i]; 
            
                wx[i] = vx[i] + (h/2.)*ax[i];
                wy[i] = vy[i] + (h/2.)*ay[i];   
            }                                  
        }
        fprintf(f4, "\n");
        
        //Inicializo en cada iteración las energías
        Ec=0.;
        Ep=0.;
        Em=0.;

        //Obtenció de la energia cinetica
        for(l=0; l<N; l++)
        {
            Ec += (1./2.)*m[l]*pow(modulo(vx[l], vy[l]), 2);
        }
        //printf("Energía cinética: %lf\n", Ec);
        
        //Obtención de la energia potencial
        for(i=0; i<N; i++)
        {
            w1[i]=0.;
            w2[i]=0.;
            for(j=0; j<N; j++)
            {
                if(i!=j)
                {
                    w1[j]=x[i]-x[j];
                    w2[j]=y[i]-y[j];
                    Ep -= m[i]*m[j]/modulo(w1[j], w2[j]);
                }           
            }   
        }
        //printf("Energía potencial: %lf\n", Ep);

        //Obtengo la energía mecánica y la copio en un  fichero
        Em=Ec+Ep;
        fprintf(f7, "%lf\n", Em);

        //Periodo       
        for(i=0; i<N; i++)
        {
            
            if(aux[i]*y[i] < 0)
            {
                vuelta[i]=vuelta[i]+0.5;;
                periodo[i]=(t*58.1)/ vuelta[i];
                printf("El periodo del planteta %i es: %lf\n", i, periodo[i]);
                fprintf(f9, "\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t \n", periodo[1], periodo[2], periodo[3], periodo[4], periodo[5], periodo[6], periodo[7], periodo[8]);
            }

        }

        //Obtengo la nueva aceleración y velocidades
        aceleracion(m, x, y, ax, ay);

        for(j=1; j<N; j++)
        {
            vx[j] = wx[j] + (h/2.)*ax[j];
            vy[j] = wy[j] + (h/2.)*ay[j];
        }

        //Aumento el tiempo y cuento una iteración
        t += + h;
        cont ++;
    }

    printf("Ha habido %i iteraciones\n", cont);
    fclose(f4);
    fclose(f7);
    fclose(f9);

    return;
}

int main()
{
    double x[N], y[N], vx[N], vy[N], a, b, e, d, m[N], M, ax[N], ay[N], h, ener[N], mec, mac;

    //Tamaño del paso de los días
    h=0.1;

    FILE *f1, *f2, *f3;

    f1=fopen("Distancia.txt", "r");
    f2=fopen("Masa.txt", "r");
    f3=fopen("Velocidades.txt", "r");

    //Copio de ficheros los datos iniciales y los reescalo
    int i=0;
    while (fscanf(f1, "%lf\t%lf\n", &a, &b) !=EOF)
    {
        x[i] = a*10e8/c;
        y[i] = b*10e8/c;

        //printf("x: %lf, y: %lf\n", x[i], y[i]);
        i++;
    }
    fclose(f1);
    
    int j=0;
    while (fscanf(f2, "%lf\n", &M) !=EOF)
    {
        m[j]=M*10e24/Ms;

        //printf("M: %lf\n", m[j]);
        j++;
    }
    fclose(f2);

    int k=0;
    while (fscanf(f3, "%lf\t%lf\n", &e, &d) !=EOF)
    {
        vx[k] = e*10e2*sqrt(c/(G*Ms));
        vy[k] = d*10e2*sqrt(c/(G*Ms));

        //printf("vx: %lf, vy: %lf\n", vx[k], vy[k]);
        k++;
    }
    fclose(f3);

    verlet(m, x, y, vx, vy, ax, ay, h);


    return 0;
}