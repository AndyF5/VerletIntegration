/*
 **	Author:			Andrew G Faulkner
 **	Username:		af303
 **	Date:			25/11/13
 **
 **	Name:			3 Body Problem
 **	
 **	Description:	Code to find how three bodies move from given coordinates velocities and mass.
 **
 **
 **
 */

#include <stdlib.h> 
#include <stdio.h> 
#include <math.h> 

#define PROGNAME "verlet" 
#define VERSION  "1.0.0" 
/*#define G 6.67384e-11*/
#define NO_OF_BODIES 3 
#define NO_OF_STEPS 10 
#define TIMESTEP 2e3 
#define G 6.67384e-11 

typedef struct body 
/*Structure to store information about the bodies*/
{ 
    double mass; /*Mass*/
    double x; /*Distance from origin in x and y direction*/
    double y;
	double z;
    double vx; /*speeds in x and y direction*/
    double vy;
	double vz;
} BODY; 

typedef struct vector 
/*Vector structure for returning x/y coordinates from functions*/
{ 
    double x; 
    double y;
	double z;
} VECTOR; 

BODY bodies[NO_OF_BODIES]; /*Array of bodies*/

double data[3*NO_OF_BODIES +1][NO_OF_STEPS]={0}; /*Array for storing information on positions of bodies*/


/*Declaring Functions*/
VECTOR func( int bodyno ); 
int output (FILE *opf, double KE[NO_OF_STEPS], double xCOM[NO_OF_STEPS], double yCOM[NO_OF_STEPS]); 
double KEng (); 


int main () 
{ 
    int i; 
    FILE *opf = fopen ( "data.txt", "w" ); 
    double KE[NO_OF_STEPS]={0}; 
    double xCOM[NO_OF_STEPS], yCOM[NO_OF_STEPS]; 
	
	
    /***********INITAL PARAMETERS************/
    data[0][0]=0.0; 
    bodies[0].mass=2e30; 
    bodies[1].mass=6e24; 
    bodies[2].mass=7.35e22; 
    bodies[0].x = data[1][0]=1e10;   /*Initial x coordinate of body 1*/
    bodies[0].y = data[2][0]=0.0;/*Initial y coordinate of body 1*/
	bodies[0].z = data[3][0]=0.0;
    bodies[1].x = data[4][0]=-1e10;  /*Initial x coordinate of body 2*/
    bodies[1].y = data[5][0]=0.0;   /*Initial y coordinate of body 2*/
	bodies[1].z = data[6][0]=0.0;
    bodies[2].x = data[7][0]=0.0; 
    bodies[2].y = data[8][0]=5e10; 
	bodies[2].z = data[9][0]=0.0;
    bodies[0].vx = 0.0;               /*Initial x velocity, body 1*/
    bodies[0].vy = 10000;             /*Initial y velocity, body 1*/
	bodies[0].vz = 0.0;
    bodies[1].vx = 0.0;               /*Initial x velocity, body 2*/
    bodies[1].vy = -10000;            /*Initial y velocity, body 2*/
	bodies[1].vz = 0.0;
    bodies[2].vx = 50000; 
    bodies[2].vy = 0.0; 
	bodies[2].vz = 0.0;
    /****************************************/
	
    KE[0] = KEng (); 
    double MASS = bodies[0].mass + bodies[1].mass + bodies[2].mass; 
    for ( i = 1; i < NO_OF_STEPS; i++) 
    { 
        VECTOR vnew1, vnew2, vnew3; 
        data[0][i] = 0 + i*TIMESTEP; 
        /*printf("%.10e    %.10e  %.10e   %.10e\n",bodies[2].vx,bodies[2].vy,bodies[1].vx,bodies[1].vy);*/
        vnew1 = func(0); 
        vnew2 = func(1); 
        vnew3 = func(2); 
		
        /*********Calculating new positions*********/
        data[1][i] = data[1][i-1] + vnew1.x*TIMESTEP; 
        data[2][i] = data[2][i-1] + vnew1.y*TIMESTEP;
		data[3][i] = data[3][i-1] + vnew1.z*TIMESTEP;
        data[4][i] = data[4][i-1] + vnew2.x*TIMESTEP; 
        data[5][i] = data[5][i-1] + vnew2.y*TIMESTEP;
		data[6][i] = data[6][i-1] + vnew2.z*TIMESTEP;
        data[7][i] = data[7][i-1] + vnew3.x*TIMESTEP; 
        data[8][i] = data[8][i-1] + vnew3.y*TIMESTEP;
		data[9][i] = data[9][i-1] + vnew3.z*TIMESTEP;
        /*******************************************/
		
        /**Changing postions and velocities of bodies to new values**/
        bodies[0].vx= vnew1.x; 
        bodies[0].vy= vnew1.y;
		bodies[0].vz= vnew1.z;
        bodies[1].vx= vnew2.x; 
        bodies[1].vy= vnew2.y; 
		bodies[1].vz= vnew2.z;
        bodies[2].vx= vnew3.x; 
        bodies[2].vy= vnew3.y;
		bodies[2].vz= vnew3.z;
        bodies[0].x = data[1][i]; 
        bodies[0].y = data[2][i];
		bodies[0].z = data[3][i];
        bodies[1].x = data[4][i]; 
        bodies[1].y = data[5][i];
		bodies[1].z = data[6][i];
        bodies[2].x = data[7][i]; 
        bodies[2].y = data[8][i];
		bodies[2].z = data[9][i];
        /************************************************************/
		
        /*****Conservation of energy and angular momentum checks*****/
        KE[i] = KEng(); 
        int f; 
        double COMx=0,COMy=0; 
        for (f=0; f<NO_OF_BODIES; f++) 
        { 
            COMx = COMx + bodies[i].mass*bodies[i].x; 
            COMy = COMy + bodies[i].mass*bodies[i].y; 
			
        } 
        /*printf("%.10e   %.10e   %.10e\n", COMx, COMy, MASS);*/
        xCOM[i] = COMx/MASS; 
        yCOM[i] = COMy/MASS; 
        /************************************************************/
    } 
    output(opf, KE, xCOM, yCOM); 
    /*print out parameters at finish so more range can be extended without crashing*/
    printf("%.10e %.10e %.10e %.10e %.10e %.10e\n %.10e %.10e %.10e %.10e %.10e %.10e",bodies[0].x,bodies[0].y,bodies[1].x,bodies[1].y,bodies[2].x,bodies[2].y,bodies[0].vx,bodies[0].vy,bodies[1].vx,bodies[1].vy,bodies[2].vx,bodies[2].vy); 
    return 0; 
} 

int output (FILE *opf, double KE[NO_OF_STEPS], double xCOM[NO_OF_STEPS], double yCOM[NO_OF_STEPS]) 
/*Outputs to display and file*/
{ 
    int i; 
    /*fprintf ( opf, "%s v%s Verlet Integration of 3-Body Problem. \n", PROGNAME, VERSION);*/
    /********Printing body positions data to file for plotting and analysis******/
    for ( i=0; i< NO_OF_STEPS; i++) 
    { 
        /*printf("%.10e   %.10e   %.10e   %.10e   %.10e\n", data[0][i], data[1][i], data[2][i], data[3][i], data[4][i]);*/
        fprintf( opf, "%.10e   %.10e   %.10e   %.10e   %.10e    %.10e   %.10e	%.10e	%.10e	%.10e\n", data[0][i], data[1][i], data[2][i], data[3][i], data[4][i], data[5][i], data[6][i], data[7][i], data[8][i], data[9][i]); 
		printf( "%.10e  \n %.10e   %.10e   %.10e   \n%.10e    %.10e   %.10e	\n%.10e	%.10e	%.10e\n", data[0][i], data[1][i], data[2][i], data[3][i], data[4][i], data[5][i], data[6][i], data[7][i], data[8][i], data[9][i]);
    } 
    fclose(opf); 
    /****************************************************************************/
	
    /******Printing checks to new file for analysis******/
    FILE *check = fopen ( "KE.txt", "w" ); 
    for (i=0; i< NO_OF_STEPS; i++) 
    { 
        fprintf(check, "%.10e %.10e   %.10e\n", KE[i], xCOM[i], yCOM[i] ); 
    } 
    /****************************************************/
    return 0; 
} 

VECTOR func ( int bodyno ) 
{ 
    double accmag[NO_OF_BODIES-1], xdisp, ydisp, zdisp, distance; 
    VECTOR accdirec[NO_OF_BODIES-1], vector; 
    int f=0, i=0; 
	
    /******Calculating accleration magnitude and direction*******/
    for (i=0;i<NO_OF_BODIES;i++) 
    { 
        if(i != bodyno ) 
        { 
            xdisp = bodies[i].x - bodies[bodyno].x; 
            ydisp = bodies[i].y - bodies[bodyno].y;
			zdisp = bodies[i].z - bodies[bodyno].z;
            distance = sqrt(xdisp*xdisp + ydisp*ydisp + zdisp*zdisp); 
            accmag[f] = ( G * bodies[i].mass )/(pow(distance+1e1,2)); 
            accdirec[f].x = xdisp/distance; 
            accdirec[f].y = ydisp/distance;
			accdirec[f].z = zdisp/distance;
            f++;
			
			
        } 
		
    } 
    /*************************************************************/
	printf("acceleration %d %.10e %.10e %.10e\n", bodyno, accmag[0]+accmag[1],accdirec[0].x+accdirec[1].x,accdirec[0].y+accdirec[1].y);
    double vnewx = bodies[bodyno].vx + accmag[0]*accdirec[0].x*TIMESTEP/2 + accmag[1]*accdirec[1].x*TIMESTEP/2; 
    double vnewy = bodies[bodyno].vy + accmag[0]*accdirec[0].y*TIMESTEP/2 + accmag[1]*accdirec[1].y*TIMESTEP/2; 
	double vnewz = bodies[bodyno].vz + accmag[0]*accdirec[0].z*TIMESTEP/2 + accmag[1]*accdirec[1].z*TIMESTEP/2;
    /*printf("%.10e   %.10e   %.10e\n", accmag[0], vnewx, vnewy);*/
    vector.x =vnewx; 
    vector.y =vnewy;
	vector.z =vnewz;
    return vector; 
} 

double KEng() 
/*Calculating kinetic energy*/
{ 
    int i; 
    double speed; 
    double KE; 
    for ( i = 0; i < NO_OF_BODIES; i++) 
    { 
        speed = sqrt ( bodies[i].vx*bodies[i].vx + bodies[i].vy*bodies[i].vy ); 
        KE = KE + 0.5*bodies[i].mass*speed*speed; 
    } 
    return KE; 
} 