/*Code running MC simulations of treadmilling filaments confined to a 2D planed (unrolled cell)*/
/*Makes use of Verlet lists*/

#include <stdbool.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/* ------- Definitions ------- */

typedef struct {
	double x,y;
} vec2D;

typedef struct {
	int id;											// molecule id
	int head;										// head particle index
	int tail;										// tail particle index
	int n;											// length
	int tnuc;										// dynamic step where it nucleated
	bool exists;									// True if present in the system, False if diluted
} FILAMENT;

typedef struct {
	int id;											// particle id
	int fil;										// filament index it belongs to
	double x, y;									// coordinates [sigma]
	int bplus;										// bound particle index (+ dir)
	int bminus;										// bound particle index (- dir)
	double angle;									// avg angle of bonds wrt x-axis
	int tnuc;										// dynamic step where it nucleated
	bool exists;									// True if present in the system, False if diluted
	int nVlist;										// Number of elements in its Verlet list
	double xV, yV;									// coordinates at time of Verlet list creation [sigma]
	int Vlist[20];									// Verlet list of this particle
} PARTICLE;

/* ------- ----------- ------- */

/* ------- Simulation parameters ------- */

int L;	 											// side length of simulation box (2D) [sigma]
int N;												// maximum # of particles
int P;												// maximum # of protofilaments
int MCSTEPS;		   								// total # of Monte Carlo steps
int RELSTEPS;										// relaxation steps
int DYNSTEP;
int DUMP;		   									// dumping interval [steps]
int SEED;											// random number generator seed
#define FILE_MAX 500								// maximum filename length
double JUMP_MAX;									// maximum jump magnitude along each direction (x & y) for translation moves

double Kbond;										// bond strength (harmonic constant) [kT/sigma^2]
double Kbend;										// angle potential strength (harmonic constant) [kT/angle^2]
double eps, reps;									// attraction strength (LJ depth) [kT]
double sigma;										// length scale
double arange;										// attraction potential range
double vrange;										// Verlet list range
double vthres;										// Threshold of displacement for recomputing Verlet lists
double Kdir;										// preferred direction potential constant
double Rmem, Rfil, Rfac;							// membrane and filament intrinsic radii of curvature plus factor between the two (Rfil/Rmem)
bool fix_Rmem;										// bool to fix Rmem independently of the box size
bool Min, rMin;										// effect of Min/Noc or other large scale localisation system
bool MinGrow;										// toggle to have modulation amplitude increase over time
bool nucleoid;										// toggle to have kon and knuc modulated as in nucleoid occlusion: start at low kon jump to high kon around midcell (gaussian added to baseline)
bool nucleoidMin;									// toggle to have kon and knuc modulated as in nucleoid occlusion + Min system: start at low kon jump to high kon around midcell and nothing else (only gaussian)
bool nucleoidMinGrow;								// toggle to have kon and knuc modulated as in nucleoid occlusion + Min system BUT DYNAMIC: start at low kon jump to high kon around midcell and nothing else (only gaussian)
double kmax;										// max growth rate inside gaussian area
double MinRate;										// modulation amplitude increase rate [fraction of max/tau_dyn]
double fracmin;										// minimum modulation value
double varMin;										// variance of Min effect gaussian
bool varMinDyn;										// variance of Min effect gaussian
bool saturate;										// set maximum monomer surface concentration
double satVal;										// value of the maximum surface concentration if set
double kon;											// polymerisation rate
double koff;										// depolymerisation rate
double knuc;										// nucleation rate
bool sizer;											// choice of hydrolisis as a sizer
int mean_size;										// target mean size for the above
bool timer;											// choice of hydrolisis as a timer (linear)
bool etimer;										// choice of hydrolisis as a timer (exponential)
int mean_time;										// target mean time for the above
double tangle;										// test filament angle with x
double polang;										// angle arc for polymerisation
int lnuctest;										// nucleus size for test run
bool flat_mem;										// toggle for flat membrane
bool tstraight;										// toggle to have test filament of flat membrane be straight
bool bands;											// initialise the system in polar bands
int nbands;											// number of bands to initialise the system with
bool olapINF;										// set overlap energy to infinity [kT]
bool thermal_noise;									// set xyMC to be thermal noise (gaussian distribution)
bool PC19;											// deactivate hydrolysis (roff = 0) 
int PC19time;										// time of hydrolysis deactivation

/* ------- --------------------- ------- */

/* ------- INITIALISE GLOBAL VARIABES AND ARRAYS ------- */

PARTICLE *particles;								// particles array (array of N particle objects)
FILAMENT *filaments;								// monomers array (array of M monomer objects)

double Etotal;
int PROGDUMP;
int step;											// steps counter
int Nparts;											// particle counter
int maxid;											// current maximum id
int Nfils;											// filament counter
int NpolTrials, NpolTrialsBand;						// polymerisation trials counters
int NpolRejectedEV, NpolRejectedEVBand;				// polymerisation rejections counters

char output_dir[FILE_MAX];							// Output directory
struct stat st = {0};
FILE *energyf;										// energy dump file (.txt)
char energy_file[FILE_MAX];
FILE *dumpf;										// movie dump file (.xyz)
char movie_file[FILE_MAX];
FILE *infof;										// info file (.txt)
char info_file[FILE_MAX];
FILE *ratesf;										// rates file (.txt)
char rates_file[FILE_MAX];
FILE *rejectf;										// rejects file (.txt)
char reject_file[FILE_MAX];

bool verbose;										// print system info
bool long_verbose;									// print more system info
bool test;											// run test on a system 3-monomer system
bool SIF;
bool progress;										// print progress
bool relax;											// bool for relaxation only
bool curvcos;										// bool for using c_sensed = c_0*cos(theta) for curv_energy
bool ICempty;										// bool for empty initial conditions

int xytrials, xyrejects;

/* ------- ------------------------------------- ------- */

/* ------- DOT PRODUCT ------- */
double dot(vec2D a, vec2D b) {
	double result = 0;
	result = a.x*b.x+a.y*b.y;
	return result;
}
/* ------- ----------- ------- */

/* ------- NORMAL DISTRIBUTION WITH VARIABLE STD DEV AND MEAN ------- */
double normal(double loc, double scale) {
	double u1, u2, g1;								// Box-Muller transform
	u1 = drand48();
	u2 = drand48();
	g1 = sqrt(-2*log(u1))*cos(2*M_PI*u2)*scale+loc;
	return g1;
}
/* ------- -------------------------------------------------- ------- */

/* ------- LENNARD-JONES POTENTIAL ------- */
double LJ(double epsilon, int sigma, double range, double r) {
	double result;
	if (r<range) {result=4*epsilon*(pow((sigma/r),12)-pow((sigma/r),6)-pow((sigma/range),12)+pow((sigma/range),6));}
	else {result=0;}
	return result;
}
/* ------- ----------------------- ------- */

/* ------- CHECK PBCs ------- */
void pbcs(int k){
	if (particles[k].x>=L/2.0) particles[k].x-=L; // printf("PBCs (right) implemented on particle %d\n", k);
	if (particles[k].y>=L/2.0) particles[k].y-=L; // printf("PBCs (top) implemented on particle %d\n", k);
	if (particles[k].x<-L/2.0) particles[k].x+=L; // printf("PBCs (left) implemented on particle %d\n", k);
	if (particles[k].y<-L/2.0) particles[k].y+=L; // printf("PBCs (bottom) implemented on particle %d\n", k);
	return;
}
/* ------- ---------- ------- */

/* ------- CHECK PBCs ------- */
vec2D vec_pbcs(vec2D a){
	if (a.x>=L/2.0) a.x-=L; //printf("PBCs (right) implemented on particle %d\n", k);
	if (a.y>=L/2.0) a.y-=L; //printf("PBCs (top) implemented on particle %d\n", k);
	if (a.x<-L/2.0) a.x+=L; //printf("PBCs (left) implemented on particle %d\n", k);
	if (a.y<-L/2.0) a.y+=L; //printf("PBCs (bottom) implemented on particle %d\n", k);
	return a;
}
/* ------- ---------- ------- */

/* ------- DISTANCE BETWEEN PARTICLES k and j ------- */
double distance_particles(int k, int j) {
	double rdist;
	double px, py, nx, ny, dx, dy;
	// Positions
	px = particles[k].x;
	py = particles[k].y;
	nx = particles[j].x;
	ny = particles[j].y;
	// Distances
	dx = fabs(px - nx);
	dy = fabs(py - ny);
	// PBCs
	if (dx > L/2.0) {dx = L - dx;}
	if (dy > L/2.0) {dy = L - dy;}
	// Radial distance
	rdist = sqrt(pow(dx,2) + pow(dy,2));
	return rdist;
}
/* ------- ---------------------------------- ------- */

/* ------- DISTANCE BETWEEN VECTORS part AND neigh ------- */
double distance_vectors(vec2D part, vec2D neigh) {
	double rdist;
	double px, py, nx, ny, dx, dy;
	// Positions
	px = part.x;
	py = part.y;
	nx = neigh.x;
	ny = neigh.y;
	// Distances
	dx = fabs(px - nx);
	dy = fabs(py - ny);
	// PBCs
	if (dx > L/2.0) {dx = L - dx;}
	if (dy > L/2.0) {dy = L - dy;}
	// Radial distance
	rdist = sqrt(pow(dx,2) + pow(dy,2));
	return rdist;
}
/* ------- ---------------------------------- ------- */

/* ------- UNITARY VECTOR BETWEEN PARTICLES k and j ------- */
vec2D uvec_particles(int k, int j) {
	double rdist;
	double px, py, nx, ny, dx, dy;
	vec2D ur;
	// Positions
	px = particles[k].x;
	py = particles[k].y;
	nx = particles[j].x;
	ny = particles[j].y;
	// PBCs
	if (px - nx < -L/2.0) {nx -= L;}
	if (py - ny < -L/2.0) {ny -= L;}
	if (px - nx > L/2.0) {nx += L;}
	if (py - ny > L/2.0) {ny += L;}
	// Distances
	dx = px - nx;
	dy = py - ny;
	// Radial distance
	rdist = sqrt(pow(dx,2) + pow(dy,2));
	ur.x = dx/rdist;
	ur.y = dy/rdist;
	return ur;
}
/* ------- ----------------------------------------- ------- */

/* ------- VERLET LIST ------- */
void new_vlist() {
	int idx, ipart;
	double dx, dy, dr;
	// Reset
	for (idx = 0; idx < N; idx++) {
		if (particles[idx].exists) {
			particles[idx].nVlist = 0;
			particles[idx].xV = particles[idx].x;
			particles[idx].yV = particles[idx].y;
			for (ipart = 0; ipart < 20; ipart++) {particles[idx].Vlist[ipart] = -1;}
		}
	}
	// Run through particles to redefine lists
	for (idx = 0; idx < N-1; idx++) {
		if (particles[idx].exists) {
			for (ipart = idx+1; ipart < N; ipart++) {
				if (particles[ipart].exists) {
					// dx = fabs(particles[idx].x - particles[ipart].x);
					// dy = fabs(particles[idx].x - particles[ipart].x);
					// if (dx > L/2.0) {dx = L - dx;}
					// if (dy > L/2.0) {dy = L - dy;}
					// dr = sqrt(dx*dx+dy*dy);
					dr = distance_particles(idx, ipart);
					if (dr <= vrange) {
						particles[idx].nVlist++;
						particles[idx].Vlist[particles[idx].nVlist-1] = ipart;
						particles[ipart].nVlist++;
						particles[ipart].Vlist[particles[ipart].nVlist-1] = idx;
						if (long_verbose) {printf("i = %d - j = %d - r = %.2f - nni = %d - nnj = %d\n", idx, ipart, dr, particles[idx].nVlist, particles[ipart].nVlist);}
					}
				}
			}
		}
	}
}
/* ------- ----------- ------- */

/* ------- SET PARTICLE ------- */
void set_particle(int idx, int id, int findx, double x, double y, int bplus, int bminus) {
	particles[idx].id = id;
	particles[idx].fil = findx;
	particles[idx].x = x;
	particles[idx].y = y;
	particles[idx].bplus = bplus;
	particles[idx].bminus = bminus;
	particles[idx].angle = 0.0;
	particles[idx].exists = true;
	particles[idx].tnuc = (int)(step/DYNSTEP);
	particles[idx].nVlist = 0;
	particles[idx].xV = x;
	particles[idx].yV = y;
	return;
}
/* ------- ----------- ------- */

/* ------- REMOVE PARTICLE ------- */
void remove_particle(int idx) {
	int ipart;
	particles[idx].id = -1;
	particles[idx].fil = -1;
	particles[idx].x = 0.0;
	particles[idx].y = 0.0;
	particles[idx].bplus = -1;
	particles[idx].bminus = -1;
	particles[idx].exists = false;
	particles[idx].nVlist = 0;
	for (ipart = 0; ipart < 20; ipart++) {particles[idx].Vlist[ipart] = -1;}
	return;
}
/* ------- --------------- ------- */

/* ------- CHECK VIABILITY OF NEW POSITION ------- */
bool check_viability(double posx, double posy, int Np) {
	double rdist;
	int ipart;
	double px, py, dx, dy;
	for (ipart = 0; ipart < Np; ipart++) {
		// Positions
		px = particles[ipart].x;
		py = particles[ipart].y;
		// Distances
		dx = fabs(px - posx);
		dy = fabs(py - posy);
		// PBCs
		if (dx > L/2.0) {dx = L - dx;}
		if (dy > L/2.0) {dy = L - dy;}
		// Radial distance
		rdist = sqrt(dx*dx + dy*dy);
		if (rdist < 1.0) {return false;}
	}
	return true;
}
/* ------- ---------------------------------- ------- */

/* ------- INITIALISE SYSTEM ------- */
void init_system(void) {
	// Define variables and arrays
	int ifil, ip, idx;
	int lfil;
	double pos0x, pos0y, ang;
	int bpidx, bmidx;
	int headidx, tailidx;
	double dx, dy;
	particles = (PARTICLE *) malloc((1000 * N + 1) * sizeof(PARTICLE));
	filaments = (FILAMENT *) malloc((1000 * P + 1) * sizeof(FILAMENT));
	// Build filaments
	ifil = 0;
	Nparts = 0;
	maxid = 0;
	Nfils = (int) floor(drand48()*5);
	if (ICempty) {Nfils = 0;}
	while (ifil < Nfils) {
		lfil = (int) floor(drand48()*20) + 2;
		pos0x = L*(drand48()-0.5);
		pos0y = L*(drand48()-0.5);
		ang = 2*M_PI*drand48();
		dx = cos(ang);
		dy = sin(ang);
		ip = 0;
		while (ip < lfil) {
			if (ip == 0) {bpidx = Nparts+1; bmidx = -1; tailidx = Nparts;}
			else if (ip == lfil-1) {bpidx = -1; bmidx = Nparts-1; headidx = Nparts;}
			else {bpidx = Nparts+1; bmidx = Nparts-1;}
			set_particle(Nparts, maxid+1, ifil, pos0x+ip*dx, pos0y+ip*dy, bpidx, bmidx);
			ip++;
			Nparts++;
			maxid++;
		}
		filaments[ifil].id = ifil+1;
		filaments[ifil].head = headidx;
		filaments[ifil].tail = tailidx;
		filaments[ifil].n = lfil;
		filaments[ifil].exists = true;
		filaments[ifil].tnuc = 0;
		ifil++;
	}
	printf("\nSystem initialised!\n");
}
/* ------- ----------------- ------- */

/* ------- INITIALISE SYSTEM ------- */
void wrong_init_system(void) {
	// Define variables and arrays
	int ifil, ip, idx, ipif;
	int lfil;
	double pos0x, pos0y, ang, posx, posy;
	int bpidx, bmidx;
	int headidx, tailidx;
	double dx, dy;
	bool accept;
	particles = (PARTICLE *) malloc((1000 * N + 1) * sizeof(PARTICLE));
	filaments = (FILAMENT *) malloc((1000 * P + 1) * sizeof(FILAMENT));
	// Build filaments
	ifil = 0;
	Nparts = 0;
	maxid = 0;
	Nfils = (int) floor(drand48()*5);
	if (ICempty) {Nfils = 0;}
	while (ifil < Nfils) {
		lfil = (int) floor(drand48()*20) + 2;
		pos0x = L*(drand48()-0.5);
		pos0y = L*(drand48()-0.5);
		ang = 2*M_PI*drand48();
		dx = cos(ang);
		dy = sin(ang);
		ip = 0;
		ipif = 0;
		accept = true;
		while (accept && ip < lfil) {
			if (ip == 0) {bpidx = Nparts+1; bmidx = -1; tailidx = Nparts;}
			else if (ip == lfil-1) {bpidx = -1; bmidx = Nparts-1; headidx = Nparts;}
			else {bpidx = Nparts+1; bmidx = Nparts-1;}
			posx = pos0x+ip*dx;
			if (posx < -(double)L/2.0) {posx += L;}
			if (posx > (double)L/2.0) {posx -= L;}
			posy = pos0y+ip*dy;
			if (posy < -(double)L/2.0) {posy += L;}
			if (posy > (double)L/2.0) {posy -= L;}
			accept = check_viability(posx, posy, Nparts);
			printf("fil %d posx %f posy %f\n", ifil, posx, posy);
			if (accept) {
				set_particle(Nparts, maxid+1, ifil, posx, posy, bpidx, bmidx);
				ip++;
				Nparts++;
				maxid++;
			}
		}
		if (ip > 0) {
			filaments[ifil].id = ifil+1;
			filaments[ifil].head = headidx;
			filaments[ifil].tail = tailidx;
			filaments[ifil].n = ip;
			filaments[ifil].exists = true;
			filaments[ifil].tnuc = 0;
			ifil++;
		}
	}
	printf("\nSystem initialised!\n");
}
/* ------- ----------------- ------- */

/* ------- INITIALISE BANDS SYSTEM ------- */
void init_bands_system(void) {
	// Define variables and arrays
	int ifil, ip, idx, irow, ifrow, nrows, nfprow, lfil, iband, wband;
	double pos0x, pos0y;
	int bpidx, bmidx;
	int headidx, tailidx;
	double dx = 1.0, dy = 0.0;
	particles = (PARTICLE *) malloc((1000 * N + 1) * sizeof(PARTICLE));
	filaments = (FILAMENT *) malloc((1000 * P + 1) * sizeof(FILAMENT));
	// Build filaments
	ifil = 0;
	Nparts = 0;
	maxid = 0;
	lfil = (int)(-kon*mean_time*log(1-kon));
	if (lfil > L) {lfil = (int)((double)L/2.0);}
	if (lfil < 2) {lfil = 2;}
	lfil = L-1;
	// nrows = L+1;
	nrows = (int)(0.8*(double)L);
	wband = nrows/nbands;
	nfprow = (int)((double)L/(double)lfil);
	Nfils = (int)(nrows*nfprow);
	for (irow = 0; irow < nrows; irow++) {
		// pos0y = (double)irow+0.5-(double)L/2.0;
		pos0y = ((double)irow*1.25)+0.625-(double)L/2.0;
		iband = floor(irow/wband);
		for (ifrow = 0; ifrow < nfprow; ifrow++) {
			if (iband%2 == 0) {
				pos0x = ((double)ifrow*lfil)+1.0-(double)L/2.0;
				dx = 1.0;
			}
			else {
				pos0x = ((double)L/2.0-(((double)ifrow*lfil)+1.0));
				dx = -1.0;
			}
			ip = 0;
			while (ip < lfil) {
				if (ip == 0) {bpidx = Nparts+1; bmidx = -1; tailidx = Nparts;}
				else if (ip == lfil-1) {bpidx = -1; bmidx = Nparts-1; headidx = Nparts;}
				else {bpidx = Nparts+1; bmidx = Nparts-1;}
				set_particle(Nparts, maxid+1, ifil, pos0x+ip*dx, pos0y+ip*dy, bpidx, bmidx);
				ip++;
				Nparts++;
				maxid++;
			}
			filaments[ifil].id = ifil+1;
			filaments[ifil].head = headidx;
			filaments[ifil].tail = tailidx;
			filaments[ifil].n = lfil;
			filaments[ifil].exists = true;
			filaments[ifil].tnuc = 0;
			ifil++;
		}
	}
	printf("\nSystem initialised!\n");
}
/* ------- ----------------- ------- */

/* ------- INITIALISE TEST SYSTEM ------- */
void init_test_system(void) {
	// Define variables and arrays
	int ifil, ip, idx;
	int lfil;
	double pos0x, pos0y, ang;
	int bpidx, bmidx;
	int headidx, tailidx;
	double dx, dy;
	particles = (PARTICLE *) malloc((100 * N + 1) * sizeof(PARTICLE));
	filaments = (FILAMENT *) malloc((100 * P + 1) * sizeof(FILAMENT));
	// Build filaments
	ifil = 0;
	Nparts = 0;
	maxid = 0;
	Nfils = 1;
	while (ifil < Nfils) {
		// lfil = 10;
		// pos0x = -5.5;
		// pos0y = 0.0;
		lfil = lnuctest;
		// pos0x = -0.5;
		// pos0y = 0.0;
		// dx = 1.0;
		// dy = 0.0;
		pos0x = -0.5*lnuctest*cos(tangle/180*M_PI);
		pos0y = -0.5*lnuctest*sin(tangle/180*M_PI);
		dx = 1.0*cos(tangle/180*M_PI);
		dy = 1.0*sin(tangle/180*M_PI);
		ip = 0;
		while (ip < lfil) {
			if (ip == 0) {bpidx = Nparts+1; bmidx = -1; tailidx = Nparts;}
			else if (ip == lfil-1) {bpidx = -1; bmidx = Nparts-1; headidx = Nparts;}
			else {bpidx = Nparts+1; bmidx = Nparts-1;}
			if (flat_mem && !tstraight) {
				dx = 1.0*cos(ip*2.0/Rfil);
				dy = 1.0*sin(ip*2.0/Rfil);
				pos0x += dx;
				pos0y += dy;
				set_particle(Nparts, maxid+1, ifil, pos0x, pos0y, bpidx, bmidx);
			}
			else {
				set_particle(Nparts, maxid+1, ifil, pos0x+ip*dx, pos0y+ip*dy, bpidx, bmidx);
			}
			// set_particle(Nparts, Nparts+1, ifil, pos0x+ip*dx, pos0y+ip*dy, bpidx, bmidx);
			ip++;
			Nparts++;
			maxid++;
		}
		filaments[ifil].id = ifil+1;
		filaments[ifil].head = headidx;
		filaments[ifil].tail = tailidx;
		filaments[ifil].n = lfil;
		filaments[ifil].exists = true;
		filaments[ifil].tnuc = 0;
		ifil++;
	}
	printf("\nSystem initialised!\n");
}
/* ------- ----------------- ------- */

/* ------- DUMP MOVIE ------- */
void dump_movie_frame(int step) {
	int ip;
	int nparts=0;
	for (ip = 0; ip < N; ip++) {
		if (particles[ip].exists) {nparts++;}
	}
	dumpf=fopen(movie_file, "a");
	fprintf(dumpf, "ITEM: TIMESTEP\n");
	fprintf(dumpf, "%d\n",step);
	fprintf(dumpf, "ITEM: NUMBER OF ATOMS\n");
	fprintf(dumpf, "%d\n",nparts);
	fprintf(dumpf, "ITEM: BOX BOUNDS pp ff ff\n");
	fprintf(dumpf, "-%.2f %.2f\n",L/2.0,L/2.0);
	fprintf(dumpf, "-%.2f %.2f\n",L/2.0,L/2.0);
	fprintf(dumpf, "-0.25 0.25\n");
	fprintf(dumpf, "ITEM: ATOMS id type mol x y z a\n");
	for (ip = 0; ip < N; ip++) {
		if (particles[ip].exists) {fprintf(dumpf, "%d %d %d %.2f %.2f %.2f %.2f\n",particles[ip].id,1,filaments[particles[ip].fil].id,particles[ip].x,particles[ip].y,0.0,particles[ip].angle);}
	}
	fclose(dumpf);
	rejectf=fopen(reject_file, "a");
	fprintf(rejectf, "%d %d %d %f\n", step, xytrials, xyrejects, (double)xyrejects/(double)xytrials);
	fclose(rejectf);
	xytrials = 0;
	xyrejects = 0;
	return;
}
/* ------- ---------- ------- */

/* ------- COMPUTE ENERGY IN THE BONDS AROUND 1 PARTICLE ------- */
double bonds_energy(int idx) {
	// Kbond = 1000
	int i, bidx;
	double r;
	double Ebs = 0;
	// Plus direction
	bidx = particles[idx].bplus;
	if (bidx >= 0) {
		r = distance_particles(idx, bidx);
		Ebs += Kbond*pow(r-1.0,2);
	}
	// Minus direction
	bidx = particles[idx].bminus;
	if (bidx >= 0) {
		r = distance_particles(idx, bidx);
		Ebs += Kbond*pow(r-1.0,2);
	}
	if (long_verbose) {printf("idx = %d -- Ebonds = %f\n", idx, Ebs);}
	return Ebs;
}
/* ------- --------------------------------------------- ------- */

/* ------- COMPUTE ENERGY IN ANGLES AROUND 1 PARTICLE ------- */
double bend_energy(int idx) {
	double Ebnd = 0;
	int idx0, idx1, idx2;
	vec2D r1, r2, mr1, ux;
	double ang, pang1, pang2, dang;
	ux.x = 1.0;
	ux.y = 0.0;
	// Self bending E
	if ((particles[idx].bplus >= 0) && (particles[idx].bminus >= 0)) {			// edge particles don't have self bending E
		idx1 = particles[idx].bplus;
		idx2 = particles[idx].bminus;
		r1 = uvec_particles(idx, idx1);
		r2 = uvec_particles(idx, idx2);
		mr1.x = -r1.x;
		mr1.y = -r1.y;
		ang = acos(dot(r1,r2)); 
		pang1 = 180/M_PI*acos(dot(r2,ux));
		if (r2.y < 0.0) {pang1 = 360-pang1;}
		pang2 = 180/M_PI*acos(dot(mr1,ux));
		if (mr1.y < 0.0) {pang2 = 360-pang2;}
		particles[idx].angle = 0.5*(pang1+pang2);
		dang = fabs(pang1-pang2);
		if (dang > 180) {particles[idx].angle -= 180;}
		if (particles[idx].angle < 0) {particles[idx].angle += 360;}
		if (particles[idx].angle >= 360) {particles[idx].angle -= 360;}
		if (long_verbose) {printf("angself = %f -- ", ang);}
		if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/2-1/Rfil,2); printf("0 -> ta: %f - a: %f - da: %f - c: %f - dc: %f - dE: %f\n", (1.0-2.0/Rfil/M_PI)*180, ang/M_PI*180, (1.0-2.0/Rfil/M_PI)*180-ang/M_PI*180, (M_PI-ang)/2, (M_PI-ang)/2-1/Rfil, 0.5*Kbend*pow((M_PI-ang)/2-1/Rfil,2));}
		else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/2,2);}
		// Ebnd += 0.5*Kbend*pow((M_PI-ang)/2,2);
	}
	else {
		if (particles[idx].bplus >= 0) {
			idx1 = particles[idx].bplus;
			r1 = uvec_particles(idx, idx1);
			mr1.x = -r1.x;
			mr1.y = -r1.y;
			pang1 = 180/M_PI*acos(dot(mr1,ux));
			if (mr1.y < 0.0) {pang1 = 360-pang1;}
			particles[idx].angle = pang1;
			if (particles[idx].angle < 0) {particles[idx].angle += 360;}
			if (particles[idx].angle >= 360) {particles[idx].angle -= 360;}
		}
		if (particles[idx].bminus>= 0) {
			idx2 = particles[idx].bminus;
			r2 = uvec_particles(idx, idx2);
			pang2 = 180/M_PI*acos(dot(r2,ux));
			if (r2.y < 0.0) {pang2 = 360-pang2;}
			particles[idx].angle = pang2;
			if (particles[idx].angle < 0) {particles[idx].angle += 360;}
			if (particles[idx].angle >= 360) {particles[idx].angle -= 360;}
		}
	}
	// Plus neigh bending E
	if (particles[idx].bplus >= 0) {
		idx0 = particles[idx].bplus;
		idx1 = particles[idx0].bplus;
		idx2 = particles[idx0].bminus;
		if ((idx1 >= 0) && (idx2 >= 0)) {
			r1 = uvec_particles(idx0, idx1);
			r2 = uvec_particles(idx0, idx2);
			ang = acos(dot(r1,r2));
			if (long_verbose) {printf("angplus = %f -- ", ang);}
			if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/2-1/Rfil,2); printf("+ -> ta: %f - a: %f - da: %f - c: %f - dc: %f - dE: %f\n", (1.0-2.0/Rfil/M_PI)*180, ang/M_PI*180, (1.0-2.0/Rfil/M_PI)*180-ang/M_PI*180, (M_PI-ang)/2, (M_PI-ang)/2-1/Rfil, 0.5*Kbend*pow((M_PI-ang)/2-1/Rfil,2));}
			else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/2,2);}
			// Ebnd += 0.5*Kbend*pow((M_PI-ang)/2,2);
		}
	}
	// Minus neigh bending E
	if (particles[idx].bminus >= 0) {
		idx0 = particles[idx].bminus;
		idx1 = particles[idx0].bplus;
		idx2 = particles[idx0].bminus;
		if ((idx1 >= 0) && (idx2 >= 0)) {
			r1 = uvec_particles(idx0, idx1);
			r2 = uvec_particles(idx0, idx2);
			ang = acos(dot(r1,r2));
			if (long_verbose) {printf("angminus = %f -- ", ang);}
			if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/2-1/Rfil,2); printf("- -> ta: %f - a: %f - da: %f - c: %f - dc: %f - dE: %f\n", (1.0-2.0/Rfil/M_PI)*180, ang/M_PI*180, (1.0-2.0/Rfil/M_PI)*180-ang/M_PI*180, (M_PI-ang)/2, (M_PI-ang)/2-1/Rfil, 0.5*Kbend*pow((M_PI-ang)/2-1/Rfil,2));}
			else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/2,2);}
			// Ebnd += 0.5*Kbend*pow((M_PI-ang)/2,2);
		}
	}
	if (long_verbose) {printf("idx = %d -- Ebend = %f\n", idx, Ebnd);}
	return Ebnd;
}
/* ------- ------------------------------------------ ------- */

/* ------- COMPUTE ENERGY IN ANGLES AROUND 1 PARTICLE ------- */
double old_bend_energy(int idx) {
	double Ebnd = 0;
	int idx0, idx1, idx2;
	vec2D r1, r2;
	double ang;
	// Self bending E
	if ((particles[idx].bplus >= 0) && (particles[idx].bminus >= 0)) {			// edge particles don't have self bending E
		idx1 = particles[idx].bplus;
		idx2 = particles[idx].bminus;
		r1 = uvec_particles(idx, idx1);
		r2 = uvec_particles(idx, idx2);
		ang = acos(dot(r1,r2));
		if (long_verbose) {printf("angself = %f -- ", ang);}
		if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/2-1/Rfil,2); printf("0 -> ta: %f - a: %f - da: %f - c: %f - dc: %f - dE: %f\n", (1.0-2.0/Rfil/M_PI)*180, ang/M_PI*180, (1.0-2.0/Rfil/M_PI)*180-ang/M_PI*180, (M_PI-ang)/2, (M_PI-ang)/2-1/Rfil, 0.5*Kbend*pow((M_PI-ang)/2-1/Rfil,2));}
		else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/2,2);}
		// Ebnd += 0.5*Kbend*pow((M_PI-ang)/2,2);
	}
	// Plus neigh bending E
	if (particles[idx].bplus >= 0) {
		idx0 = particles[idx].bplus;
		idx1 = particles[idx0].bplus;
		idx2 = particles[idx0].bminus;
		if ((idx1 >= 0) && (idx2 >= 0)) {
			r1 = uvec_particles(idx0, idx1);
			r2 = uvec_particles(idx0, idx2);
			ang = acos(dot(r1,r2));
			if (long_verbose) {printf("angplus = %f -- ", ang);}
			if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/2-1/Rfil,2); printf("+ -> ta: %f - a: %f - da: %f - c: %f - dc: %f - dE: %f\n", (1.0-2.0/Rfil/M_PI)*180, ang/M_PI*180, (1.0-2.0/Rfil/M_PI)*180-ang/M_PI*180, (M_PI-ang)/2, (M_PI-ang)/2-1/Rfil, 0.5*Kbend*pow((M_PI-ang)/2-1/Rfil,2));}
			else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/2,2);}
			// Ebnd += 0.5*Kbend*pow((M_PI-ang)/2,2);
		}
	}
	// Minus neigh bending E
	if (particles[idx].bminus >= 0) {
		idx0 = particles[idx].bminus;
		idx1 = particles[idx0].bplus;
		idx2 = particles[idx0].bminus;
		if ((idx1 >= 0) && (idx2 >= 0)) {
			r1 = uvec_particles(idx0, idx1);
			r2 = uvec_particles(idx0, idx2);
			ang = acos(dot(r1,r2));
			if (long_verbose) {printf("angminus = %f -- ", ang);}
			if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/2-1/Rfil,2); printf("- -> ta: %f - a: %f - da: %f - c: %f - dc: %f - dE: %f\n", (1.0-2.0/Rfil/M_PI)*180, ang/M_PI*180, (1.0-2.0/Rfil/M_PI)*180-ang/M_PI*180, (M_PI-ang)/2, (M_PI-ang)/2-1/Rfil, 0.5*Kbend*pow((M_PI-ang)/2-1/Rfil,2));}
			else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/2,2);}
			// Ebnd += 0.5*Kbend*pow((M_PI-ang)/2,2);
		}
	}
	if (long_verbose) {printf("idx = %d -- Ebend = %f\n", idx, Ebnd);}
	return Ebnd;
}
/* ------- ------------------------------------------ ------- */

/* ------- COMPUTE ENERGY IN ANGLES AROUND 1 PARTICLE ------- */
double confused_bend_energy(int idx) {
	// Kbend = 0.5
	double Ebnd = 0;
	int idx0, idx1, idx2, idx11, idx22, idx111, idx222;
	vec2D r1, r2, r11, r22, r111, r222;
	double ang;
	// Self bending E
	if ((particles[idx].bplus >= 0) && (particles[idx].bminus >= 0)) {			// edge particles don't have self bending E
		idx1 = particles[idx].bplus;
		idx2 = particles[idx].bminus;
		r1 = uvec_particles(idx, idx1);
		r2 = uvec_particles(idx, idx2);
		// ang = acos(dot(r1,r2))/M_PI*180;
		ang = acos(dot(r1,r2));
		if (long_verbose) {printf("angself = %f -- ", ang);}
		// Ebnd += Kbend*pow(ang-180,2);
		// Ebnd += Kbend*pow(ang,2);
		// Ebnd += Kbend*pow(1+dot(r1,r2),2);
		if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/2-1/Rfil,2); printf("0 -> ta: %f - a: %f - da: %f - c: %f - dc: %f - dE: %f\n", (1.0-2.0/Rfil/M_PI)*180, ang/M_PI*180, (1.0-2.0/Rfil/M_PI)*180-ang/M_PI*180, (M_PI-ang)/2, (M_PI-ang)/2-1/Rfil, 0.5*Kbend*pow((M_PI-ang)/2-1/Rfil,2));}
		else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/2,2);}
		// Ebnd += 0.5*Kbend*pow((M_PI-ang)/2,2);
		// Extend to 2nd neighs!!
		if (particles[idx1].bplus >= 0) {
			idx11 = particles[idx1].bplus;
			r11 = uvec_particles(idx1, idx11);
			ang = acos(dot(r11,r2));
			if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/4-1/Rfil,2); printf("0- -> ta: %f - a: %f - da: %f - c: %f - dc: %f - dE: %f\n", (1.0-4.0/Rfil/M_PI)*180, ang/M_PI*180, (1.0-4.0/Rfil/M_PI)*180-ang/M_PI*180, (M_PI-ang)/4, (M_PI-ang)/4-1/Rfil, 0.5*Kbend*pow((M_PI-ang)/4-1/Rfil,2));}
			else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/4,2);}
			// Extend to 3rd neighs!!
			if (particles[idx11].bplus >= 0) {
				idx111 = particles[idx11].bplus;
				r111 = uvec_particles(idx11, idx111);
				ang = acos(dot(r111,r2));
				if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/6-1/Rfil,2);}
				else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/6,2);}
			}
		}
		if (particles[idx2].bminus >= 0) {
			idx22 = particles[idx2].bminus;
			r22 = uvec_particles(idx2, idx22);
			ang = acos(dot(r1,r22));
			if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/4-1/Rfil,2); printf("0+ -> ta: %f - a: %f - da: %f - c: %f - dc: %f - dE: %f\n", (1.0-4.0/Rfil/M_PI)*180, ang/M_PI*180, (1.0-4.0/Rfil/M_PI)*180-ang/M_PI*180, (M_PI-ang)/4, (M_PI-ang)/4-1/Rfil, 0.5*Kbend*pow((M_PI-ang)/4-1/Rfil,2));}
			else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/4,2);}
			// Extend to 3rd neighs!!
			if (particles[idx22].bminus >= 0) {
				idx222 = particles[idx22].bminus;
				r222 = uvec_particles(idx22, idx222);
				ang = acos(dot(r1,r222));
				if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/6-1/Rfil,2);}
				else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/6,2);}
			}
		}
	}
	// Plus neigh bending E
	if (particles[idx].bplus >= 0) {
		idx0 = particles[idx].bplus;
		idx1 = particles[idx0].bplus;
		idx2 = particles[idx0].bminus;
		if ((idx1 >= 0) && (idx2 >= 0)) {
			r1 = uvec_particles(idx0, idx1);
			r2 = uvec_particles(idx0, idx2);
			// ang = acos(dot(r1,r2))/M_PI*180;
			ang = acos(dot(r1,r2));
			if (long_verbose) {printf("angplus = %f -- ", ang);}
			// Ebnd += Kbend*pow(ang-180,2);
			// Ebnd += Kbend*pow(ang,2);
			// Ebnd += Kbend*pow(1+dot(r1,r2),2);
			if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/2-1/Rfil,2); printf("+ -> ta: %f - a: %f - da: %f - c: %f - dc: %f - dE: %f\n", (1.0-2.0/Rfil/M_PI)*180, ang/M_PI*180, (1.0-2.0/Rfil/M_PI)*180-ang/M_PI*180, (M_PI-ang)/2, (M_PI-ang)/2-1/Rfil, 0.5*Kbend*pow((M_PI-ang)/2-1/Rfil,2));}
			else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/2,2);}
			// Ebnd += 0.5*Kbend*pow((M_PI-ang)/2,2);
			// Extend to 2nd neighs!!
			if (particles[idx1].bplus >= 0) {
				idx11 = particles[idx1].bplus;
				r11 = uvec_particles(idx1, idx11);
				ang = acos(dot(r11,r2));
				if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/4-1/Rfil,2); printf("+- -> ta: %f - a: %f - da: %f - c: %f - dc: %f - dE: %f\n", (1.0-4.0/Rfil/M_PI)*180, ang/M_PI*180, (1.0-4.0/Rfil/M_PI)*180-ang/M_PI*180, (M_PI-ang)/4, (M_PI-ang)/4-1/Rfil, 0.5*Kbend*pow((M_PI-ang)/4-1/Rfil,2));}
				else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/4,2);}
				// Extend to 3rd neighs!!
				if (particles[idx11].bplus >= 0) {
					idx111 = particles[idx11].bplus;
					r111 = uvec_particles(idx11, idx111);
					ang = acos(dot(r111,r2));
					if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/6-1/Rfil,2);}
					else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/6,2);}
				}
			}
			if (particles[idx2].bminus >= 0) {
				idx22 = particles[idx2].bminus;
				r22 = uvec_particles(idx2, idx22);
				ang = acos(dot(r1,r22));
				if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/4-1/Rfil,2); printf("++ -> ta: %f - a: %f - da: %f - c: %f - dc: %f - dE: %f\n", (1.0-4.0/Rfil/M_PI)*180, ang/M_PI*180, (1.0-4.0/Rfil/M_PI)*180-ang/M_PI*180, (M_PI-ang)/4, (M_PI-ang)/4-1/Rfil, 0.5*Kbend*pow((M_PI-ang)/4-1/Rfil,2));}
				else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/4,2);}
				// Extend to 3rd neighs!!
				if (particles[idx22].bminus >= 0) {
					idx222 = particles[idx22].bminus;
					r222 = uvec_particles(idx22, idx222);
					ang = acos(dot(r1,r222));
					if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/6-1/Rfil,2);}
					else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/6,2);}
				}
			}
		}
	}
	// Minus neigh bending E
	if (particles[idx].bminus >= 0) {
		idx0 = particles[idx].bminus;
		idx1 = particles[idx0].bplus;
		idx2 = particles[idx0].bminus;
		if ((idx1 >= 0) && (idx2 >= 0)) {
			r1 = uvec_particles(idx0, idx1);
			r2 = uvec_particles(idx0, idx2);
			// ang = acos(dot(r1,r2))/M_PI*180;
			ang = acos(dot(r1,r2));
			if (long_verbose) {printf("angminus = %f -- ", ang);}
			// Ebnd += Kbend*pow(ang-180,2);
			// Ebnd += Kbend*pow(ang,2);
			// Ebnd += Kbend*pow(1+dot(r1,r2),2);
			if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/2-1/Rfil,2); printf("- -> ta: %f - a: %f - da: %f - c: %f - dc: %f - dE: %f\n", (1.0-2.0/Rfil/M_PI)*180, ang/M_PI*180, (1.0-2.0/Rfil/M_PI)*180-ang/M_PI*180, (M_PI-ang)/2, (M_PI-ang)/2-1/Rfil, 0.5*Kbend*pow((M_PI-ang)/2-1/Rfil,2));}
			else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/2,2);}
			// Ebnd += 0.5*Kbend*pow((M_PI-ang)/2,2);
			// Extend to 2nd neighs!!
			if (particles[idx1].bplus >= 0) {
				idx11 = particles[idx1].bplus;
				r11 = uvec_particles(idx1, idx11);
				ang = acos(dot(r11,r2));
				if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/4-1/Rfil,2); printf("-- -> ta: %f - a: %f - da: %f - c: %f - dc: %f - dE: %f\n", (1.0-4.0/Rfil/M_PI)*180, ang/M_PI*180, (1.0-4.0/Rfil/M_PI)*180-ang/M_PI*180, (M_PI-ang)/4, (M_PI-ang)/4-1/Rfil, 0.5*Kbend*pow((M_PI-ang)/4-1/Rfil,2));}
				else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/4,2);}
				// Extend to 3rd neighs!!
				if (particles[idx11].bplus >= 0) {
					idx111 = particles[idx11].bplus;
					r111 = uvec_particles(idx11, idx111);
					ang = acos(dot(r111,r2));
					if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/6-1/Rfil,2);}
					else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/6,2);}
				}
			}
			if (particles[idx2].bminus >= 0) {
				idx22 = particles[idx2].bminus;
				r22 = uvec_particles(idx2, idx22);
				ang = acos(dot(r1,r22));
				if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/4-1/Rfil,2); printf("-+ -> ta: %f - a: %f - da: %f - c: %f - dc: %f - dE: %f\n", (1.0-4.0/Rfil/M_PI)*180, ang/M_PI*180, (1.0-4.0/Rfil/M_PI)*180-ang/M_PI*180, (M_PI-ang)/4, (M_PI-ang)/4-1/Rfil, 0.5*Kbend*pow((M_PI-ang)/4-1/Rfil,2));}
				else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/4,2);}
				// Extend to 3rd neighs!!
				if (particles[idx22].bminus >= 0) {
					idx222 = particles[idx22].bminus;
					r222 = uvec_particles(idx22, idx222);
					ang = acos(dot(r1,r222));
					if (flat_mem) {Ebnd += 0.5*Kbend*pow((M_PI-ang)/6-1/Rfil,2);}
					else {Ebnd += 0.5*Kbend*pow((M_PI-ang)/6,2);}
				}
			}
		}
	}
	if (long_verbose) {printf("idx = %d -- Ebend = %f\n", idx, Ebnd);}
	return Ebnd;
}
/* ------- ------------------------------------------ ------- */

/* ------- COMPUTE ENERGY IN BINDING AND OVERLAPS FOR 1 PARTICLE ------- */
double bind_energy(int idx) {
	int nidx;
	double Ebind = 0;
	double r;
	// for (nidx = 0; nidx < N; nidx++) {
	// 	if (nidx == idx) {continue;}
	// 	if (particles[nidx].exists) {
	// 		if (filaments[particles[nidx].fil].id == filaments[particles[idx].fil].id) {continue;}
	// 		r = distance_particles(idx, nidx);
	// 		if (long_verbose) {printf("f0 = %d -- f1 = %d -- r = %f\n", filaments[particles[nidx].fil].id, filaments[particles[idx].fil].id, r);}
	// 		// if (step == 120000 && (particles[idx].id == 273 || particles[idx].id == 275)) {printf("f0 = %d -- f1 = %d -- r = %f\n", filaments[particles[nidx].fil].id, filaments[particles[idx].fil].id, r);}
	// 		if (r < 1.0) {Ebind += 1e6;}
	// 		else {Ebind += LJ(eps, sigma, arange, r);}
	// 	}
	// }
	for (nidx = 0; nidx < N; nidx++) {
		if (nidx == idx) {continue;}
		if (particles[nidx].exists) {
			if (particles[idx].bplus == nidx || particles[idx].bminus == nidx) {continue;}
			r = distance_particles(idx, nidx);
			if (r < 1.0) {
                if (olapINF) {Ebind += INFINITY;}
                else {Ebind += 1e6;}
            }
			else {
                if (filaments[particles[nidx].fil].id == filaments[particles[idx].fil].id) {
                    Ebind += 0.0;
                }
                Ebind += LJ(eps, sigma, arange, r);
            }
			if (long_verbose) {printf("f0 = %d -- f1 = %d -- r = %f\n", filaments[particles[nidx].fil].id, filaments[particles[idx].fil].id, r);}
		}
	}
	if (long_verbose) {printf("idx = %d -- Ebind = %f\n", idx, Ebind);}
	// if (step == 120000 && (particles[idx].id == 273 || particles[idx].id == 275)) {printf("idx = %d -- Ebind = %f\n", idx, Ebind);}
	return Ebind;
}
/* ------- ----------------------------------------------------- ------- */

/* ------- COMPUTE ENERGY IN BINDING AND OVERLAPS FOR 1 PARTICLE USING VERLET LISTS ------- */
double Vbind_energy(int idx) {
	int vidx, nidx;
	double Ebind = 0;
	double r;
	if (long_verbose) {printf("idx = %d -- neighs = %d\n", idx, particles[idx].nVlist);}
	for (vidx = 0; vidx < particles[idx].nVlist; vidx++) {
		nidx = particles[idx].Vlist[vidx];
		if (nidx == idx) {continue;}
		if (particles[nidx].exists) {
			if (particles[idx].bplus == nidx || particles[idx].bminus == nidx) {continue;}
			r = distance_particles(idx, nidx);
			if (r < 1.0) {
                if (olapINF) {Ebind += INFINITY;}
                else {Ebind += 1e6;}
            }
			else {
                if (filaments[particles[nidx].fil].id == filaments[particles[idx].fil].id) {
                    Ebind += 0.0;
                }
                Ebind += LJ(eps, sigma, arange, r);
            }
			if (long_verbose) {printf("i = %d - j = %d - fi = %d - fj = %d - r = %f - vi = %d - jj = %d\n", idx, nidx, filaments[particles[nidx].fil].id, filaments[particles[idx].fil].id, r, vidx, particles[idx].Vlist[vidx]);}
		}
	}
	if (long_verbose) {printf("idx = %d -- Ebind = %f\n", idx, Ebind);}
	return Ebind;
}
/* ------- ------------------------------------------------------------------------ ------- */

/* ------- COMPUTE ENERGY FROM MISALIGNMENT WITH PREFERRED DIRECTION ------- */
double direc_energy(int idx) {
	// Kdir = 10
	if (Kdir == 0.0) {return 0.0;}
	if (flat_mem) {return 0.0;}
	double Edir = 0;
	int fillen = filaments[particles[idx].fil].n;
	int idx1 = particles[idx].bplus;
	int idx2 = particles[idx].bminus;
	double proj, ang, xval, rloc, dc;
	vec2D r, ux;
	ux.x = 1.0;
	ux.y = 0.0;
	// Plus neigh bending E
	if (idx1 >= 0) {
		r = uvec_particles(idx, idx1);
		proj = fabs(dot(r,ux));
		ang = acos(proj);
		if (curvcos) {
			// dc = 1/Rfil*(1-cos(ang));
			rloc = Rmem/cos(ang);
		}
		else {
			xval = tan(ang)*sin(proj/Rmem);
			rloc = Rmem*sqrt(1+xval*xval);
		}
		dc = 1/Rfil - 1/rloc;
		if (flat_mem) {dc = 0.0;}
		Edir += 0.5*Kbend*dc*dc;
		// Edir += 0.5*Kdir*dc*dc*fillen;
		// Edir += Kdir*fillen*pow(1-proj,2);
		// if (long_verbose) {printf("idx = %d - idx1 = %d -- proj = %f\n", idx, idx1, proj);}
	}
	// Minus neigh bending E
	if (idx2 >= 0) {
		r = uvec_particles(idx, idx2);
		proj = fabs(dot(r,ux));
		ang = acos(proj);
		if (curvcos) {
			// dc = 1/Rfil*(1-cos(ang));
			rloc = Rmem/cos(ang);
		}
		else {
			xval = tan(ang)*sin(proj/Rmem);
			rloc = Rmem*sqrt(1+xval*xval);
		}
		dc = 1/Rfil - 1/rloc;
		if (flat_mem) {dc = 0.0;}
		Edir += 0.5*Kbend*dc*dc;
		// Edir += 0.5*Kdir*dc*dc*fillen;
		// Edir += Kdir*fillen*pow(1-proj,2);
		// if (long_verbose) {printf("idx = %d - idx2 = %d -- proj = %f\n", idx, idx2, proj);}
	}
	if (long_verbose) {printf("idx = %d -- Edir = %f\n", idx, Edir);}
	return Edir;
}
/* ------- --------------------------------------------------------- ------- */

/* ------- COMPUTE ENERGY FROM CURVATURE MISALIGNMENT ------- */
double curv_energy(int idx) {
	// Kdir = 10
	if (Kdir == 0.0) {return 0.0;}
	if (flat_mem) {return 0.0;}
	double Edir = 0;
	int filid = particles[idx].fil;
	int idx0, idx1, idx2;
	double proj, ang, xval, rloc, dc;
	double lfil = filaments[filid].n;
	double xfil;
	vec2D r, ux;
	ux.x = 1.0;
	ux.y = 0.0;
	for (idx0 = 0; idx0 < N; idx0++) {
		if (particles[idx0].fil != filid) {continue;}
		idx1 = particles[idx0].bplus;
		idx2 = particles[idx0].bminus;
		// xfil = distance_particles(idx0, filaments[filid].tail) - lfil/2.0;
		// Plus neigh bending E
		if (idx1 >= 0) {
			r = uvec_particles(idx0, idx1);
			proj = fabs(dot(r,ux));
			ang = acos(proj);
			xval = tan(ang)*sin(proj/Rmem);
			rloc = Rmem*sqrt(1+xval*xval);
			dc = 1/Rfil - 1/rloc;
			if (flat_mem) {dc = 0;}
			Edir += 0.5*Kbend*dc*dc;
			// Edir += Kdir*fillen*pow(1-proj,2);
			// if (long_verbose) {printf("idx = %d - idx1 = %d -- proj = %f\n", idx, idx1, proj);}
		}
		// Minus neigh bending E
		if (idx2 >= 0) {
			r = uvec_particles(idx0, idx2);
			proj = fabs(dot(r,ux));
			ang = acos(proj);
			xval = tan(ang)*sin(proj/Rmem);
			rloc = Rmem*sqrt(1+xval*xval);
			dc = 1/Rfil - 1/rloc;
			if (flat_mem) {dc = 0;}
			Edir += 0.5*Kbend*dc*dc;
			// Edir += Kdir*fillen*pow(1-proj,2);
			// if (long_verbose) {printf("idx = %d - idx2 = %d -- proj = %f\n", idx, idx2, proj);}
		}
	}
	if (long_verbose) {printf("idx = %d -- Edir = %f\n", idx, Edir);}
	return Edir;
}
/* ------- --------------------------------------------------------- ------- */

/* ------- TRANSLATION MOVE FOR 1 PARTICLE ------- */
void xyMC(int idx) {
	double oldE, newE;
	double dx, dy, dr;
	vec2D opos;
	double rand;
	xytrials++;
	// Check Verlet list validity
	dx = particles[idx].x - particles[idx].xV;
	dy = particles[idx].y - particles[idx].yV;
	dr = sqrt(dx*dx + dy*dy);
	if (dr > vthres) {
		new_vlist();
	}
	// Compute old E
	if (long_verbose) {printf("\nOLD E:\n");}
	// oldE = bonds_energy(idx) + bend_energy(idx) + bind_energy(idx) + direc_energy(idx);
	oldE = bonds_energy(idx) + bend_energy(idx) + Vbind_energy(idx) + direc_energy(idx);
	if (thermal_noise) {
		dx = sqrt(1.0/(double)DYNSTEP)*(drand48()-0.5);
		dy = sqrt(1.0/(double)DYNSTEP)*(drand48()-0.5);
	}
	else {
		dx = JUMP_MAX*(drand48()-0.5);
		dy = JUMP_MAX*(drand48()-0.5);
	}
	opos.x = particles[idx].x;
	opos.y = particles[idx].y;
	particles[idx].x += dx;
	particles[idx].y += dy;
	pbcs(idx);
	if (long_verbose) {printf("idx = %d -- dx = %f -- dy = %f -- npos = [%f, %f] -- opos = [%f, %f]\n", idx, dx, dy, particles[idx].x, particles[idx].y, opos.x, opos.y);}
	// Compute new E
	if (long_verbose) {printf("NEW E:\n");}
	// newE = bonds_energy(idx) + bend_energy(idx) + bind_energy(idx) + direc_energy(idx);
	newE = bonds_energy(idx) + bend_energy(idx) + Vbind_energy(idx) + direc_energy(idx);
	if (long_verbose) {printf("idx = %d -- oldE = %f -- newE = %f\n", idx, oldE, newE);}
	if (newE > oldE) {
		rand = drand48();
		if (rand > exp(-(newE-oldE))) {								// reject
			particles[idx].x = opos.x;
			particles[idx].y = opos.y;
			xyrejects++;
		}
	}
	return;
}
/* ------- ------------------------------- ------- */

/* ------- POLYMERISATION MOVE ------- */
void polMC(int fidx) {
	// kon = 0.6
	int hidx, nidx, idx;
	vec2D hpos, npos, r;
	double angbase, ang;
	bool findpos;
	int ntriallocal;
	int noverlaps;
	double dx, dy;
	double konl = kon;
	double rand;
	double rnorm;
	double frac;
	if (! filaments[fidx].exists) {return;}
	// Find head position
	hidx = filaments[fidx].head;
	hpos.x = particles[hidx].x;
	hpos.y = particles[hidx].y;
	nidx = particles[hidx].bminus;
	npos.x = particles[nidx].x;
	npos.y = particles[nidx].y;
		// PBCs
	if (hpos.x - npos.x < -L/2.0) {npos.x -= L;}
	if (hpos.y - npos.y < -L/2.0) {npos.y -= L;}
	if (hpos.x - npos.x > L/2.0) {npos.x += L;}
	if (hpos.y - npos.y > L/2.0) {npos.y += L;}
	r.x = hpos.x - npos.x;
	r.y = hpos.y - npos.y;
	angbase = atan(r.y/r.x);
	if (r.x < 0) {angbase += M_PI;}
	findpos = true;
	ntriallocal = 0;
	while (findpos) {
		ntriallocal++;
		// ang = angbase + (drand48()-0.5)*M_PI/4.0;
		ang = angbase + (drand48()-0.5)*M_PI/180*polang;
		// ang = angbase;
		npos.x = hpos.x + cos(ang);
		npos.y = hpos.y + sin(ang);
		npos = vec_pbcs(npos);
		noverlaps = 0;
		for (nidx = 0; nidx < N; nidx++) {
			if (particles[nidx].exists) {
				// Distances
				dx = fabs(npos.x - particles[nidx].x);
				dy = fabs(npos.y - particles[nidx].y);
				// PBCs
				if (dx > L/2.0) {dx = L - dx;}
				if (dy > L/2.0) {dy = L - dy;}
				rnorm = sqrt(pow(dx,2)+pow(dy,2));
				if (rnorm < 1.0) {noverlaps++;}
				// if (rnorm < 0.9) {noverlaps++;}
				// if (step == 187500 && rnorm < 1.0) {printf("fil id = %d -- nid = %d -- r = %f -- noverlaps = %d -- ntrials = %d\n", filaments[fidx].id, particles[nidx].id, rnorm, noverlaps, ntriallocal);}
			}
		}
		if (noverlaps == 0 || ntriallocal > 50) {findpos = false;}
	}
	NpolTrials++;
	if (fabs(npos.y) <= 4) {NpolTrialsBand++;}
	// if (step == 187500 && filaments[fidx].id == 49) {printf("noverlaps = %d -- ntriallocal = %d\n", noverlaps, ntriallocal);}
	if (noverlaps > 0) {NpolRejectedEV++; if (fabs(npos.y) <= 4) {NpolRejectedEVBand++;} return;}
	if (Min) {konl = konl*exp(-pow(npos.y,2)/(varMin));}
	if (nucleoid) {
		if (nucleoidMinGrow) {
			frac = MinRate*(kmax-kon)*(float)(step-RELSTEPS)/(float)DYNSTEP;
			if (frac > (kmax-kon)) {frac = kmax-kon;}
			if (frac < 0.0) {frac = 0.0;}
			konl = kon+frac*exp(-pow(npos.y,2)/(varMin));
		}
		else {
			if (step > RELSTEPS) {konl = kon+(kmax-kon)*exp(-pow(npos.y,2)/(varMin));}
		}
	}
	if (nucleoidMin) {
		if (nucleoidMinGrow) {
			double fracminl = kon;
			if (step > RELSTEPS) {fracminl = kon-MinRate*kon*(float)(step-RELSTEPS)/(float)DYNSTEP;}
			if (fracminl < 0.0) {fracminl = 0.0;}
			// frac = MinRate*kmax*(float)(step-RELSTEPS)/(float)DYNSTEP;
			// if (frac > kmax) {frac = kmax;}
			// if (frac < 0.0) {frac = 0.0;}
			// konl = fracminl+(frac-fracminl)*exp(-pow(npos.y,2)/(varMin));
			frac = MinRate*(kmax-kon)*(float)(step-RELSTEPS)/(float)DYNSTEP;
			if (frac > (kmax-fracminl)) {frac = kmax-fracminl;}
			if (frac < 0.0) {frac = 0.0;}
			konl = fracminl+frac*exp(-pow(npos.y,2)/(varMin));
		}
		else {
			if (step > RELSTEPS) {konl = kmax*exp(-pow(npos.y,2)/(varMin));}
		}
	}
	if (varMinDyn) {
		if (step < MCSTEPS-MCSTEPS/4.0) {konl = exp(-pow(npos.y,2)/(varMin))*step/(float)(MCSTEPS-MCSTEPS/4.0);}
		else {konl = exp(-pow(npos.y,2)/(varMin));}
	}
	if (MinGrow) {
		frac = MinRate*(float)(step-RELSTEPS)/(float)DYNSTEP;
		if (frac > 1.0) {frac = 1.0;}
		if (frac < 0.0) {frac = 0.0;}
		konl = kon*(1-fracmin)*frac*exp(-pow(npos.y,2)/(varMin))+fracmin;
	}
	// if (step%PROGDUMP == 0) {printf("kon_local = %f\n", konl);}
	rand = drand48();
	// if (step == 187500 && filaments[fidx].id == 49) {printf("rand = %f -- kon = %f\n", rand, konl);}
	if (rand < konl) {
		// Find free index in particles array
		for (idx = 0; idx < N; idx++) {
			if (! particles[idx].exists) {break;}
		}
		// if (step == 187500 && filaments[fidx].id == 49) {printf("adding particle\n");}
		// Add new particle
		Nparts++;
		maxid++;
		set_particle(idx, maxid, fidx, npos.x, npos.y, -1, filaments[fidx].head);
		particles[filaments[fidx].head].bplus = idx;
		filaments[fidx].head = idx;
		filaments[fidx].n++;
	}
	return;
}
/* ------- ------------------- ------- */

/* ------- DEPOLYMERISATION MOVE ------- */
void depolMC(int fidx) {
	// koff = 0.4
	double rand;
	int tidx, hidx, ntidx;
	double koffl = koff;
	if (PC19 && step-RELSTEPS > PC19time) {return;}
	if (filaments[fidx].n <= 2) {
		rand = drand48();
		if (sizer) {koffl = koff*(double)filaments[fidx].n;}
		if (timer) {koffl = koff*(double)(step/DYNSTEP-particles[tidx].tnuc);}
		if (etimer) {koffl = koff*(1.0-exp(-(double)(step/DYNSTEP-particles[tidx].tnuc)/(double)mean_time));}
		if (rand < koffl) {
			tidx = filaments[fidx].tail;
			hidx = filaments[fidx].head;
			remove_particle(tidx);
			Nparts--;
			remove_particle(hidx);
			Nparts--;
			filaments[fidx].exists = false;
		}
		return;
	}
	tidx = filaments[fidx].tail;
	rand = drand48();
	if (sizer) {koffl = koff*(double)filaments[fidx].n;}
	// if (timer) {koffl = koff*(double)(step/DYNSTEP-filaments[fidx].tnuc);}
	if (timer) {koffl = koff*(double)(step/DYNSTEP-particles[tidx].tnuc);}
	if (etimer) {koffl = koff*(1.0-exp(-(double)(step/DYNSTEP-particles[tidx].tnuc)/(double)mean_time));}
	if (rand < koffl) {
		ntidx = particles[tidx].bplus;
		filaments[fidx].tail = ntidx;
		filaments[fidx].n--;
		particles[ntidx].bminus = -1;
		remove_particle(tidx);
		Nparts--;
	}
	return;
}
/* ------- --------------------- ------- */

/* ------- NUCLEATION MOVE ------- */
void nucMC() {
	// knuc = 0.6
	int nidx, idx1, idx2, fidx;
	int attempts, noverlaps, nattempts, nnoverlaps;
	double dx, dy, r, angle, rand, knucl;
	bool findpos, nfindpos;
	double frac;
	vec2D postail, poshead;
	findpos = true;
	attempts = 0;
	knucl = knuc;
	// Find nucleating position
	while (findpos && attempts < 1000) {
		if (Min) {
			postail.x = L*(drand48()-0.5);
			postail.y = normal(0.0, sqrt(varMin));
			// postail.y = L*(drand48()-0.5);
		}
		else {
			postail.x = L*(drand48()-0.5);
			postail.y = L*(drand48()-0.5);
		}
		// Find overlaps
		noverlaps = 0;
		for (nidx = 0; nidx < N; nidx++) {
			if (particles[nidx].exists) {
				// Distances
				dx = fabs(postail.x - particles[nidx].x);
				dy = fabs(postail.y - particles[nidx].y);
				// PBCs
				if (dx > L/2.0) {dx = L - dx;}
				if (dy > L/2.0) {dy = L - dy;}
				r = sqrt(pow(dx,2)+pow(dy,2));
				if (r < 1.0) {noverlaps++;}
			}
		}
		if (noverlaps == 0) {
			nfindpos = true;
			nattempts = 0;
			while (nfindpos && nattempts < 1000) {
				angle = 2*M_PI*drand48();
				poshead.x = postail.x + cos(angle);
				poshead.y = postail.y + sin(angle);
				// Find overlaps
				nnoverlaps = 0;
				for (nidx = 0; nidx < N; nidx++) {
					if (particles[nidx].exists) {
						// Distances
						dx = fabs(poshead.x - particles[nidx].x);
						dy = fabs(poshead.y - particles[nidx].y);
						// PBCs
						if (dx > L/2.0) {dx = L - dx;}
						if (dy > L/2.0) {dy = L - dy;}
						r = sqrt(pow(dx,2)+pow(dy,2));
						if (r < 1.0) {nnoverlaps++;}
					}
				}
				if (nnoverlaps == 0) {nfindpos = false;}
				nattempts++;
			}
			findpos = false;
		}
		attempts++;
	}
	if (noverlaps > 0 || nnoverlaps > 0) {return;}
	// if (varMinDyn) {
	// 	if (step < MCSTEPS-MCSTEPS/4.0) {knucl = knuc*step/(float)(MCSTEPS-MCSTEPS/4.0);}
	// 	else {knucl = knuc;}
	// }
	if (Min) {knucl = knuc*exp(-pow(postail.y,2)/(varMin));}
	if (varMinDyn) {
		if (step < MCSTEPS-MCSTEPS/4.0) {knucl = exp(-pow(postail.y,2)/(varMin))*step/(float)(MCSTEPS-MCSTEPS/4.0);}
		else {knucl = exp(-pow(postail.y,2)/(varMin));}
	}
	if (nucleoid) {
		if (step > RELSTEPS) {knucl = knucl+(kmax/kon*knuc-knuc)*exp(-pow(postail.y,2)/(varMin));}
	}
	if (nucleoidMin) {
		if (step > RELSTEPS) {knucl = kmax/kon*knuc*exp(-pow(postail.y,2)/(varMin));}
	}
	if (MinGrow) {
		frac = MinRate*(float)(step-RELSTEPS)/(float)DYNSTEP;
		if (frac > 1.0) {frac = 1.0;}
		if (frac < 0.0) {frac = 0.0;}
		knucl = knuc*(1-fracmin)*frac*exp(-pow(postail.y,2)/(varMin))+fracmin;
	}
	// Accept according to knuc
	rand = drand48();
	if (rand < knucl) {
		// Find free index in filaments array
		for (fidx = 0; fidx < P; fidx++) {
			if (! filaments[fidx].exists) {break;}
		}
		// Add new filament
		Nfils++;
		filaments[fidx].exists = true;
		filaments[fidx].n = 2;
		filaments[fidx].id = Nfils+1;
		filaments[fidx].tnuc = (int)(step/DYNSTEP);
		// Find 1st free index in particles array
		for (idx1 = 0; idx1 < N; idx1++) {
			if (! particles[idx1].exists) {break;}
		}
		filaments[fidx].tail = idx1;
		particles[idx1].exists = true;
		// Find 2nd free index in particles array
		for (idx2 = 0; idx2 < N; idx2++) {
			if (! particles[idx2].exists) {break;}
		}
		filaments[fidx].head = idx2;
		particles[idx2].exists = true;
		// Add tail particle
		Nparts++;
		maxid++;
		set_particle(idx1, maxid, fidx, postail.x, postail.y, filaments[fidx].head, -1);
		// Add head particle
		Nparts++;
		maxid++;
		set_particle(idx2, maxid, fidx, poshead.x, poshead.y, -1, filaments[fidx].tail);
	}
	return;
}
/* ------- --------------- ------- */

/* ------- COMPUTE SYSTEM ENERGY ------- */
void compute_system_energy(int step) {
	double Ebonds = 0, Ebend = 0, Ebind = 0, Edir = 0, Etotal;
	int Nparts = 0;
	int idx;
	for (idx = 0; idx < N; idx++) {
		if (! particles[idx].exists) {continue;}
		Ebonds += bonds_energy(idx);
		Ebend += bend_energy(idx);
		// Ebind += bind_energy(idx);
		Ebind += Vbind_energy(idx);
		Edir += direc_energy(idx);
		Nparts += 1;
	}
	Etotal = Ebonds+Ebend+Ebind+Edir;
	energyf=fopen(energy_file, "a");
	fprintf(energyf, "%d %lf %lf %lf %lf %lf %d\n", step, Etotal, Ebonds, Ebend, Ebind, Edir, Nparts);
	fclose(energyf);
}
/* ------- --------------------- ------- */

/* ------- MAIN ------- */
int main(int argc, char *argv[]){

	// Timer variables
	clock_t start_t, end_t;
	double total_t;
	start_t = clock();

	// Default argument values
	sprintf(output_dir, "/Users/christian/");
	L=320;
	N=102400; //1600;
	P=51200; //800;
	MCSTEPS=3000000; //2000000; //MCSTEPS=500000; // 10 minutes
	RELSTEPS=600000; //100000; // 2 minutes
	// MCSTEPS=2100000; //2000000; //MCSTEPS=500000; // 7 minutes
	// RELSTEPS=150000; //100000; // 30 seconds
	DUMP=500;
	DYNSTEP=500;
	// SEED=161803398;
	SEED=123456;
	JUMP_MAX=0.01; //0.1
	Kbond = 1000;
	Kbend = 100000.0; //0.5;
	eps = 24.0;
	sigma = 1.0;
	arange = 1.5;
	vrange = 2.0;
	vthres = vrange-arange;
	Kdir = 100000.0;
	Rmem = 50;
	Rfac = 0.8;
	fix_Rmem = false;
	curvcos = false;
	ICempty = false;
	flat_mem = false;
	Min = false;
	MinGrow = false;
	MinRate = 0.0;
	fracmin = 0.2;
	varMin = 400.0;
	varMinDyn = false;
	nucleoid = false;
	nucleoidMin = false;
	nucleoidMinGrow = false;
	kmax = 0.9;
	saturate = false;
	satVal = 1.0;
	kon = 0.6;
	koff = 1.0;
	knuc = 0.2;
	sizer = false;
	mean_size = 20;
	timer = false;
	etimer = false;
	mean_time = 75;
	tangle = 0.0;
	lnuctest = 2;
	tstraight = false;
	bands = false;
	nbands = 0;
	olapINF = false;
	thermal_noise = false;
	PC19 = false;
	PC19time = 1500000;
	polang = 45;
	// Default control values
	verbose=false;
	long_verbose=false;
	test=false;
	relax=false;
	progress=false;

	// Read arguments
	int argcount;
	for (argcount = 1; argcount<argc; argcount++) {
		if (!strcmp(argv[argcount],"-p")) sprintf(output_dir, "%s/",argv[++argcount]);
		if (!strcmp(argv[argcount],"-L")) L=atoi(argv[++argcount]);
		if (!strcmp(argv[argcount],"-relsteps")) RELSTEPS=atoi(argv[++argcount]);
		if (!strcmp(argv[argcount],"-steps")) MCSTEPS=atoi(argv[++argcount]);
		if (!strcmp(argv[argcount],"-dump")) DUMP=atoi(argv[++argcount]);
		if (!strcmp(argv[argcount],"-dynstep")) DYNSTEP=atoi(argv[++argcount]);
		if (!strcmp(argv[argcount],"-seed")) SEED=atoi(argv[++argcount]);
		if (!strcmp(argv[argcount],"-JUMP_MAX")) JUMP_MAX=atof(argv[++argcount]);
		if (!strcmp(argv[argcount],"-olapINF")) olapINF=true;
		//
		if (!strcmp(argv[argcount],"-Kbond")) Kbond=atof(argv[++argcount]);
		if (!strcmp(argv[argcount],"-Kbend")) Kbend=atof(argv[++argcount]);
		if (!strcmp(argv[argcount],"-eps")) eps=atof(argv[++argcount]);
		if (!strcmp(argv[argcount],"-arange")) arange=atof(argv[++argcount]);
		if (!strcmp(argv[argcount],"-Kdir")) Kdir=atof(argv[++argcount]);
		if (!strcmp(argv[argcount],"-Rfac")) Rfac=atof(argv[++argcount]);
		if (!strcmp(argv[argcount],"-Rmem")) Rmem=atof(argv[++argcount]);
		if (!strcmp(argv[argcount],"-ICempty")) ICempty=true;
		if (!strcmp(argv[argcount],"-fix_Rmem")) fix_Rmem=true;
		if (!strcmp(argv[argcount],"-curvcos")) curvcos=true;
		if (!strcmp(argv[argcount],"-flat_mem")) flat_mem=true;
		if (!strcmp(argv[argcount],"-tstraight")) tstraight=true;
		if (!strcmp(argv[argcount],"-Min")) Min=true; //printf("Modulation is activated!\n");
		if (!strcmp(argv[argcount],"-MinGrow")) {MinGrow=true; MinRate=atof(argv[++argcount]);} // printf("Modulation growth is activated at rate %f\n", MinRate);
		if (!strcmp(argv[argcount],"-nucleoid")) {nucleoid=true; kmax=atof(argv[++argcount]);}
		if (!strcmp(argv[argcount],"-nucleoidMin")) {nucleoidMin=true; kmax=atof(argv[++argcount]);}
		if (!strcmp(argv[argcount],"-nucleoidMinGrow")) {nucleoidMinGrow=true; MinRate=atof(argv[++argcount]);}
		if (!strcmp(argv[argcount],"-fracmin")) {fracmin=atof(argv[++argcount]);}
		if (!strcmp(argv[argcount],"-varMin")) varMin=atof(argv[++argcount]);
		if (!strcmp(argv[argcount],"-saturate")) saturate=true;
		if (!strcmp(argv[argcount],"-satVal")) satVal=atof(argv[++argcount]);
		if (!strcmp(argv[argcount],"-varMinDyn")) varMinDyn=true;
		if (!strcmp(argv[argcount],"-kon")) kon=atof(argv[++argcount]);							//growth rate [monomers per dynamic step]
		if (!strcmp(argv[argcount],"-koff")) koff=atof(argv[++argcount]);
		if (!strcmp(argv[argcount],"-knuc")) knuc=atof(argv[++argcount]);
		if (!strcmp(argv[argcount],"-sizer")) {mean_size=atoi(argv[++argcount]); sizer=true;}
		if (!strcmp(argv[argcount],"-timer")) {mean_time=atoi(argv[++argcount]); timer=true;}
		if (!strcmp(argv[argcount],"-etimer")) {mean_time=atoi(argv[++argcount]); etimer=true;}	// hydrolisis time [dynamic steps]
		if (!strcmp(argv[argcount],"-tangle")) {tangle=atof(argv[++argcount]);}
		if (!strcmp(argv[argcount],"-polang")) {polang=atof(argv[++argcount]);}
		if (!strcmp(argv[argcount],"-lnuctest")) {lnuctest=atoi(argv[++argcount]);}
		if (!strcmp(argv[argcount],"-thermal_noise")) {thermal_noise = true;}
		if (!strcmp(argv[argcount],"-nothing")) {eps = 0.0; Kdir = 0.0; Min = false;}
		if (!strcmp(argv[argcount],"-bands")) {bands = true; nbands=atoi(argv[++argcount]);}
		if (!strcmp(argv[argcount],"-PC19")) {PC19 = true; PC19time=atoi(argv[++argcount]);}
		//
		if (!strcmp(argv[argcount],"-v")) verbose=true;
		if (!strcmp(argv[argcount],"-vv")) long_verbose=true;
		if (!strcmp(argv[argcount],"-test")) test=true;
		if (!strcmp(argv[argcount],"-relax")) relax=true;
		if (!strcmp(argv[argcount],"-SIF")) SIF=true;
		if (!strcmp(argv[argcount],"-progress")) progress=true;
	}
	// printf("Relaxation steps: %d", RELSTEPS);
	N = L*L;
	P = (int)(N/(double)2);
	if (test) {
		MCSTEPS = 1000000;
		eps = 0.0;
		// Kdir = 0.0;
		Min = false;
	}
	if (SIF) {
		MCSTEPS = 1000000;
	}
	if (relax) {MCSTEPS = 500;}
	PROGDUMP = (int)(MCSTEPS/(double)10);
	if (sizer) {koff = kon/(double)mean_size;}
	// if (timer) {koff = kon/(double)mean_time;}
	if (timer) {koff = 1.0/(double)mean_time;}
	if (!fix_Rmem) {
		Rmem = L/(2*M_PI);
	}
	Rfil = Rfac*Rmem;

	// Initialise & open files
	if (stat(output_dir, &st) == -1) {mkdir(output_dir, S_IRWXU);}
	printf("\nSimulation results stored in:\n%s\n", output_dir);
	sprintf(movie_file,"%s/output.xyz",output_dir);
	dumpf=fopen(movie_file, "w");
	fclose(dumpf);
	sprintf(energy_file,"%s/energy.txt",output_dir);
	energyf=fopen(energy_file, "w");
	fprintf(energyf, "step total bonds bend bind dir Nparts\n");
	fclose(energyf);
	sprintf(info_file,"%s/info.txt",output_dir);
	infof=fopen(info_file, "w");
	fprintf(infof, "\n----------------\nINFORMATION FILE\n----------------\n\n");
	fprintf(infof, "path: %s\n", output_dir);
	fprintf(infof, "L: %d\n", L);
	fprintf(infof, "relaxation steps: %d\n", RELSTEPS);
	fprintf(infof, "steps: %d\n", MCSTEPS);
	fprintf(infof, "dump: %d\n", DUMP);
	fprintf(infof, "dynstep: %d\n", DYNSTEP);
	fprintf(infof, "seed: %d\n", SEED);
	fprintf(infof, "JUMP_MAX: %f\n\n", JUMP_MAX);
	fprintf(infof, "Kbond: %f\n", Kbond);
	fprintf(infof, "Kbend: %f\n", Kbend);
	fprintf(infof, "epsilon: %f\n", eps);
	fprintf(infof, "arange: %f\n", arange);
	fprintf(infof, "Kdir: %f\n", Kdir);
	fprintf(infof, "Rmem: %f\n", Rmem);
	fprintf(infof, "Rfil: %f\n", Rfil);
	fprintf(infof, "Rfil/Rmem: %f\n", Rfac);
	if (ICempty) {fprintf(infof, "Initial conditions: EMPTY SYSTEM\n");}
	else  {fprintf(infof, "Initial conditions: RANDOM SYSTEM\n");}
	if (curvcos) {fprintf(infof, "Type of curvature sensing protocol: cosine\n");}
	else {fprintf(infof, "Type of curvature sensing protocol: precise\n\n");}
	if (Min) {
		fprintf(infof, "Modulation: True\n");
		if (MinGrow) {
			fprintf(infof, "Growing modulation amplitude: True\n");
			fprintf(infof, "\tModulation amplitude growth rate [fraction of max / tau_dyn]: %f\n\n", MinRate);
		}
		else {
			if (varMinDyn) {
				fprintf(infof, "Dynamic modulation profile: True\n");
				fprintf(infof, "\tvarMin: %f\n\n", varMin);
			 }
			else {fprintf(infof, "Dynamic modulation profile: False\n\n");}
		}
	}
	else {fprintf(infof, "Modulation: False\n\n");}
	fprintf(infof, "kon: %f\n", kon);
	if (nucleoid || nucleoidMin) {
		fprintf(infof, "kmax: %f\n", kmax);
	}
	fprintf(infof, "koff: %f\n", koff);
	fprintf(infof, "knuc: %f\n", knuc);
	if (sizer) {fprintf(infof, "nmean: %d\n", mean_size);}
	if (timer) {fprintf(infof, "thyd (linear distribution): %d\n", mean_time);}
	if (etimer) {fprintf(infof, "thyd (exponential distribution): %d\n", mean_time);}
	fprintf(infof, "\n");
	if (test) {
		fprintf(infof, "TEST SIMULATION\n");
		fprintf(infof, "test filament length: %d\n", lnuctest);
		fprintf(infof, "test filament theta angle: %f [degrees]\n", tangle);
	}
	if (SIF) {fprintf(infof, "SINGLE FILAMENT SIMULATION\n\n");}
	fprintf(infof, "\n");
	fclose(infof);
	sprintf(rates_file,"%s/rates.txt",output_dir);
	ratesf=fopen(rates_file, "w");
	fprintf(ratesf, "step full_trials full_rejects band_trials band_rejects\n");
	fclose(ratesf);
	sprintf(reject_file,"%s/rejects.txt",output_dir);
	rejectf=fopen(reject_file, "w");
	fprintf(rejectf, "step xy_trials xy_rejects reject_rate\n");
	fclose(rejectf);

	// Initialise variables
	srand48(SEED);
	step = 0;
	NpolTrials = 0;
	NpolTrialsBand = 0;
	NpolRejectedEV = 0;
	NpolRejectedEVBand = 0;

	xytrials = 0;
	xyrejects = 0;

	// Initialise the system
	if (test) {init_test_system();}
	else if (SIF) {init_test_system();}
	else if (bands) {init_bands_system();}
	else {init_system();}
	new_vlist();

	// Save IC
	compute_system_energy(step);
	dump_movie_frame(step);

	// Relax IC
	rMin = Min;
	reps = eps;
	// eps = 0.0;						//uncomment this line to switch off the attraction in the rel-steps too!
	Min = false;
	if (progress) {printf("\n");}
	while (step < RELSTEPS) {
		int fid, k, pid;
		double flipcoin;
		// Dynamic steps
		if (step%DYNSTEP == 0) {
			if (verbose) {printf("Dynamic step...\n");}
			// if (verbose) {printf("Reject rate: %f\n", (double)xyrejects/(double)xytrials);}
			for (fid = 0; fid < P; fid++) {
				k = (int)floor(drand48()*P);
				if (! filaments[k].exists) {continue;}
				flipcoin = drand48();
				if (flipcoin < 0.5) {
					if (saturate) {
						if (verbose) {printf("Polymerisation... current #: %d - would be #: %d - max #: %d\n", Nparts, Nparts+1, (int)(L*L*satVal));}
						if ((Nparts+1) <= (int)(L*L*satVal)) {
							polMC(k);
						}
					}
					else {
						polMC(k);
					}
					// polMC(k);
					depolMC(k);
				}
				else {
					depolMC(k);
					if (saturate) {
						if (verbose) {printf("Polymerisation... current #: %d - would be #: %d - max #: %d\n", Nparts, Nparts+1, (int)(L*L*satVal));}
						if ((Nparts+1) <= (int)(L*L*satVal)) {
							polMC(k);
						}
					}
					else {
						polMC(k);
					}
					// polMC(k);
				}
			}
			if (saturate) {
				if (verbose) {printf("Nucleation... current #: %d - would be #: %d - max #: %d\n", Nparts, Nparts+2, (int)(L*L*satVal));}
				if ((Nparts+2) <= (int)(L*L*satVal)) {
					nucMC();
				}
			}
			else {
				nucMC();
			}
			// nucMC();
		}
		// Translation relaxation steps
		new_vlist();
		for (pid = 0; pid < N; pid++) {
			k = (int)floor(drand48()*N);
			if (! particles[k].exists) {continue;}
			xyMC(k);
		}
		// Increase step
		step ++;
		// Dump & print
		if (step%DUMP == 0) {
			if (verbose) {printf("Dumping...\n");}
			compute_system_energy(step);
			dump_movie_frame(step);
			NpolTrials = 0;
			NpolRejectedEV = 0;
			NpolTrialsBand = 0;
			NpolRejectedEVBand = 0;
		}
	}
	eps = reps;
	Min = rMin;
	printf("Relaxation done!\n");

	// Evolve
	if (progress) {printf("\n");}
	while (step < MCSTEPS+RELSTEPS) {
		int fid, k, pid;
		double flipcoin;
		// Dynamic steps
		if (step%DYNSTEP == 0) {
			if (verbose) {printf("Dynamic step...\n");}
			// if (verbose) {printf("Reject rate: %f\n", (double)xyrejects/(double)xytrials);}
			if (varMinDyn) {
				if (step < MCSTEPS-MCSTEPS/4.0) {varMin = (L-(L-4.0)*step/(float)(MCSTEPS-MCSTEPS/4.0))*(L-(L-4.0)*step/(float)(MCSTEPS-MCSTEPS/4.0));}
				else {varMin = 16.0;}
				// printf("step = %d -- sqrt(varMin) = %f\n", step, sqrt(varMin));
			}
			for (fid = 0; fid < P; fid++) {
				k = (int)floor(drand48()*P);
				if (! filaments[k].exists) {continue;}
				flipcoin = drand48();
				if (flipcoin < 0.5) {
					if (saturate) {
						if (verbose) {printf("Polymerisation... current #: %d - would be #: %d - max #: %d\n", Nparts, Nparts+1, (int)(L*L*satVal));}
						if ((Nparts+1) <= (int)(L*L*satVal)) {
							polMC(k);
						}
					}
					else {
						polMC(k);
					}
					// polMC(k);
					depolMC(k);
				}
				else {
					depolMC(k);
					if (saturate) {
						if (verbose) {printf("Polymerisation... current #: %d - would be #: %d - max #: %d\n", Nparts, Nparts+1, (int)(L*L*satVal));}
						if ((Nparts+1) <= (int)(L*L*satVal)) {
							polMC(k);
						}
					}
					else {
						polMC(k);
					}
					// polMC(k);
				}
			}
			if (saturate) {
				if (verbose) {printf("Nucleation... current #: %d - would be #: %d - max #: %d\n", Nparts, Nparts+2, (int)(L*L*satVal));}
				if ((Nparts+2) <= (int)(L*L*satVal)) {
					nucMC();
				}
			}
			else {
				nucMC();
			}
			// nucMC();
		}
		// Translation relaxation steps
		new_vlist();
		if (long_verbose) {printf("Relaxation step...\n");}
		for (pid = 0; pid < N; pid++) {
			k = (int)floor(drand48()*N);
			if (! particles[k].exists) {continue;}
			xyMC(k);
		}
		// Increase step
		step ++;
		// Dump & print
		if (step%DUMP == 0) {
			if (verbose) {printf("Dumping...\n");}
			compute_system_energy(step);
			dump_movie_frame(step);
			ratesf=fopen(rates_file, "a");
			fprintf(ratesf, "%d %d %d %d %d\n", step, NpolTrials, NpolRejectedEV, NpolTrialsBand, NpolRejectedEVBand);
			fclose(ratesf);
			NpolTrials = 0;
			NpolRejectedEV = 0;
			NpolTrialsBand = 0;
			NpolRejectedEVBand = 0;
		}
		// Progress
		if (progress && ((step-RELSTEPS)%PROGDUMP == 0)) {printf("%d %% done\n", (step-RELSTEPS)/PROGDUMP*10);}
	}

	end_t = clock();
	total_t = ((double)end_t - (double)start_t) / (double)CLOCKS_PER_SEC;
	int hrs, mins;
	double secs;
	mins = (int)floor(total_t/60.0);
	secs = (total_t-mins*60.0);
	hrs = (int)floor((float)mins/60.0);
	mins = (int)(mins-hrs*60.0);
	printf("\nWall time: %d hrs %d mins %f secs\n",hrs,mins,secs);

	printf("\n");
	return 0;
}														// end of MAIN
/* ------- ---- ------- */
