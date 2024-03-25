
// � A J S Hamilton 2001
// � A Yang 2009 - modified for histogram system

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "vector.h"
#include "globals.h"
// Parse arguments.

const char *optstr = "b:t:n:e:q:rv";

// Global parameters
vec1 gp_eps; // epsilon
vec1 gp_dt; // base timestep
int gp_burst; // burst steps per frame
int gp_ntgt; // target number of particles
// Termination conditions
vec1 gp_term_tc; // crossing time
unsigned int gf_verbose; // verbose mode
unsigned int gp_term_ns; // max steps
unsigned int gp_term_mode; // termination mode - 0: manual, 1: auto, bit 2: set, bit 3: condition

void parse_args(int argc, char *argv[]) {
	char opt;
	char tc;
	vec1 tmpn;
	int iscan;

	gp_eps = 0.0001;
	gp_dt = 0.01;
	gp_ntgt = 1024;
	gp_burst = 10;
	gp_term_mode = 0;
	gf_verbose = 0;

	/* turn off getopt complaints */
	opterr = 0;
	
	/* parse arguments */
	while ((opt = getopt(argc, argv, optstr)) != -1) {
		switch (opt) {
			case 'v': // verbose mode
				gf_verbose = 1;
				break;
			case 'r': // batch mode
				gp_term_mode = gp_term_mode | 1;
				break;
			case 'q': // termination conditions
#ifdef PREC_DOUBLE
				iscan = sscanf(optarg, "%c%lf", &tc, &tmpn);
#else
				iscan = sscanf(optarg, "%c%f", &tc, &tmpn);
#endif
				if (iscan != 2) {
					fprintf(stderr, "-%c%s: expecting [n|t][number]\n", opt, optarg);
					die(1);
				}
				if (tc == 'n') {
					gp_term_ns = tmpn;
					gp_term_mode = gp_term_mode | (3 << 1);
				}
				else if (tc == 't') {
					gp_term_tc = tmpn;
					gp_term_mode = gp_term_mode | (1 << 1);
					printf("run time interval: %f\n", gp_term_tc);
				}
				else {
					printf("Unrecognised termination mode: %c\n", tc);
					die(1);
				}
				break;
			case 'b': // burst count
				iscan = sscanf(optarg, "%d", &gp_burst);
				if (iscan != 1) {
					fprintf(stderr, "-%c%s: expecting number\n", opt, optarg);
					die(1);
				}
				break;
			case 't': // timestep
#ifdef PREC_DOUBLE
				iscan = sscanf(optarg, "%lf", &gp_dt);
#else
				iscan = sscanf(optarg, "%f", &gp_dt);
#endif
				if (iscan != 1) {
					fprintf(stderr, "-%c%s: expecting number\n", opt, optarg);
					die(1);
				}
				break;
			case 'n': // galaxies
				iscan = sscanf(optarg, "%d", &gp_ntgt);
				if (iscan != 1) {
					fprintf(stderr, "-%c%s: expecting number\n", opt, optarg);
					die(1);
				}
				break;
			case 'e': // epsilon
#ifdef PREC_DOUBLE
				iscan = sscanf(optarg, "%lf", &gp_eps);
#else
				iscan = sscanf(optarg, "%f", &gp_eps);
#endif
				if (iscan != 1) {
					fprintf(stderr, "-%c%s: expecting number\n", opt, optarg);
					die(1);
				}
				break;
			case ':':
			case '?':
			default:
				break;
		}
	}
}
