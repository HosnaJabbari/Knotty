
// a simple driver for the CCJ
// Hosna: base pair maximization version, Feb 9, 2014

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>

// include the simfold header files
#include "simfold.h"
#include "externs.h"
#include "h_globals.h"
#include "h_externs.h"
#include "constants.h"
#include "params.h"

#include "cmd_line_options.h"

#include "CCJ.h"

int main (int argc, char *argv[])
{
    char sequence[MAXSLEN];
    char structure[MAXSLEN];
   // char restricted[MAXSLEN];
    double energy;
    char structures[MAXSUBSTR][MAXSLEN];
    double energies[MAXSUBSTR];

    bool cmd_line_error = false;
    bool w = false;

    if (argc < 2)
        cmd_line_error = true;
    else
    {
        strcpy (sequence, argv[1]);
	// important that this is before set_shape_file
        cmd_line_options.set_sequence_length(strlen(sequence));

        if (argc > 2) {
            for (int i = 2; i < argc; ++i) {
                char * arg = argv[i];

                if (!strcmp(arg, "-ns"))
                    cmd_line_options.set_use_sparse(false);
                else
                if (!strcmp(arg, "-ngc"))
                    cmd_line_options.set_use_garbage_collection(false);
                else
                if (!strcmp(arg, "-w"))
                    w = true;
                else
                if (!strcmp(arg, "-pta"))
                    cmd_line_options.set_print_trace_arrow_info(1);
                else
                if (!strcmp(arg, "-pta-v"))
                    cmd_line_options.set_print_trace_arrow_info(2);
                else
                if (!strcmp(arg, "-pcl"))
                    cmd_line_options.set_print_candidate_list_info(1);
                else
                if (!strcmp(arg, "-pcl-v"))
                    cmd_line_options.set_print_candidate_list_info(2);
                else
                if (!strncmp(arg, "-shape=", 7))
                    cmd_line_options.set_shape_file(std::string(arg));
                else
                if (!strncmp(arg, "-b=",3)) {
                    std::string str = std::string(arg);
                    str = str.substr(3,str.length()-2);
                    if (shape.is_number(str))
                        shape.set_b(atof(str.c_str()));
                    else
                        cmd_line_error = true;
                }
                else
                if (!strncmp(arg, "-m=",3)) {
                    std::string str = std::string(arg);
                    str = str.substr(3,str.length()-2);
                    if (shape.is_number(str))
                        shape.set_m(atof(str.c_str()));
                    else
                        cmd_line_error = true;
                }
                else
                    cmd_line_error = true;
            }
        }
    }


    if (cmd_line_error) {
        printf ("Usage: %s <sequence> <arguments>\n", argv[0]);
        printf ("Valid agruments include: \n");
        printf ("-ns to use non-sparse or \"Modifed CCJ\" version\n");
        printf ("-ngc to not use garbage collection \n \n");

        //printf ("-shape=\"filename\" to specify a file for shape data\n");
        //printf ("-b=number to specify an intercept for the shape data (default is %f)\n",shape.b());
        //printf ("-m=number to specify a slope for the shape data (default is %f)\n\n",shape.m());

        printf ("-w to print only the result and energy\n");
        printf ("-pta to print information on the number of trace arrows\n");
        printf ("-pta-v to print verbose trace arrow information\n");
        printf ("-pcl to print information on the candidate lists\n");
        printf ("-pcl-v to print verbose candidate list information\n");
        printf ("Example: %s GCAACGAUGACAUACAUCGCUAGUCGACGC \n", argv[0]);
        return -1;
    }
    cmd_line_options.set_done();

    //strcpy (restricted, argv[2]);

    // Before calling any function in the library, you have to initialize config_file, dna_or_rna, temperature
    //     and to call the function init_data, which loads the thermodynamic parameters into memory

    // configuration file, the path should be relative to the location of this executable
    char config_file[200];
    strcpy (config_file, SIMFOLD_HOME "/params/multirnafold.conf");

    // what to fold: RNA or DNA
    int dna_or_rna;
    dna_or_rna = RNA;

    // temperature: any integer or real number between 0 and 100
    // represents degrees Celsius
    double temperature = 37.0;


    // initialize the thermodynamic parameters
    // call init_data only once for the same dna_or_rna and same temperature
    // if one of them changes, call init_data again  
    init_data (argv[0], config_file, dna_or_rna, temperature);

	// Hosna, July 18, 2012
	// In simfold we have the following for RNA && temp=37
	fill_data_structures_with_new_parameters (SIMFOLD_HOME "/params/turner_parameters_fm363_constrdangles.txt");

	// Hosna, July 25, 2012
	// in HotKnots and ComputeEnergy package the most up-to-date parameters set is DP09.txt
	// so we add it here
	fill_data_structures_with_new_parameters (SIMFOLD_HOME "/params/parameters_DP09.txt");

	energy = ccj(sequence, structure);

    if (w) {
        printf ("%s %.2lf\n", structure, energy);
    } else {
        printf ("Seq: %s\n", sequence);
        printf ("RES: %s  %.2lf\n", structure, energy);
    }

    return 0;
}



