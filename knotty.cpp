
// a simple driver for Knotty
// Hosna: base pair maximization version, Feb 9, 2014

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>

// include the simfold header files
#include "simfold.h"
#include "externs.h"
#include "h_externs.h"
#include "constants.h"
#include "params.h"

#include "cmd_line_options.h"

#include "knotty.h"

int main (int argc, char *argv[])
{
    char sequence[MAXSLEN];
    char structure[MAXSLEN];
    double energy;

    bool cmd_line_error = false;
    bool w = false;

    // reading arguments
    if (argc < 2)
        cmd_line_error = true; // need at least 2 arguments (knotty and sequence)
    else
    {
	// get sequence
        strcpy (sequence, argv[1]);
	// important that this is before set_shape_file
        cmd_line_options.set_sequence_length(strlen(sequence));

	// addtional arguments
        if (argc > 2) {
            for (int i = 2; i < argc; ++i) {
                char * arg = argv[i];

		// -ns uses non-sparse version
                if (!strcmp(arg, "-ns"))
                    cmd_line_options.set_use_sparse(false);
                else

		// ** The following arguments are primarily debug options **
		// -ngc does not use garbage collection (sparse version only)
                if (!strcmp(arg, "-ngc"))
                    cmd_line_options.set_use_garbage_collection(false);
                else
		// -w for web printing (used for web server)
                if (!strcmp(arg, "-w"))
                    w = true;
                else
		// -pta prints extra trace arrow info
                if (!strcmp(arg, "-pta"))
                    cmd_line_options.set_print_trace_arrow_info(1);
                else
		// -pta-v prints even more verbose trace arrow info
                if (!strcmp(arg, "-pta-v"))
                    cmd_line_options.set_print_trace_arrow_info(2);
                else
		// -pcl print extra candidate list info
                if (!strcmp(arg, "-pcl"))
                    cmd_line_options.set_print_candidate_list_info(1);
                else
		// -pcl-v prints even more verbose candidate list info
                if (!strcmp(arg, "-pcl-v"))
                    cmd_line_options.set_print_candidate_list_info(2);
                else
                    cmd_line_error = true; // argument is invalid
            }
        }
    }


    if (cmd_line_error) {
        printf ("Usage: %s <sequence> <arguments>\n\n", argv[0]);
        printf ("Arguments are largely for debugging.\n");
        printf ("Valid arguments include: \n");
        printf ("-ns to use non-sparse version\n");
        printf ("-ngc to not use garbage collection \n \n");

        printf ("-w to print only the result and energy\n");
        printf ("-pta to print information on the number of trace arrows\n");
        printf ("-pta-v to print verbose trace arrow information\n");
        printf ("-pcl to print information on the candidate lists\n");
        printf ("-pcl-v to print verbose candidate list information\n\n");

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

	energy = knotty(sequence, structure);

    if (w) {
        printf ("%s %.2lf\n", structure, energy);
    } else {
        printf ("Seq: %s\n", sequence);
        printf ("RES: %s  %.2lf\n", structure, energy);
    }

    return 0;
}
