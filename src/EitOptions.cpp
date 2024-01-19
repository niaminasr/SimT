#include "ParEITConfig.hpp"

void
eitGetOptions(int argc, char* argv[], float* radius, int* elec, unsigned int* dims, int* solver)
{
    int opt;

    /**
     * default values
     *
     */
    {
        *radius = 1;
        *elec = 4;
        *dims = 10;
        *solver = 1;
    }

    while ((opt = getopt(argc, argv, "r:n:d:s")) != -1) {
        switch (opt) {
            case 'r':
                *radius = atof(optarg);
                break;
            case 'n':
                *elec = atoi(optarg);
                break;
            case 'd':
                *dims = atoi(optarg);
                break;
            case 's':
                *solver = atoi(optarg);
                break;
            default:
                fprintf(stderr, "Usage: -r radius of the form\n"\
                                "       -n number of electrodes\n"\
                                "       -d grid dimensions\n"\
                                "       -s chosen solver\n");
                exit(EXIT_FAILURE);
        }
    }

}