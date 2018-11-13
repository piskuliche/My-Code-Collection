#include "main.h"

void pbc_dist(double& dr,double& L)
{
    dr = dr - L*nearbyint(dr/L);
}

void norm_comp(double& ex, double& ey, double& ez)
{
    double norm = sqrt(pow(ex,2) + pow(ey,2) + pow(ez,2));
    ex /= norm, ey /= norm, ez /= norm;
}

double P2(double& edote)
{
    double value = 0.5*(3*pow(edote,2) - 1.0);
    return value;
}



int main(int argc, char* argv[])
{
    /*Initializations*/
    int n_corr = 400, start_config = 0, end_config = 0, nevery = 20, nframes = 0;
    int bin_length = 0, nblocks = 0;
    int num_mols, natoms = 0, atoms_per_mol=0, type=0;
    int l_per_frame = 0, start_skip = 0, end_skip = 0;
    int msdindex = 0, c2index1 = 0, c2index2 = 0;
    double Lx = 0.0, Ly = 0.0, Lz = 0.0;
    double dump_freq = 0.0;
    vector<double> mass;
    vector<double> X,Y,Z;   // Holds masses
    vector<double> exo, eyo, ezo;
    vector<double> MSD, C1, C2; // Holds Mean Squared Displacement
    vector< vector<double>> bl_MSD,bl_C1, bl_C2;
    string molfile  = "molfile";
    string filename = "default.dat";
    string compound = "test";
    string buffer;
    string output_string;
    
    // Read Command Line Arguments
    if(argc < 8)
    {
        printf("Usage: diff.exe molecfile trajfile start_config end_config num_mols nblocks dumpfreq [nevery] [n_corr] [startskip] [endskip]\n");
        printf("Note: for optional arguments you should include all arguments that precede the optional arguments you need.\n I.E. if you need startskip, you must include nevery and n_corr.\n");
        printf("Defaults for Optional Arguments: nevery=20, n_corr = 400, startskip = 0, endskip = 0\n");
        printf("Note: you need not include startskip = 2 for an xyz file - the top two lines are already accounted for\n");
        exit(0);
    }
    for(int counter=1; counter<argc;counter++)
    {
        char argname[12][15] = {"exe","molfile","traj_file", "start_config", "end_config", "num_mols","nblocks", "dump_freq","nevery", "n_corr","start_skip", "end_skip"};
        printf("Argument %d: %s = %s\n", counter, argname[counter], argv[counter]);
    }
    molfile = argv[1], filename = argv[2], start_config = strtol(argv[3],nullptr, 0), end_config = strtol(argv[4],nullptr, 0), num_mols = strtol(argv[5],nullptr, 0), nblocks = strtol(argv[6],nullptr, 0), dump_freq = strtod(argv[7],nullptr);

    if(argc == 9)
    {
        nevery = strtol(argv[8], nullptr, 0);
    }
    else if (argc == 10)
    {
        nevery = strtol(argv[8], nullptr, 0);
        n_corr = strtol(argv[9], nullptr, 0);
    }
    else if(argc == 11)
    {
        nevery = strtol(argv[8], nullptr, 0);
        n_corr = strtol(argv[9], nullptr, 0);
        start_skip = strtol(argv[10], nullptr, 0);
    }
    else
    {
        nevery = strtol(argv[8], nullptr, 0);
        n_corr = strtol(argv[9], nullptr, 0);
        start_skip = strtol(argv[10],nullptr,0);
        end_skip   = strtol(argv[11],nullptr,0);
    }

    //read molfile
    ifstream mfile(molfile);
    if(mfile.good() == false)
    {
        printf("Error: molfile does not exist\n");
        exit(1);
    }
    mfile >> atoms_per_mol;
    mfile >> msdindex;
    mfile >> c2index1 >> c2index2;
    for(int i=0;i<atoms_per_mol;i++)
    {
        int m = 0;
        mfile >> m;
        mass.push_back(m);
    }
    mfile.close();

    //read box dimensions
    ifstream boxfile("box.info");
    if(boxfile.good() == false)
    {
        printf("Error: boxfile does not exist\n");
        exit(1);
    }
    double xlo = 0.0, xhi = 0.0, ylo = 0.0, yhi = 0.0, zlo = 0.0, zhi = 0.0;
    boxfile >> xlo >> xhi;
    boxfile >> ylo >> yhi;
    boxfile >> zlo >> zhi;
    Lx = xhi - xlo;
    Ly = yhi - ylo;
    Lz = zhi - zlo;
    boxfile.close();


    //read trajectory
    ifstream tin(filename);
    if(tin.good() == false)
    {
        printf("Error: trajfile does not exist\n");
        exit(1);
    }
    l_per_frame = 2+atoms_per_mol*num_mols+start_skip+end_skip;
    int count = 0;
    printf("Starting to Read Trajectory\n");
    while(getline(tin,buffer)){count++;}
    nframes = count/l_per_frame;
    printf("There are %d lines, %d frames and %d lines per frame\n", count, nframes,l_per_frame);
    tin.close();
    tin.open(filename);
    double x=0, y=0, z=0;
    for(int f=0; f<nframes; f++)
    {
        if((f+1) % 1000 == 0) {printf("Reading Frame %d of %d\n", f+1,nframes);}
        getline(tin,buffer);
        getline(tin,buffer);
        if(start_skip != 0)
        {
            for(int i=0; i<start_skip;i++)
            {
                getline(tin,buffer);
            }
        }
        for(int m=0; m<num_mols; m++)
        {
            for(int a=0; a<atoms_per_mol; a++)
            {
                getline(tin,buffer);
                istringstream stream(buffer);
                string tmptype, tmpx, tmpy, tmpz;
                stream >> tmptype >> tmpx >> tmpy >> tmpz;

                x = stod(tmpx);
                y = stod(tmpy);
                z = stod(tmpz);
                X.push_back(x);
                Y.push_back(y);
                Z.push_back(z);
            }
        }
        if(end_skip != 0)
        {
            for(int i=0; i<end_skip;i++)
            {
                getline(tin,buffer);
            }
        }
    }
    printf("Trajectory read completed.\n");
    // Create MSD, C1, C2 Vectors
    for(int i=0; i<n_corr; i++)
    {
        MSD.push_back(0.0);
        C1.push_back(0.0);
        C2.push_back(0.0);
    }
    for(int i=0; i<num_mols; i++)
    {
        exo.push_back(0.0);
        eyo.push_back(0.0);
        ezo.push_back(0.0);
    }


    printf("Calculating Correlation Functions\n");
    for(int block = 0; block<nblocks;block++)
    {
        int sep_bl   = floor((end_config-start_config-n_corr)/(double) nblocks);
        int start_bl = start_config+block*sep_bl;
        int end_bl   = start_config+(block+1)*sep_bl;
        printf("-------------\nBLOCK: %d\nSTART: %d\nEND: %d\n-------------\n", block+1, start_bl, end_bl);
        // For loop over time
        double dx = 0.0, dy = 0.0, dz = 0.0, dr2 = 0.0, dr2o = 0.0;
        double ex = 0.0, ey = 0.0, ez = 0.0;
        //Zero MSD C1 C2 Vectors
        for(int i=0; i<n_corr; i++)
        {
            MSD[i] = 0.0;
            C1[i]  = 0.0;
            C2[i]  = 0.0;
        }
        //Loop over time origins
        int count = 0;
        for(int to = start_bl; to<end_bl; to+=nevery)
        {
            count += 1;
            // Loop over times
            dr2o = 0.0;
            for(int t = to+1; t < (to+n_corr); t++)
            {
                dr2 = 0.0;
                // Loop over molecules
                for(int j=0; j<num_mols; j++)
                {
                    // Set indices for the first frame (fo), previous frame(ao) and the current frame(a)
                    int fo = to*num_mols*atoms_per_mol + j*atoms_per_mol;
                    int ao = (t-1)*num_mols*atoms_per_mol+j*atoms_per_mol;
                    int a  = t*num_mols*atoms_per_mol + j*atoms_per_mol;
                    // Calculates distance r(t)-r(t-1)
                    dx = X[a]-X[ao];
                    dy = Y[a]-Y[ao];
                    dz = Z[a]-Z[ao];
                    // Periodic boundary conditions (minimum image distance)
                    pbc_dist(dx,Lx);
                    pbc_dist(dy,Ly);
                    pbc_dist(dz,Lz);
                    // Increment dr2 by the sq. of the distance as calculated
                    dr2 += pow(dx,2.) + pow(dy,2.) + pow(dz,2.);
                    // This section will include the calculation for the C2 correlation function
                    if(atoms_per_mol > 1)
                    {
                        // Calculates bond vector between atoms for the original time
                        exo[j] = X[fo+c2index1] - X[fo+c2index2];
                        eyo[j] = Y[fo+c2index1] - Y[fo+c2index2];
                        ezo[j] = Z[fo+c2index1] - Z[fo+c2index2];
                        // Periodic boundary conditions, minimum image distance for the original vector
                        pbc_dist(exo[j],Lx);
                        pbc_dist(eyo[j],Ly);
                        pbc_dist(ezo[j],Lz);
                        // Normalizes components of the e(t) vector
                        norm_comp(exo[j], eyo[j], ezo[j]);
                        // Calculates the new bond vector
                        ex = X[a+c2index1] - X[a+c2index2];
                        ey = Y[a+c2index1] - Y[a+c2index2];
                        ez = Z[a+c2index1] - Z[a+c2index2];
                        // Minimum image distance for the vector
                        pbc_dist(ex,Lx);
                        pbc_dist(ey,Ly);
                        pbc_dist(ez,Lz);
                        // Calculates the magnitude of the vector and normalizes
                        norm_comp(ex,ey,ez);
                        // Calculate e(t)•e(0) and then increments C1(t) and C2(t)
                        double edot = ex*exo[j]+ey*eyo[j]+ez*ezo[j];
                        C1[t-to] += edot;
                        C2[t-to] += P2(edot);
                    }
                    
                }
                // Divide by number of mols
                dr2    /= (double) num_mols;
                // Increment MSD(t)
                MSD[t-to]  += dr2o+dr2;
                // Divide C1(t) and C2(t) by number of mols
                if(atoms_per_mol > 1)
                {
                    C1[t-to]  /= num_mols;
                    C2[t-to]  /= num_mols;
                }
                // increment dr2o by dr2
                dr2o += dr2;
            }
        }
        cout << "There are " << count << " total tos per block" << endl;
        // This piece of code divides the total by the number of time origins per block -> average
        for(int i = 0; i < n_corr; i++)
        {
            MSD[i] /= ((sep_bl/(double) nevery));
            if(atoms_per_mol > 1)
            {
                C1[i]  /= (sep_bl/(double) nevery);
                C2[i]  /= (sep_bl/(double) nevery);
            }
        }
        // Stores the MSD, C1, C2 vectors in an array.
        bl_MSD.push_back(MSD);
        if(atoms_per_mol > 1)
        {
            bl_C1.push_back(C1);
            bl_C2.push_back(C2);
        }
    }
    // Zeros the MSD, C1, C2 vectors
    for(int i = 0; i<n_corr; i++)
    {
        MSD[i] = 0;
        C1[i]  = 0;
        C2[i]  = 0;
    }
    // Calculates the total correlation functions
    for(int block = 0; block<nblocks;block++)
    {
        for(int i = 0; i<n_corr; i++)
        {
            MSD[i] += bl_MSD[block][i];
            if(atoms_per_mol > 1)
            {
                C1[i] += bl_C1[block][i];
                C2[i] += bl_C2[block][i];
            }
        }
    }
    // Outputs the files to a file
    ofstream msdout("msd.dat");
    for(int i =0; i<n_corr; i++)
    {
        MSD[i] /= nblocks;
        double time = i*dump_freq;
        msdout << time << " " <<MSD[i] << " ";
        for(int block = 0; block<nblocks; block++)
        {
            msdout << bl_MSD[block][i] << " ";
        }
        msdout << "\n";
    }
}
    



