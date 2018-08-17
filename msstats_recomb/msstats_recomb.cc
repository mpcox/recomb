/*
 msstats_recomb - read data from ms via stdin, calculate n-tuples and recombination summary statistics
 
 Original code version by Kevin Thornton
 2002
 University of California Irvine
 
 Modified extensively by Murray Cox <m.p.cox@massey.ac.nz>
 March 2013
 Massey University, New Zealand
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 
 Compile: g++ -O3 -o msstats_recomb msstats_recomb.cc -lsequence
 */

#include <iostream>
#include <vector>
#include <cstdio>
#include <numeric>
#include <utility>
#include <math.h>
#include <Sequence/SimParams.hpp>
#include <Sequence/SimData.hpp>
#include <Sequence/PolySIM.hpp>
#include <Sequence/SeqConstants.hpp>
#include <Sequence/PolyTableFunctions.hpp>
#include <Sequence/Recombination.hpp>

using namespace std;
using namespace Sequence;

/*
 Returns a simple bound on the minumum number of recombination events, calculated according to equation 4 of Myers & Griffiths paper
 */
unsigned Rm_MG( const Sequence::SimData & matrix, unsigned segsites, unsigned nhaps ){
    
    bool anc = ( std::find(matrix.begin(),matrix.end(),std::string(segsites,'0')) != matrix.end() );
    return (nhaps > segsites) ? ( (anc) ? nhaps-segsites-1 : nhaps-segsites ) : 0;
}

/* MPC implementation of min_d */
signed int min_d(const SimData & matrix){
	
	unsigned int n = matrix.size();
	unsigned int sites = matrix.numsites();
	
	signed int current_min = sites;
	
	if( sites == 0 ){
		return current_min = -1;
	}
	
	for ( unsigned i = 0; i < n; ++i ){				//iterate over all pairs individuals
		for ( unsigned j = i + 1; j < n; ++j ){
			
			signed int this_dist = 0;
			
			for ( unsigned a = 0; a < sites; ++a ){	//iterate over SNPs
				
				if( matrix[i][a] != matrix[j][a] ){
					++this_dist;
				}
			}
			
			if( this_dist < current_min ){
				current_min = this_dist;
			}
			
		}
	}
    
	return current_min;
}

int main(int argc, char *argv[]){
    
    // global variables
    std::vector<int> config;
    unsigned int n_runs, n_inds, n_sites;
    int mincount = 1;
    
    // command line variables
    unsigned int n_subsamples = 0, n_permutations = 1;
    signed int flag_c;
	
    // command line usage
    char* usage = "usage: msstats_recomb [-q sequences_in_subsample -p number_of_permutations -h]\n";
	
    // read values from command line
    for(flag_c = 1; flag_c < argc; ++flag_c ) {
          if(argv[flag_c][0] == '-') {
                
             switch( argv[flag_c][1] ) {
                        
                  // set number of individuals in subsample
                 case 'q':
                       ++flag_c;
                     if(argv[flag_c] == NULL){
                           fprintf(stderr, "error: no values following -q flag\n");
                          exit(1);
                     }
                     n_subsamples = atoi(argv[flag_c]);
                     break;
                        
                 // set number of permutations
                 case 'p':
                    ++flag_c;
                      if(argv[flag_c] == NULL){
                           fprintf(stderr, "error: no values following -p flag\n");
                           exit(1);
                     }
                      n_permutations = atoi(argv[flag_c]);
                      break;
                        
                 // print help information
                 case 'h':
                    fprintf(stderr, "%s", usage);
					exit(1);
                        
                      default:
                     fprintf(stderr, "error: option %s unknown\n", argv[flag_c]);
                     exit(1);
			}
				
		}else ++flag_c;
	}
	
    SimParams p;
    p.fromfile(stdin);
    
    n_runs = p.runs();     // number of ms replicates
    n_inds = p.totsam();   // the total sample size (number of sequences)
    
    if( n_subsamples == 0 ){
        n_subsamples = n_inds;
    }
    
    SimData d;
    
    #if __GNUG__ && __GNUC__ >= 3
    std::ios_base::sync_with_stdio(true);
    #endif
    
    // input request test
    if( n_subsamples > n_inds ){
         std::cerr << "error: requested subsample size greater than number of simulated sequences" << endl;
         exit(1);
     }
        
    // print header line
    std::cout << "S\t"
    << "min_d\t"
    << "thetaW\t"
    << "ThetaPi\t"
    << "Rmin\t"
    << "rmmg\t"
    << "nhaps\t"
    << "hapdiv\t"
    << "wallsb\t"
    << "wallsq\t"
    << "hudsonsc\t"
    << "zns\t"
    << endl;
    
    int rv;
    while( (rv=d.fromfile(stdin)) != EOF ){
        
        n_sites = d.numsites();
        
        // iteration loop
        for( unsigned int i = 0; i < n_permutations; ++i ){
            
            // shuffle rows in original data object
            random_shuffle(d.begin(),d.end());
            
            // make new, reduced data object
            SimData d2;
            d2.assign(&*d.pbegin(),d.numsites(),&d[0],n_subsamples);
            
            // remove sites lacking mutations
            RemoveInvariantColumns(&d2);
            
            // call summaries
            PolySIM P(&d2);
            
            // write summaries
            std::cout << P.NumPoly() << "\t";
            
            signed int mind = min_d(d2);
            if( mind < 0 ){
                std::cout << "nan" << '\t';
            }else{
                std::cout << mind << '\t';
            }
            
            std::cout << P.ThetaW() << "\t"
            << P.ThetaPi() << "\t";
            unsigned rm = P.Minrec(), nhaps = P.DandVK();
            if( rm != SEQMAXUNSIGNED ){
                std::cout << rm << "\t";
            }else{
                std::cout << "NAN" << "\t";
            }
            
            std::cout << Rm_MG(d2,P.NumPoly(),nhaps) << "\t"
            << nhaps << "\t"
            << P.DandVH() << "\t"
            << P.WallsB() << "\t"
            << P.WallsQ() << "\t"
            << P.HudsonsC() << "\t";
            if(d2.numsites() > 1){
            
                unsigned site1=0,site2=1;
                vector<double> LDSTATS(6);
                bool allskipped=true;
                double zns=0.;
                unsigned npairsLD=0;
                while( Recombination::Disequilibrium(&d2,LDSTATS,&site1,&site2,false,0,mincount)){
                    if( !LDSTATS[LDSTATS.size()-1] ){  //if site pair NOT skipped!
                    	
                        allskipped=false;
                        zns += LDSTATS[2];
                        ++npairsLD;
                    }
                }
                zns /= double(npairsLD);
                if( allskipped ){
                    cout << "nan" << endl;
                }else{
                    cout << zns << endl;
                }
            }
            else{
                cout << "nan" << endl;
            }
              
        } 
    }
}
