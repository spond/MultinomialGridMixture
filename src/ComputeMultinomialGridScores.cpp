
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <vector>
#include <time.h>
#include <float.h>
#include <math.h>

using namespace std;

vector<double> log_factorials;
vector<double> multinomial_coefficients_for_sites; 

double compute_log_factorial (unsigned long N) {

  unsigned long already_done = log_factorials.size();

  if (N >= already_done) {
    if (already_done == 0) {
      log_factorials.push_back(0.0);
      already_done++;
    } 
    
    double current_factorial = log_factorials[already_done-1],
           float_N = already_done;

    for (unsigned long counter = already_done; counter <= N; counter ++, float_N+=1.) {
      current_factorial += log (float_N);
      log_factorials.push_back(current_factorial);
    }
  }

  return log_factorials[N];

}

double compute_multinomial_coefficient (const double countA, const double countC, const double countG, const double countT) {
  return compute_log_factorial (countA + countC + countG + countT) - 
         compute_log_factorial (countA)- 
         compute_log_factorial (countC)- 
         compute_log_factorial (countG)- 
         compute_log_factorial (countT); 
}


static char Usage[] = "ngs_cmgs multimomial_grid_file [see README for a format specification]"
                      "\n\t(input data to stdin as "
                      "\n\tcount of A at site 1, count of C at site 1, count of G at site 1, count of T as site 1"
                      "\n\tcount of A at site 2, count of C at site 2, count of G at site 2, count of T as site 2"
                      "\n\t...)";
                      

//---------------------------------------------------------------

int main (int argc, const char * argv[])
{
    if (argc != 2) {
        cerr << "Usage is `" << Usage << "'." << endl;
        return 1;
    }
    
    ifstream rate_file(argv[1]);

    if (rate_file.fail())
    {
        cerr << "Cannot open file `" << argv[1] << "'." << endl;
        return 1;
    }

    
    istream_iterator<unsigned long> end_of_input,
                                    number_eater (cin);

    rate_file.close();                                  
                     
    vector<double> data,
                   rates;
                     
    rate_file.seekg (0, rate_file.end);
    cerr << rate_file.tellg() << endl;
    rate_file.seekg (0, rate_file.beg);
    cerr << rate_file.get() << endl;
 
    long numbers_read = 0;
    double last3 [3];
                     
    while (rate_file >> last3[numbers_read]) {
         rates.push_back(last3[numbers_read]);
         numbers_read ++; 
         cerr << last3 [0] << " " << numbers_read << endl;
         if (numbers_read == 3) {
          double sum = last3[0] + last3[1] + last3 [2];
          if (last3[0] <= 0.0 || last3[1] <= 0.0 || last3[2] <= 0.0 || sum >= 1.) {
            cerr << "Invalid multinomial rate specificaion " << last3[0] << " " << last3[1] << " " << last3[2] << endl;
            return 1;
          }
          numbers_read = 0;
         } 
    }
    
    cerr << rate_file.eof() << endl;
    cerr << rate_file.fail() << endl;
    cerr << rate_file.bad() << endl;
    

    while (number_eater != end_of_input) {
         data.insert (data.end(), *number_eater++);
    }
    
    if (data.size() % 4 != 0 && data.size () < 4) {
        cerr << "Expected an even number of values (at least 2" << endl;
        return 1;
    }
    
    unsigned long data_points = data.size()/4;
    
    // populate the vector of multinomial (log space) coefficients
    
    for (unsigned long site = 0; site < data_points; site ++) {
      multinomial_coefficients_for_sites.push_back (
        compute_multinomial_coefficient (
          data[site<<2],
          data[(site<<2)+1],
          data[(site<<2)+2],
          data[(site<<2)+3]
        )
      );
           
      cout << site << ":" << multinomial_coefficients_for_sites[multinomial_coefficients_for_sites.size()-1] << endl;
    }
    
    
    
    /*srand (time(NULL)); rand ();

    double best_aic = DBL_MAX;
    
    double   * coverage  = new double [data_points],
             * majority  = new double [data_points];

    for (unsigned long i = 0; i < data_points; i++) {
        if (data.at (2*i) < data.at (2*i+1)) {
            cerr << "Coverage values must not be less than majority counts" << endl;
            return 1;
        }
        coverage[i] = data.at (2*i);
        majority[i] = data.at (2*i+1);
    }
    
    bool go_on = true;
    
    double *previous_rates = NULL,
           *previous_weights = NULL;
           
    cout << "{" << endl;       
    
    while (go_on) {
    
        double * rates       =  new double [rate_class_count],
               * weights     =  new double [rate_class_count],
               * try_rates   =  new double [rate_class_count],
               * try_weights =  new double [rate_class_count];
                       
        //cerr << "Fitting " << data_points << " observations to a binomial mixture with " << rate_class_count << " components" << endl;
        
        random_starting_values (rates, weights, rate_class_count,previous_rates,previous_weights);
        double final_logL = runEM(coverage, majority, data_points, rates, weights, rate_class_count);
        
        for (unsigned long repl = 1; repl < 50*rate_class_count; repl ++) {
             if (repl < 10) {
                random_starting_values (try_rates, try_weights, rate_class_count, previous_rates,previous_weights);
             } else {
                random_starting_values (try_rates, try_weights, rate_class_count, NULL, NULL);                
             }
             double try_logL = runEM(coverage, majority, data_points, try_rates, try_weights, rate_class_count);
             if (try_logL > final_logL) {
                final_logL = try_logL;
                for (unsigned long comp = 0;  comp < rate_class_count; comp += 1) {
                    rates [comp]  = try_rates[comp];
                    weights[comp] = try_weights[comp];
                }
              }
        }
        
        
        double parameter_count = 2*rate_class_count-1.,
               c_aic           = 2*(parameter_count*data_points/(data_points - parameter_count -1) - final_logL);
        
        if (c_aic < best_aic) {
            cout << (rate_class_count > 1? ',' : ' ') << endl << '"' << rate_class_count << '"' << " : ";
            report_model ( final_logL,  c_aic,  rates,  weights, rate_class_count);
            if (previous_rates)   delete [] previous_rates;
            if (previous_weights) delete [] previous_weights;
            previous_rates   = rates;
            previous_weights = weights;
            best_aic = c_aic;
        } else {
            go_on = false;
        }
        
        if (go_on == false && rate_class_count > 2) {
        
            double  * pij          = new double [data_points*rate_class_count],
                    min_rate       = 1.,
                    min_rate_index = 0;
            
            rate_class_count --;
            for (unsigned long rclass = 0; rclass < rate_class_count; rclass ++) {
                    if (previous_rates[rclass] < min_rate) {
                        min_rate = previous_rates[rclass];
                        min_rate_index = rclass;
                    }
            }
            
            cout << ",\n\"background_rate\" : " << min_rate;

            compute_pij(pij, coverage, majority, data_points, previous_rates, previous_weights, rate_class_count);
            
            cout << ",\n\"prob_above_background\" : [";
            
            for (unsigned long i = 0; i < data_points; i++) {
                double sum = 0.;
                for (unsigned long j = 0; j < rate_class_count; j++) {
                    if (j != min_rate_index) {
                        sum += pij[i*rate_class_count + j];
                    }
                }
                cout << (i>0 ? ',' : ' ') << sum;
            }
            
            cout << "]\n";
            delete [] pij;
        }
        
        rate_class_count ++;
                
        if (!go_on) {
            delete [] rates;
            delete [] weights;
        }
        delete [] try_rates;
        delete [] try_weights;
    }
    
    cout << "}" << endl;       
    delete [] coverage;
    delete [] majority;
    if (previous_rates)   delete [] previous_rates;
    if (previous_weights) delete [] previous_weights;*/
    
    return 0;

}

