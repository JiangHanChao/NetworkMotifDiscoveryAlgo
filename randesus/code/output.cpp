#include "output.hpp"
using std::setw;

double calc_mean(double* array, int num){

   if (num == 0) return std::sqrt(-1.0); //Return NaN
    double ret = 0;
    for (int i=1; i <= num; ++i){
        ret += array[i];
    }
    return ret / num;
}

double calc_deviation(double* array, int num, double mean){

    if (num <= 1) return std::sqrt(-1.0); //Return NaN
    double ret = 0;
    for (int i=1; i <= num; ++i){
        ret += std::pow(array[i] - mean, 2);
    }
    ret /= (num-1);
    return std::sqrt(ret);
}

// Converts an int into string
string int_to_str(const int & i)
{
   std::stringstream s;
   s<<i;
   return s.str();
}

// Get the number of digits an uint64 has.
int num_digits(uint64 i)
{
    int ret=1;
    while (i > 9) {
        i /= 10;
        ++ret;
    }
    return ret;
}

void
pretty_output(hash_map < graph64, uint64* > & res_graphs, short G_N, uint64* count_subgr, int num_r_nets, std::ofstream & outfile)
{
    vector < subgr_res > result_vec(res_graphs.size());
    subgr_res gr;
    int total_num_nets = num_r_nets + 1; // Random nets plus the original
    double *concentration = new double[total_num_nets];
    uint64 max_ID = 0; //For the result table: max_ID has to fit in first column

    int idx = 0;
    for (hash_map < graph64, uint64* >::const_iterator iter =
	 res_graphs.begin(); iter != res_graphs.end(); ++iter) {
		gr.id = iter->first;
        // Get the largest ID
        if (gr.id > max_ID) max_ID = gr.id;

        // Calculate frequency in the original graph
		gr.freq = ((double) (iter->second)[0] / (double) count_subgr[0]);
        concentration[0] = gr.freq;

        if (num_r_nets > 0) { // at least one random network was sampled
            gr.p_value = 0;
            for (int i=1; i < total_num_nets; ++i){
                // Build the concentration array for this subgraph
                concentration[i] = ((double) (iter->second)[i] / (double) count_subgr[i]);
                // Calculate the p-value of the subgraph
                if (concentration[i] > gr.freq) ++gr.p_value;
            }
            gr.p_value /= double(num_r_nets);
            // Calculate the Z-Score via mean and standard deviation
            gr.rand_mean = calc_mean(concentration, num_r_nets);
            gr.rand_sd = calc_deviation(concentration, num_r_nets, gr.rand_mean);
            if (gr.rand_sd==0) gr.z_score = std::sqrt(-1.0); // Return NaN
            else gr.z_score = (gr.freq - gr.rand_mean) / gr.rand_sd;
        }
		result_vec[idx++] = gr;
    }
    sort(result_vec.begin(), result_vec.end(), compare);
    const uint64 bitmask = 1ULL << (G_N * G_N - 1);

    // Width constants
    const int id_width = num_digits(max_ID), adj_width = G_N+1, freq_w = 10,
              mfreq_w = 13, sd_width = 14, z_width = 11, p_width = 9;
    outfile.precision(5);
    outfile.flags(std::ofstream::right);

    outfile << "Result overview:\n" << endl;
    // First header line of the result table
    outfile << setw(id_width) << "ID" << setw(adj_width) << "Adj"
            << setw(freq_w+2) << "Frequency";
    if (num_r_nets > 0)
        outfile << setw(mfreq_w+1) << "Mean-Freq" << setw(sd_width) << "Standard-Dev"
                << setw(z_width) << "Z-Score" << setw(p_width) << "p-Value";
    outfile << endl;

    // Second header line of the result table
    outfile << setw(id_width) << ' ' << setw(adj_width) << ' '
            << setw(freq_w+2) << "[Original]";
    if (num_r_nets > 0)
        outfile << setw(mfreq_w+1) << "[Random]" << setw(sd_width) << "[Random]";
    outfile << endl << endl;

    // Datalines of the result table
    for (vector < subgr_res >::const_iterator iter = result_vec.begin();
         iter != result_vec.end(); ++iter) {
        graph64 g = iter->id;
  		outfile << setw(id_width) << iter->id << ' ';
        // Output the first line of the adj-matrix
        for (int j = 0; j != G_N; ++j) {
        outfile << (((g & bitmask) == bitmask) ? '1' : '0');
        g <<= 1;
        }
        outfile << ' ';
		outfile	<< setw(freq_w) << iter->freq*100 << '%';
        if (num_r_nets > 0) // Only if random nets have been sampled
            outfile << setw(mfreq_w) << iter->rand_mean*100 << '%'
                    << setw(sd_width) << iter->rand_sd
                    << setw(z_width) << iter->z_score
                    << setw(p_width) << iter->p_value;
        outfile << endl;
        // Output the remaining lines of the adj-matrix
		for (int i = 1; i != G_N; ++i) {
            outfile << setw(id_width+1) << ' ';
			for (int j = 0; j != G_N; ++j) {
			outfile << (((g & bitmask) == bitmask) ? '1' : '0');
			g <<= 1;
			}
			outfile << endl;
		}
		outfile << endl;
    }

    delete[] concentration;
}
