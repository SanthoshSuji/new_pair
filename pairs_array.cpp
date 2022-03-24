#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <complex>
#include <valarray>
#include <numeric>
#include <ctime>
#include <cmath>
 
using namespace std;

//parameters
double spread_ma_period = 3;
double bb_sma_period = 20;
double std_range_1 = 2;
double std_range_2 = 2;
double std_increment = 0.5;
double allocation = 400000;
string oscillator = "ma_spread";
bool george_values[] = {0};
// bool george = false;
string ratio_type = "Equally_weighted";
string today = "2022_03_11";
int freq = 10;
string resolution = "minute";
string new_file = "CL_BZ_"+to_string(freq)+"_"+resolution+"_2022-03-08.csv";
string contract = "Futures";

vector<string> date1, date2, data1, data2;
string source = "close";
int data_avl = 1;

vector<double> standard_deviation_range;

void get_std_range(){
	double v1 = std_range_1, v2 = std_range_2;    
    double current =  std::min(v1, v2);
    double max = std::max(v1, v2);
	double increment_range = std_increment;
    for(int m=0; m<200; m++){ 
		if(current <= max){
			if (m==0){
				standard_deviation_range.push_back(current);
			}else{ 
				standard_deviation_range.push_back(current += increment_range);
			}
		}
    }  
}

double stand_deviation(vector<double> prices, double bb_sma_period)
{	
	double mean, sum=0, SD=0;
	for(int m=0; m<bb_sma_period; m++){
		sum += prices[m];
	}
	mean = sum/bb_sma_period;

	for(int n=0; n<bb_sma_period; n++)
	SD += pow(prices[n]-mean, 2);
	return sqrt(SD/bb_sma_period);
}

double get_sma(vector<double> values, double bb_sma_period)
{	
	double mean, sum=0, SD=0;
	for(int m=0; m<bb_sma_period; m++){
		sum += values[m];
	}
	mean = sum/bb_sma_period;
	return mean;
}

int main(){
    const clock_t start_time = clock();

    string fname = "./pairs_symbols.csv";
	string data_file = "./data_2022_03_08/" + new_file;
	vector<vector<string> > content, d_content;
	vector<string> row, d_row;
	string line, word, d_line, d_word, symbol1, symbol2;

	int no_data, row_count = 1;

    void read_csv_files(string fname, string data_file, vector<string> row, vector<string> d_row, string line, string word, string d_line, string d_word, vector<vector<string> > content, vector<vector<string> > d_content, string symbol1, string symbol2);
    read_csv_files(fname, data_file, row, d_row, line, word, d_line, d_word, content, d_content, symbol1, symbol2);

    return 0;
}

void read_csv_files(string fname, string data_file, vector<string> row, vector<string> d_row, string line, string word, string d_line, string d_word, vector<vector<string> > content, vector<vector<string> > d_content, string symbol1, string symbol2){
    // Read csv file
	fstream file (fname, ios::in);
	if(file.is_open())
	{
		while(getline(file, line))
		{
			row.clear();
            // cout<<"row 1\n";
 
			stringstream str(line);
 
			while(getline(str, word, ','))
                // cout<<"row in 1\n";
				row.push_back(word);
			content.push_back(row);
		}
	}
	else
		cout<<"Could not open the file\n";
    
    // Read data file
	fstream d_file (data_file, ios::in);
	if(d_file.is_open())
	{
		while(getline(d_file, d_line))
		{
			d_row.clear();
            // cout<<"row 1\n";
			stringstream str(d_line);
 
			while(getline(str, d_word, ','))
                // cout<<"row in 1\n";
				d_row.push_back(d_word);
			d_content.push_back(d_row);
		}
	}
	else{
		cout<<"Could not open the file\n";
	}

    void get_data(vector<vector<string> > content, vector<vector<string> > d_content, string symbol1, string symbol2);
    get_data(content, d_content, symbol1, symbol2);

//    cout<<d_content.size();
}

void get_data(vector<vector<string> > content, vector<vector<string> > d_content, string symbol1, string symbol2){
    for(int z=0; z<1; z++){ 
		cout <<"--------- George -----------"<<george_values[z]<<"\n";
		bool george = george_values[z];

        // call standard deviation func
        get_std_range();

        for(int m=0; m<standard_deviation_range.size()-1; m++){ 
			cout <<"-------- standard_deviation_range ----------"<<standard_deviation_range[m] <<"\n";
			cout <<"------- freq ----------" <<freq<<"\n";
			double std_arange = standard_deviation_range[m];

            for(int j=1;j<content.size();j++){
				symbol1 = content[j][1];
				symbol2 = content[j][2];

                try{
					cout<<"symbols"<< " : " << symbol1 << ", " << symbol2<< "\n";
					// get_data(d_content, content[j][1], content[j][2], george, std_arange);
                    for(int i=0;i<d_content[0].size();i++){
                        if((d_content[0][i] == symbol1) or (d_content[0][i] == symbol2))
                        {
                            for(int k=i;k<d_content.size();k++){
                                if(d_content[0][i] == symbol1){
                                    if (d_content[1][i] == source){
                                        // cout<<" --symbol1-- "<<"\n";
                                        date1.push_back(d_content[k][0]);
                                        data1.push_back(d_content[k][i]);
                                    }
                                }else if(d_content[0][i] == symbol2){
                                    if (d_content[1][i] == source){
                                        // cout<<" --symbol2-- "<<"\n";
                                        date2.push_back(d_content[k][0]);
                                        data2.push_back(d_content[k][i]);
                                    }
                                }
                            }
                        }
                        // cout<<"data size : " << data1.size() <<" , " << data2.size() <<"\n";
                        if ((data1.size() & data2.size()) > 0){
                            cout<<"real pairs : " <<symbol1<<" , "<<symbol2 <<"\n";
                            void data_collection(string, string, vector<string>, vector<string>, bool, double);
                            data_collection(symbol1, symbol2, data1, data2, george, std_arange);
                        }
                        else{
                            cout <<"No data : " << data1.size() <<" , " << data2.size()<<"\n";
                        }
                    }
				}
				catch(exception& e){
					cout<<"No data available : " << symbol1 <<" , " << symbol2 <<"\n";
					cout << e.what();
				}
            }
        }
    }
}

void data_collection(string symbol1, string symbol2, vector<string> data1, vector<string>data2, bool george, double std_arange){

    int size_items;
    if(data1.size()<data2.size()){
        size_items = data1.size();
    }else{
        size_items = data2.size();
    }

    // cout<<"size "<< size_items<<"\n";
    // spread calculation
	vector<double> spread;

	if(george){
		cout<<"--------------GeorgeT--------------"<<"\n";
		for (int m=0;m<size_items;m++){
			double spd1 = stof(data1[m].c_str()); //convert str to double
			double spd2 = stof(data2[m].c_str());
			double spread1 = spd2 - spd1;
			spread.push_back(spread1);
		}
	}else{
		cout<<"--------------GeorgeF--------------"<<"\n";
		for (int m=0;m<size_items;m++){
			double spd1 = stof(data1[m].c_str()); //convert str to double
			double spd2 = stof(data2[m].c_str());
			double spread1 = spd1 - spd2;
			spread.push_back(spread1);
		}
	}

	//ma_spread
	double three_dim_runningTotal = 0.0;
	int three_dim_windowSize = spread_ma_period;

	vector<double> ma_spread;
	vector<double> ma_spread1;

	if(contract == "Futures"){
		for(int a = 0; a < spread.size(); a++)
		{
			three_dim_runningTotal += spread[a];   // add
			if(a >= three_dim_windowSize)
				three_dim_runningTotal -= spread[a - three_dim_windowSize];   // subtract
			else
				ma_spread1.push_back(0);
			if(a >= (three_dim_windowSize - 1))  // output moving average
				// cout << "Your SMA: " << runningTotal / (double)windowSize;
				ma_spread1.push_back(three_dim_runningTotal / (double)three_dim_windowSize);
			
			if(a>0){
				ma_spread.push_back(ma_spread1[a]);
			}
		}
	}else{
		for(int a = 0; a < spread.size(); a++)
		{
			three_dim_runningTotal += spread[a];   // add
			if(a >= three_dim_windowSize)
				three_dim_runningTotal -= spread[a - three_dim_windowSize];   // subtract
			else
				ma_spread.push_back(0);
			if(a >= (three_dim_windowSize - 1))  // output moving average
				// cout << "Your SMA: " << runningTotal / (double)windowSize;
				ma_spread.push_back(three_dim_runningTotal / (double)three_dim_windowSize);
		}
		
	}

	vector<double> oscillator;
	oscillator = ma_spread;

	//standard deviation
	vector<double> stand2_dev;
	vector<double> oscil;
	int index_val = 20;
	int first_val = 0;
	
	for(int a = 0; a < oscillator.size(); a++){
		for(int b = first_val; b < index_val; b++){
			oscil.push_back(oscillator[b]);
		}
		stand2_dev.push_back(stand_deviation(oscil, bb_sma_period));
		first_val += 1;
		index_val += 1;
		oscil.clear();
	}

	vector<double> stand_dev;
	for(int a = 0; a < oscillator.size(); a++){
		if(a>bb_sma_period-spread_ma_period+1){
			for(int b = 0; b < oscillator.size(); b++){
				stand_dev.push_back(stand2_dev[b]);
			}
		}else{
			stand_dev.push_back(0);
		}
	}

	//sma
	double sma_runningTotal = 0.0;
	vector<double> bb2_sma;
	vector<double> sma_value;

	double mean_sum = 0;

	int n_val = 20;
	int f_val = 0;
	for(int a = 0; a < oscillator.size(); a++){
		for(int b = f_val; b < n_val; b++){
			sma_value.push_back(oscillator[b]);
		}
		bb2_sma.push_back(get_sma(sma_value, bb_sma_period));

		f_val += 1;
		n_val += 1;
		sma_value.clear();
	}

	vector<double> bb_sma, bb_up, bb_down, entry_signal, entry_count, entry_position, entry_position_count;
	
	int sum = 1;
	int entry_sum = 1;
	int position_count = 1;

	for(int a = 0; a < oscillator.size(); a++){
		if(a>bb_sma_period-spread_ma_period+1){
			for(int b = 0; b < oscillator.size(); b++){
				bb_sma.push_back(bb2_sma[b]);
			}
		}else{
			bb_sma.push_back(0);
		}

		double bollinger_up = bb_sma[a] + (stand_dev[a] * std_arange);
		double bollinger_down = bb_sma[a] - (stand_dev[a] * std_arange);
		bb_up.push_back(bollinger_up);
		bb_down.push_back(bollinger_down);

		if(a>=bb_sma_period){
			if(oscillator[a]<=bb_down[a]){
				entry_signal.push_back(1);
			}else if(oscillator[a]>=bb_up[a]){
				entry_signal.push_back(-1);
			}else{
				entry_signal.push_back(0);
			}
		}else{
			entry_signal.push_back(0);
		}

		if(a>=bb_sma_period){
			if(entry_signal[a]==1 or entry_signal[a]==-1){
				entry_count.push_back(sum);
				sum += 1;
			}else{
				entry_count.push_back(0);
			}
		}else{
			entry_count.push_back(0);
		}

		if(a>=bb_sma_period){
			if(entry_count[a]!=0){
				entry_position.push_back(entry_sum);
				entry_sum += 1;
			}else{
				entry_position.push_back(0);
			}
		}else{
			entry_position.push_back(0);
		}

		if(a>=bb_sma_period){
			if(entry_count[a]!=0){
				entry_position_count.push_back(position_count);
				position_count += 1;
			}
		}
	}

	//running pnl
	vector<double> running_pnl, initial_price_sym1_vec, initial_price_sym2_vec, qty_sym1_vec, qty_sym2_vec;

	double initial_price_sym1, initial_price_sym2, qty_sym1, qty_sym2 = 0;
	int start_signal = 0;

	for(int k=0;k<entry_position_count.size();k++){
		for(int c = 0; c < oscillator.size()-1; c++){
			if(c>=bb_sma_period){
				if(entry_position[c] == entry_position_count[k]){
					cout<<"entry_position --------------------- "<<entry_position_count[k]<<" -- " << entry_position[c] <<"\n";
					start_signal = int(entry_signal[c]);
					initial_price_sym1 = stof(data1[c].c_str());
					initial_price_sym2 = stof(data2[c].c_str());

					initial_price_sym1_vec.push_back(initial_price_sym1);
					initial_price_sym2_vec.push_back(initial_price_sym2);

					if(contract=="Futures"){
						qty_sym1 = round((allocation/2)/(1000*initial_price_sym1));
						qty_sym2 = round((allocation/2)/(1000*initial_price_sym2));
						qty_sym1_vec.push_back(qty_sym1);
						qty_sym2_vec.push_back(qty_sym2);
					}else{
						qty_sym1 = round((allocation/2)/(initial_price_sym1));
						qty_sym2 = round((allocation/2)/(initial_price_sym2));
						qty_sym1_vec.push_back(qty_sym1);
						qty_sym2_vec.push_back(qty_sym2);
					}

					for(int a = c; a < oscillator.size()-1; a++){
						if(start_signal == 1){
							running_pnl.push_back((-initial_price_sym1 + stof(data1[a+1].c_str()) * qty_sym1) + (initial_price_sym2 - stof(data2[a+1].c_str()) * qty_sym2));
						}else if(start_signal == -1){
							running_pnl.push_back((initial_price_sym1 - stof(data1[a+1].c_str()) * qty_sym1) + (-initial_price_sym2 + stof(data2[a+1].c_str()) * qty_sym2));
						}
					}
				}else{
					initial_price_sym1_vec.push_back(0);
					initial_price_sym2_vec.push_back(0);
					qty_sym1_vec.push_back(0);
					qty_sym2_vec.push_back(0);
					running_pnl.push_back(0);
				}
			}else{
				initial_price_sym1_vec.push_back(0);
				initial_price_sym2_vec.push_back(0);
				qty_sym1_vec.push_back(0);
				qty_sym2_vec.push_back(0);
				running_pnl.push_back(0);
			}
		}

		vector<double> running_pnl2;
		for(int a = 0; a < running_pnl.size(); a++){
			if(a>0){
				for(int b = 0; b < oscillator.size(); b++){
					running_pnl2.push_back(running_pnl[b]);
				}
			}else{
				running_pnl2.push_back(0);
			}
		}

		string georgeStr;

		if(george){
			georgeStr = "True";
		}else{
			georgeStr = "False";
		}

		std::ofstream myfile;
		string filename = "./ind_results_what_if/what_if_loop_" + symbol1 + "_" + symbol2 + "_" + to_string(std_arange) +  "_" + to_string(freq)+ "_" +resolution+  "_" + georgeStr + "_" + to_string(entry_position_count[k]) + ".csv";
		// std::ofstream myfile( filename, std::ios::app ) ;

		myfile.open(filename);

		myfile <<"date"<<","<<"ratio_type"<<","<<"freq"<<","<<"resolution"<<","<<"stand_dev_period"<<","<<symbol1<<","<<symbol2<<","<<"allocation"<<","<<"george"<<","<<"price_spread"<<","<<"ma_spread"<<","<<"oscillator"<<","<<"bb_sma"<<","<<"bb_up"<<","<<"bb_down"<<","<<"std"<<","<<"entry_signal"<<","<<"entry_count"<<","<<"entry_position"<<","<<"entry_pr_"+symbol1<<","<<"entry_pr_"+symbol2<<","<<"qty_"+symbol1<<","<<"qty_"+symbol2<<","<<"running_pnl"<<"\n";
		for(int z=bb_sma_period; z<oscillator.size()-1; z++){
			myfile <<date1[z+1]<<","<<ratio_type<<","<<freq<<","<<resolution<<","<<std_arange<<","<<data1[z]<<","<<data2[z]<<","<<allocation<<","<<georgeStr<<","<<spread[z]<<","<<ma_spread[z]<<","<<oscillator[z]<<","<<bb_sma[z]<<","<<bb_up[z]<<","<<bb_down[z]<<","<<stand_dev[z]<<","<<entry_signal[z]<<","<<entry_count[z]<<","<<entry_position[z]<<","<<initial_price_sym1_vec[z]<<","<<initial_price_sym2_vec[z]<<","<<qty_sym1_vec[z]<<","<<qty_sym2_vec[z]<<","<<running_pnl2[z]<<"\n";
		}
		myfile.close();

		initial_price_sym1_vec.clear();
		initial_price_sym2_vec.clear();
		qty_sym1_vec.clear();
		qty_sym2_vec.clear();
		running_pnl.clear();
		running_pnl2.clear();
	}
}