//
// Copyright 2010-2011,2014 Ettus Research LLC
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <uhd/types/tune_request.hpp>
#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/exception.hpp>

#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/thread.hpp>
#include <iostream>
#include <fstream>
#include <csignal>
#include <complex>
#include <queue>
#include <cmath>
#include <vector>

//The following libraries have been added only because they were
//added in the bpsk_continuous test code.

#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>

#define BUFFER_SIZE 0x200 //This is used for writing the data
			  //onto the FIFO in packets of 512

namespace po = boost::program_options;

static bool stop_signal_called = false;
void sig_int_handler(int){stop_signal_called = true;}

float mysign(float x){
/*	if (x < 0.0f)
		return -1.0f;
	else
		return 1.0f;*/
	if(x < 0.0f) return -1.0f;
	if(x > 0.0f) return 1.0f;
	return 0.0f;
} 

void bpsk(
    uhd::usrp::multi_usrp::sptr usrp,
    const std::string &cpu_format,
    const std::string &wire_format,
    const std::string &file,
    size_t samps_per_buff,
    unsigned long long num_requested_samples,
    double time_requested = 0.0,
    bool bw_summary = false,
    bool stats = false,
    bool null = false,
    bool enable_size_map = false,
    bool continue_on_bad_packet = false,
	double fpga_freq = 100000
){
	//variables defined by Ettus

	size_t sample_number;
	unsigned long long num_total_samps = 0;

	//variables for BPSK demod 

	// matched filter
	float DF[20] = {   0.010378066969709,
					   0.023688987704657,
					   0.009767134822858,
					  -0.027017804469398,
					  -0.089762303133391,
					  -0.110346523809347,
					  -0.051853991233850,
					   0.154921158891652,
					   0.568943123186263,
					   0.792392871766106,
					   0.792392871766106,
					   0.568943123186263,
					   0.154921158891652,
					  -0.051853991233850,
					  -0.110346523809347,
					  -0.089762303133391,
					  -0.027017804469398,
					   0.009767134822858,
					   0.023688987704657,
					   0.010378066969709};	
	const int FILTERLENGTH = sizeof(DF)/sizeof(DF[0]);//gives the number of elements in DF
    
	//std::vector<std::complex<float> > S4Di, S4Dq;
	std::vector<float>   S4Di, S4Dq;
	std::vector<std::complex<float> > sample_packet;
	std::complex<float> processing_sample, processing_sample1;
	float FX1, FX2, FX3, FX4, FX5, FX6, FX7, FX8, FX9, FY1, FY2, FY3, FY4, FY5, FY6, FY7, FY8, YI2, XI3;
	float EP1, VP1, ET1, VT1, MU, NCO, THETA, CTHETA, STHETA;
	bool STROBE;
	std::complex<float> ctemp;
	float x, y, xr, yr, xi, yi, xi1, yi1, xi2, ri, rq, ri1, rq1, tempFx, tempFy, v1, v2;
	float et, ep, vt, w, vp, temp;
	int idx;
	bool d0, d1 = {false};
	bool DBIT1 = {false};
	
	//PLL loop filter gains 	
	float k1t = -0.001285244042460;
	float k2t = -5.140976169840373e-06;
	float k1p = 0.001435532746273;
	float k2p = 1.914068136249203e-05;

	float b0t, b1t, b0p, b1p;

	uint32_t FIFO_out = open("/dev/xillybus_write_32", O_WRONLY);//Opens access to the DMA
	if(FIFO_out < 0){
		perror("Failed to open devfile");
		std::cout<<"This is in case you don't see the perror error."<<std::endl<<"Your FIFO_out didn't open correctly"<<std::endl;
		exit(1);
	}

	//convert fpga_freq to hold_cycles (which is the number of clock cycles to hold the bert clk high and low)
	const double fpga_clk_freq = 125000000; //125M Hz as defined by the fpga clk frequency
	int hold_cycles = (int)(fpga_clk_freq / fpga_freq / 2);
	uint32_t fpga_out_freq = open("/dev/uio0", O_RDWR);
	if(fpga_out_freq < 0){
		perror("Failed to open devFile");
		exit(1);
	}
		
	void *map_addr;	
	int size = 4;
	map_addr = mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_SHARED, fpga_out_freq, 0);
	if(map_addr == MAP_FAILED){
		perror("Failed to mmap");
		exit(1);
	}
	
	volatile unsigned int *pointer = (volatile unsigned int *)map_addr;
	*pointer = hold_cycles;

	//write to fpga register 
	
			
	uint16_t write_buffer_ptr = 0x0000;//points to the address the the buffer. Initialized at the beginning	
	uint32_t write_buffer[BUFFER_SIZE] = {0};//creates the size of the buffer
	uint32_t buffer32 = 0;//deals with creating the 32 bit packets (buffer packet of 32 bits)
	uint8_t shifter = 0;//deals with creating the 32 bit packets (how much to shift in order to create 32 bit buffer packet)

	// compute loop filter coefficients for streamlined computations
	b0t = k1t + k2t;
	b1t = -k1t;
	b0p = k1p + k2p;
	b1p = -k1p;

    //initialize the state variables

    ctemp.real(0.0f);
    ctemp.imag(0.0f);
    for (idx=0; idx< FILTERLENGTH; idx++){
    	S4Di.push_back(ctemp.real());
		S4Dq.push_back(ctemp.imag());
    }
    
    FX1 = 0.0f;
	FX2 = 0.0f;
	FX3 = 0.0f;
	FX4 = 0.0f;
	FX5 = 0.0f;
	FX6 = 0.0f;
	FX7 = 0.0f;
	FX8 = 0.0f;
	FX9 = 0.0f;

	FY1 = 0.0f;
	FY2 = 0.0f;
	FY3 = 0.0f;
	FY4 = 0.0f;
	FY5 = 0.0f;
	FY6 = 0.0f;
	FY7 = 0.0f;
	FY8 = 0.0f;

	XI3 = 0.0f;
	YI2 = 0.0f;

    VP1 = 0.0f;
    EP1 = 0.0f;
    VT1 = 0.0f;
    ET1 = 0.0f;

    MU = 0.0f;
    NCO = 0.0f;
    THETA = 0.0f;
    CTHETA = 1.0f;
    STHETA = 0.0f;
    STROBE = false;

    //create a receive streamer
    uhd::stream_args_t stream_args(cpu_format,wire_format);
    uhd::rx_streamer::sptr rx_stream = usrp->get_rx_stream(stream_args);

    uhd::rx_metadata_t md;
    std::vector<std::complex<float> > buff(samps_per_buff);
    std::ofstream outfile ;
    bool overflow_message = true;


    //setup streaming
    uhd::stream_cmd_t stream_cmd((num_requested_samples == 0)?
        uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS:
        uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE
    );
    stream_cmd.num_samps = num_requested_samples;
    stream_cmd.stream_now = true;
    stream_cmd.time_spec = uhd::time_spec_t();
    rx_stream->issue_stream_cmd(stream_cmd);
    
    boost::system_time start = boost::get_system_time();
    unsigned long long ticks_requested = (long)(time_requested * (double)boost::posix_time::time_duration::ticks_per_second());
    boost::posix_time::time_duration ticks_diff;
    boost::system_time last_update = start;
    unsigned long long last_update_samps = 0;
    
    typedef std::map<size_t,size_t> SizeMap;
    SizeMap mapSizes;

	std::queue<std::vector<std::complex<float> > > processing_queue;

    while(not stop_signal_called and (num_requested_samples != num_total_samps or num_requested_samples == 0)){
		boost::system_time now = boost::get_system_time();

        size_t num_rx_samps = rx_stream->recv(&buff.front(), buff.size(), md, 3.0, enable_size_map);

		for (sample_number=0; sample_number<samps_per_buff; sample_number += 2){
			processing_sample = 2000.0f*buff.at(sample_number);
			processing_sample1 = 2000.0f*buff.at(sample_number + 1);
			//compute current Matched Filter output
			
			ri1 = processing_sample1.real();
			rq1 = processing_sample1.imag();
			ri = processing_sample.real();
			rq = processing_sample.imag();	
			
			x =   DF[0] * (ri1 + S4Di[17]) 
				+ DF[1] * (ri + S4Di[16]) 
				+ DF[2] * (S4Di[0] + S4Di[15]) 
				+ DF[3] * (S4Di[1] + S4Di[14]) 
				+ DF[4] * (S4Di[2] + S4Di[13]) 
				+ DF[5] * (S4Di[3] + S4Di[12]) 
				+ DF[6] * (S4Di[4] + S4Di[11]) 
				+ DF[7] * (S4Di[5] + S4Di[10])
				+ DF[8] * (S4Di[6] + S4Di[9]) 
				+ DF[9] * (S4Di[7] + S4Di[8]);


			y =   DF[0] * (rq1 + S4Dq[17]) 
				+ DF[1] * (rq + S4Dq[16])
				+ DF[2] * (S4Dq[0] + S4Dq[15])
				+ DF[3] * (S4Dq[1] + S4Dq[14])
				+ DF[4] * (S4Dq[2] + S4Dq[13])
				+ DF[5] * (S4Dq[3] + S4Dq[12])
				+ DF[6] * (S4Dq[4] + S4Dq[11])
				+ DF[7] * (S4Dq[5] + S4Dq[10])
				+ DF[8] * (S4Dq[6] + S4Dq[9])
				+ DF[9] * (S4Dq[7] + S4Dq[8]);
    

			//rotate Matched Filter output
			xr = x * CTHETA + y * STHETA;
			yr = -(x * STHETA) + y * CTHETA;

			//update error signals and bit decisions
			if(!STROBE){
					et = 0.0f;
					ep = 0.0f;
			} else{
				//interpolate rotated Matched Filter output to compute xi and yi

				/* taken out to match matlab code				
				tempFx = -0.5f * xr;
				v2 = -tempFx + FX1 + FX2 - FX3;
				v1 = tempFx - FX1 + FX5 + FX2 + FX3;
				xi = (v2 * MU + v1) * MU + FX6;
				*/
	
				// compute yi 				
				tempFy = -0.5f * yr;
				v2 = -tempFy + FY1 + FY2 - FY3;
				v1 = tempFy - FY1 + FY6 + FY2 + FY3;
				yi = (v2 * MU + v1) * MU + FY7;

				//compute interpolants xi1 and yi1
				v2 = -FX1 + FX2 + FX3 - FX4;
				v1 = FX1 - FX2 + FX7 + FX3 + FX4;
				xi1 = (v2 * MU + v1) * MU + FX8;

				v2 = -FY1 + FY2 + FY3 - FY4;
				v1 = FY1 - FY2 + FY7 + FY3 + FY4;
				yi1 = (v2 * MU + v1) * MU + FY8;

				//compute interpolant xi2 
		        v2 = -FX2 + FX3 + FX4 - FX5;
				v1 = FX2 - FX3 + FX8 + FX4 + FX5;
				xi2 = (v2 * MU + v1) * MU + FX9;
				
				// compute et
				et = mysign(xi2) * (xi1 - XI3) + mysign(yi1) * (yi - YI2);
				
				// compute ep
				ep = mysign(xi2) * YI2 - mysign(yi1) * xi1;

				// output
				if (xi2 > 0.0f){
					d0 = true;
				}
				else  {
					d0 = false;
				}

				if (yi1 > 0.0f){
					d1 = true;
				}
				else  {
					d1 = false;
				}

				//////////////////// code from working bpsk test ////////////////;
				buffer32 += (!(d0^ DBIT1) << shifter); //31 - shifter; 
				++shifter;
				buffer32 += ((d1^ d0) << shifter);
				++shifter;
				if (shifter == 32){
					shifter = 0;
					write_buffer[write_buffer_ptr] = buffer32;
					buffer32 = 0;
					++write_buffer_ptr;
				}
				
				//0x200buffer
				if (write_buffer_ptr >= BUFFER_SIZE){	
					int32_t error;
					while (1){ //writing
						//error = write(FIFO_out, &write_buffer[0], 0x1C14);
						error = write(FIFO_out, &write_buffer[0], sizeof(write_buffer));
						//returns number of bytes written
						//cout << "bytes written: " << error << endl;
						if ((error < 0) && (errno == EINTR)) {
							continue;
							//interrupted try again
					
							std::cout << "write error: " << error;
							//close(FIFO_out);
							//return -1;
						}
						if (error < 0){
							perror("write() failed");
						//std::cout<<"Write failed"<<std::endl;
							break;
						}
						if (error == 0){
							fprintf(stderr, "Reached write EOF (?!) \n");
							//std::cout<<"Reached Write EOF"<<std::endl;
							break;
						}
						//do something with "error" bytes of data
						//std::cout << "Error=" << error << std::endl;
						if(error == sizeof(write_buffer)){
							break;
						}		
					}
					write_buffer_ptr = 0x0000;
				}
			}

			// compute timing loop filter output
			vt = VT1 + b0t * et + b1t * ET1;

			// compute phase loop filter output;
			vp = VP1 + b0p * ep + b1p * EP1;

			// compute NCO input
			w =  vt + 0.5f;

			//update Matched Filter buffer
			S4Di[17] = S4Di[15];
			S4Di[16] = S4Di[14];
			S4Di[15] = S4Di[13];
			S4Di[14] = S4Di[12];
			S4Di[13] = S4Di[11];
			S4Di[12] = S4Di[10];
			S4Di[11] = S4Di[9];
			S4Di[10] = S4Di[8];
			S4Di[9] = S4Di[7];
			S4Di[8] = S4Di[6];
			S4Di[7] = S4Di[5];
			S4Di[6] = S4Di[4];
			S4Di[5] = S4Di[3];
			S4Di[4] = S4Di[2];
			S4Di[3] = S4Di[1];
			S4Di[2] = S4Di[0];
			S4Di[1] = ri;
			S4Di[0] = ri1;
		
			S4Dq[17] = S4Dq[16];
			S4Dq[16] = S4Dq[14];
			S4Dq[15] = S4Dq[13];
			S4Dq[14] = S4Dq[12];
			S4Dq[13] = S4Dq[11];
			S4Dq[12] = S4Dq[10];
			S4Dq[11] = S4Dq[9];
			S4Dq[10] = S4Dq[8];
			S4Dq[9] = S4Dq[7];
			S4Dq[8] = S4Dq[6];
			S4Dq[7] = S4Dq[5];
			S4Dq[6] = S4Dq[4];
			S4Dq[5] = S4Dq[3];
			S4Dq[4] = S4Dq[2];
			S4Dq[3] = S4Dq[1];
			S4Dq[2] = S4Dq[0];
			S4Dq[1] = rq;
			S4Dq[0] = rq1;

			// update states

    
			FX5 = FX4;
			FX4 = FX3;
			FX3 = FX2;
			FX2 = FX1;
			FX1 = -0.5*xr;
			FX9 = FX8;
			FX8 = FX7;
			FX7 = FX6;
			FX6 = xr;
	
			FY4 = FY3;
			FY3 = FY2;
			FY2 = FY1;
			FY1 = -0.5*yr;
			FY8 = FY7;
			FY7 = FY6;
			FY6 = yr;

			if (STROBE){
				XI3 = xi1;
				YI2 = yi;
				DBIT1 = d1;
			}
			
			VP1 = vp;
			EP1 = ep;
			VT1 = vt;
			ET1 = et;
				
			temp = NCO - w;
			if (temp < 0){
				STROBE = true;
				MU = NCO / w;
				NCO = 1 + temp;
			} else{
				STROBE = false;
				NCO = temp; 
			}
			
			THETA = THETA + vp;
			if(THETA > 6.283185307179586f){
				THETA = THETA - 6.283185307179586f;
			}
			if(THETA < 0.0f){
				THETA = THETA + 6.283185307179586f;
			}

			CTHETA = (float)cosf(THETA);
			STHETA = (float)sinf(THETA);
		} // end for loop
		
		
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT) {
            std::cout << boost::format("Timeout while streaming") << std::endl;
            break;
        }
        if (md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW){
            if (overflow_message){
                overflow_message = false;
                std::cerr << boost::format(
                    "Got an overflow indication. Please consider the following:\n"
                    "  Your write medium must sustain a rate of %fMB/s.\n"
                    "  Dropped samples will not be written to the file.\n"
                    "  Please modify this example for your purposes.\n"
                    "  This message will not appear again.\n"
                ) % (usrp->get_rx_rate()*sizeof(std::complex<float>)/1e6);
            }
            continue;
        }
        if (md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE){
            std::string error = str(boost::format("Receiver error: %s") % md.strerror());
            if (continue_on_bad_packet){
                std::cerr << error << std::endl;
                continue;
            }
            else
                throw std::runtime_error(error);
        }

        if (enable_size_map){
			SizeMap::iterator it = mapSizes.find(num_rx_samps);
			if (it == mapSizes.end())
				mapSizes[num_rx_samps] = 0;
			mapSizes[num_rx_samps] += 1;
		}

        num_total_samps += num_rx_samps;
         //std::cout << num_total_samps << " samples processed" << std::endl;

		if (bw_summary){
			last_update_samps += num_rx_samps;
			boost::posix_time::time_duration update_diff = now - last_update;
			if (update_diff.ticks() > boost::posix_time::time_duration::ticks_per_second()) {
				double t = (double)update_diff.ticks() / (double)boost::posix_time::time_duration::ticks_per_second();
				double r = (double)last_update_samps / t;
				std::cout << boost::format("\t%f Msps") % (r/1e6) << std::endl;
				last_update_samps = 0;
				last_update = now;
			}
		}
        
        ticks_diff = now - start;
		if (ticks_requested > 0){
			if ((unsigned long long)ticks_diff.ticks() > ticks_requested)
				break;
		}
    } // end while loop
    
    if (stats){
		std::cout << std::endl;

		double t = (double)ticks_diff.ticks() / (double)boost::posix_time::time_duration::ticks_per_second();
		std::cout << boost::format("Received %d samples in %f seconds") % num_total_samps % t << std::endl;
		double r = (double)num_total_samps / t;
		std::cout << boost::format("%f Msps") % (r/1e6) << std::endl;
		
		if (enable_size_map) {
			std::cout << std::endl;
			std::cout << "Packet size map (bytes: count)" << std::endl;
			for (SizeMap::iterator it = mapSizes.begin(); it != mapSizes.end(); it++)
				std::cout << it->first << ":\t" << it->second << std::endl;
		}
	}	
} // close bpsk function

typedef boost::function<uhd::sensor_value_t (const std::string&)> get_sensor_fn_t;

bool check_locked_sensor(std::vector<std::string> sensor_names, const char* sensor_name, get_sensor_fn_t get_sensor_fn, double setup_time){
	if (std::find(sensor_names.begin(), sensor_names.end(), sensor_name) == sensor_names.end())
		return false;
	
	boost::system_time start = boost::get_system_time();
	boost::system_time first_lock_time;
	
	std::cout << boost::format("Waiting for \"%s\": ") % sensor_name;
	std::cout.flush();
	
	while (true){
		if ((not first_lock_time.is_not_a_date_time()) and
			(boost::get_system_time() > (first_lock_time + boost::posix_time::seconds(setup_time))))
		{
			std::cout << " locked!" << std::endl;
			break;
		}
		
		if (get_sensor_fn(sensor_name).to_bool()){
			if (first_lock_time.is_not_a_date_time())
				first_lock_time = boost::get_system_time();
			std::cout << "+";
			std::cout.flush();
		}
		else{
			first_lock_time = boost::system_time();	//reset to 'not a date time'
			
			if (boost::get_system_time() > (start + boost::posix_time::seconds(setup_time))){
				std::cout << std::endl;
				throw std::runtime_error(str(boost::format("timed out waiting for consecutive locks on sensor \"%s\"") % sensor_name));
			}
			
			std::cout << "_";
			std::cout.flush();
		}
		
		boost::this_thread::sleep(boost::posix_time::milliseconds(100));
	}
	
	std::cout << std::endl;
	return true;
}

int UHD_SAFE_MAIN(int argc, char *argv[]){
    uhd::set_thread_priority_safe();

    //variables to be set by po
    std::string args, file, type, ant, subdev, ref, wirefmt;
    size_t total_num_samps, spb;
    double rate, freq, gain, bw, total_time, setup_time, fpga_freq;

    //setup the program options
    po::options_description desc("Allowed options.");
    desc.add_options()
        ("help", "help message")
        ("args", po::value<std::string>(&args)->default_value(""), "multi uhd device address args")
        ("file", po::value<std::string>(&file)->default_value("usrp_samples.dat"), "name of the file to write binary samples to")
        ("type", po::value<std::string>(&type)->default_value("short"), "sample type: double, float, or short")
        ("nsamps", po::value<size_t>(&total_num_samps)->default_value(0), "total number of samples to receive")
        ("time", po::value<double>(&total_time)->default_value(0), "total number of seconds to receive")
        ("spb", po::value<size_t>(&spb)->default_value(10000), "samples per buffer")
        ("rate", po::value<double>(&rate)->default_value(1e6), "rate of incoming samples")
        ("freq", po::value<double>(&freq)->default_value(0.0), "RF center frequency in Hz")
		("fpga_freq", po::value<double>(&fpga_freq)->default_value(100000), "fpga ouput frequency in Hz")
        ("gain", po::value<double>(&gain), "gain for the RF chain")
        ("ant", po::value<std::string>(&ant), "daughterboard antenna selection")
        ("subdev", po::value<std::string>(&subdev), "daughterboard subdevice specification")
        ("bw", po::value<double>(&bw), "daughterboard IF filter bandwidth in Hz")
        ("ref", po::value<std::string>(&ref)->default_value("internal"), "waveform type (internal, external, mimo)")
        ("wirefmt", po::value<std::string>(&wirefmt)->default_value("sc16"), "wire format (sc8 or sc16)")
        ("setup", po::value<double>(&setup_time)->default_value(1.0), "seconds of setup time")
        ("progress", "periodically display short-term bandwidth")
        ("stats", "show average bandwidth on exit")
        ("sizemap", "track packet size and display breakdown on exit")
        ("null", "run without writing to file")
        ("continue", "don't abort on a bad packet")
        ("skip-lo", "skip checking LO lock status")
        ("int-n", "tune USRP with integer-N tuning")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    //print the help message
    if (vm.count("help")){
        std::cout << boost::format("UHD RX samples to file %s") % desc << std::endl;
        return ~0;
    }
    
	// check for float -- added by mdr
	if (type != "float"){
		std::cout << "type was not chosen to be float ... d'oh!" << std::endl;
		exit(-1);
	}

	// check for positive value --added by Chris
	if(fpga_freq < 0){
		std::cout << "fpga_freq must be greater than 0 ... :)" << std::endl;
		exit(-1);
	}
    bool bw_summary = vm.count("progress") > 0;
    bool stats = vm.count("stats") > 0;
    bool null = vm.count("null") > 0;
    bool enable_size_map = vm.count("sizemap") > 0;
    bool continue_on_bad_packet = vm.count("continue") > 0;
    
    if (enable_size_map)
		std::cout << "Packet size tracking enabled - will only recv one packet at a time!" << std::endl;

    //create a usrp device
    std::cout << std::endl;
    std::cout << boost::format("Creating the usrp device with: %s...") % args << std::endl;
    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(args);

    //Lock mboard clocks
    usrp->set_clock_source(ref);

    //always select the subdevice first, the channel mapping affects the other settings
    if (vm.count("subdev")) usrp->set_rx_subdev_spec(subdev);

    std::cout << boost::format("Using Device: %s") % usrp->get_pp_string() << std::endl;

    //set the sample rate
    if (rate <= 0.0){
        std::cerr << "Please specify a valid sample rate" << std::endl;
        return ~0;
    }
    std::cout << boost::format("Setting RX Rate: %f Msps...") % (rate/1e6) << std::endl;
    usrp->set_rx_rate(rate);
    std::cout << boost::format("Actual RX Rate: %f Msps...") % (usrp->get_rx_rate()/1e6) << std::endl << std::endl;

    //set the center frequency
    if (vm.count("freq")){	//with default of 0.0 this will always be true
		std::cout << boost::format("Setting RX Freq: %f MHz...") % (freq/1e6) << std::endl;
        uhd::tune_request_t tune_request(freq);
        if(vm.count("int-n")) tune_request.args = uhd::device_addr_t("mode_n=integer");
		usrp->set_rx_freq(tune_request);
		std::cout << boost::format("Actual RX Freq: %f MHz...") % (usrp->get_rx_freq()/1e6) << std::endl << std::endl;
	}

    //set the rf gain
    if (vm.count("gain")){
        std::cout << boost::format("Setting RX Gain: %f dB...") % gain << std::endl;
        usrp->set_rx_gain(gain);
        std::cout << boost::format("Actual RX Gain: %f dB...") % usrp->get_rx_gain() << std::endl << std::endl;
    }

    //set the IF filter bandwidth
    if (vm.count("bw")){
        std::cout << boost::format("Setting RX Bandwidth: %f MHz...") % bw << std::endl;
        usrp->set_rx_bandwidth(bw);
        std::cout << boost::format("Actual RX Bandwidth: %f MHz...") % usrp->get_rx_bandwidth() << std::endl << std::endl;
    }

    //set the antenna
    if (vm.count("ant")) usrp->set_rx_antenna(ant);

    boost::this_thread::sleep(boost::posix_time::seconds(setup_time)); //allow for some setup time

    //check Ref and LO Lock detect
    if (not vm.count("skip-lo")){
		check_locked_sensor(usrp->get_rx_sensor_names(0), "lo_locked", boost::bind(&uhd::usrp::multi_usrp::get_rx_sensor, usrp, _1, 0), setup_time);
		if (ref == "mimo")
			check_locked_sensor(usrp->get_mboard_sensor_names(0), "mimo_locked", boost::bind(&uhd::usrp::multi_usrp::get_mboard_sensor, usrp, _1, 0), setup_time);
		if (ref == "external")
			check_locked_sensor(usrp->get_mboard_sensor_names(0), "ref_locked", boost::bind(&uhd::usrp::multi_usrp::get_mboard_sensor, usrp, _1, 0), setup_time);
	}

    if (total_num_samps == 0){
        std::signal(SIGINT, &sig_int_handler);
        std::cout << "Press Ctrl + C to stop streaming..." << std::endl;
    }

#define bpsk_continuous_args(format) \
	(usrp, format, wirefmt, file, spb, total_num_samps, total_time, bw_summary, stats, null, enable_size_map, continue_on_bad_packet, fpga_freq)
    //recv to file
//    if (type == "double") recv_to_file<std::complex<double> >recv_to_file_args("fc64");
//    else if (type == "float") recv_to_file<std::complex<float> >recv_to_file_args("fc32");
//    else if (type == "short") recv_to_file<std::complex<short> >recv_to_file_args("sc16");
//    else throw std::runtime_error("Unknown type " + type);

	bpsk bpsk_continuous_args("fc32");

    //finished
    std::cout << std::endl << "Done!" << std::endl << std::endl;

    return EXIT_SUCCESS;
}
