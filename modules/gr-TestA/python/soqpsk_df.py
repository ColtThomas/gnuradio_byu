#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
# Copyright 2016 <+YOU OR YOUR COMPANY+>.
# 
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# 
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
# 

import numpy
from gnuradio import gr
import math
import csv
class soqpsk_df(gr.basic_block):
    """
    docstring for block soqpsk_df
    """
    samples_per_bit = 2
    state = 0
    DF = [
         0.010378066969709,
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
         0.010378066969709
         ]
    S4Di = [0.0] * 18
    S4Dq = [0.0] * 18
    
    def __init__(self):
        gr.basic_block.__init__(self,
            name="soqpsk_df",
            in_sig=[numpy.complex64],
            out_sig=[numpy.complex64])
    def clear(self):
        S4Di = [0.0] * 18
        S4Dq = [0.0] * 18
    def forecast(self, noutput_items, ninput_items_required):
        #setup size of input_items[i] for work call
        for i in range(len(ninput_items_required)):
           ninput_items_required[i] = 2*noutput_items

    def general_work(self, input_items, output_items):
        print '---------------------------'
        print 'Calling After Filtering Function:'
        print '---------------------------'
        print 'inputs: ',input_items[0].size,' outputs: ',output_items[0].size
        samples_per_buffer = 2*len(output_items[0])
        
        #pk = 0;
        #bk = 1;
        #n = 1;
        output_items[0] = [0]*len(output_items[0])
        debug=[0]*(len(output_items[0]))
        for sample_idx in range(0,samples_per_buffer-2,2):
            # compute outputs
            #print 'index: ',sample_idx
            # you filter this and then downsample at the same time
            ## get the next two samples
            ri1 = input_items[0][sample_idx+1].real
            rq1 = input_items[0][sample_idx+1].imag
            ri = input_items[0][sample_idx].real
            rq = input_items[0][sample_idx].imag
            
            
            # compute detection filter outputs
            # are we missing DF values?
            x = self.DF[0] * (ri1 + self.S4Di[17]) \
                + self.DF[1] * (ri + self.S4Di[16]) \
                +self.DF[2] * (self.S4Di[0] + self.S4Di[15]) \
                +self.DF[3] * (self.S4Di[1] + self.S4Di[14]) \
                +self.DF[4] * (self.S4Di[2] + self.S4Di[13]) \
                +self.DF[5] * (self.S4Di[3] + self.S4Di[12]) \
                +self.DF[6] * (self.S4Di[4] + self.S4Di[11]) \
                +self.DF[7] * (self.S4Di[5] + self.S4Di[10]) \
                +self.DF[8] * (self.S4Di[6] + self.S4Di[9]) \
                +self.DF[9] * (self.S4Di[7] + self.S4Di[8])
                
            y =self.DF[0] * (rq1 + self.S4Dq[17]) \
                +self.DF[1] * (rq + self.S4Dq[16]) \
                +self.DF[2] * (self.S4Dq[0] + self.S4Dq[15]) \
                +self.DF[3] * (self.S4Dq[1] + self.S4Dq[14]) \
                +self.DF[4] * (self.S4Dq[2] + self.S4Dq[13]) \
                +self.DF[5] * (self.S4Dq[3] + self.S4Dq[12]) \
                +self.DF[6] * (self.S4Dq[4] + self.S4Dq[11]) \
                +self.DF[7] * (self.S4Dq[5] + self.S4Dq[10]) \
                +self.DF[8] * (self.S4Dq[6] + self.S4Dq[9]) \
                +self.DF[9] * (self.S4Dq[7] + self.S4Dq[8])
            debug[sample_idx/2] = (x-y*1j)
            if sample_idx<2:
                print 'input',ri,' ',rq, ' and ',ri1,' ',rq1
                print (x), ' iteration: ',sample_idx,' of ' ,samples_per_buffer
            
            #Test phase 1
            output_items[0][(sample_idx/2)] = (x-y*1j) #NOTE! YOU NEED TO PUT X+Y*1J
            #pk = pk + 1
            
            #update the states
            self.S4Di[17] = self.S4Di[15];
            self.S4Di[16] = self.S4Di[14];
            self.S4Di[15] = self.S4Di[13];
            self.S4Di[14] = self.S4Di[12];
            self.S4Di[13] = self.S4Di[11];
            self.S4Di[12] = self.S4Di[10];
            self.S4Di[11] = self.S4Di[9];
            self.S4Di[10] = self.S4Di[8];
            self.S4Di[9] = self.S4Di[7];
            self.S4Di[8] = self.S4Di[6];
            self.S4Di[7] = self.S4Di[5];
            self.S4Di[6] = self.S4Di[4];
            self.S4Di[5] = self.S4Di[3];
            self.S4Di[4] = self.S4Di[2];
            self.S4Di[3] = self.S4Di[1];
            self.S4Di[2] = self.S4Di[0];
            self.S4Di[1] = ri;
            self.S4Di[0] = ri1;
            
            self.S4Dq[17] = self.S4Dq[15];
            self.S4Dq[16] = self.S4Dq[14];
            self.S4Dq[15] = self.S4Dq[13];
            self.S4Dq[14] = self.S4Dq[12];
            self.S4Dq[13] = self.S4Dq[11];
            self.S4Dq[12] = self.S4Dq[10];
            self.S4Dq[11] = self.S4Dq[9];
            self.S4Dq[10] = self.S4Dq[8];
            self.S4Dq[9] = self.S4Dq[7];
            self.S4Dq[8] = self.S4Dq[6];
            self.S4Dq[7] = self.S4Dq[5];
            self.S4Dq[6] = self.S4Dq[4];
            self.S4Dq[5] = self.S4Dq[3];
            self.S4Dq[4] = self.S4Dq[2];
            self.S4Dq[3] = self.S4Dq[1];
            self.S4Dq[2] = self.S4Dq[0];
            self.S4Dq[1] = rq;
            self.S4Dq[0] = rq1;
        
        
        filename = 'debugDF'+str(self.state)+'.csv'
        with open(filename, "a") as output:
            writer = csv.writer(output, lineterminator='\n')
            for val in debug:
                writer.writerow([val]) 
        self.state = self.state + 1
        self.consume(0, len(output_items[0]))
        return len(output_items[0])
