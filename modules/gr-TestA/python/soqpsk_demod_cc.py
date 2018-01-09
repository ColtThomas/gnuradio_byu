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
class soqpsk_demod_cc(gr.basic_block):
    """
    docstring for block soqpsk_demod_cc
    """
    samples_per_bit = 2
    #state variable test
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
    
    FX1 = 0;
    FX2 = 0;
    FX3 = 0;
    FX4 = 0;
    FX5 = 0;
    FX6 = 0;
    FX7 = 0;
    FX8 = 0;
    FX9 = 0;

    FY1 = 0;
    FY2 = 0;
    FY3 = 0;
    FY4 = 0;
    FY6 = 0;
    FY7 = 0;
    FY8 = 0;
    
    XI3 = 0;
    YI2 = 0;
    
    VP1 = 0;
    EP1 = 0;
    VT1 = 0;
    ET1 = 0;
    
    STROBE = 0;
    MU = 0;
    NCO = 0;
    CTHETA = 1;
    STHETA = 0;
    THETA = 0;
    def __init__(self):
        gr.basic_block.__init__(self,
            name="soqpsk_demod_cc",
            in_sig=[numpy.complex64],
            out_sig=[numpy.complex64])

    def forecast(self, noutput_items, ninput_items_required):
        for i in range(len(ninput_items_required)):
           ninput_items_required[i] = 4*noutput_items
           
    def sign(self, number):
        if (number>0.0): 
            return 1.0
        else:
            return -1.0
    def general_work(self, input_items, output_items):
        print '---------------------------'
        print 'Calling Filtering Function:'
        print '---------------------------'
        print 'inputs: ',input_items[0].size,' outputs: ',output_items[0].size
        #print 'states: '
        #print 'fx: ',self.FX1,self.FX2,self.FX3,self.FX4,self.FX5,self.FX6,self.FX7,self.FX8,self.FX9
        #print 'fy: ',self.FY1,self.FY2,self.FY3,self.FY4,self.FY6,self.FY7,self.FY8
        #print 'strobe: ', self.STROBE
        #print 'First input: ',input_items[0][0]
        #Timing PLL variables
        BnTst = 0.01;
        zetat = 1;
        N = 2;
        kpt = 12.35;
        k0t = -1;
        temp = BnTst/(zetat + 0.25/zetat);
        denom = 1 + 2*zetat/N*temp + temp*temp/(N*N);
        k0kpk1t = 4*zetat/N*temp/denom;
        k0kpk2t = 4*temp*temp/(N*N*denom);
        k1t = k0kpk1t/(kpt*k0t);
        k2t = k0kpk2t/(kpt*k0t);

        b0t = k1t + k2t;
        b1t = -k1t;
        
        
        #Phase PLL variables
        BnTsp = 0.02;
        zetap = 0.7071;
        N = 2;
        kpp = 18.33;
        k0p = 1;
        temp = BnTsp/(zetap + 0.25/zetap);
        denom = 1 + 2*zetap/N*temp + temp*temp/(N*N);
        k0kpk1p = 4*zetap/N*temp/denom;
        k0kpk2p = 4*temp*temp/(N*N*denom);
        k1p = k0kpk1p/(kpp*k0p);
        k2p = k0kpk2p/(kpp*k0p);
        
        b0p = k1p + k2p;
        b1p = -k1p;
        
        
        
        
        
        
        
        
        
        
        samples_per_buffer = 4*len(output_items[0])
        
        output_items[0] = [0]*len(output_items[0])
        debug=[0]*(len(output_items[0]))
        df_out=[0]*(len(output_items[0]))
        
        pk = 0;
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
            #df_out[sample_idx/2] = (x-y*1j) #NOTE! YOU NEED TO PUT X+Y*1J
            
            xr = x*self.CTHETA + y*self.STHETA;
            yr = -x*self.STHETA + y*self.CTHETA;
            
            #debug[sample_idx/2] = (xr -yr*1j); # DEBUG ----------
            
            if (self.STROBE == 0):
                et = 0;
                ep = 0;
            else:
                tempFx = -0.5*xr;
                tempFy = -0.5*yr;
                
                #compute interpolant yi from rotated DF outputs
                v2 = -tempFy + self.FY1 + self.FY2 - self.FY3;
                v1 = tempFy - self.FY1 + self.FY6 + self.FY2 + self.FY3;
                yi = (v2 * self.MU + v1) * self.MU + self.FY7;   
                
                #compute interpolants xi1 and yi1 rotated DF outputs
                v2 = -self.FX1 + self.FX2 + self.FX3 - self.FX4;
                v1 = self.FX1 - self.FX2 + self.FX7 + self.FX3 + self.FX4;
                xi1 = (v2 * self.MU + v1) * self.MU + self.FX8;
                
                v2 = -self.FY1 + self.FY2 + self.FY3 - self.FY4;
                v1 = self.FY1 - self.FY2 + self.FY7 + self.FY3 + self.FY4;
                yi1 = (v2 * self.MU + v1) * self.MU + self.FY8;   
                
                #compute interpolant xi2 from rotated DF outputs        
                v2 = -self.FX2 + self.FX3 + self.FX4 - self.FX5;
                v1 = self.FX2 - self.FX3 + self.FX8 + self.FX4 + self.FX5;
                xi2 = (v2 * self.MU + v1) * self.MU + self.FX9;  
                
                # compute et
                et = self.sign(xi2) * (xi1 - self.XI3) + self.sign(yi1) * (yi - self.YI2);
                
                # compute ep
                ep = self.sign(xi2)*self.YI2 - self.sign(yi1)*xi1;
                
                output_items[0][pk] = (xi2+1j*yi1) #invert on launch!!
                debug[pk] = (xi2+1j*yi1) #invert on launch!! (x+yj)
                #print 'output: ',debug[pk]
                pk = pk+1;
            #compute timing loop filter output
            vt = self.VT1 + b0t*et + b1t*self.ET1;
        
            #compute phase loop filter output
            vp = self.VP1 + b0p*ep + b1p*self.EP1;
            
            #compute NCO input
            w = vt + 0.5;
            
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
        
            self.FX5 = self.FX4;
            self.FX4 = self.FX3;
            self.FX3 = self.FX2;
            self.FX2 = self.FX1;
            self.FX1 = -0.5*xr;
            self.FX9 = self.FX8;
            self.FX8 = self.FX7;
            self.FX7 = self.FX6;
            self.FX6 = xr;
            
            self.FY4 = self.FY3;
            self.FY3 = self.FY2;
            self.FY2 = self.FY1;
            self.FY1 = -0.5*yr;
            self.FY8 = self.FY7;
            self.FY7 = self.FY6;
            self.FY6 = yr;
        
            if self.STROBE:
                self.XI3 = xi1;
                self.YI2 = yi;
            
            self.VP1 = vp;
            self.EP1 = ep;
            self.VT1 = vt;
            self.ET1 = et;     
            
            temp = self.NCO - w;
            if (temp < 0):
                self.STROBE = 1;
                self.MU = self.NCO/w;
                self.NCO = 1 + temp;
            else:
                self.STROBE = 0;
                self.NCO = temp;
                
                
            self.THETA = self.THETA + vp;
            self.CTHETA = math.cos(self.THETA);
            self.STHETA = math.sin(self.THETA);   
        filename = 'debugDF.csv'
        with open(filename, "a") as output:
            writer = csv.writer(output, lineterminator='\n')
            writer.writerow('filter called')
            for val in output_items[0]:
                writer.writerow([val]) 
        self.state = self.state + 1
        self.consume(0, len(output_items[0]))
        print 'done'
        return len(output_items[0])
