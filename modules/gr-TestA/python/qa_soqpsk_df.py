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

from gnuradio import gr, gr_unittest
from gnuradio import blocks
from soqpsk_df import soqpsk_df
import numpy as np
import csv
class qa_soqpsk_df (gr_unittest.TestCase):

    def setUp (self):
        self.tb = gr.top_block ()

    def tearDown (self):
        self.tb = None

    def test_001_t (self):
        print 'calling test'
        # extract from the csv files
        with open('r_data.csv','rb') as f:
            reader1 = csv.reader(f)
            r_data = list(reader1)[0]
        with open('filter_test.csv','rb') as g:
            reader2 = csv.reader(g)
            filter_data = list(reader2)[0]
        
        # testing with just zeros
        #src_data = [0]*10000
        expected_results = [0] * 10000
        
        # testing with actual data
        src_data = [complex(i) for i in r_data]
        print src_data[0],' is the first piece'
        #expected_results = [complex(i) for i in filter_data]
        
        
        print 'data length: ',len(src_data)
        print 'length of expected: ',len(expected_results)
        
        
        src = blocks.vector_source_c (src_data) # declare source
        demod = soqpsk_df () # declare the demod
        demod.clear()
        dst = blocks.vector_sink_c () # declare the sink
        
        # Connections
        self.tb.connect (src, demod) 
        self.tb.connect (demod, dst)
        
        # Test run and check results
        self.tb.run ()
        result_data = dst.data ()
        print 'length of results: ',len(result_data)
        self.assertFloatTuplesAlmostEqual (expected_results[0:9999], result_data, 6)



if __name__ == '__main__':
    gr_unittest.run(qa_soqpsk_df, "qa_soqpsk_df.xml")
