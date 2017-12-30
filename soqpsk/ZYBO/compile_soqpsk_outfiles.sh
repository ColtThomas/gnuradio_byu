#!/bin/bash

g++ -Ofast -g3 -Wall -L/usr/local/lib -L/usr/local/include/uhd -L/usr/include/boost -o soqpsk_outfiles soqpsk_outfiles.cpp -luhd -lboost_thread -lboost_program_options
