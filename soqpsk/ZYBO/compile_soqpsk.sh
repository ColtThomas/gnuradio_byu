#!/bin/bash

g++ -Ofast -g3 -Wall -L/usr/local/lib -L/usr/local/include/uhd -L/usr/include/boost -o soqpsk soqpsk.cpp -luhd -lboost_thread -lboost_program_options
