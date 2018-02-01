# gnuradio_byu

This repository is used to contain the material used for the Software Defined Radio group at Brigham Young University. The following page contains details of the project:

http://gnu-radio.groups.et.byu.net/dokuwiki/doku.php?id=start


## Notes pertaining to modules

During December, I went through most of the modules that I worked on over the last two years and wrote some notes down about each of them. This is copied straight from my handwritten notes. They may be useful if you are trying to track down a particular example or implementation

### Mar3_2016 (file)
First attemt to write demod in c++. Most likely broken

### Mar27_2016 (file)
Old modules: gr-dirt,gr-research, and gr-TestA. gr-dirt is worth deleting.

### SOQPSK (file)
Code written by Chris Nash and Hogstrom that could be useful. Also Matlab code of SOQPSK implementation given by Dr. Rice.

### gr-simple
examples of certain concepts and block types implemented in python:
-history functionality
-interpolation
-synchronous
-general

### gr-TestA
Contains another qpsk instance of the same tutorial, and a block called "cleanslate". Cleanslate was a sanity check that passes samples through. Not that useful.

Contains more pythin stuff that is useful:
-derp.py passes samples (not useful)
-lookup_table.py (not done)
-matfile.py - an attempt to load .mat files for testing. Didn't work.
-old_soqpsk_demod.py - I was trying to debug the matched filter... might be useful
-pnSequence.py (incomplete)
-qa_soqpsk.....py - some better written qa code
-other qa scripts
-several .csv files that were generated/used in testing. I remember a matched filter working, and a demodulator close to working...

### gr-lime
contains the first blocks used for the pn correlator, and the matched filter test for soqpsk. Run the qa code and see what happens...

### gr-kiwi
I messed this module up by accidentally installing gnuradio inside of it... may be worth trying some qa code.
