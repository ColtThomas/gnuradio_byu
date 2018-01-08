%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   QPSK matlab to python conversion
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Purpose of this is to prepare to convert a simple QPSK
% modulator/demodulator into python or C for GNURadio. This script breaks
% down Matlab functions into for loops for easy conversion.


%% Modulator

% Data to send
M = 4; % for qpsk
statement = 'A'; 
data = asciitobin(statement,M); % you still need to convert this to symbols

%%%%%%%%%%%%%%%%%%%%%%%
% Look-up Tables
%%%%%%%%%%%%%%%%%%%%%%%

% I just copied the 3/sqrt(2) from a 485 lab. I think it represents the 
% average energy or something. I'll have to check

LUT_I = [1 -1 1 -1]*3/sqrt(2); 
LUT_Q = [1 1 -1 -1]*3/sqrt(2);

data_i = zeros(1,length(data));
data_q = zeros(1,length(data));


% With graycoding in our QPSK constellation

% ex 3 -> (1,1) -> i=1 q=1 and 2-> (1,-1) -> i=1 q=-1

% Ready for Python
for idx = 1:length(data)    
    data_i(idx) = LUT_I(data(idx) + 1); % +1 because MATLAB is sacreligious and starts arrays from 1
    data_q(idx) = LUT_Q(data(idx) + 1);
end

%%%%%%%%%%%%%%%%%%%%%%%
% Upsampling
%%%%%%%%%%%%%%%%%%%%%%%
N = 8;

% TODO: Implement with for loop
upsampledData_I = upsample(data_i,N,N-1); % upsampling by 8, shifted by 7
upsampledData_Q = upsample(data_q,N,N-1);


%%%%%%%%%%%%%%%%%%%%%%%
% Matched Filter
%%%%%%%%%%%%%%%%%%%%%%%

% TODO: create array of coefficients instead
srrc = rcosdesign(0.5,12,8); % rolloff of 0.5, 8 samples per symbol, 12?
srrc = wrev(srrc); % makes convolving easier

% use a discrete filter function instead of a convolution function
% TODO: Use for loop to convolve. With downsampling?
mf_out_i = conv(upsampledData_I,srrc(end:-1:1)); % You will have to manually do this with a for loop in C
mf_out_q = conv(upsampledData_Q,srrc(end:-1:1));

% if you plot these you will see the srrc pulse shapes for each +1 and -1 

%%%%%%%%%%%%%%%%%%%%%%%
% Carrier Frequency Mixing
%%%%%%%%%%%%%%%%%%%%%%%

% purpose: When transmitting a signal wirelessly, transmit it at an 
% expected frequency

phase = pi/2; % separate sin and cos
freq = 2*pi/4; % I think this can be whatever radian you want it at


carrier_i = zeros(1,length(mf_out_i));
carrier_q = zeros(1,length(mf_out_q));

t = -pi:0.001:pi;

% Our orthonormal basis functions (pg 238 eq 5.50) are created by
% multiplying our pulse shape by a sinusoid
% TODO: implement with convolution
for idx = 1:length(mf_out_i)
    carrier_i(idx) = mf_out_i(idx) * sqrt(2)*sin(freq*idx+phase); 
    carrier_q(idx) = mf_out_q(idx) *-sqrt(2)*sin(freq*idx);
end

%%%%%%%%%%%%%%%%%%%%%%%
% Send Modulated Signal
%%%%%%%%%%%%%%%%%%%%%%%
% TODO: implement with convolution
transmitOut = carrier_i+carrier_q;


%% Receiver

% notes: GNU Radio will hand you a buffer with all input samples from the 
% USRP. The buffer size will vary; the faster your computer is, the larger
% the buffers. If the data to operate on is finite (ie. you aren't running
% real time) the buffer size will diminish to a size of 1 as you compute
% your last input sample. Buffer sizes are always in powers of 2.

receive_i = zeros(1,length(transmitOut));
receive_q = zeros(1,length(transmitOut));
% bring back to baseband
% TODO: combine with matched filter
for idx = 1:length(transmitOut)
   receive_i(idx) = transmitOut(idx)*sqrt(2)*sin(freq*idx+phase); 
   receive_q(idx) = transmitOut(idx)*-sqrt(2)*sin(freq*idx);
end


%%%%%%%%%%%%%%%%%%
% Matched Filter
%%%%%%%%%%%%%%%%%%

% Since GNU Radio has a realtime buffer of a variable size, we need to
% immitate the buffer mechanics. Suppose we obtain 1000 samples from the
% USRP. GNU Radio will vary the buffer sizes between computational blocks
% (ie filtering, downsampling,...) to allow real time computation. This 
% prevents excess accumulation of samples between blocks. In your
% convolution block, you will notice the buffer changes size as follows:

% Buf. Size | Samples Remaining
%------------------------------
%   512     |   1000 - 2^9 = 488
%   256     |   488 - 2^8 = 232
%   128     |   232 - 2^7 = 104
%   64      |   104 - 2^6 = 40
%   32      |   40 - 2^5 = 8
%   8       |   8 - 2^3 = 0
%
%   Done! Buffer size always a power of 2


% Below is an implementation of the buffer mechanism (implemented by 
% default in GNU Radio). Since GNU Radio handles all the buffer sizing for 
% you, just focus on the 'GNU Radio Computation' section.

% IGNORE: buffer setup
    k = 1;
    while length(transmitOut) - 2^k > 0
       k = k+1; % Find initial buffer size
    end
    n = k-1; % initial power of 2
    k = 0; % recycling variable for next while loop
    buffer = zeros(1,2^n);
    buffer_i = zeros(1,2^n);
    buffer_q = zeros(1,2^n);
    placehold= 0; % placeholder
    initLen = length(buffer);
    bufIdx = 1;
    tmp = 0;
    done = 0; % indicates end of data
    bufLen = 2^n;
    history = zeros(1,length(srrc)-1);
% END IGNORE

% Relevant Initializations
offset = 1; % placeholder for convolutions
mf_recd_i_prog = zeros(1,bufLen+length(srrc)-1);
mf_recd_q_prog = zeros(1,bufLen+length(srrc)-1);

out_i = [];
out_q = [];
while done==0   
    
    % IGNORE: fill up current buffer
        for idx=1+placehold:initLen+tmp
            idx;
%             buffer(bufIdx) = transmitOut(idx);
            buffer_i(bufIdx) = receive_i(idx);
            buffer_q(bufIdx) = receive_q(idx);
            bufIdx = bufIdx+1;
        end

        disp('samples processed: ');
        disp(1+placehold);
        disp(initLen+tmp);
    % END IGNORE
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % GNU Radio Computation
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % You are given access to the buffer (buffer), an output vector 
    % (mf_out), sample history (history), and the current buffer length 
    % (bufLen). In GNU Radio, you can get the buffer length from the input
    % vector itself.
    
    
    % history buffer needs to be the length of the filter minus 1.
    % convolve to the last data point
    
    tmp_buffer_i= [history buffer_i]; % append filter length-1 prior samples
    tmp_buffer_q= [history buffer_q];
%     disp('length of current buffer w/ history: ')
%     length(tmp_buffer)
    
    for idx_a=1:length(tmp_buffer_i)-length(srrc)+1
        mf_recd_i_prog(idx_a) = 0;
        mf_recd_q_prog(idx_a) = 0;
        for idx_b = 1:length(srrc)
%             mf_recd_i_prog(idx_a) = mf_recd_i_prog(idx_a)+srrc(idx_b)*tmp_buffer(idx_a+idx_b-1)*sqrt(2)*sin(freq*(idx_a+idx_b-1)+phase);
%             mf_recd_q_prog(idx_a) = mf_recd_q_prog(idx_a)+srrc(idx_b)*tmp_buffer(idx_a+idx_b-1)*-sqrt(2)*sin(freq*(idx_a+idx_b-1));
            mf_recd_i_prog(idx_a) = mf_recd_i_prog(idx_a)+srrc(idx_b)*tmp_buffer_i(idx_a+idx_b-1);
            mf_recd_q_prog(idx_a) = mf_recd_q_prog(idx_a)+srrc(idx_b)*tmp_buffer_q(idx_a+idx_b-1);
        end
    end
   
    
    
%     history = 
    % In gnuradio, you will assign the given output vector the values of
    % mf_out (can be a real double or complex double, see tutorials). You
    % will also use the 'consume' function to report the number of inputs
    % you used. At this point, you will be finished with the filtering 
    % operation. For the sake of matlab, we will just append mf_out to an 
    % external vector.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % GNU Radio End
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % IGNORE: Check to see if we are done
        if initLen+tmp>=length(transmitOut)
            done=1;

        else
            bufIdx = 1; % start from the beginning
            placehold = initLen+tmp
            sampRemain = length(transmitOut)-placehold % remaining samples
            k = 0;
            while sampRemain-2^k>=0
                k=k+1;            
            end
            tmp = tmp+2^(k-1)
            disp('processing range: ' );
            disp(placehold+1)
            disp(initLen+tmp);
            bufLen = tmp; % used for GNU Radio portion
        end
    % END IGNORE
end

disp('finished');


%% Downsample

% purpose: discard all the extra samples between the symbols. At this point
% we want to take the peaks of the signal and translate them to 1's and 0's

% the sign function takes the positive or negative magnitude to +1's and
% -1's

% the +1 divide by two changes the array from +1 and -1's to 0's and 1's

downsampledData_I = (sign(downsample(mf_recd_i_prog,8,7) ));
downsampledData_Q = (sign(downsample(mf_recd_q_prog,8,7)) ); % I'm cheating here and letting the Q be a more significant digit so I don't need a LUT


downsampledData = [downsampledData_I;downsampledData_Q];


%% Detection

constellation_i = 3*[1 -1 1 -1] /sqrt(2);
constellation_q = 3*[1 1 -1 -1] /sqrt(2);

decision_prog = zeros(1,length(downsampledData_I));

for idx=1:length(downsampledData_I)
    i_bit = (constellation_i-downsampledData_I(idx)).^2;
    q_bit = (constellation_q-downsampledData_Q(idx)).^2;
    
    [temp, decision_prog(idx)] = min(i_bit+q_bit);
end

decision_prog = decision_prog - 1;
ascii_prog = M2ascii(decision_prog,M)
