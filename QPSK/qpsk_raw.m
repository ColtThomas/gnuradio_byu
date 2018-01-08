%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   QPSK example in pure matlab
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Purpose of this is to prepare to convert a simple QPSK
% modulator/demodulator into python or C for GNURadio. This script
% specifically uses Matlab functions for readability.

%% Data to send
M = 4; % for qpsk
statement = 'Totally Working';
data = asciitobin(statement,M); % custom functions used
%% Look-up Tables

% purpose: this represents our designed constellation

% I just copied the 3/sqrt(2) from a 485 lab. I think it represents the 
% average energy or something. I'll have to check

LUT_I = [1 -1 1 -1]*3/sqrt(2); 
LUT_Q = [1 1 -1 -1]*3/sqrt(2);

data_i = zeros(1,length(data));
data_q = zeros(1,length(data));


% With graycoding in our QPSK constellation

% ex 3 -> (1,1) -> i=1 q=1 and 2-> (1,-1) -> i=1 q=-1
for idx = 1:length(data)    
    data_i(idx) = LUT_I(data(idx) + 1); % +1 because MATLAB is sacreligious and starts arrays from 1
    data_q(idx) = LUT_Q(data(idx) + 1);
end

%% Upsampling

% purpose: Avoid Inter-Symbol Interference; associates with the symbol span

N = 8;

upsampledData_I = upsample(data_i,N,N-1); % upsampling by 8, shifted by 7
upsampledData_Q = upsample(data_q,N,N-1);

%% Matched Filter

% purpose: Helps to identify waveforms at the receiver

% rcosine(1,8,'fir/sqrt',0.5,6) % use this, but figure out how to do it in
% python or C

srrc = rcosdesign(0.5,12,8); % rolloff of 0.5, 8 samples per symbol, 12?
% srrc = rcosine(1,8,'fir/sqrt',0.5,6);  % !!Problem!! We are getting errors because of the rcosdesign parameters



% use a discrete filter function instead of a convolution function
mf_out_i = conv(upsampledData_I,srrc(end:-1:1)); % You will have to manually do this with a for loop in C
mf_out_q = conv(upsampledData_Q,srrc(end:-1:1));

% if you plot these you will see the srrc pulse shapes for each +1 and -1 

%% Carrier Frequency Mixing

% purpose: When transmitting a signal wirelessly, transmit it at an 
% expected frequency

phase = pi/2; % separate sin and cos
freq = 2*pi/4; % I think this can be whatever radian you want it at


carrier_i = zeros(1,length(mf_out_i));
carrier_q = zeros(1,length(mf_out_q));

t = -pi:0.001:pi;

% Our orthonormal basis functions (pg 238 eq 5.50) are created by
% multiplying our pulse shape by a sinusoid
for idx = 1:length(mf_out_i)
    carrier_i(idx) = mf_out_i(idx) * sqrt(2)*sin(freq*idx+phase); 
    carrier_q(idx) = mf_out_q(idx) *-sqrt(2)*sin(freq*idx);
end

%% Send (add)

% purpose: combine the inphase and quadrature to form the signal being
% transmitted off the antenna

transmitOut = carrier_i+carrier_q;
%% Receive

% purpose: bring the signal back to baseband
receive_i = zeros(1,length(transmitOut));
receive_q = zeros(1,length(transmitOut));


for idx = 1:length(transmitOut)
   receive_i(idx) = transmitOut(idx)*sqrt(2)*sin(freq*idx+phase); 
   receive_q(idx) = transmitOut(idx)*-sqrt(2)*sin(freq*idx);
end

%% Matched Filter

% purpose: make sure that the signal we are reading is actually what we
% want vs. some garbage/noise sent at our carrier frequency

mf_recd_i = conv(receive_i,srrc(end:-1:1));
mf_recd_q = conv(receive_q,srrc(end:-1:1));
mf_recd = [mf_recd_i;mf_recd_q];

%% Downsample

% purpose: discard all the extra samples between the symbols. At this point
% we want to take the peaks of the signal and translate them to 1's and 0's

% the sign function takes the positive or negative magnitude to +1's and
% -1's

% the +1 divide by two changes the array from +1 and -1's to 0's and 1's

downsampledData_I = (sign(downsample(mf_recd_i,8,7) ));
downsampledData_Q = (sign(downsample(mf_recd_q,8,7)) ); % I'm cheating here and letting the Q be a more significant digit so I don't need a LUT


downsampledData = [downsampledData_I;downsampledData_Q];
%% Detection

constellation_i = 3*[1 -1 1 -1] /sqrt(2);
constellation_q = 3*[1 1 -1 -1] /sqrt(2);

decision = zeros(1,length(downsampledData_I));

for idx=1:length(downsampledData_I)
    i_bit = (constellation_i-downsampledData_I(idx)).^2;
    q_bit = (constellation_q-downsampledData_Q(idx)).^2;
    
    [temp, decision(idx)] = min(i_bit+q_bit);
end

decision = decision - 1;
ascii = M2ascii(decision,M,statement(1));
%% Plots

% figure(1); % constellation
% % plot(real(aout(1:pk-1)),imag(aout(1:pk-1)),'.');
% plot(mf_recd_i,mf_recd_q,'.');
% grid on;
% grid on;
% axis square;
% axis(4*[-1 1 -1 1]);
% xlabel('inphase');
% ylabel('quadrature');
