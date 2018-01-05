%%% Demo4Colt

pn = GeneratePN7; % this is a length 2^7 - 1 = 127 seqeunce
% pn = GeneratePN3;

apn = 2*pn - 1; % converts 0 -> -1, 1 -> +1

%%%
%%% demonstrate the cool correlation properties

bits = [randi([0 1],1000,1); pn; randi([0 1],900,1)]; % sandwich your pn code in random bits
abits = 2*bits - 1; % converts 0 -> -1, 1 -> +1

colt = zeros(length(bits)+length(pn)-1,1);

% doing an autocorrelation
for idx = 1:length(pn)-1
    colt(idx) = sum(abits(1:idx).*apn(end-idx+1:end));
end
for idx = length(pn):length(bits)
    colt(idx) = sum(abits(idx-length(pn)+1:idx).*apn);
end
for idx = length(bits)+1:length(bits)+length(pn)-1
    colt(idx) = sum(abits(idx-length(pn)+1:end).*apn(1:end+length(abits)-idx));
end
% figure(1);
plot(colt); grid on;

%%% by construction the pn seqeunce starts at index 1001
%%% the correlation will occur at n_peak = 1001 + length(pn) - 1
%%% so, the plot should have a peak at 1001 + 127 - 1 = 1127.


%%% Correlation is convolution excecpt the second one is not flipped before
%%% sliding. So, I can use

colt2 = conv(abits,apn(end:-1:1));
figure(2); plot(colt2); grid on;


%%% suppse the bit stream is 10 copies of the pn sequence
%%% what does the correlation look like?

abits = repmat(apn,10,1);
colt3 = conv(abits,apn(end:-1:1));
% figure(3); plot(colt3); grid on;

%%% you see 10 peaks, this tells you where the pn sequence
%%% occurs in bits (or abits).

%%% let's put some errors in there
abits(50) = -abits(50);
% abits(30) =  -abits(30);
% abits(31) =  -abits(31);
% abits(32) =  -abits(32);
% abits(33) =  -abits(33);
% abits(34) =  -abits(34);
% abits(35) =  -abits(35);
% abits(36) =  -abits(36);
% abits(37) =  -abits(37);
% abits(38) =  -abits(38);
% abits(39) =  -abits(39);
% abits(40) =  -abits(40);
% abits(800) =  -abits(800);
colt4 = conv(abits,apn(end:-1:1));
figure(4); plot(colt4); grid on;

% figure(5); plot(colt4); grid on;
axis([120 400 120 128]);
bits = (abits+1)/2;
%%% the correlation peak corresponding to the copy of the pn sequence where
%%% the error occurred is smaller by two.
%%% In general, the number of errors in an occurrence of one cycle of the
%%% pn sequence is (length(pn) - correlation_peak_value)_/2. 
%%% So, in this window you would update your bit error rate as follows
%%%
%%% First peak:
%%% bit_errors = bit_errors + 0
%%% total_bits = totbal_bits + 127
%%%
%%% Second peak:
%%% bit_errors = bit_errors + (127-125)/2;
%%% total_bits = total_bits + 127
%%%
%%% Third peak:
%%% bit_errors = bit_errors + 0
%%% total_bits = totbal_bits + 127
 filestuff = fopen('testDataPN.bin','w')
fwrite(filestuff,colt4','float32','n')
fclose(filestuff);