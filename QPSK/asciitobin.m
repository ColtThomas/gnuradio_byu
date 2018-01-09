function [ y ] = asciitobin(x,M)
% Takes a string and converts it to M-ary symbols (8 bit ascii)
% x - input vector
% M - bits per symbol
bits = '';
for idx=1:length(x)
    bits = [bits,dec2bin(x(idx),8)]; % was 7, now 8 for simplicity
end


BPS = log2(M); % bits per symbol

y = zeros(1,ceil(length(x)/BPS)); % if we have leftover bits for whatever reason, we will add zeros to the end of the data
% if mod(length(bits),2)==1
%    display('odd');
%    bits = [bits, 0];
% end


idx_n = 1;
for idx=1:BPS:length(bits)-1
   y(idx_n) = bin2dec(bits(idx:idx+1));
   idx_n = idx_n + 1;
end

end
