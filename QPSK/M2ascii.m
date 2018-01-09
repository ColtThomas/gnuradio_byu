function [ y ] = M2ascii( x,M,pream)
% M2ascii- Converts M-ary symbols to 8 bit ascii characters (instead of 7 like 485 labs)
%   x - double array - binary array
%   M - int - log2(M) bits per symbol
%   pream - string - preamble to search in message
%   shift - bool - shifts bits by one for offset
if any(x) > M
    error('input vector contains values larger than M-1')
elseif any(x) < 0
    error('no negative numbers permitted')
else
Nbit = log2(M);

a = dec2bin(x,Nbit).';  % you have to convert the decimal value to binary, but let Nbit tell how many bits to represent the data (avoid 0000000101, etc)
b = reshape(a,1,[]) - '0';
NB = 8*floor(length(b)/8); % used to be 7, now redundant
c=b(1:NB);

switch nargin
    case 3
        temp = bintoascii(c);
        idx = strfind(temp,pream);
        if(length(idx)>0)
            y = temp(idx(1):end);    
        else
           error('unable to detect preamble'); 
        end
    otherwise
        y = bintoascii(c);
end

end
