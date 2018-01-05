function pn = GeneratePN3

% Generates the maximal-length shift register
% sequence based on X^7 + X^6 + 1
% See "Linear Feedback Shift Register" in Wikipedia


%  this is for a 7 bit lfsr with a period of 127
% pn = zeros(127,1);
% STATE = zeros(7,1);
% STATE(1) = 1;

% for the 7 bit lfsr
% % for idx = 1:127
% %     
% %     %%% generate output
% %     pn(idx) = STATE(1);
% %     
% %     %%% update states
% %     STATE(1) = STATE(2);
% %     STATE(2) = STATE(3);
% %     STATE(3) = STATE(4);
% %     STATE(4) = STATE(5);
% %     STATE(5) = STATE(6);
% %     STATE(6) = mod(STATE(7)+pn(idx),2); % basically acts as an xor
% %     STATE(7) = pn(idx);
% % end

pn = zeros(7,1);
STATE = zeros(3,1);
STATE(1) = 1;
for idx = 1:7
    pn(idx) = STATE(1);
    
    %% update states
     STATE(1) = STATE(2);
    STATE(2) = mod(STATE(3) + pn(idx),2); % Xor. Between states 3 and 2 for x^3+x^2+1.
    STATE(3) =pn(idx);
    
end


return

a = 2*pn-1;
b = [a(end); a(1:end-1)];
x = zeros(126+126+1,1);

x=zeros(6+6+1,1);

for idx = 1:length(x)
    x(idx) = b'*a;
    b = [b(end); b(1:end-1)];
end
figure(99);
% stem(-126:126,x); grid on;
stem(-6:6,x);grid on;
