%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SOQPSK_demod.m
%%% Created for the USRP Radio team
%%% The data are sampled at 2 samples/bit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The easiest way to keep track of all the variables and their
%%% relationships to each other is to use the state-space formulation
%%%
%%% Notation: 
%%% X_k = the vector of states at time index k
%%% U_k = the vector inputs at time k
%%% Y_K = the vector of outputs at time k
%%% Update equations: Y_{k} = g(X_k,U_k)
%%%                   X_{k+1} = f(X_k,U_k)
%%% where f() is the function defining the relationship between the next
%%% states and the current states and inputs and g() is the function
%%% defining the relationship between the output and the current states and
%%% inputs.
%%%
%%% The states correspond to the contents of any memory elements in the
%%% simulik model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% initializations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Phase PLL

BnTsp = 0.02;
zetap = 0.7071;
N = 2;
kpp = 18.33;
k0p = 1;
temp = BnTsp/(zetap + 0.25/zetap);
denom = 1 + 2*zetap/N*temp + temp*temp/(N*N);
k0kpk1p = 4*zetap/N*temp/denom;
k0kpk2p = 4*temp*temp/(N*N*denom);
k1p = k0kpk1p/(kpp*k0p);
k2p = k0kpk2p/(kpp*k0p);

b0p = k1p + k2p;
b1p = -k1p;

%%% Timing PLL

BnTst = 0.01;
zetat = 1;
N = 2;
kpt = 12.35;
k0t = -1;
temp = BnTst/(zetat + 0.25/zetat);
denom = 1 + 2*zetat/N*temp + temp*temp/(N*N);
k0kpk1t = 4*zetat/N*temp/denom;
k0kpk2t = 4*temp*temp/(N*N*denom);
k1t = k0kpk1t/(kpt*k0t);
k2t = k0kpk2t/(kpt*k0t);

b0t = k1t + k2t;
b1t = -k1t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define the detection filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DF = [ 0.010378066969709;
       0.023688987704657;
       0.009767134822858;
      -0.027017804469398;
      -0.089762303133391;
      -0.110346523809347;
      -0.051853991233850;
       0.154921158891652;
       0.568943123186263;
       0.792392871766106;
       0.792392871766106;
       0.568943123186263;
       0.154921158891652;
      -0.051853991233850;
      -0.110346523809347;
      -0.089762303133391;
      -0.027017804469398;
       0.009767134822858;
       0.023688987704657;
       0.010378066969709];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% load the I/Q baseband data from the USRP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%load all captured data
fp = fopen('ctheta.txt','r');
ctheta_zybo = fscanf(fp,'%f');
fclose(fp);

fp = fopen('stheta.txt','r');
stheta_zybo = fscanf(fp,'%f');
fclose(fp);

fp = fopen('THETA.txt','r');
theta_zybo = fscanf(fp,'%f');
fclose(fp);

fp = fopen('decodedBits.txt','r');
decodedBits_zybo = fscanf(fp,'%f');
fclose(fp);

fp = fopen('ep.txt','r'); 
ep_zybo = fscanf(fp,'%f');
fclose(fp);

fp = fopen('et.txt','r'); 
et_zybo = fscanf(fp,'%f');
fclose(fp);

fp = fopen('mfout.txt','r');
mfout_zybo = fscanf(fp,'%f');
fclose(fp);

fp = fopen('mu.txt','r'); 
mu_zybo = fscanf(fp,'%f');
fclose(fp);

fp = fopen('strobe.txt','r');
strobe_zybo = fscanf(fp,'%f');
fclose(fp);

fp = fopen('xi.txt','r');
xi_zybo = fscanf(fp,'%f');
fclose(fp);

fp = fopen('yi.txt','r');
yi_zybo = fscanf(fp,'%f');
fclose(fp);

fp = fopen('xr.txt','r');
xr_zybo = fscanf(fp,'%f');
fclose(fp);

fp = fopen('yr.txt','r');
yr_zybo = fscanf(fp,'%f');
fclose(fp);

% fp = fopen('PN11.txt','r');
% pn11 = fscanf(fp,'%f');
% fclose(fp);

%end load captured data

scale = 1;
fp = fopen('mfin.txt','r'); %in the c++ code we already scaled mfin by 2000
r_temp = fscanf(fp,'(%f,%f)\n');
processing_samples = r_temp(1:2:end) + 1i*r_temp(2:2:end);
processing_samples = scale*processing_samples;
fclose(fp);
figure(10)
plot(processing_samples, '.b' ); axis square; grid on;
% plot(r(:,1), r(:,2), '.b'); axis square;
samples_per_buffer = length(processing_samples);
samples_per_buffer2 = samples_per_buffer/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize the states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S4Di = zeros(18,1);
S4Dq = zeros(18,1);

FX1 = 0;
FX2 = 0;
FX3 = 0;
FX4 = 0;
FX5 = 0;
FX6 = 0;
FX7 = 0;
FX8 = 0;
FX9 = 0;

FY1 = 0;
FY2 = 0;
FY3 = 0;
FY4 = 0;
FY6 = 0;
FY7 = 0;
FY8 = 0;

XI3 = 0;
YI2 = 0;

DBIT1 = 0;

VP1 = 0;
EP1 = 0;
VT1 = 0;
ET1 = 0;

STROBE = 0;
MU = 0;
NCO = 0;
THETA = 0;
CTHETA = 1;
STHETA = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% End initialize the states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% variables for plotting
%%% only for demonstration purposes

plotstrobe = NaN*zeros(samples_per_buffer2,1);
plotmu = NaN*zeros(samples_per_buffer2,1);
plotet = zeros(samples_per_buffer2,1);
plotep = zeros(samples_per_buffer2,1);
plotctheta = NaN*zeros(samples_per_buffer2,1);
plotstheta = NaN*zeros(samples_per_buffer2,1);

%%% end initializing the plot variables

aout = zeros(fix(samples_per_buffer2/2),1);
bits = zeros(fix(samples_per_buffer2),1);
pk = 1;
bk = 1;
n = 1;
for sample_idx = 1:2:samples_per_buffer
    
    %%% compute the outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% get next two samples
    ri1 = real(processing_samples(sample_idx+1));
    rq1 = imag(processing_samples(sample_idx+1));
    ri = real(processing_samples(sample_idx));
    rq = imag(processing_samples(sample_idx));
    
    %%% compute detection filter outputs
    x = DF(1) * (ri1 + S4Di(18)) ...
        + DF(2) * (ri + S4Di(17)) ...
        + DF(3) * (S4Di(1) + S4Di(16)) ...
        + DF(4) * (S4Di(2) + S4Di(15)) ...
        + DF(5) * (S4Di(3) + S4Di(14)) ...
        + DF(6) * (S4Di(4) + S4Di(13)) ...
        + DF(7) * (S4Di(5) + S4Di(12)) ...
        + DF(8) * (S4Di(6) + S4Di(11)) ...
        + DF(9) * (S4Di(7) + S4Di(10)) ...
        + DF(10) * (S4Di(8) + S4Di(9));
    
    y = DF(1) * (rq1 + S4Dq(18)) ...
        + DF(2) * (rq + S4Dq(17)) ...
        + DF(3) * (S4Dq(1) + S4Dq(16)) ...
        + DF(4) * (S4Dq(2) + S4Dq(15)) ...
        + DF(5) * (S4Dq(3) + S4Dq(14)) ...
        + DF(6) * (S4Dq(4) + S4Dq(13)) ...
        + DF(7) * (S4Dq(5) + S4Dq(12)) ...
        + DF(8) * (S4Dq(6) + S4Dq(11)) ...
        + DF(9) * (S4Dq(7) + S4Dq(10)) ...
        + DF(10) * (S4Dq(8) + S4Dq(9));
    
    %%% rotate DF outputs
    xr = x*CTHETA + y*STHETA;
    yr = -x*STHETA + y*CTHETA;
    
    %%% if STROBE make decisions and compute timing and phase errors
    if STROBE == 0
        et = 0;
        ep = 0;
    else
        tempFx = -0.5*xr;
        tempFy = -0.5*yr;
        
        %%% compute interpolant yi from rotated DF outputs
        v2 = -tempFy + FY1 + FY2 - FY3;
        v1 = tempFy - FY1 + FY6 + FY2 + FY3;
        yi = (v2 * MU + v1) * MU + FY7;
        
        %%% compute interpolants xi1 and yi1 rotated DF outputs
        v2 = -FX1 + FX2 + FX3 - FX4;
        v1 = FX1 - FX2 + FX7 + FX3 + FX4;
        xi1 = (v2 * MU + v1) * MU + FX8;
        
        v2 = -FY1 + FY2 + FY3 - FY4;
        v1 = FY1 - FY2 + FY7 + FY3 + FY4;
        yi1 = (v2 * MU + v1) * MU + FY8;
        
        %%% compute interpolant xi2 from rotated DF outputs
        
        v2 = -FX2 + FX3 + FX4 - FX5;
        v1 = FX2 - FX3 + FX8 + FX4 + FX5;
        xi2 = (v2 * MU + v1) * MU + FX9;
        
        %%% compute et
        et = sign(xi2) * (xi1 - XI3) + sign(yi1) * (yi - YI2);
        
        %%% compute ep
        ep = sign(xi2)*YI2 - sign(yi1)*xi1;
        
        %%% output
        aout(pk) = xi2 + 1i*yi1;
        pk = pk + 1;
        if xi2 > 0
            d0 = 1;
        else
            d0 = 0;
        end
        if yi1 > 0
            d1 = 1;
        else
            d1 = 0;
        end
        bits(bk) = ~xor(d0,DBIT1);
        bits(bk+1) = xor(d1,d0);
        bk = bk+2;
    end
    
    %%% compute timing loop filter output
    vt = VT1 + b0t*et + b1t*ET1;
    
    %%% compute phase loop filter output
    vp = VP1 + b0p*ep + b1p*EP1;

    %%% compute NCO input
    w = vt + 0.5;
        
    %%% end compute the outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% update the plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    plotstrobe(n) = STROBE;
    plotmu(n) = MU;
    plotctheta(n) = CTHETA;
    plotstheta(n) = STHETA;
    plotep(n) = ep;
    plotet(n) = et;
    n = n+1;
    
    %%% end update the plot data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% compute the next states %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    S4Di(18) = S4Di(16);
    S4Di(17) = S4Di(15);
    S4Di(16) = S4Di(14);
    S4Di(15) = S4Di(13);
    S4Di(14) = S4Di(12);
    S4Di(13) = S4Di(11);
    S4Di(12) = S4Di(10);
    S4Di(11) = S4Di(9);
    S4Di(10) = S4Di(8);
    S4Di(9) = S4Di(7);
    S4Di(8) = S4Di(6);
    S4Di(7) = S4Di(5);
    S4Di(6) = S4Di(4);
    S4Di(5) = S4Di(3);
    S4Di(4) = S4Di(2);
    S4Di(3) = S4Di(1);
    S4Di(2) = ri;
    S4Di(1) = ri1;
    
    S4Dq(18) = S4Dq(16);
    S4Dq(17) = S4Dq(15);
    S4Dq(16) = S4Dq(14);
    S4Dq(15) = S4Dq(13);
    S4Dq(14) = S4Dq(12);
    S4Dq(13) = S4Dq(11);
    S4Dq(12) = S4Dq(10);
    S4Dq(11) = S4Dq(9);
    S4Dq(10) = S4Dq(8);
    S4Dq(9) = S4Dq(7);
    S4Dq(8) = S4Dq(6);
    S4Dq(7) = S4Dq(5);
    S4Dq(6) = S4Dq(4);
    S4Dq(5) = S4Dq(3);
    S4Dq(4) = S4Dq(2);
    S4Dq(3) = S4Dq(1);
    S4Dq(2) = rq;
    S4Dq(1) = rq1;
    
    
    FX5 = FX4;
    FX4 = FX3;
    FX3 = FX2;
    FX2 = FX1;
    FX1 = -0.5*xr;
    FX9 = FX8;
    FX8 = FX7;
    FX7 = FX6;
    FX6 = xr;
    
    FY4 = FY3;
    FY3 = FY2;
    FY2 = FY1;
    FY1 = -0.5*yr;
    FY8 = FY7;
    FY7 = FY6;
    FY6 = yr;
    
    if STROBE
        XI3 = xi1;
        YI2 = yi;
        DBIT1 = d1;
    end
    
    VP1 = vp;
    EP1 = ep;
    VT1 = vt;
    ET1 = et;
    
    temp = NCO - w;
    if temp < 0
        STROBE = 1;
        MU = NCO/w;
        NCO = 1 + temp;
    else
        STROBE = 0;
        NCO = temp;
    end
    
    THETA = THETA + vp;
    CTHETA = cos(THETA);
    STHETA = sin(THETA);
    
    %%% end compute the next states %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%% plot stuff
figure(1);
mycolors = get(gca,'colororder');

subplot(211);
plot(1:n-1,plotctheta,'-',1:n-1,plotstheta,'--');
grid on;
ylabel('cos/sin PHASE');
axis([1 n-1 -1.2 1.2]);

subplot(212);
plot(1:n-1,plotmu,'-');
grid on;
ylabel('\mu');
axis([1 n-1 -0.2 1.2]);

figure(2);
% plot(real(aout(1:pk-1)),imag(aout(1:pk-1)),'.');
plot(real(aout(1500:pk-1)),imag(aout(1500:pk-1)),'.');
grid on;
grid on;
axis square;
axis(4*[-1 1 -1 1]);
xlabel('inphase');
ylabel('quadrature');

figure(3);

subplot(211);
plot(1:n-1,plotep);
grid on;
axis([1 n-1 -5 5]);
ylabel('e_p');

subplot(212);
plot(1:n-1,plotet);
grid on;
axis([1 n-1 -10 10]);
ylabel('e_t');

figure(6);
load pn11.mat;
% maggie = xcorr(2*pn11-1,2*decodedBits_zybo(1:bk-3)-1);
maggie = xcorr(2*pn11-1,2*bits(1:bk-3)-1);
plot(1:bk-1,maggie(1:bk-1));
grid on;
