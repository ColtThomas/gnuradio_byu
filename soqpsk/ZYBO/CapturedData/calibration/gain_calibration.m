i = 1;
scale = 1;
fp = fopen('1_13at400M.txt','r'); 
r_temp = fscanf(fp,'(%f,%f)\n');
processing_samples = r_temp(1:2:end) + 1i*r_temp(2:2:end);
processing_samples = scale*processing_samples;
fclose(fp);
figure(i)
plot(processing_samples, '.b' ); axis square; grid on;
i = i + 1;

for n = 1 : 10
    absolute(n) = abs(processing_samples(n));
end
alpha = 1; 
beta = .5;
for n = 1 :10
    if(abs(real(processing_samples(n))) > abs(imag(processing_samples(n))))
        temp(n) = alpha*abs(real(processing_samples(n))) + beta*abs(imag(processing_samples(n))) ;
    else
        temp(n) = alpha*abs(imag(processing_samples(n))) + beta*abs(real(processing_samples(n))) ;
    end
end
return %delete when done

fp = fopen('1_13at400M0G.txt','r'); 
r_temp = fscanf(fp,'(%f,%f)\n');
processing_samples = r_temp(1:2:end) + 1i*r_temp(2:2:end);
processing_samples = scale*processing_samples;
fclose(fp);
figure(i)
plot(processing_samples, '.b' ); axis square; grid on;
i = i + 1;

fp = fopen('1_13at1485M0G.txt','r'); 
r_temp = fscanf(fp,'(%f,%f)\n');
processing_samples = r_temp(1:2:end) + 1i*r_temp(2:2:end);
processing_samples = scale*processing_samples;
fclose(fp);
figure(i)
plot(processing_samples, '.b' ); axis square; grid on;
i = i + 1;

fp = fopen('1_13at1485M20G.txt','r'); 
r_temp = fscanf(fp,'(%f,%f)\n');
processing_samples = r_temp(1:2:end) + 1i*r_temp(2:2:end);
processing_samples = scale*processing_samples;
fclose(fp);
figure(i)
plot(processing_samples, '.b' ); axis square; grid on;
i = i + 1;

fp = fopen('1_13at 1485M38G.txt','r'); 
r_temp = fscanf(fp,'(%f,%f)\n');
processing_samples = r_temp(1:2:end) + 1i*r_temp(2:2:end);
processing_samples = scale*processing_samples;
fclose(fp);
figure(i)
plot(processing_samples, '.b' ); axis square; grid on;
i = i + 1;
