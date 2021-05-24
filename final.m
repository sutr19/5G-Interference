clc;
close all;
clear all;

d12=1;
d23=3;
d13=4;
tofmsg=input('enter type of message 1:URLLC;;2:MMTC;;3:HARQ');

if(tofmsg==2||tofmsg==3)
       disp('higher priority') 
priority=1;
else
    disp('low priority')
    priority=0;
end
cqi=input('enter channel quality indicator(0-1)');  


    
[filename pathname] = uigetfile({'*.txt'},'Select mixed   Text File');
if(filename ~=0)
    fid=fopen([pathname filename]);
    c=fscanf(fid','%f',[100,inf]);
end
inputdata=c;

txantennas=3;
rxantennas=3;

figure(1)
stem(inputdata), grid on;
title('INPUT DATA');
ld=length(inputdata);
axis([ 0 ld 0 1.5]);

b=c;

ln=length(b);

% Converting bit 0 to -1
for i=1:ln
    if b(i)==0
        b(i)=-1;
    end
end   

k=1;
for i=1:ln
   for j=1:127
       bb(k)=b(i);
       j=j+1;
       k=k+1;
   end
    i=i+1;
end
len=length(bb);
figure(2)
subplot(2,1,1);
stairs(bb); axis([0 len -2 3]);
title('input bit seq');

pr_sig=round(rand(1,len));
for i=1:len
    if pr_sig(i)==0
        pr_sig(i)=-1;
    end
end    
subplot(2,1,2);
stairs(pr_sig,'linewidth',2); axis([0 len -2 3]);
title('PSEUDORANDOM BIT SEQUENCE');

for i=1:len
   bbs(i)=bb(i).*pr_sig(i);
end

dsss=[];
t=0:1/10:2*pi;   
c1=cos(t);
c2=cos(t+pi);

for k=1:len
    if bbs(1,k)==-1
        dsss=[dsss c1];
    else
        dsss=[dsss c2];
    end
end
figure(3)
subplot(2,1,1);
stairs(bbs); axis([0 len -2 3]);
title('MULTIPLIER OUTPUT SEQUENCE');
subplot(2,1,2);
plot(dsss);
title('spreaded signal');


inputdata_NZR=2*inputdata-1;

serial2parallel=reshape(inputdata_NZR,2,length(inputdata)/2);  



br=input('enter data rate');     %   9.6*10^3    ,9600

f=br;
T=1/br; 
t=T/99:T/99:T; 



y=[];
y_in=[];
y_qd=[];
for(i=1:length(inputdata)/2)
    y1=serial2parallel(1,i)*cos(2*pi*f*t);
    y2=serial2parallel(2,i)*sin(2*pi*f*t) 
    
    y_in=[y_in y1];
    
    y_qd=[y_qd y2];
    
    y=[y y1+y2]; 
end
Tx_sig=y; %signal after modulation
tt=T/99:T/99:(T*length(inputdata))/2;

figure(4)

plot(tt,Tx_sig), grid on;
title('QPSK modulated signal i phase+q phase)');
xlabel('time');
ylabel(' amplitude');


Txdataawgn=awgn(Tx_sig,0.5);



slope = 30e6/(2e-3);
fb = 1/1.2288e6;    %chirp frequency
r = beat2range(fb,slope);

 Win=fft(Txdataawgn);
 
 M2=(fft(Win)); %First FFT for range information
M3=fftshift(fft(M2')); %   %interference power


[My,Ny]=size(M3);


idir=max(max(M3)/(M3-max(M3)));

noipow=abs(real(idir))+abs(imag(idir));


interf1=input('enter A1 interference for DIR');

interf2=input('enter A2 interference for DIR');

interf3=input('enter A3 interference for DIR');

dir1=abs(10*log(abs(real(idir))+abs(imag(idir))))+interf1;
dir2=abs(10*log(abs(real(idir))+abs(imag(idir))))+interf2;
dir3=abs(10*log(abs(real(idir))+abs(imag(idir))))+interf3;


if(txantennas>1)
   beta1=(noipow*d12/3)/(1+noipow);
end
if(txantennas==1)
beta1=1;
end

if(txantennas>1)
   beta2=(noipow*d23/3)/(1+noipow);
end
if(txantennas==1)
beta2=1;
end

if(txantennas>1)
   beta3=(noipow*d13/3)/(1+noipow);
end
if(txantennas==1)
beta3=1;
end



syms gamma1
val1=abs(solve((beta1+((noipow/3)/(1+noipow/(d12*gamma1)))==0),gamma1))

syms gamma2
val2=abs(solve((beta2+((noipow/3)/(1+noipow/(d12*gamma2)))==0),gamma2))

syms gamma3
val3=abs(solve((beta3+((noipow/3)/(1+noipow/(d12*gamma3)))==0),gamma3))


if(cqi>0.5 && priority==1 && dir1>2 && dir1<6)
antenna1=1;
else
    antenna1=0;
end
if(cqi>0.5 && priority==1 && dir2>2 && dir2<6)
antenna2=1;
else
antenna2=0;
end

if(cqi>0.5 && priority==1 && dir3>2 && dir3<6)
antenna3=1;
else
antenna3=0;
end

if(antenna1==1)
if(val1<val2 && val1<val3)
    disp('antenna 1 selected');
end
end
if(antenna2==1)
if(val2<val1 && val2<val3)
    disp('antenna 2 selected');
end
end
if(antenna3==1)
if(val3<val2 && val3<val1)
    disp('antenna 3 selected');
end
end


Fs = 12000;
figure(5)
%subplot(2,1,2)
plot(Txdataawgn)
xlabel('time');ylabel('Amplitude');title('qpsk noise signal')



hlpf = fdesign.lowpass('Fp,Fst,Ap,Ast',4.0e3,5.5e3,01,50,Fs);
D = design(hlpf);
yyr = filter(D,Txdataawgn);
figure(6)
plot(y)
xlabel('time');ylabel('Amplitude');title('despreaded signal and white gaussian removed')





Rx_inputdata=[];

Rx_sig=yyr;
for(i=1:1:length(inputdata)/2)

    Z_in=Rx_sig((i-1)*length(t)+1:i*length(t)).*cos(2*pi*f*t); 
    
    
    Z_in_intg=(trapz(t,Z_in))*(2/T);
    if(Z_in_intg>0)
        Rx_in_inputdata=1;
    else
       Rx_in_inputdata=0; 
    end
    
    Z_qd=Rx_sig((i-1)*length(t)+1:i*length(t)).*sin(2*pi*f*t);
    
    Z_qd_intg=(trapz(t,Z_qd))*(2/T);
        if (Z_qd_intg>0)
        Rx_qd_inputdata=1;
        else
       Rx_qd_inputdata=0; 
        end
        
        
        Rx_inputdata=[Rx_inputdata  Rx_in_inputdata  Rx_qd_inputdata];
end



figure(7)
stem(Rx_inputdata) 
title('Retreived data');
axis([ 0 ld 0 1.5]), grid on;

mse=10*sum(sum(Rx_inputdata))/length(Rx_inputdata);
display('psnr value is')
psnr=100-(sum(10*log10((10^2)./mse)));

ber=mse;
display('bit error rate')
display(ber);
