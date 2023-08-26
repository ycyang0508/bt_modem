clear, clf,
A=[-3  10   -4  -6   -3.5  20  -4  -6    -4   20       % power gain (dB)
   3    3    4   6.5  3.5    6   4   6.5   4    6        % Noise figure (dB)
   60  20  60   0   60    20  60  10    60  20       % Ip3o (dBm)
   80  20  80  10   80    30  80  20    80  30 ];   % IP2o (dBm)

Dr=1;        % data rate (Mbps)
N=100;       % N/2 : span (MHz)
Nl=100;      % resolution in points/MHz
Nip=-110;    % input noise floor (dBm/MHz)
Nop=-114;    % equvelant theormal noise floor (dBm/MHz)
Fr=[5];      % center frequency of the carrier (MHz)
Frm=[-20 ];  % power of the carrier (dBm)

h=.2;      % modulation factor
bt=.5;     % relative bandwidth for Gaussian shaping

Bp1=[0 4 6 0 3 60 5]; % [* fl3 fh3 * * UATT order]
%%% prn11 generator

n=11; m=2^n-1;
a=ones(1,n); a(1:11)=[1 1 1 0 1 0 1 0 1 1 1]; 
for i=1:m-n
    ni=n+i; a(ni)=rem(a(ni-2)+a(ni-11),2);
end
a=[a a]; e=a*2-1;  e(1:2)=zeros(1,2);

%%% source
lf=N*Nl;  % total points
q2=sqrt(2); p2=2*pi; qN=sqrt(N); qk=lf/1000;
qf=sqrt(lf);ql=sqrt(Nl); rand('seed',56);
    
F=0:lf-1; Fh=F(1:lf/2); Fu=F/Nl;
FT=F/N;  % time axis
FM=Fh/Nl;  % frequency axis


nip=10^((Nip)/20);  
Noise=randn(1,lf)*nip/ql; 

dr=N/Dr;  % number of points/bit
ds=lf/dr;  % number of bits
dh=round(dr/2);
ons=ones(1,dr);
for i=1:ds
    Di((i-1)*dr+1:i*dr)=e(i)*ons;
end


Dip=Di; bd=Dr*bt; io=bd*Nl; c=.347/io^2;
gau(1)=1;
for i=2:lf/2
    gau(i)=exp(-c*(i-1)^2); gau(lf-i+2)=gau(i);      
end
gau(lf/2+1)=gau(lf/2);

uf=fft(Di); uff=uf.*gau; Di=real(ifft(uff));
Dio=Di; Di=Di/dr*pi*h; U0=0; 


%%% GFSK generator      
for i=1:lf
    U0=U0+Di(i); Ui(i)=U0;
end  
Fe=10^(Frm/20)*cos(p2*Fr*FT+Ui)/qf*q2; Fa=Fe;

%%% AWGN channel  
so=Fa+Noise;    Faf=fft(so)/qf; 
Na=20*log10(abs(Faf));
Nb(1)=Na(1);  Nb(2:lf/2)=Na(2:lf/2)+3;

%%% 2nd AMP1
nt=2;G=A(1,nt); NF=A(2,nt); Pi3=A(3,nt); Pi2=A(4,nt);
%csamp
%%% 1st stage: BPF1
Bp=Bp1; nt=1; G=A(1,nt); NF=A(2,nt); 
%csbpfbut

%%% GFSK demodulator  
     to=so.*exp(-j*p2*Fr*FT)/q2;
     tfo=fft(to)/qf; 
     h0=1; sf=Dr*h0*2*Nl; tfz=zeros(1,lf); 
     tfz(1:sf)=tfo(1:sf); tfz(lf-sf:lf)=tfo(lf-sf:lf);
     txo=ifft(tfz)*qf; 

ang=angle(txo);
angtem=ang; 
     
for i=1:lf-1
    if ang(i+1)-ang(i)>5.6, ang(i+1:lf)=ang(i+1:lf)-p2;
    elseif ang(i+1)-ang(i)<-5.6, ang(i+1:lf)=ang(i+1:lf)+p2;
    end
end
anh=ang*dr/(pi*h);

fui=diff([anh(1) anh]);
Fui=fui(dh+1:dr:lf); 
Dui=Dio(dh+1:dr:lf); 

ss=Fui-Dui; 
lss=length(ss);
ss=ss(2:lss);

ferr=sqrt(ss*ss')/(lss-1);


