function [X] = snapshot_gen_sto(doas, snapshot, ncov, scov)
%SNAPSHOT_GEN_STO Generates snapshots for the stochastic model.
%Inputs:
%   doas - DOA vector. For 2D DOAs, each column represents a DOA pair.
%   ncov - Covariance matrix of the additive complex circular-symmetric
%          Gaussian noise. Can be a scalar, vector (for uncorrelated noise
%          with different powers), or a matrix.
%   scov - Covariance matrix of the source signals. Can be a scalar, vector
%          (for uncorrelated sources with different powers), or a matrix.
%Outputs:
%   X - Snapshots, where each columns is a single snapshot.
% if nargin <= 4
%     scov = 1;
% end
% if nargin <= 3
%     ncov = 1;
% end
% if nargin <= 2
%     snapshot = 1;
% end
c = 340; 
fl=2000; 
f0=3000; 
fh=4000; 
Fs=10000; 
bw=fh-fl;
Nfft=128;%选取每段的时域数据长度，即dft长度 
N=Nfft*snapshot;%产生时域信号的总点数 
Pn=ncov;  %因为噪声肯定一样，所以要令其固定，来确定不同信号源的功率 
SNR = 10*log10(scov / ncov);
Ps=Pn*10 .^(SNR/10); 

Nsignal = length(doas);
%下面产生高斯白噪声，并且低通滤波得到宽带信号，存放在st中 
ss=randn(N,Nsignal)*diag(sqrt(Ps)); 
st=[]; 
[bb,aa]=butter(10,[fl/(Fs/2),fh/(Fs/2)]); 
for r=1:Nsignal 
    temp=filter(bb,aa,ss(:,r)); 
    st=[st,temp]; 
end
%产生信号的功率谱 
window=[];%求功率谱加的窗 
noverlap=[]; 
[Pxx,f]=psd(ss(:,1),Nfft,Fs,window,noverlap); 
% figure,plot(f,Pxx); 
[Pyy,f]=psd(st(:,1),Nfft,Fs,window,noverlap); 
% figure,plot(f,Pyy); 

J1=round(Nfft*fl/Fs+1);%最DFT的第J1低频率fl对应个点 
J2=round(Nfft*fh/Fs+1);  
J=J2-J1+1;%选用的频率点数目 
fll=(J1-1)*Fs/Nfft;%第J1个点对应的模拟频率 
fhh=(J2-1)*Fs/Nfft;%第J2个点对应的模拟频率 
subband=linspace(fll,fhh,J);%各个子带频点 

wavelength = c/subband(1);
d = wavelength / 2;
design = design_array_1d('coprime', [3 4], d);
Nsensor = design.element_count;
Xfallsnapshot=zeros(Nsensor,J,snapshot);%用于存放多次快拍频域数据，Xfallsnapshot为3维矩阵

for number=1:snapshot%快拍循环开始 
     
    Sf=fft(st((number-1)*Nfft+1:number*Nfft,:)); 
    %fft是按列进行的，每一列都变换，对应各个信号 
    %Sf的规模为（Nfft，Nsignal） 
     
 
    Sf=Sf(J1:J2,:);%将Sf频率分量为零的部分去掉，保留规模为（J，Nsignal） 
     
    %产生出复噪声 
    randn('state',sum(100*clock)); 
    nn=sqrt(Pn/2)*randn(Nsensor,J); 
    randn('state',sum(100*clock)); 
    nn=nn+j*sqrt(Pn/2)*randn(Nsensor,J); 
    Nf=fft(conj(nn'));%因为fft对列进行，而这里需要对nt的各行进行fft，所以先对nt进行转置，fft后再转置回来 
    Nf=conj(Nf');%Nf的规模为(Nsensor,J) 
     
    %下面建立模型得到频域接受数据表达式，存放在Afsf中 
    Afsf=zeros(Nsensor,J); 

    for m=1:Nsignal 
        a=steerMultiFre(Nsensor, doas(m),c,subband);%a为规模为（Nsensor，J）
        Sf1=repmat(conj(Sf(:,m))',Nsensor,1);%将每一个信号变成(Nsensor,J)的规模 
        b = a.';
        as1=b.*Sf1; 
        Afsf=Afsf+as1; %各个信号累加得到频域信号数据   
    end 
    Xf=Afsf+Nf;%接收到的频域数据 
    Xfallsnapshot(:,:,number)=Xf; 
     
end%快拍循环完毕 

X = Xfallsnapshot;

end

%多频率点阵列流形形成函数

function s=steerMultiFre(Nsensor, doas,c,subband)

A = zeros(length(subband), Nsensor);
for ii = 1:length(subband)
    wavelength = c/subband(ii);
    d = wavelength / 2;
    design = design_array_1d('coprime', [3 4], d);
       A(ii, :) = exp(2j*pi/wavelength*(design.element_positions' * sin(doas)));
end
% 
s = A;
    
end
