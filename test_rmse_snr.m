%% test for underdetermined wideband DOA estimation using coprime arrays v.s. snr
clear(); close all;


c = 340;
fl=2000; 
f0=3000; 
fh=4000; 
Fs=10000; 
bw=fh-fl;

Nfft=128;%选取每段的时域数据长度，即dft长度 
J1=round(Nfft*fl/Fs+1);%最DFT的第J1低频率fl对应个点 
J2=round(Nfft*fh/Fs+1);  
J=J2-J1+1;%选用的频率点数目 
fll=(J1-1)*Fs/Nfft;%第J1个点对应的模拟频率 
fhh=(J2-1)*Fs/Nfft;%第J2个点对应的模拟频率 
subband=linspace(fll,fhh,J);%各个子带频点

% but we have 12 sources here
doas = linspace(-pi/3, pi/3, 12);  

power_source = 1;
snapshot_count = 200;
source_count = length(doas);

grid_size = 360;
grid_size1 = 360
% we vary SNR from -20dB to 20dB
n_param_power_noise = 20;
param_power_noise = 10.^(-linspace(-20, 20, n_param_power_noise)/10);
n_repeat = 1% 200 

crb = zeros(n_param_power_noise, 1);
snr = zeros(n_param_power_noise, 1);
mse_SBL = zeros(n_param_power_noise, 1);
mse_SOMP_LS = zeros(n_param_power_noise, 1);
mse_SOMP_TLS = zeros(n_param_power_noise, 1);
mse_SBI = zeros(n_param_power_noise, 1);

% YHQin 2017/12/20
rmse_crb = zeros(n_param_power_noise, 1);
rmse_SBL = zeros(n_param_power_noise, 1);
rmse_SOMP_LS = zeros(n_param_power_noise, 1);
rmse_SOMP_TLS = zeros(n_param_power_noise, 1);
rmse_SBI = zeros(n_param_power_noise, 1);
% end YHQin

for ii = 1:n_param_power_noise
    power_noise = param_power_noise(ii);
    snr(ii) = 10*log10(power_source / power_noise);
   
    cur_SBL = 0;
    cur_SOMP_LS = 0;
    cur_SOMP_TLS = 0;
    cur_SBI = 0;
    
    for rr = 1:n_repeat
        [X] = snapshot_gen_sto(doas, snapshot_count, power_noise, power_source);
        for n=1:J 
            Mx=X(:,n,:); 
            Mx=squeeze(Mx);%此函数将维数为1的去掉，得到单一频率点的多次快拍频域接收数据矩阵 
            [Rowx,Colx]=size(Mx);
             R=(Mx*Mx')/snapshot_count;
        
             ff = subband(n);
             wavelength = c/ff;
            d = wavelength / 2;
            design_cp = design_array_1d('coprime', [3 4], d);
                    
            sp_SBL = sparse_SBL_1d(R, source_count, design_cp, wavelength, grid_size1, 9);
            if(length(sp_SBL.x_est) == length(doas))
                cur_SBL = cur_SBL + sum((sp_SBL.x_est - doas).^2);
            end
            
            sp_SOMP_LS = sparse_SOMP_LS_1d(R, source_count, design_cp, wavelength, grid_size1, 9);
            if(length(sp_SOMP_LS.x_est) == length(doas))
                cur_SOMP_LS = cur_SOMP_LS + sum((sp_SOMP_LS.x_est - doas).^2);
            end
            
            sp_SOMP_TLS = sparse_SOMP_TLS_1d(R, source_count, design_cp, wavelength, grid_size1, 9);
            if(length(sp_SOMP_TLS.x_est) == length(doas))
                cur_SOMP_TLS = cur_SOMP_TLS + sum((sp_SOMP_TLS.x_est - doas).^2);
            end
            
            sp_SBI = sparse_SBI_1d(R, source_count, design_cp, wavelength, grid_size1, 9);
            if(length(sp_SBI.x_est) == length(doas))
                cur_SBI = cur_SBI + sum((sp_SBI.x_est - doas).^2);
            end
        end
    end

    mse_SBL(ii) = cur_SBL / (source_count * n_repeat * J);
    mse_SOMP_LS(ii) = cur_SOMP_LS / (source_count * n_repeat * J);
    mse_SOMP_TLS(ii) = cur_SOMP_TLS / (source_count * n_repeat * J);
    mse_SBI(ii) = cur_SBI / (source_count * n_repeat * J);
    
    % YHQin 2017/12/20

    rmse_SBL(ii) = sqrt(mse_SBL(ii));
    rmse_SOMP_LS(ii) = sqrt(mse_SOMP_LS(ii));
    rmse_SOMP_TLS(ii) = sqrt(mse_SOMP_TLS(ii));
    rmse_SBI(ii) = sqrt(mse_SBI(ii));
      
    mse_an_fre = 0
    crb_fre = 0;
    for nn = 1:J
        ff = subband(nn);
             wavelength = c/ff;
            d = wavelength / 2;
            design = design_array_1d('coprime', [3 4], d);
            crb_fre = crb_fre + mean(diag(crb_uc_sto_1d(design, wavelength, doas, power_source, power_noise, snapshot_count)));
    end
    crb(ii) = crb_fre / J; 
 
    rmse_crb(ii) = sqrt(crb(ii));
end
fprintf('\n');

% plot
semilogy(snr, crb, '--', snr, mse_SBL, '-s', snr, mse_SOMP_LS, '-*', snr, mse_SOMP_TLS, '->', snr, mse_SBI, '-<');
xlabel('SNR (dB)'); ylabel('MSE / rad^2'); grid on;
legend( 'CRLB', 'SBL MSE', 'SOMP_LS MSE', 'SOMP_TLS MSE', 'OGSBI MSE' );


figure
% plot
semilogy(  snr, rmse_crb, '--', snr, rmse_SBL, '-s', snr, rmse_SOMP_LS, '-*', snr, rmse_SOMP_TLS, '->', snr, rmse_SBI, '-<');
xlabel('SNR (dB)'); ylabel('RMSE / rad'); grid on;
legend(  'CRLB',  'SBL', 'SOMP-LS', 'SOMP-TLS', 'OGSBI' );


