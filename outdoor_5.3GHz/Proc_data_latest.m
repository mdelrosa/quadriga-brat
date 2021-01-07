%%
clear all
CCM_UL =  [];
CCM_DL = [];
 count = 0;
for iii = 1:10 % recommend to process 10 files a time, otherwise, it will be out of memory.
    name = ['matnew',num2str(iii),'.mat'];
    load(name);
    for idx = 1:1000
        CCM_UL = [CCM_UL; CCM{1,idx}];
        CCM_DL = [CCM_DL; CCM{2,idx}];
    end
    iii
end
clear CCM 
CCM_DL_final = zeros(count, N*M,  No_subcarrier, 1);
CCM_UL_final = zeros(count, N*M,  No_subcarrier, 1);
for iii = 1:count
    CCM_DL_final(iii,:, :) = CCM_DL((iii-1)*N*M+1:iii*N*M, :);
    CCM_UL_final(iii,:, :) = CCM_UL((iii-1)*N*M+1:iii*N*M, :);
end
%  angular fft
MAP_2D_UL = (fft(CCM_DL_final, M, 2));
MAP_2D_DL =  (fft(CCM_UL_final, M, 2));

% delay fft
MAP_3D_UL = fftshift(fft(MAP_2D_UL, No_subcarrier, 3));
MAP_3D_DL = fftshift(fft(MAP_2D_DL, No_subcarrier, 3));

% discard reduntdunt information
list_remove = [[1:No_subcarrier/2-32] [No_subcarrier/2+33:No_subcarrier]];
MAP_3D_UL(:,:,list_remove) = [];
MAP_3D_DL(:,:,list_remove) = [];


X_AE = MAP_3D_DL;
X_AE_UL = MAP_3D_UL;

% % % iiii = 1;
% % % no_sc = 1024;
% % % VVVV(:,:) = X_AE(iiii,:,:);
% % % Map = abs(VVVV);
% % % AA = angle(VVVV);
% % % RR = real(VVVV);
% % % II = imag(VVVV);
% % % figure(112);
% % % subplot(411)
% % % imagesc(Map.');
% % % subplot(412)
% % % imagesc(AA.');
% % % subplot(413)
% % % imagesc(RR.');
% % % subplot(414)
% % % imagesc(II.');
% X_UL_norm = permute(X_UL_norm,[2 1 3]);
% X_AE = permute(X_AE,[2 1 3]);
% Y_AE = permute(Y_AE,[2 1 3]);

% iii =600;
% XX(:,:) = X_AE(iii,:,:); YY(:,:) = X_UL_norm(iii,:,:);
% figure(1);subplot(211);imagesc(abs(XX));subplot(212);imagesc(abs(YY))

save('all_data10.mat','X_AE','X_AE_UL')