clear all
H_ul = [];
H_dl = [];
for iii = 1:10
    name = ['all_data',num2str(iii),'.mat'];
    load(name);
    H_ul = cat(1, H_ul, X_AE_UL);
    H_dl = cat(1, H_dl, X_AE);
end
% save('H_QuaDRiGa','H_dl','H_ul')
% % MM1(:,:) = H_ul(1050,:,:);
% % MM2(:,:) = H_dl(1050,:,:);
% % figure(1)
% % subplot(211)
% % imagesc(abs(MM1))
% % subplot(212)
% % imagesc(abs(MM2))