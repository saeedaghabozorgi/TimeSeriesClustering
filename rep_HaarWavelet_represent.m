function HaarWavelet_represent()

load leleccum; 
 s = leleccum;
 a=s(1:365);


a = (a-mean(a))/std(a);         % z-normalization

 plot(a); 
maxlevels = wmaxlev(length(a),'haar');
[Ca, La] = wavedec(a,maxlevels,'haar');

% Plot coefficients and MRA
for level = 1:maxlevels
    cla;
    subplot(3,1,1);
    plot(detcoef(Ca,La,level)); axis tight;
    title(sprintf('Wavelet coefficients (detail) – Level %d',level));
    subplot(3,1,2);
    plot(wrcoef('d',Ca,La,'haar',level)); axis tight;
    title(sprintf('MRA-  reconstructed coefficients – Level %d',level));
    subplot(3,1,3);
    plot(appcoef(Ca,La,'haar',level)); axis tight;
    title(sprintf('Approx coefficients of signal – Level %d',level));
    
    pause;
end

% Top-20 coefficient reconstruction
[Ca_sorted, Ca_sortind] = sort(Ca);
Ca_top20 = Ca; Ca_top20(Ca_sortind(1:end-19)) = 0;
a_top20 = waverec(Ca_top20,La,'haar');
figure; hold on;
plot(a, 'b'); plot(a_top20, 'r');
end