function DimensionReduction
 load leleccum; 
 s = leleccum;
 s=s(1:365);

% Decompose the signal s at level 5 using the wavelet db3. 
w = 'haar'; 
ml=wmaxlev(size(s),'haar');
[c,l] = wavedec(s,ml,w);
% ca3 = appcoef(c,l,w,8);
% tt = 1:length(ca3); 
% t=1:365;
% plot(t,s(t),'r'); 
%plot(tt,ca3(tt),'r'); 
% 
% D(1,:)=wrcoef('d',c,l,w,20);
% tt = 1+100:length(s)-100; 
% 
% subplot(3,1,1); 
% plot(tt,s(tt),'r'); 
% 
% subplot(3,1,2);
% plot(tt,D(1,tt),'g'); 




% Avoid edge effects by suppressing edge values and plot.
tt = 1:length(s); 
subplot(9,2,1); 
plot(tt,s(tt),'r'); 
title('Orginal Time Series'); 
for i = 1:ml, 
    % Reconstruct the details using the decomposition structure. 
    D(i,:) = wrcoef('a',c,l,w,i);
    subplot(9,2,2*i+1);
    plot(tt,D(i,tt),'g');
    combinedStr = strcat('reconstructed coefficients, level:',num2str(i));
    title(combinedStr);
    
    ca5 = appcoef(c,l,w,i);
    t=1:length(ca5);
    subplot(9,2,2*i+2);
    plot(t,ca5,'g');
    combinedStr = strcat('approximation coefficients, level:',num2str(i));
    title(combinedStr);  
end