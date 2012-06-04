function dis_DWT
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
subplot(5,2,1); 
plot(tt,s(tt),'r'); 
title('Orginal Time Series'); 
    xlabel('Day')
    ylabel('Balance')
for i = 1:4, 
    % Reconstruct the details using the decomposition structure. 
    D(2*i,:) = wrcoef('a',c,l,w,2*i);
    subplot(5,2,2*i+1);
    plot(tt,D(2*i,tt),'g');
    combinedStr = strcat('Reconstructed coefficients, level:',num2str(2*i));
    title(combinedStr);

    
    ca5 = appcoef(c,l,w,2*i);
    t=1:length(ca5);
    subplot(5,2,2*i+2);
    plot(t,ca5,'g');
    combinedStr = strcat('Approximation coefficients, level:',num2str(2*i));
    title(combinedStr);  
end