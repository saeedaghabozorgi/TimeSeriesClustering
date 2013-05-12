function Fourier_representation()
load leleccum; 
 a = leleccum;
 a=a(1:365);
a = (a-mean(a))/std(a);         % z-normalization

fa = fft(a);

maxInd = ceil(length(a)/2);     % until the middle
N = length(a);                  

energy = zeros(maxInd-1, 1);
E = sum(a.^2);                  % energy of a

for ind=2:maxInd,
    
    fa_N = fa;                  % copy fourier
    fa_N(ind+1:N-ind+1) = 0;    % zero out unused
    r = real(ifft(fa_N));       % reconstruction
    
    plot(r, 'r','LineWidth',2); hold on;
    plot(a,'k');
    title(['Reconstruction using ' num2str(ind-1) 'coefficients']);
    set(gca,'plotboxaspectratio', [3 1 1]);
    axis tight
    pause;		      % wait for key
    cla;	         	      % clear axis
end
