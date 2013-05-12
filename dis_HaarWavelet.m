function HaarWavelet_distance()

load leleccum; 
 s = leleccum;
 a=s(1:365);
 b=s(101:465);


figure;
hold off;
a = (a-mean(a))/std(a);         % z-normalization
b = (b-mean(b))/std(b); 
 plot(a); 
 hold on;
 plot(b);
 
 
maxlevels = wmaxlev(length(a),'haar');
[Ca, La] = wavedec(a,maxlevels,'haar');
[Cb, Lb] = wavedec(b,maxlevels,'haar');


% Top-20 coefficient reconstruction
[Ca_sorted, Ca_sortind] = sort(Ca);
Ca_top20 = Ca; 
Ca_top20(Ca_sortind(21:end-20)) = 0;
a_top20 = waverec(Ca_top20,La,'haar');

 [Ca_sorted, Cb_sortind] = sort(Cb);
 Cb_top20 = Cb; 
 Cb_top20(Cb_sortind(21:end-20)) = 0;
 b_top20 = waverec(Cb_top20,Lb,'haar');

 
 d1=dis_euclidean(a,b);
 d2=dis_euclidean(a_top20,b_top20);
%  
%  x1=Cb(Cb_sortind(1:20));
%  x2=Cb(Cb_sortind(end-20:end));
%  x=[x1, x2]
%  y1=Ca(Ca_sortind(1:20));
%  y2= Ca(Ca_sortind(end-20:end));
%  y=[y1,y2];

%  d3=dis_euclidean(x,y);
 
 hold on;
plot(a_top20, 'r');
 hold on;
 plot(b_top20, 'r');
end