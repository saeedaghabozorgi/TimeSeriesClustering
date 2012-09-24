function clustering_ds_Bank
% for demo of similar time series, I choose 86,125,904

ab=[];

k=30;
fprintf('Importing data .. \n\r');
nor_traj=[];
[nor_traj,t_traj]=Import_Data_Bank(1,200,'sequence',1,360);

% 
% 
%  hold off;
%            str=demo_sax(nor_traj{73},alphabet_size,compression_ratio);
%  hold off;
%            str=demo_sax(nor_traj{253},alphabet_size,compression_ratio);
%  hold off;
%            str=demo_sax(nor_traj{352},alphabet_size,compression_ratio);

 details=clustering_Hybrid_3Level(nor_traj,10,[],[]);

end
