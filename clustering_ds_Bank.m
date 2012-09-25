function clustering_ds_Bank
% for demo of similar time series, I choose 86,125,904

ab=[];

k=30;
fprintf('Importing data .. \n\r');
nor_traj=[];
[nor_traj,t_traj]=Import_Data_Bank(1,700,'sequence',1,360);


  dist_mtx_DTW_file=['..\dataset Bank\Bank_1000_dismat_DTW.mat'];
    load(dist_mtx_DTW_file,'dismat');
    dist_mtx_DTW=dismat;
    
% 
% 
%  hold off;
%            str=demo_sax(nor_traj{73},alphabet_size,compression_ratio);
%  hold off;
%            str=demo_sax(nor_traj{253},alphabet_size,compression_ratio);
%  hold off;
%            str=demo_sax(nor_traj{352},alphabet_size,compression_ratio);

 details=clustering_Hybrid_3Level_bank(nor_traj,20,[],dist_mtx_DTW,[]);

end
