function [SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality]= do_Evaluate(p,c,nor_traj,class_center,cluster_center)

if cluster_center~=-1 SSEP=Calculate_Cluster_SSE(p,nor_traj,class_center); else SSEP=0; end
if class_center~=-1 SSEC=Calculate_Cluster_SSE(c,nor_traj,cluster_center); else SSEC=0; end
RI= Calculate_Cluster_Rand_statistic( c,p);
purity= Calculate_Cluster_Purity( c,p,2);
[BCubed,~,~]=Calculate_Cluster_BCubed_precision( c,p);
entropy=Calculate_cluster_entropy( c,p);
ConEntropy=1-entropy;
f_measure=Calculate_Cluster_F_measure( c,p);
jacard=Calculate_Cluster_Jaccard(c,p);
FM=Calculate_Cluster_FM(c,p);
quality=(purity+BCubed+ConEntropy+f_measure+RI)/5;

end