function [SSEP,SSEC,RI,ARI,purity,ConEntropy,f_measure,jacard,FM,NMI,CSM,quality]= do_Evaluate(p,c,nor_traj,class_center,cluster_center)
SSEP=0;
SSEC=0;
RI=0;
ARI=0;
purity=0;
ConEntropy=0;
f_measure=0;
jacard=0;
FM=0;
NMI=0;
CSM=0;
if isempty(p)
    %    SSEP=Calculate_Cluster_SSE(p,nor_traj,class_center);
    if  ~isempty (cluster_center)  
        SSEC=Calculate_Cluster_SSE(c,nor_traj,cluster_center);
    end
else
    RI= Calculate_Cluster_Rand_statistic(c,p);
    ARI=Calculate_Cluster_AdjustedRandIndex(p,c);
    purity= Calculate_Cluster_Purity( c,p,2);
    %[BCubed,~,~]=Calculate_Cluster_BCubed_precision( c,p);
    entropy=Calculate_Cluster_entropy( c,p);
    ConEntropy=1-entropy;
    f_measure=Calculate_Cluster_F_measure( c,p);
    jacard=Calculate_Cluster_Jaccard(c,p);
    FM=Calculate_Cluster_FM(c,p);
    NMI=Calculate_Cluster_NMI(p,c);
    CSM=Calculate_Cluster_CSM(c,p);
end
quality=(RI+ARI+purity+ConEntropy+f_measure+jacard+FM+NMI+CSM)/9;
end