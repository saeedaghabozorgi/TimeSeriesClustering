function [c,details]=clustering_anytime(D,A,k,p,nor_traj,alg)

start=length(D)-1;
for i=start:-1:1
    for j=i+1:length(D)
        if A(i,j)~=1
            D(i,j)=dis_dtw3(nor_traj{i},nor_traj{j},length(nor_traj{i}));
            D(j,i)=D(i,j);
        end
    end
    Dist=D;
    Dist(eye(size(Dist))~=0)=0;
    Dist=squareform(Dist);
    [c,~]= do_kMedoids_keogh(k,Dist);
    %evaluation
    [SSEP,SSEC,RI,ARI,purity,BCubed,ConEntropy,f_measure,jacard,FM,NMI,quality]= do_Evaluate(p,c,[],[],[]);
    details=[SSEP,SSEC,RI,ARI,purity,BCubed,ConEntropy,f_measure,jacard,FM,NMI,quality];
    disp(['  --> i:', num2str(i), ' | ','incremental quality:',num2str(quality)]);
end
end