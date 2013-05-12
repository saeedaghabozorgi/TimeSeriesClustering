function s=Calculate_intra_cluster(c,clust,nortraj)
mems=nortraj((c==clust));
n=length(mems);
s=0;
if n==1
    return;
else
    
    dis=Mtx_Euclid_Distance(mems,mems,'same');
    
    s=(sum(sum(dis)))/(n*(n-1)/2);
end
end