function [nor_traj,t_traj]=NormalizeTS(data)

nor_traj=[]; % normalized orginal traj
t_traj=[];
for z=1:size(data,1);
    a=data(z,:);
    nor_traj{z}=(a-mean(a))/std(a);
    t_traj{z}=(1:1:size(data,2));
end