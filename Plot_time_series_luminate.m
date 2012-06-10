function Plot_time_series_luminate(wplot,hplot,c,p,centroid,nor_traj,t_traj,clusterCount,showMems,luminate,fig_no)

if (wplot==0 || hplot==0)
    wplot=ceil(nthroot(clusterCount,2));
    if (wplot*(wplot-1))>=clusterCount
        hplot= wplot-1;
    else
        hplot=wplot;
    end
end
%--incremental Plot--------------------------------------------------------------
if clusterCount> wplot*hplot;
    rrr=wplot*hplot;
else
    rrr=clusterCount;
end

if isempty(t_traj)
    for i=1:length(nor_traj)
    t_traj{i}=[1:1:length(nor_traj{1})];
    end
end

clus=cell(1,clusterCount);
figure(fig_no);
clf(fig_no);
hold off;
cc=hsv(max(p));
for j=1:rrr
    memStr=[];
    clus{j}=find(c(:,1)==j);
    
    %combinedStr=strcat('(cluster:',num2str(j),')');
    ax = subplot(hplot,wplot,j);
    cla(ax);
    for i=1:length(clus{j})
        ind=clus{j}(i,1);
        if showMems==1
            memStr = strcat(memStr,',',num2str(ind));
        elseif showMems==2
            memStr = strcat(memStr,',',num2str(p(ind)));
        end
        color=cc(p(ind),:);
        %plot(t_traj{ind},nor_traj{ind},'color',color,'LineWidth',2)
        h = patch([0,t_traj{ind},t_traj{ind}(1,end)+1],[0,nor_traj{ind},0] ,'r');
        
        set(h,'EdgeColor',color);
        set(h,'facealpha',0);
        set(h,'edgealpha',luminate);
        axis([0 length(nor_traj{1}) -4 4]);
        hold on
    end
    if ~isempty(centroid)
        ti=(1:length(centroid{j}));
        xx=p(clus{j});
        centroid_color=cc(mode(xx),:);
        plot(ti,centroid{j},'color',centroid_color,'LineWidth',1)
    end
    combinedStr = strcat('Clus',num2str(j),'(',num2str(length(clus{j})),')');
    
    if showMems>0
        combinedStr = strcat(combinedStr,'-',memStr);
    end
    title(combinedStr);
end