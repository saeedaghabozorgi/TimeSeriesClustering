function data=data_normalize(data,method)

%   method     description
%   'var'      Variance is normalized to one (linear operation).
%   'range'    Values are normalized between [0,1] (linear operation).

data_Xold=data;
if strcmp(method,'range')
     data_min=min(data);
     data_max=max(data);
     data=(data-repmat(min(data),size(data,1),1))./(repmat(max(data),...
         size(data,1),1)-repmat(min(data),size(data,1),1));
 elseif strcmp(method,'var')
     data=(data-repmat(mean(data),size(data,1),1))./(repmat(std(data),size(data,1),1));
     data_mean=mean(data);
     data_std=std(data);
 else
     error('Unknown method given')
end
