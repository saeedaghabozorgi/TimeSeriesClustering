function dataset=CBF_Generator(N)
    %N=200;
    c=CBF('CYLINDER',N);
    c=[ones(1,N)',c];

    
    b=CBF('BELL',N);
    b=[ones(1,N)'*2,b];
 
    
    f=CBF('FUNNEL',N);
    f=[ones(1,N)'*3,f];
    dataset=[c;b;f];
%         dlmwrite('CBF_TRAIN.txt',c ,'-append','delimiter', '\t','newline','pc');
%            dlmwrite('CBF_TRAIN.txt',b ,'-append','delimiter', '\t','newline','pc');
%     dlmwrite('CBF_TRAIN.txt',f ,'-append','delimiter', '\t','newline','pc');
end

function TS=CBF(shape,N)

a= randi([16,32],N,1);
b= randi([32,96],N,1)+ a;
n=rand(N,1);

for j=1:N
    for t=1:128
        x=xab(a(j,1),b(j,1),t);
        if strmatch(shape,'CYLINDER')
            TS(j,t)=(6+rand)*x+rand;
        elseif strmatch(shape,'BELL')
            TS(j,t)=(6+rand)*x*(t - a(j,1)) ./ (b(j,1) - a(j,1))+rand;
        elseif strmatch(shape,'FUNNEL')
            TS(j,t)=(6+rand)*x*(b(j,1) - t) ./ (b(j,1) - a(j,1))+rand;
        end
    end
end
end

function x=xab(a,b,t)
if ( (t>=a) && (t <=b) )
    x= 1.0;
else
    x= 0.0;
end
end
