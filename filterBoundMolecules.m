function ID_bound = filterBoundMolecules(x,N,Alimit)

%clear;

% N=30;
% dx=randn(1,1e6);
% %x=cumsum(dx);
% x=dx;

xN=zeros(N,length(x)-N+1);
for i=1:N
    xN(i,:)=x(i:end-N+i);
end

A=std(xN,1,1)./std(diff(xN,1,1),1,1);

ID_bound=find(A<Alimit);
if length(ID_bound)>0
ID_bound=ID_bound+[0:N-1]'*ones(1,length(ID_bound));
ID_bound=unique(reshape(ID_bound,1,[]));
end


% return
% clear;
% A=[1,10,100];
% N=[1e3,1e4,1e5,1e6];
% for i=1:3
%     for j=1:4
%         a=A(i)*randn(1,N(j));
%         STDa(i,j)=std(a);
%         STDda(i,j)=std(diff(a));
%     end
% end
% 
% 
% return


