function [Tk,pvals]=sp_rMANOVA1(data,Label,maxRepeat,tm,drawplot)
% return test statistics and p-values for k group
%   robustified one-way MANOVA obtained 
%   through permutation test
%
% data: n by p by N
%    n: pooled sample size
%    p: dimension of r.v.
%    N: # of subjects
% maxRepeat: maximum permutation times
% tm: 'trace' or 'median'
% Label: in a cell structure
% drawplot: true or false
%
% e.g.,
%   data=randn(4,3,1000);
%   Label={'1';'1';'2';'2'};
%   [Tk,pvals]=sp_rMANOVA1(data,Label,200,'trace',true);
%

% initial
[sampleSize mNumProbe mNumProbeSet]=size(data);
NumChip=sampleSize;
pvals=[];
% grouping
Label0=Label;
clear Label;
for i=1:length(Label0)
    if isnumeric(Label0{i})
        Label(i,1)=cellstr(num2str(Label0{i}));
    else
        Label(i,1)=Label0(i);
    end
end %  here to deal with numerical labels
GroupLabel=unique(Label);
NumGroup=length(GroupLabel);
for i=1:NumGroup;
    GroupDataId{i,1}=strmatch(GroupLabel(i),Label);
    NumRep(i)=length(GroupDataId{i,1});
end % GroupId in cell, # of Reps. each group

% compute test stat.
Tk=rManovaStat1(data,GroupDataId,tm);

% estimate null dist
CombId=multinomiald(1:sampleSize,NumRep,maxRepeat);  % in cell 
cTk=[];
lenSimu=length(CombId);
display(['# of permutation = ',num2str(lenSimu)])
for iComb=1:lenSimu %!!
   mTk=rManovaStat1(data,CombId{iComb},tm);
   cTk=[cTk mTk];
end

% p-val after quantile normalization
Tkhat=mean(sort(cTk),2);
for i=1:mNumProbeSet
    pvals(i,1)=mean(Tkhat>Tk(i));
end

% draw figure
if drawplot
figure
subplot(2,2,1)
hist(Tk,100)
colormap([1 1 1])        
title('G in full range')
ax=axis;
text(ax(1)+.9*(ax(2)-ax(1)),ax(3)+.9*(ax(4)-ax(3)),'a','FontSize',14,'Color','r');
subplot(2,2,2)
hist(Tkhat,100)
colormap([1 1 1])        
mx=axis;
title('G^*')
ax=axis;
text(ax(1)+.9*(ax(2)-ax(1)),ax(3)+.9*(ax(4)-ax(3)),'b','FontSize',14,'Color','r');
subplot(2,2,4)
qqplot(Tk,Tkhat)
xlabel('quantile of G');
ylabel('quantile of G^*');
title('Q-Q plot')
ax=axis;
text(ax(1)+.9*(ax(2)-ax(1)),ax(3)+.9*(ax(4)-ax(3)),'d','FontSize',14,'Color','r');
subplot(2,2,3)
hist(Tk(Tk<mx(2)),100)
colormap([1 1 1])        
axis(mx);
title('G in range as G^*')
ax=axis;
text(ax(1)+.9*(ax(2)-ax(1)),ax(3)+.9*(ax(4)-ax(3)),'c','FontSize',14,'Color','r');
drawnow
end

   
% an inner function computing test stat. for k group comparison
function mTk=rManovaStat1(Data, GroupDataId,tm)
mNumGroup=length(GroupDataId);
[s1 mNumProbe mNumProbeSet]=size(Data); % n by p by N
for i=1:mNumGroup;
    mNumRep(i)=length(GroupDataId{i,1});
end
% Nomi = tr(B)
grandMedian=median(Data,1); % 1 by p by N
for i=1:mNumGroup
   GroupMedian(i,1:mNumProbe,1:mNumProbeSet)=median(Data(GroupDataId{i,1},:,:),1);
   GroupMedian(i,:,:)=[GroupMedian(i,:,:)-grandMedian].^2*mNumRep(i);
end % use the same notation for computing the squared error k by p by N
switch tm
    case 'trace'
  Nomi=reshape(sum(sum(GroupMedian,1),2),1,mNumProbeSet); % 1 by N
    case 'median'
  Nomi=reshape(median(sum(GroupMedian,1),2),1,mNumProbeSet); % 1 by N
end
% Deno
trW=0;
for i=1:mNumGroup
    if mNumRep(i)>1
        iGroupMad=mad(Data(GroupDataId{i,1},:,:),1,1); % 1 by p by N
        trW=trW+mNumRep(i)*iGroupMad.^2; % 1 by p by N
    end
end
if trW==0
    Deno=1;
else
    
  switch tm
    case 'trace'
      Deno=reshape(sum(trW,2),1,mNumProbeSet); 
    case 'median'
      Deno=reshape(median(trW,2),1,mNumProbeSet); 
   end
end
% ratio
mTk=[Nomi./Deno]';

% an inner function computing distinct multinomial combinations  
function [mM,NumEle]=multinomiald(vec,pattern,N)
% return distinct multinomial combinations in cells
%   with number of outcome less than N,
%   NumEle is the number of possible distinct combinations
% e.g.,
%   vec=[1:12];
%   pattern=[3 4 5];
%   N=10000;
%   [mM,NumEle]=multinomiald(vec,pattern,N);


% validity check
n=length(vec);
if n~=sum(pattern)
    error('Mismatched pattern!')
end

% compute possible number of combinations
if nargout==2
    NumEle=factorial(n)/prod(factorial(pattern));
    for i=1:floor(n/2)
        iNumTie=sum(pattern==i);
        NumEle=NumEle/factorial(iNumTie);
    end
end

% main
mM=multinomiald0(vec,pattern,N);
r=length(mM); % if more than N, cut randomly
if r>N
  rPerm=randperm(r);
  mM=mM(rPerm(1:N),:);
end

function mM=multinomiald0(vec,pattern,N)
% return distinct multinomial combinations in cells
%   with number of outcome less than N+10% more

% initial
n=length(vec);
vec=vec(randperm(n)); % permute original vector
pattern=sort(pattern); % order pattern 
k=length(pattern);
mM=[];

% main 
switch k
    
  case 1 % just the whole vector 
    mM={vec};
    
  case 2 % binomial partition
    subpattern=pattern(pattern>1);
    if isempty(subpattern) % patttern=[1 1]
        mM={[{vec(1)};{vec(2)}]};
    else % pattern ~= [1 1]
      if pattern(1)~=pattern(2) % unbalanced
          sub1=nchoosekb(vec,pattern(1),N);
      else % pattern(1)==pattern(2), balanced 
          sub1=nchoosekb(vec(2:end),pattern(1)-1,N);
          [numrow numcol]=size(sub1);
          sub1=[repmat(vec(1),numrow,1) sub1];
      end
      [numrow numcol]=size(sub1); % numrow<=N
      for i=1:numrow
         mM{i,1}=[{sub1(i,:)};...
                 {setdiff(vec,sub1(i,:))}];
      end
    end
    
  otherwise % k>=3,4,5,... multinomial cases

    subpattern=pattern(pattern>1); % non-1 terms
    subk=length(subpattern); % # of non-1 terms
    
    if subk==0 % pattern=[1 1 ... 1]
      
      for i=1:length(vec)
        mM{1,1}{i,1}=vec(i);
      end
      
    else %subk>0 pattern not all 1s
        
      numOnes=length(find(pattern==1));
      if numOnes>0 % some 1s
        sub1=nchoosekb(vec,numOnes,N);
        % grow with 1s
        [numrow numcol]=size(sub1);
        for i=1:numrow
          remain=setdiff(vec,sub1(i,:));
          sub2=multinomiald0(remain,subpattern,ceil(N/numrow)); 
          for j=1:length(sub2)
            subOnes=num2cell(sub1(i,:)');
            extra{j,1}=[subOnes;...
                        sub2{j,1}];
          end
          mM=[mM;extra];
        end % of grow with 1s
        
      else % no 1s
        
        % deal with >1 terms
        
        [up cp]=fcount(subpattern);
        
        singleton=find(cp==1);
        
        if ~isempty(singleton) % exist singleton pattern, e.g., [2 3 4]
           sub1=nchoosekb(vec,subpattern(up(singleton(1))),N);
           % grow with singleton pattern
          [numrow numcol]=size(sub1); % numrow<=N
          for i=1:numrow
            remain=setdiff(vec,sub1(i,:));
            sub2=multinomiald0(remain,subpattern(subpattern~=up(singleton(1))),ceil(N/numrow));
            for j=1:length(sub2)
              extra{j,1}=[{sub1(i,:)};...
                          sub2{j,1}];
            end
            mM=[mM;extra];
          end % of grow without repeated pattern
          
        else % no singleton, all repeated pattern, e.g., [2 2 3 3]
            
          % two cases: balanced or unbalanced
          if length(up)==1 % balanced
            % first balanced subgroup   
            sub1=nchoosekb(vec(2:end),subpattern(1)-1,N); 
            [numrow numcol]=size(sub1); % numrow<=N
            sub1=[repmat(vec(1),numrow,1) sub1];
            % grow with first group
            for i=1:numrow
              remain=setdiff(vec,sub1(i,:));
              sub2=multinomiald0(remain,subpattern(2:subk),ceil(N/numrow));
              for j=1:length(sub2)
                extra{j,1}=[{sub1(i,:)};sub2{j,1}];
              end
              mM=[mM;extra];
            end % of grow with first balanced subgroup  
          else % unbalanced              
          % start with the first repeated pattern  
            subv1=nchoosekb(vec,up(1)*cp(1),N);  
            [numSubv1, temp]=size(subv1);
            for i=1:numSubv1
              sub1=multinomiald0(subv1(i,:),repmat(up(1),1,cp(1)),ceil(N/numSubv1));
              [numsub1,numcol]=size(sub1); % actual length of sub1
              remain=setdiff(vec,subv1(i,:));
              remainsubpattern=subpattern(subpattern~=up(1));
              sub2=multinomiald0(remain,remainsubpattern,ceil(N/numSubv1/numsub1));
              [numsub2,numcol]=size(sub2); % actual length of sub2
              cj=1;
              for j1=1:numsub1
                for j2=1:numsub2
                   extra{cj,1}=[sub1{j1,1};sub2{j2,1}];
                   cj=cj+1;
                end
              end
              mM=[mM;extra];
            end % of grow with repeated pattern
          end % balanced or unbalanced
                    
        end % of non-singleton case
                         
      end % of some of 1s and no 1s 
      
    end % of if pattern=[1 1 ... 1]

end % of switch k

% an inner function computing bounded nchoosek
function Y=nchoosekb(vecn,k,b)
% return a matrix with distinct rows of k numbers
%   from vector of 1 by n without replacement
%   the number of rows is not great than b

Y=nchoosekb0(vecn,k,b);
[r c]=size(Y);
if r>b
  rPerm=randperm(r);
  Y=Y(rPerm(1:b),:);
end

% an inner function of nchoosekb
function Y=nchoosekb0(vecn,k,b)
Y=[];
n=length(vecn);
vecn=vecn(randperm(n));
B=nchoosek(n,k);
if B<=b
    Y=nchoosek(vecn,k);
else % B>b
    if n<=15 || k<=2 % can be handled directly
        Y=nchoosek(vecn,k);
        rPerm=randperm(B);
        Y=Y(rPerm(1:b),:);
    else % n>15 and k>2 may not be handled directly
        Y0=nchoosek(vecn,2);
        lenY0=length(Y0);
        if lenY0<=b
          for i=1:lenY0
            Y1=nchoosekb0(setdiff(vecn,Y0(i,:)),k-2,ceil(b/lenY0));
            [r c]=size(Y1);
            Y=[Y;[repmat(Y0(i,:),r,1) Y1]]; % attach
          end
        else % lenY0>b
          rPerm=randperm(lenY0);
          Y0=Y0(rPerm(1:b),:);  
          for i=1:b
            Y1=setdiff(vecn,Y0(i,:));
            rPerm=randperm(n-2);
            Y(i,:)=[Y0(i,:) Y1(rPerm(1:k-2))]; % attach
          end
        end
    end
end
    
% an inner function of frequency count
function [u c]=fcount(x)
% return unique elements and corresponding frequencies
u=unique(x);
c=[];
for i=1:length(u)
    c(i)=sum(x==u(i));
end                