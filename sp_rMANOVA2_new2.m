function [Tk,pvals]=sp_rMANOVA2(data,GroupLabel,maxRepeat,drawplot)
% return test statistics and p-values for k group
%   robustified two-way MANOVA obtained 
%   through permutation test
%
% data: n by p by N
%    n: pooled sample size
%    p: dimension of r.v.
%    N: # of subjects
% maxRepeat: maximum permutation times
% GroupLabel: in a cell, 1 by 2, {n by 1}, 
% drawplot: true or false
%
% e.g.,
%   data=randn(8,3,1000);
%   GroupLabel={{'A';'A';'A';'A';'B';'B';'B';'B'},{'1';'1';'2';'2';'1';'1';'2';'2'}};
%   [Tk,pvals]=sp_rMANOVA2(data,GroupLabel,200,true);
%

% initial
[n p N]=size(data);
pvals=[];

% grouping
rGroupLabel=unique(GroupLabel{1});
row=length(rGroupLabel);
cGroupLabel=unique(GroupLabel{2});
col=length(cGroupLabel);
for i=1:row
  for j=1:col
    GroupDataId{i,j}=intersect(strmatch(rGroupLabel{i},GroupLabel{1}),...
                               strmatch(cGroupLabel{j},GroupLabel{2}));
    rcNumRep(i,j)=length(GroupDataId{i,j});                           
  end
end

% compute test stat.
Tk=rManovaStat2(data,GroupDataId);

% estimate null dist
CombId=getComb(rcNumRep,maxRepeat); % comb matrix nested in cell 
cTk{1,3}=[];
lenSimu=length(CombId);
display(['# of permutation = ',num2str(lenSimu)])
for iComb=1:lenSimu %!!
   mTk=rManovaStat2(data,CombId{iComb}); % 1 by 3 cell
   for iC=1:3
     cTk{iC}=[cTk{iC} mTk{iC}];
   end
end
% p-val after quantile normalization
for iC=1:3
  Tkhat{iC}=mean(sort(cTk{iC}),2);
  for i=1:N
    pvals{iC}(i,1)=mean(Tkhat{iC}>Tk{iC}(i));
  end
end 

% draw figure
for iC=1:3
if drawplot
figure
subplot(2,2,1)
hist(Tk{iC},100)
title('G in full range')
ax=axis;
text(ax(1)+.9*(ax(2)-ax(1)),ax(3)+.9*(ax(4)-ax(3)),'a','FontSize',14,'Color','r');
subplot(2,2,2)
qqplot(Tk{iC},Tkhat{iC})
xlabel('quantile of G');
ylabel('quantile of G^*');
title('Q-Q plot')
ax=axis;
text(ax(1)+.9*(ax(2)-ax(1)),ax(3)+.9*(ax(4)-ax(3)),'b','FontSize',14,'Color','r');
subplot(2,2,4)
hist(Tkhat{iC},100)
mx=axis;
title('G^*')
ax=axis;
text(ax(1)+.9*(ax(2)-ax(1)),ax(3)+.9*(ax(4)-ax(3)),'d','FontSize',14,'Color','r');
subplot(2,2,3)
hist(Tk{iC}(Tk{iC}<mx(2)),100)
axis(mx);
title('G in range as G^*')
ax=axis;
text(ax(1)+.9*(ax(2)-ax(1)),ax(3)+.9*(ax(4)-ax(3)),'c','FontSize',14,'Color','r');
drawnow
end
end

% an inner function computing all possible combination
%   for subgroup testing
function mCombId=getComb(rcNumRep,maxRepeat) 
[r,c]=size(rcNumRep);
mNumRep=reshape(rcNumRep,1,r*c);
if ~any(mNumRep>1)
    error('Not enough replicates!')
end
mN=sum(mNumRep);
mCombId=multinomiald([1:mN],mNumRep,maxRepeat);
lenCombId=length(mCombId);
display([num2str(lenCombId) ' available permutations...'])
for i=1:length(mCombId)
    mCombId{i}=reshape(mCombId{i},r,c);
end
   
% an inner function computing test stat. for k group comparison
function mTk=rManovaStat2(Data, GroupDataId)
[row col]=size(GroupDataId); %
[n p N]=size(Data); % n by p by N

for i=1:row
  rGroupDataId{i}=unionc(GroupDataId(i,:));
  rNumRep(i)=length(rGroupDataId{i});
  for j=1:col
    rcNumRep(i,j)=length(GroupDataId{i,j});
  end
end
for j=1:col
    cGroupDataId{j}=unionc(GroupDataId(:,j));
    cNumRep(j)=length(cGroupDataId{j});
end                % summarize number of replicates 

% Nomis = [tr(Ba) tr(Bb) tr(Bab)]
grandMedian=median(Data,1); % 1 by p by N
for i=1:row
  rGroupMedian(i,1:p,1:N)=median(Data(rGroupDataId{i},:,:),1); % 1 by p by N
  rGroupMedianS(i,1:p,1:N)=[rGroupMedian(i,:,:)-grandMedian].^2*rNumRep(i);
end % use the same notation for computing the squared error
Nomi{1}=reshape(sum(sum(rGroupMedianS,1),2),1,N); % 1 by N trBa
clear rGroupMedianS;

for j=1:col
  cGroupMedian(j,1:p,1:N)=median(Data(cGroupDataId{j},:,:),1); % 1 by p by N
  cGroupMedianS(j,1:p,1:N)=[cGroupMedian(j,:,:)-grandMedian].^2*cNumRep(j);
end % use the same notation for computing the squared error
Nomi{2}=reshape(sum(sum(cGroupMedianS,1),2),1,N); % 1 by N trBb
clear cGroupMedianS;

ct=0;
for i=1:row
  for j=1:col
    ct=ct+1;
    rcGroupMedian(ct,1:p,1:N)=median(Data(GroupDataId{i,j},:,:),1); % 1 by p by N
    rcGroupMedian(ct,:,:)=[rcGroupMedian(ct,:,:)-rGroupMedian(i,:,:)-...
                           cGroupMedian(j,:,:)+grandMedian].^2*rcNumRep(i,j); % use same notation!
  end
end
Nomi{3}=reshape(sum(sum(rcGroupMedian,1),2),1,N); % 1 by N trBab

% Deno
trW=0;
for i=1:row
  for j=1:col
    if rcNumRep(i,j)>1
        ijGroupMad=mad(Data(GroupDataId{i,j},:,:),1,1); % 1 by p by N
        trW=trW+rcNumRep(i,j)*ijGroupMad.^2; % 1 by p by N
    end
  end
end
if trW==0
    Deno=1;
else
    Deno=reshape(sum(trW,2),1,N); 
end

% stats
for iC=1:3
  mTk{iC}=[Nomi{iC}./Deno]';
end

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