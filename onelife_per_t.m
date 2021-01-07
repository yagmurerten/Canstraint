function [deathcause,mat_time,addedbudget,new_bs,N]=onelife_per_t(strategy,mort,...
    cancerdanger,maturitycriterion,withextmort,celldeath,N,mat_time,age,nofonco)

% this is pretty similar to the onelife.m that we used previously for optimising the
% ontogenetic strategy, but now evaluated one step at a time for every
% individual

% strategy: 
% 1 - probability of asymmetric division; P

% 2 - conditional probability of differentiation GIVEN there was a
%   symmetric division (= given the option of asymmetric division was 
%   not used); Q

% 3 - telomere length: how many telomere layers does the organism
%   use, dim(1) of the N matrix; Hayflick limits, H

% 4: apoptosis threshold; A

% 5: apoptosis percent; S

% 6: number of differentiation levels; T

% 7: division propensities: relative use of different 
%   cells for division, compared with the previous level; X

%%%%%%%%%

% mort(1): per-cell death rate (v)
% mort(2) whole-organism death rate, both expressed per time step (mu)

%%%%%%%%%

% maturitycriterion is the target number of fully differentiated cells

%%%%%%%%%

% N will contain:
% rows = H-telomeres (1 = longest, strategy(3) = so short you won't divide 
%   any more)
% columns = K-damage level (1 = no damage, nof_oncosteps = so damaged you're
%   recorded dead because of cancer [assuming 1 such cell is considered 
%   bad enough])
% layers = T-differentiation level (1 = stem cell, T = fully differentiated 
%   cell)


apoptosis_thr=strategy(4); % A
apoptosis_percent=strategy(5); % S

% this implements noise around damage detection
noisesd=0.5;
k=0.05; % parameter for the growth trajectory
minpercenttissue=0.8; 
% the organism needs to have at least 80% of its target terminally differentiated cells
% after maturation - related to death cause #5
% intituitvely, we would not hit to this limit as before this, the organism
% would die due to lack of stem cells
% unless the apoptosis rate or baseline cell death rate is perhaps really high    T=focal_str(6); % differentiation layers
T=strategy(6); % differentiation layers
H=strategy(3);
K=nofonco;
divisionpropensities=strategy(7).^linspace(1,T-1,T-1);
% update your info on dividibles
dividibles=[N(1,1,1) zeros(1,T-2)]; % these cells can divide

for i=1:length(dividibles)
    dividibles(i)=sum(sum(N(1:H-1,1:K,i)));
end

extmort=mort(2);

% we check maturity with the previous information, instead of the body
% size, because the organism might have been mature but then lost some of
% its cells and went back slightly down to a smaller body size
mature=mat_time>0;alive=1;
deathcause=0;

if withextmort % if we want extrinsic mortality to be applied stochastically througout sim. time
    if rand(1)<extmort
        % whole organism dies of extrinsic causes
        alive=0; deathcause=1; 
    end
end
if ~mature && alive % not mature yet, let's grow
     % this kind of leads to vanB. growth        
    tissue=sum(sum(N(:,:,end)));
    % how many cells are needed
    total_divisions=ceil(k*(maturitycriterion-tissue));
    % if the stem cells are depleted before maturation, it is
    % dead --> it won't mature anyway
    if sum(dividibles)==0 
        alive=0; deathcause=2;
    end
    if alive && (total_divisions > 0) 
        N=divisionsandmutations(N,strategy,total_divisions,...
            dividibles,divisionpropensities,cancerdanger,H,K);
    end
    % now check if maturity has been reached
    tissue=sum(sum(N(:,:,end)));
    % maturitycriterion=target number of differentiated cells that
    % means the organism has matured
    if tissue>=maturitycriterion
        % now it is mature!
        mature=1;mat_time=age;
        % update your info on total body size
    end
% what to do when mature, slightly changed from the previous onelife 
% as here one can die from ext mort stochastically
elseif mature && alive 
    tissue=sum(sum(N(:,:,end)));
    % counteract the loss of cells that may have happened
    total_divisions=ceil(k*(maturitycriterion-tissue));     
    if sum(dividibles)==0 
        alive=0; deathcause=8;
    elseif total_divisions > 0 % do this part only if you need divisions
        N=divisionsandmutations(N,strategy,total_divisions,...
            dividibles,divisionpropensities,cancerdanger,H,K);
    end
end

% if you managed to survive the first two death causes:

% now some cells also die randomly (baseline cell death, sort of like 
% cell extrinsic mortality)
if celldeath && alive
    deadcells=binornd(N,mort(1)); N=N-deadcells;
    N(H,:,:)=0;
    killed=apoptosis(N,noisesd,apoptosis_thr,apoptosis_percent,H,K,T);
    N=N-killed;
    N(N<0)=0;
end

% now check if you've got cancer, even 1 cell in
% the O (maximum onco-state) category kills you!
cancercells=sum(sum(N(:,end,:)));
if cancercells>0 && alive
    alive=0; deathcause=3;
end

% update your info on total body size and tissue size
bodysize=sum(sum(sum(N,3)));

% death causes related to cell numbers = losing all cells, 
% imbalance between differentiated cells and stem cells

if bodysize==0 && alive
    alive=0; deathcause=4; 
end
% the organism actually vanished (this cause should only be applied if not mature yet)

% Then once you mature and reach your 'final' body size, 
% you cannot go below a certain percent of that; i.e. you need a
% certain number of tissue cells to maintain your functions

if alive && mature && (tissue < minpercenttissue*maturitycriterion)
    alive=0; deathcause=5;
end  

bodysize=sum(sum(sum(N,3)));tissue=sum(sum(N(:,:,end)));new_bs=tissue;

% calculating the budget for reproduction
if mature
    addedbudget=tissue/bodysize;
else
    addedbudget=0;
end