function population=createpop(initial_pop_str,populationsize,targetsize,init_var,evol_bounds,w_bs,w_cancer)
    length_strategy=length(initial_pop_str);
    initial_pop_str_dip=ones(1,length_strategy*6);
    
    index=1;    
    for i=1:length_strategy
      initial_pop_str_dip(1,index:index+5)=initial_pop_str(i);      
      index=index+6;
    end
    % now we need to add the other alleles 
    init_pop=ones(populationsize,length_strategy*6).*initial_pop_str_dip;

    % if traits are evolving add some variation around the initial set of
    % strategies
    if ~isempty(evol_bounds)
        sz=size(init_pop);
        probs=evol_bounds(:,1)+rand(length_strategy*6,populationsize).*(evol_bounds(:,2)-evol_bounds(:,1));
        % how far the resident strategy is from the new random strategy
        difference=probs'-init_pop;
        pd=makedist('Normal','mu',0,'sigma',init_var);pd_t=truncate(pd,0,1);
        % move by a random percent
        percentchange=random(pd_t,sz);
        newpop=init_pop+difference.*percentchange; 
    else
        newpop=init_pop;
    end
    % body size always evolves in our similations (so far), so we start
    % with some standing variation in that as well
    targetbodysize=exp(log(targetsize)+normrnd(0,init_var,[populationsize,6])); 
    pop=horzcat(newpop,targetbodysize);
    % 0s alleles from the father, 1s alleles from the mother
    indices=repmat([0 0 0 1 1 1],1,8);
    mean_parent_val=(pop(:,indices==0)+pop(:,indices==1))/2;    
    sex=vertcat(zeros(populationsize/2,1),ones(populationsize/2,1));
    % male-female-shared
    phenotype_cancer=zeros(populationsize,length_strategy);
    phenotype_cancer(sex==0,:)=mean_parent_val(sex==0,1:3:19)*w_cancer(2)+mean_parent_val(sex==0,3:3:21)*w_cancer(1);
    phenotype_cancer(sex==1,:)=mean_parent_val(sex==1,2:3:20)*w_cancer(2)+mean_parent_val(sex==1,3:3:21)*w_cancer(1);
    
    phenotype_bs(sex==0,:)=mean_parent_val(sex==0,22)*w_bs(2)+mean_parent_val(sex==0,24)*w_bs(1);
    phenotype_bs(sex==1,:)=mean_parent_val(sex==1,23)*w_bs(2)+mean_parent_val(sex==1,24)*w_bs(1);
    phenotype=horzcat(phenotype_cancer,phenotype_bs);
    % round the things that need to be integers (Hayflick limit and
    % differentiation levels)
    % trait 3=Hayflick
    % trait 6=dif layers
    phenotype(:,3)=round(phenotype(:,3));
    phenotype(:,6)=round(phenotype(:,6));    
    % at the beginning all these things are zero (everyone is a zygote)
    bs_age_budget_mattime_noffspring=zeros(populationsize,5);
    population=horzcat(pop,phenotype,sex,bs_age_budget_mattime_noffspring);
    population =  population(randperm(end),:);
end