function [population,Ns,offspringproduced,parents]=reproduce(population,type_f,type_m,Ns,nofonco,init_var,evol_bounds,...
    index_sex,index_budget,index_csize,index_mat,index_offspringnum,index_age,w_bs,w_cancer)

    emptyspots=find(population(:,1)==0); % if there are spaces to make offspring
    totalspots=length(emptyspots);
    parents=zeros(totalspots*2,4);
    offspringproduced=0;
    if type_f==string('budget') 
        % the relevant trait for competition is budget
        femcompetitiveness=population(:,index_budget).*(population(:,index_sex)==1 & population(:,index_mat)>0);
    elseif type_f==string('size') 
        % the relevant trait for competition is size
        femcompetitiveness=population(:,index_csize).*(population(:,index_sex)==1 & population(:,index_mat)>0);
    end
    % then do the same for males
    if type_m==string('budget') 
        % the relevant trait for competition is budget
        malecompetitiveness=population(:,index_budget).*(population(:,index_sex)==0 & population(:,index_mat)>0);
    elseif type_m==string('size') 
        % the relevant trait for competition is size
        malecompetitiveness=population(:,index_csize).*(population(:,index_sex)==0 & population(:,index_mat)>0);
    end

    if sum(femcompetitiveness)>0 && sum(malecompetitiveness)>0 && totalspots>0
        mothers=ddists(femcompetitiveness,totalspots);
        sires=ddists(malecompetitiveness,totalspots);
        offspringproduced=totalspots;
        index=1;
        for j=1:totalspots
            mother=population(mothers(j),:); 
            father=population(sires(j),:);
            % gamete traits initialized
            offspring=zeros(1,length(mother));
            indices=repmat([0 0 0 1 1 1],1,8);
            
            % the gamete randomly from mother's mother or mother's father
            mother_genotype_1=mother(indices==0);
            mother_genotype_2=mother(indices==1);
            mother_gamete=zeros(1,24);ind_gamete=binornd(1,0.5,[1,24]);
            mother_gamete(ind_gamete==0)=mother_genotype_1(ind_gamete==0);
            mother_gamete(ind_gamete==1)=mother_genotype_2(ind_gamete==1);

            % same for the father
            father_genotype_1=father(indices==0);
            father_genotype_2=father(indices==1);
            father_gamete=zeros(1,24);ind_gamete=binornd(1,0.5,[1,24]);
            father_gamete(ind_gamete==0)=father_genotype_1(ind_gamete==0);
            father_gamete(ind_gamete==1)=father_genotype_2(ind_gamete==1);            
            
            rep_var=init_var*10;
            % mutation steps normally distributed around 0
            mutations=normrnd(0,rep_var,[1,24]);
            % only one parent gets mutations for each locus
            whichparent=binornd(1,0.5,[1,24]);
            % ontogeny trait mutations only when evol_bounds are defined
            % otherwise only body size mutations
            if isempty(evol_bounds)
                mutations=[zeros(1,21),1,1,1].*mutations;
            end
            mother_gamete(whichparent==1)=exp(log(mother_gamete(whichparent==1))+mutations(whichparent==1));
            father_gamete(whichparent==0)=exp(log(father_gamete(whichparent==0))+mutations(whichparent==0));
            % so that the trait values do not exceed the pre-defined boundaries
            if ~isempty(evol_bounds)
                mother_gamete(mother_gamete(1:21)>evol_bounds(:,2)')=evol_bounds(mother_gamete(1:21)>evol_bounds(:,2)',2)';
                mother_gamete(mother_gamete(1:21)<evol_bounds(:,1)')=evol_bounds(mother_gamete(1:21)<evol_bounds(:,1)',1)';
                father_gamete(father_gamete(1:21)>evol_bounds(:,2)')=evol_bounds(father_gamete(1:21)>evol_bounds(:,2)',2)';
                father_gamete(father_gamete(1:21)<evol_bounds(:,1)')=evol_bounds(father_gamete(1:21)<evol_bounds(:,1)',1)';
            end
            
            offspring(indices==0)=father_gamete;
            offspring(indices==1)=mother_gamete;         

            mean_parent_val=(offspring(:,indices==0)+offspring(:,indices==1))/2;    
            sex=rand(1)>0.5;
            offspring(index_sex)=sex;
            % calculate the phenotype based on the sex-specificity of
            % expression
            % each mid-parent value is indexed as: male-female-shared
            if sex==0
                phenotype_cancer=mean_parent_val(1:3:19)*w_cancer(2)+mean_parent_val(3:3:21)*w_cancer(1);
                phenotype_bs=mean_parent_val(22)*w_bs(2)+mean_parent_val(24)*w_bs(1);
            else
                phenotype_cancer=mean_parent_val(2:3:20)*w_cancer(2)+mean_parent_val(3:3:21)*w_cancer(1);
                phenotype_bs=mean_parent_val(23)*w_bs(2)+mean_parent_val(24)*w_bs(1);
            end
            phenotype=horzcat(phenotype_cancer,phenotype_bs);
            % round the things that need to be integers (Hayflick limit and
            % differentiation levels)
            % trait 3=Hayflick
            % trait 6=dif layers
            phenotype(3)=round(phenotype(3));
            phenotype(6)=round(phenotype(6));  
            offspring(49:56)=phenotype;

            offsping_ind=emptyspots(j);
            population(offsping_ind,:)=offspring(:);
            % create the new 'individual'
            H=phenotype(3); % Hayflick limit
            K=nofonco; % how many steps until cancer
            T=phenotype(6); % differentiation layers
            N_temp=zeros([H K T]);N_temp(1,1,1)=1;
            Ns{offsping_ind}=N_temp;

            % update the population
            population(mothers(j),index_budget)=0;
            population(sires(j),index_budget)=0;
            population(mothers(j),index_offspringnum)=population(mothers(j),index_offspringnum)+1;
            population(sires(j),index_offspringnum)=population(sires(j),index_offspringnum)+1;
            parents(index,:)=[mothers(j), mother(index_sex), mother(index_age), mother(index_csize)];
            parents(index+1,:)=[sires(j), father(index_sex), father(index_age), father(index_csize)];
            index=index+2;
        end
    end
end