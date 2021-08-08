function [obj,index]=Infill_Weighted_VB(x,search_x,kriging_model,type)

switch type
    
    case 'parallel'
        % calculate the values of x
        num_obj = length(kriging_model);
        u_g = zeros(size(x,1), num_obj);
        mse_g = zeros(size(x,1), num_obj);
        for ii = 1: num_obj
            [u_g(:, ii), mse_g(:, ii)] = predictor(x, kriging_model{ii});
             Max_error_Comp(:,ii)=bootstrap_stop_single_fidelity(search_x,kriging_model{ii},0.05);
        end
        s_g=sqrt(max(0,mse_g));
       %
        VB= Gaussian_CDF(-u_g./s_g).*(1-Gaussian_CDF(-u_g./s_g));
        
        %% identify the two situations
        index_fail=find(sum(u_g<=0,2)==num_obj);
        index_safe=find(sum(u_g<=0,2)< num_obj);
        
        if size(index_fail,1)~=0
            [obj(index_fail,:),index(index_fail,:)]=min(VB(index_fail,:),[],2);
        end
        if size(index_safe,1)~=0
            flag=u_g(index_safe,:)>0;
            [obj(index_safe,:),index(index_safe,:)]=max(flag.*VB(index_safe,:) ,[],2);
        end       
    case'series'
        % calculate the values of x
        num_obj = length(kriging_model);
        u_g = zeros(size(x,1), num_obj);
        mse_g = zeros(size(x,1), num_obj);
        for ii = 1: num_obj
            [u_g(:, ii), mse_g(:, ii)] = predictor(x, kriging_model{ii});
             Max_error_Comp(:,ii)=bootstrap_stop_single_fidelity(search_x,kriging_model{ii},0.05);
        end
        s_g=sqrt(max(0,mse_g));
        VB= Gaussian_CDF(-abs(u_g)./s_g).*(1-Gaussian_CDF(-abs(u_g)./s_g));       
        %%
        index_safe=find(sum(u_g>0,2)==num_obj);
        index_fail=find(sum(u_g>0,2)< num_obj);
        if sum(index_safe)>0
            [obj(index_safe,:),index(index_safe,:)]=max(VB(index_safe,:),[],2);
        end
        if  sum(index_fail)>0
            flag=u_g(index_fail,:)>0;
            [obj(index_fail,:),index(index_fail,:)]=min(flag.*VB(index_fail,:) ,[],2);
        end
        
        
        
end

end
