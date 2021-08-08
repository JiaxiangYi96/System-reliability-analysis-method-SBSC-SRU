function [obj,index]=Infill_SYSi_RU_without_distance_threshold(x,kriging_model,type)

switch type
    
    case 'parallel'
        % calculate the values of x
        num_obj = length(kriging_model);
        u_g = zeros(size(x,1), num_obj);
        mse_g = zeros(size(x,1), num_obj);
        for ii = 1: num_obj
            [u_g(:, ii), mse_g(:, ii)] = predictor(x, kriging_model{ii});
            
        end
        s_g=sqrt(max(0,mse_g));
        %
        U_comp=abs(u_g)./(Gaussian_PDF(abs(u_g)./s_g).*s_g);
        %% identify the two situations
        index_fail=find(sum(u_g<=0,2)==num_obj);
        index_safe=find(sum(u_g<=0,2)< num_obj);
        
        if size(index_fail,1)~=0
            [obj(index_fail,:),index(index_fail,:)]=min(U_comp(index_fail,:),[],2);
        end
        if size(index_safe,1)~=0
            flag=u_g(index_safe,:)>0;
            [obj(index_safe,:),index(index_safe,:)]=max(flag.*U_comp(index_safe,:) ,[],2);
        end       
    case'series'
        % calculate the values of x
        num_obj = length(kriging_model);
        u_g = zeros(size(x,1), num_obj);
        mse_g = zeros(size(x,1), num_obj);
        for ii = 1: num_obj
            [u_g(:, ii), mse_g(:, ii)] = predictor(x, kriging_model{ii});
        end
        s_g=sqrt(max(0,mse_g));
        U_comp=abs(u_g)./(Gaussian_PDF(abs(u_g)./s_g).*s_g);        
        %%
        index_safe=find(sum(u_g>0,2)==num_obj);
        index_fail=find(sum(u_g>0,2)< num_obj);
        if sum(index_safe)>0
            [obj(index_safe,:),index(index_safe,:)]=min(U_comp(index_safe,:),[],2);
        end
        if  sum(index_fail)>0
            flag=u_g(index_fail,:)>0;
            [obj(index_fail,:),index(index_fail,:)]=max(flag.*U_comp(index_fail,:) ,[],2);
        end
        
        
        
end

end
