function [obj,index]=Infill_eek_criterion(search_x,kriging_model,type)
% output :
% obj:  the candidate location 
%  index : the index of lsf for updating 
switch type
    case 'parallel'
        % calculate the values of x
        num_obj = length(kriging_model);
        u_g = zeros(size(search_x,1), num_obj);
        mse_g = zeros(size(search_x,1), num_obj);
        for ii = 1: num_obj
            [u_g(:, ii), mse_g(:, ii)] = predictor(search_x, kriging_model{ii});
        end
        s_g=sqrt(max(0,mse_g));
        U_comp=abs(u_g)./s_g;
        %%
        [~,p_safe_fail,p_fail_safe]=eek_stopping_criterion(search_x,kriging_model,type);
        [Prob_safe_wrong_max,location_Prob_safe_wrong_max]=max(p_safe_fail);
        [Prob_failure_wrong_max,location_Prob_failure_wrong_max]=max(p_fail_safe);
        if Prob_safe_wrong_max>Prob_failure_wrong_max   %selection from safe point
            location_safe_predict=find(max(u_g,[],2)>=0);
            location=location_safe_predict(location_Prob_safe_wrong_max);
            U_safe_safe=U_comp(location,u_g(location,:)>=0);
            [~,location_p]=max(U_safe_safe,[],2);
            location_pp=find(u_g(location,:)>=0);
            index=location_pp(location_p);
            obj=search_x(location,:);
        else   %selection from failure point
            location_failure_predict=find(max(u_g,[],2)<0);
            location=location_failure_predict(location_Prob_failure_wrong_max);
            [~,index]=min(U_comp(location,:));
            obj=search_x(location,:);
        end
        
        
    case'series'
        % calculate the values of x
        num_obj = length(kriging_model);
        u_g = zeros(size(search_x,1), num_obj);
        mse_g = zeros(size(search_x,1), num_obj);
        for ii = 1: num_obj
            [u_g(:, ii), mse_g(:, ii)] = predictor(search_x, kriging_model{ii});
        end
        s_g=sqrt(max(0,mse_g));
        U_comp=abs(u_g)./s_g;
       %%
        
        [~,p_safe_fail,p_fail_safe]=eek_stopping_criterion(search_x,kriging_model,type);
        
        [Prob_safe_wrong_max,location_Prob_safe_wrong_max]=max(p_safe_fail);
        [Prob_failure_wrong_max,location_Prob_failure_wrong_max]=max(p_fail_safe);
        if Prob_safe_wrong_max>Prob_failure_wrong_max%selection from safe point
            location_safe_predict=find(min(u_g,[],2)>=0);
            location=location_safe_predict(location_Prob_safe_wrong_max);
            [~,index]=min(U_comp(location,:));
             obj=search_x(location,:);
        else%selection from failure point
            location_failure_predict=find(min(u_g,[],2)<0);
            location=location_failure_predict(location_Prob_failure_wrong_max);
            U_failure_failure=U_comp(location,u_g(location,:)<0);
            [~,location_p]=max(U_failure_failure,[],2);
            location_pp=find(u_g(location,:)<0);
            index=location_pp(location_p);
            obj=search_x(location,:);
        end
end

end