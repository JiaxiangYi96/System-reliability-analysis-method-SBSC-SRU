function  [num_vari,num_obj,num_initial_sample,mu,sigma,design_space,type,stopping_thresholds]=test_function_system_reliability(fun_name)


switch fun_name
    case'multi_modal_function'
        num_vari=2;   num_obj=3; num_initial_sample=10; mu=[1.5,2.5];sigma=[1,1];   design_space=[mu-5*sigma;mu+5*sigma]; type='parallel';
        stopping_thresholds=[0.05 0.04 0.03 0.02 0.01];
    case'parallel_three_modes'
        num_vari=2;   num_obj=3; num_initial_sample=10; mu=[0,0];sigma=[1,1];   design_space=[mu-5*sigma;mu+5*sigma]; type='parallel';
        stopping_thresholds=[0.05 0.04 0.03 0.02 0.01];
    case'AK_SYSi_parallel_three_modes'
        num_vari=2;   num_obj=3; num_initial_sample=12; mu=[0,0];sigma=[0.3,0.3];   design_space=[mu-3*sigma;mu+3*sigma]; type='parallel';
        stopping_thresholds=[0.05 0.04 0.03 0.02 0.01];
    case'four_branches_system'
        num_vari=2;   num_obj=3; num_initial_sample=10; mu=[0,0];sigma=[1,1];   design_space=[mu-5*sigma;mu+5*sigma]; type='series';
        stopping_thresholds=[0.05 0.04 0.03 0.02 0.01];
    case 'roof_truss_system'
        num_vari=8; num_obj=3;num_initial_sample=12;mu=[20000,12,9.82e-4,0.04,2e11,3e10,3.35e8,1.34e7];
        sigma=mu.*[0.07,0.01,0.06,0.12,0.06,0.06,0.12,0.18]; design_space=[mu-3*sigma;mu+3*sigma];    type='series';
        stopping_thresholds=[0.05 0.04 0.03 0.02 0.01];
        
    otherwise
        sprintf('%s function is not defined', test_function);
end


end