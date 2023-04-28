function [modes, lambdas, phases, wave_id] = separate_waves(modes, lambdas,...
                                                            phases, trunc_e)
%% Separate evanescent waves 
lambda_x = lambdas(:,1);
lambda_y = lambdas(:,2);

eL_id = find(real(log(lambda_x))+1e-2*imag(log(lambda_x))<0);
lambdas_eL = lambdas(eL_id,:);
modes_eL = modes(:,eL_id);
phases_eL = phases(eL_id,:);

eR_id = find(real(log(lambda_x))+1e-2*imag(log(lambda_x))>0);
lambdas_eR = lambdas(eR_id,:);
modes_eR = modes(:,eR_id);
phases_eR = phases(eR_id,:);

eB_id = find(real(log(lambda_y))+1e-2*imag(log(lambda_y))<0);
lambdas_eB = lambdas(eB_id,:);
modes_eB = modes(:,eB_id);
phases_eB = phases(eB_id,:);

eT_id = find(real(log(lambda_y))+1e-2*imag(log(lambda_y))>0);
lambdas_eT = lambdas(eT_id,:);
modes_eT = modes(:,eT_id);
phases_eT = phases(eT_id,:);

p_id = setdiff([1:numel(lambda_x)]',[eL_id;eR_id;eB_id;eT_id]);
lambdas_p = lambdas(p_id,:);
modes_p = modes(:,p_id);
phases_p = phases(p_id,:);

%% re-order the waves
modes = [modes_p,modes_eL,modes_eR,modes_eB,modes_eT];
lambdas = [lambdas_p;lambdas_eL;lambdas_eR;lambdas_eB;lambdas_eT];
phases = [phases_p;phases_eL;phases_eR;phases_eB;phases_eT];
wave_id.p = [1:numel(p_id)];
wave_id.eL = [1:numel(eL_id)];
wave_id.eR = [1:numel(eR_id)];
wave_id.eB = [1:numel(eB_id)];
wave_id.eT = [1:numel(eT_id)];



                                                        


