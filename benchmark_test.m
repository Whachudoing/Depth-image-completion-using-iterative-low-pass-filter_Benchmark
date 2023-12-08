count_RMSE_2 = 0;
count_RMSE_1 = 0;
count_RMSE_0 = 1;
% % compare only matrix index >1000 elements.
criteria_m = 0;


check_RMSE_test = {};
check_MAE_test = {};
check_iRMSE_test = {};
check_iMAE_test = {};
check_Rel_test = {};

check_RMSE_test_0 = {};
check_MAE_test_0 = {};
check_iRMSE_test_0 = {};
check_iMAE_test_0 = {};
check_Rel_test_0 = {};

check_RMSE_all = {};
check_MAE_all = {};
check_iRMSE_all = {};
check_iMAE_all = {};
check_Rel_all = {};


check_RMSE ={};
check_MAE ={};
check_iRMSE ={};
check_iMAE ={};
check_Rel ={};

sum_RMSE = [];
sum_MAE = [];
sum_Rel = [];
sum_iRMSE = [];
sum_iMAE = [];
sum_iRel = [];
count_RMSE_2 = 0;
count_RMSE_1 = 0;
count_RMSE_0 = 0;
count_MAE_2 = 0;
count_MAE_1 = 0;
count_MAE_0 = 0;
count_iRMSE_2 = 0;
count_iRMSE_1 = 0;
count_iRMSE_0 = 0;
count_iMAE_2 = 0;
count_iMAE_1 = 0;
count_iMAE_0 = 0;
count_Rel_2 = 0;
count_Rel_1 = 0;
count_Rel_0 = 0;

% %  select gradient
% normal_grad = {};
% 
% for q = 1:length (I_gradient)
%     if length(I_gradient{q})>0
%         I_gradient_p = I_gradient{q};
%         normal_grad = 
%     



% 
% % % 
% for j = 1:length(RMSE_our)
% 
% %     [RMSE_our_tmp tf_rmse_our] = rmmissing(RMSE_our{1,j});
% %     [RMSE_NYU_tmp tf_rmse_NYU] = rmmissing(RMSE_NYU{1,j});
% %     [RMSE_CSPN_tmp tf_rmse_CSPN] = rmmissing(RMSE_CSPN{1,j});
% %     [RMSE_data_tmp tf_rmse_data] = rmmissing(RMSE_data{1,j});
% %     if length(tf_rmse_our)>0
% %         RMSE_our{1,j} = RMSE_our_tmp;
% %         RMSE_NYU{1,j}= RMSE_NYU{1,j}(~tf_rmse_our)
% 
%     [RMSE_our_tmp tf_rmse_our] = rmmissing(RMSE_our{1,j});
%     [RMSE_NYU_tmp tf_rmse_NYU] = rmmissing(RMSE_NYU{1,j});
%     [RMSE_CSPN_tmp tf_rmse_CSPN] = rmmissing(RMSE_CSPN{1,j});
%     [RMSE_data_tmp tf_rmse_data] = rmmissing(RMSE_data{1,j});
% 
% 
% 
% 
% % % 
m_count_tot = [];
for u = 1:length(MAE_CSPN)

    for kk = 1:length(MAE_CSPN{1, u} )
        if (m_count{1, u}{1,kk})>criteria_m 
    % %         RMSE
            if RMSE_our{1,u}{1,kk}<RMSE_CSPN{1,u}{1,kk} &&  RMSE_our{1,u}{1,kk}<RMSE_NYU{1,u}{1,kk}
                check_RMSE{u}(kk) = 2;
                check_RMSE_test = [check_RMSE_test;u kk 2 m_count{1, u}(kk) RMSE_data{1,u}{1,kk}];
                count_RMSE_2 = count_RMSE_2+1;
                sum_RMSE = [sum_RMSE;RMSE_CSPN{1,u}{1,kk}*m_count{1, u}{1,kk} RMSE_NYU{1,u}{1,kk}*m_count{1, u}{1,kk} RMSE_our{1,u}{1,kk}*m_count{1, u}{1,kk} RMSE_data{1,u}{1,kk}*m_count{1, u}{1,kk}];
                 m_count_tot = [m_count_tot;m_count{1, u}{1,kk}];
            elseif RMSE_our{1,u}{1,kk}<RMSE_CSPN{1,u}{1,kk} || RMSE_our{1,u}{1,kk}<RMSE_NYU{1,u}{1,kk}
                check_RMSE{u}(kk) = 1;
                check_RMSE_test = [check_RMSE_test;u kk 1 m_count{1, u}(kk) RMSE_data{1,u}{1,kk}];
                count_RMSE_1 = count_RMSE_1+1;
                sum_RMSE = [sum_RMSE;RMSE_CSPN{1,u}{1,kk}*m_count{1, u}{1,kk} RMSE_NYU{1,u}{1,kk}*m_count{1, u}{1,kk} RMSE_our{1,u}{1,kk}*m_count{1, u}{1,kk} RMSE_data{1,u}{1,kk}*m_count{1, u}{1,kk}];
                m_count_tot = [m_count_tot;m_count{1, u}{1,kk}];
            else
                check_RMSE{u}(kk) = 0;
%                 check_RMSE_test = [check_RMSE_test;u kk 0 m_count{1, u}(kk)];
                check_RMSE_test_0 = [check_RMSE_test_0;u kk 0 m_count{1, u}(kk) RMSE_data{1,u}{1,kk}];
                count_RMSE_0 = count_RMSE_0+1;
            end
    % %         MAE
            if MAE_our{1,u}{1,kk}<MAE_CSPN{1,u}{1,kk} &&  MAE_our{1,u}{1,kk}<MAE_NYU{1,u}{1,kk}
                check_MAE{u}(kk) = 2;
                check_MAE_test = [check_MAE_test;u kk 2 m_count{1, u}(kk) MAE_data{1,u}{1,kk}];
                count_MAE_2 = count_MAE_2+1;
                sum_MAE = [sum_MAE; MAE_CSPN{1,u}{1,kk}*m_count{1, u}{1,kk}  MAE_NYU{1,u}{1,kk}*m_count{1, u}{1,kk}  MAE_our{1,u}{1,kk}*m_count{1, u}{1,kk} ];
            elseif MAE_our{1,u}{1,kk}<MAE_CSPN{1,u}{1,kk} || MAE_our{1,u}{1,kk}<MAE_NYU{1,u}{1,kk}
                check_MAE{u}(kk) = 1;
                check_MAE_test = [check_MAE_test;u kk 1 m_count{1, u}(kk) MAE_data{1,u}{1,kk}];
                count_MAE_1 = count_MAE_1+1;
                sum_MAE = [sum_MAE; MAE_CSPN{1,u}{1,kk}*m_count{1, u}{1,kk}  MAE_NYU{1,u}{1,kk}*m_count{1, u}{1,kk}  MAE_our{1,u}{1,kk}*m_count{1, u}{1,kk} ];
            else
                check_MAE{u}(kk) = 0;
                check_MAE_test_0 = [check_MAE_test_0; u kk 0 m_count{1, u}(kk) MAE_data{1,u}{1,kk}];
                count_MAE_0 = count_MAE_0+1;
            end     
    % %         iRMSE
            if iRMSE_our{1,u}{1,kk}<iRMSE_CSPN{1,u}{1,kk} &&  iRMSE_our{1,u}{1,kk}<iRMSE_NYU{1,u}{1,kk}
                check_iRMSE{u}(kk) = 2;
                check_iRMSE_test = [check_iRMSE_test; u kk 2 m_count{1, u}(kk) iRMSE_data{1, u}(kk)];
                count_iRMSE_2 = count_iRMSE_2+1;
                sum_iRMSE = [sum_iRMSE;iRMSE_CSPN{1,u}{1,kk}*m_count{1, u}{1,kk}  iRMSE_NYU{1,u}{1,kk}*m_count{1, u}{1,kk}  iRMSE_our{1,u}{1,kk}*m_count{1, u}{1,kk} ];
            elseif iRMSE_our{1,u}{1,kk}<iRMSE_CSPN{1,u}{1,kk} || iRMSE_our{1,u}{1,kk}<iRMSE_NYU{1,u}{1,kk}
                check_iRMSE{u}(kk) = 1;
                check_iRMSE_test = [check_iRMSE_test; u kk 1 m_count{1, u}(kk) iRMSE_data{1, u}(kk)];
                count_iRMSE_1 = count_iRMSE_1+1;
                sum_iRMSE = [sum_iRMSE;iRMSE_CSPN{1,u}{1,kk}*m_count{1, u}{1,kk}  iRMSE_NYU{1,u}{1,kk}*m_count{1, u}{1,kk}  iRMSE_our{1,u}{1,kk}*m_count{1, u}{1,kk} ];
            else
                check_iRMSE{u}(kk) = 0;
                check_iRMSE_test_0 = [check_iRMSE_test_0; u kk 0 m_count{1, u}(kk) iRMSE_data{1, u}(kk)];
                count_iRMSE_0 = count_iRMSE_0+1;
            end 
    % %         iMAE
            if iMAE_our{1,u}{1,kk}<iMAE_CSPN{1,u}{1,kk} &&  iMAE_our{1,u}{1,kk}<iMAE_NYU{1,u}{1,kk}
                check_iMAE{u}(kk) = 2;
                check_iMAE_test = [check_iMAE_test; u kk 2 m_count{1, u}(kk) iMAE_data{1,u}{1,kk}];
                count_iMAE_2 = count_iMAE_2+1;
                sum_iMAE = [sum_iMAE; iMAE_CSPN{1,u}{1,kk}*m_count{1, u}{1,kk}  iMAE_NYU{1,u}{1,kk}*m_count{1, u}{1,kk}  iMAE_our{1,u}{1,kk}*m_count{1, u}{1,kk}  iMAE_data{1,u}{1,kk}*m_count{1, u}{1,kk} ];
            elseif iMAE_our{1,u}{1,kk}<iMAE_CSPN{1,u}{1,kk} || iMAE_our{1,u}{1,kk}<iMAE_NYU{1,u}{1,kk}
                check_iMAE{u}(kk) = 1;
                check_iMAE_test = [check_iMAE_test; u kk 1 m_count{1, u}(kk) iMAE_data{1,u}{1,kk}];
                count_iMAE_1 = count_iMAE_1+1;
                sum_iMAE = [sum_iMAE; iMAE_CSPN{1,u}{1,kk}*m_count{1, u}{1,kk}  iMAE_NYU{1,u}{1,kk}*m_count{1, u}{1,kk}  iMAE_our{1,u}{1,kk}*m_count{1, u}{1,kk}  iMAE_data{1,u}{1,kk}*m_count{1, u}{1,kk} ];
            else
                check_iMAE{u}(kk) = 0;
                check_iMAE_test_0 = [check_iMAE_test_0; u kk 0 m_count{1, u}(kk) iMAE_data{1,u}{1,kk}];
                count_iMAE_0 = count_iMAE_0+1;
            end         
    % %         Rel
            if Rel_our{1,u}{1,kk}<Rel_CSPN{1,u}{1,kk} &&  Rel_our{1,u}{1,kk}<Rel_NYU{1,u}{1,kk}
                check_Rel{u}(kk) = 2;
                check_Rel_test = [check_Rel_test; u kk 2 m_count{1, u}(kk) Rel_data{1,u}{1,kk}];
                count_Rel_2 = count_Rel_2+1;
                sum_Rel = [sum_Rel; Rel_CSPN{1,u}{1,kk}*m_count{1, u}{1,kk}   Rel_NYU{1,u}{1,kk}*m_count{1, u}{1,kk}   Rel_our{1,u}{1,kk} Rel_data{1,u}{1,kk}*m_count{1, u}{1,kk} ];
            elseif Rel_our{1,u}{1,kk}<Rel_CSPN{1,u}{1,kk} || Rel_our{1,u}{1,kk}<Rel_NYU{1,u}{1,kk}
                check_Rel{u}(kk) = 1;
                check_Rel_test = [check_Rel_test; u kk 1 m_count{1, u}(kk) Rel_data{1,u}{1,kk}];
                count_Rel_1 = count_Rel_1+1;
                sum_Rel = [sum_Rel; Rel_CSPN{1,u}{1,kk}*m_count{1, u}{1,kk}   Rel_NYU{1,u}{1,kk}*m_count{1, u}{1,kk}   Rel_our{1,u}{1,kk}*m_count{1, u}{1,kk}  Rel_data{1,u}{1,kk}*m_count{1, u}{1,kk} ];
            else
                check_Rel{u}(kk) = 0;
                check_Rel_test_0 = [check_Rel_test_0; u kk 0 m_count{1, u}(kk) Rel_data{1,u}{1,kk}];
                count_Rel_0 = count_Rel_0+1;
            end 
        
        end
    end

end


check_RMSE_all = [check_RMSE_test;check_RMSE_test_0];
check_MAE_all = [check_MAE_test;check_MAE_test_0];
check_iRMSE_all = [check_iRMSE_test;check_iRMSE_test_0];
check_iMAE_all = [check_iMAE_test;check_iMAE_test_0];
check_Rel_all = [check_Rel_test;check_Rel_test_0];


avg_RMSE = sum(sum_RMSE)/sum(m_count_tot);
avg_MAE = sum(sum_MAE)/sum(m_count_tot);
avg_Rel = sum(sum_Rel)/sum(m_count_tot);
avg_iRMSE = sum(sum_iRMSE)/sum(m_count_tot);
avg_iMAE = sum(sum_iMAE)/sum(m_count_tot);


