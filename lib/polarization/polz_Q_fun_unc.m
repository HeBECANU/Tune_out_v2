function [polz_q_val,polz_q_unc]=polz_Q_fun_unc(cont_val,theta_val,phi_val,...
                                cont_ci,theta_ci,phi_unc)
        polz_q_val=polz_Q_fun(cont_val,theta_val,phi_val);
        if nargout>1
            % propagate uncertianty as tolerances
            % this isnt a good way to do this but it is quick to implement
            % i can check this with numerical experiments later
            % want all possible combinations of 3 values drawn from [-1,0,1] with replacements
            signs=unique(nchoosek(repelem([-1,0,1],3), 3),'rows'); % this is a inefficeint way to get this, but all i can find with inbuilts
            ans_range=zeros(size(polz_q_val,1),size(signs,1));
            for ii=1:size(signs,1)
                if signs(ii,1)==1
                    cont_delt=cont_ci(:,2);
                elseif signs(ii,1)==-1
                    cont_delt=cont_ci(:,1);
                else 
                    cont_delt=0;
                end
                
                if signs(ii,2)==1
                    theta_delt=theta_ci(:,2);
                elseif signs(ii,2)==-1
                    theta_delt=theta_ci(:,1);
                else 
                    theta_delt=0;
                end
                
                phi_delt=signs(ii,3)*phi_unc;

                ans_range(:,ii)=polz_Q_fun(cont_val+cont_delt,theta_val+theta_delt,phi_val+phi_delt);
            end
            ans_ci=[nanmin(ans_range,[],2),nanmax(ans_range,[],2)];
            ans_ci=ans_ci-repmat(polz_q_val,[1,2]);
            %[ans_ci,polz_q_val]
            polz_q_unc=ans_ci;
        end
end