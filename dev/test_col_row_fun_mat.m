%test_col_row_fun_mat
test_mat=rand(5,10);
logic_str = {'FAIL', 'pass'};

%% simple tests of suming and slecting parts of a matrix
fprintf('testing row sum\n')
col_row_fun_result=col_row_fun_mat(@sum,test_mat,2);
inbuilt_result=sum(test_mat,2);
fprintf('TEST:Results equal          : %s \n',logic_str{isequal(inbuilt_result,col_row_fun_result)+1})
%%

fprintf('testing col sum\n')
col_row_fun_result=col_row_fun_mat(@sum,test_mat,1);
inbuilt_result=sum(test_mat,1);
fprintf('TEST:Results equal          : %s \n',logic_str{isequal(inbuilt_result,col_row_fun_result)+1})

%%
fprintf('testing row subset\n')
col_row_fun_result=col_row_fun_mat(@(x) x(1:2),test_mat,2);
inbuilt_result=test_mat(:,1:2);
fprintf('TEST:Results equal          : %s \n',logic_str{isequal(inbuilt_result,col_row_fun_result)+1})
%%
fprintf('testing col subset\n')
col_row_fun_result=col_row_fun_mat(@(x) x(1:2),test_mat,1);
inbuilt_result=test_mat(1:2,:);
fprintf('TEST:Results equal          : %s \n',logic_str{isequal(inbuilt_result,col_row_fun_result)+1})

%% Test speed compared to inbuilt sum method
mat_size=1e3;
fprintf('testing speed of sum rows\n')
test_mat=rand(mat_size,mat_size+10)-0.5;
timer_handle=tic;
col_row_fun_result=col_row_fun_mat(@sum,test_mat,2);
time_row_col_fun=toc(timer_handle);

timer_handle=tic;
inbuilt_result=sum(test_mat,2);
time_inbuilt=toc(timer_handle);

fprintf('TEST:Results equal          : %s \n',logic_str{isequal(inbuilt_result,col_row_fun_result)+1})
fprintf('INFO:col_row time           : %f \n',time_row_col_fun)
fprintf('INFO:inbuilt time           : %f \n',time_inbuilt)
fprintf('INFO:rel time factor        : %f \n',time_row_col_fun/time_inbuilt)
%unsurprisingly the matlab inbuilt functions are faster


%% Test speed of for loop vs anonymous function
fprintf('testing speed compared to explicit looping\n')


core_function=@(x) norm((sin(x*3)+1).^2);
% the old way
timer_handle=tic;
% explicitly writing out the loop is pretty messy and has a lot of places you can fail
iimax=size(test_mat,1);
reslut_old=zeros(iimax,1);
for ii=1:iimax
    reslut_old(ii)=core_function(test_mat(ii,:));
end

time_old_mehod=toc(timer_handle);


timer_handle=tic;
% the new way is a single line solution that is nearly impossible to break
% and is far easier to maintain
col_row_fun_result=col_row_fun_mat(core_function,test_mat,2);
time_row_col_fun=toc(timer_handle);

fprintf('TEST:Results equal          : %s \n',logic_str{isequal(reslut_old,col_row_fun_result)+1})
fprintf('INFO:col_row time           : %f \n',time_row_col_fun)
fprintf('INFO:old school time        : %f \n',time_old_mehod)
fprintf('INFO:speedup factor        : %f \n',time_old_mehod/time_row_col_fun)

%and the speed is generaly a little bit faster with col_row_fun













