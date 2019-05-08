%test_row_col_fun_mat

a=rand(5,10);
row_col_fun_mat(@sum,a,1)

%%
row_col_fun_mat(@(x) x(1:2),a,1)

%%
row_col_fun_mat(@(x) x(1:2),a,2)
