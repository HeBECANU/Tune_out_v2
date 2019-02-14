% cell horzcat
function c_hcat = cell_horzcat(c)
c_hcat=cell(size(c,1),1);
for ii=1:size(c,1)
    c_hcat{ii}=horzcat(c{ii,:});
end
end