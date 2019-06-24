% cell vertcat
function C = cell_vertcat(c)
C=cell(1,size(c,2));

for ii=1:size(c,2)
    C{ii}=vertcat(c{:,ii});
end

end