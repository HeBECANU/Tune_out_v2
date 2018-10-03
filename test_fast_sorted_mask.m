addpath('colormaps')

%% compare the performance
data=sort(rand(1e7,1));
min_val=0.9;
max_val=0.91;
tic; mask=data<max_val & data>min_val;
subdata1=sort(data(mask)); toc

tic;mask_idx=fast_sorted_mask(data,min_val,max_val);
subdata2=data(mask_idx(1):mask_idx(2)); toc
isequal(subdata1,subdata2)

%% show how it breaks when not sorted
data=rand(1e2,1);
min_val=0.9;
max_val=0.92;
tic; mask=data<max_val & data>min_val;
subdata1=sort(data(mask)); toc

tic;mask_idx=fast_sorted_mask(data,min_val,max_val);
subdata2=sort(data(mask_idx(1):mask_idx(2))); toc
isequal(subdata1,subdata2)


%%
points_list=round(logspace(2,8.5,100));
iimax=numel(points_list);
time_unsorted_brute=nan(1,iimax);
time_presorted_brute=nan(1,iimax); %these should be the same
time_sort=nan(1,iimax);
time_presorted_search=time_unsorted_brute;

last_update=posixtime(datetime('now'));

tic
sfigure(1);
clf
set(gcf,'color','w')
plot_colors=viridis(5);

min_val=0.89;
max_val=0.9;
m_val=1*10^(2); %compare the scaled speed if you repeat the search m times after a single sort
% min,max=[0.1,0.2]
fprintf('  \n%03u',0) 
for ii=1:iimax
fprintf('\b\b\b%03u',ii)  

unsorted_data=rand(points_list(ii),1);

tic; 
mask=unsorted_data<max_val & unsorted_data>min_val;
subdata1=unsorted_data(mask);
time_unsorted_brute(ii)=toc;
subdata1=sort(subdata1);%need to sort this for comparing later

tic;
sorted_data=sort(unsorted_data);
time_sort(ii)=toc;


tic;
mask_idx=fast_sorted_mask(sorted_data,min_val,max_val);
subdata2=sorted_data(mask_idx(1):mask_idx(2));
time_presorted_search(ii)=toc;


tic; 
mask=sorted_data<max_val & sorted_data>min_val;
subdata3=sorted_data(mask);
time_presorted_brute(ii)=toc;



if ~isequal(subdata1,subdata2,subdata3)
    error('not the same')
end
ptime=posixtime(datetime('now'));
if ptime-last_update>2 || ii==iimax
    sfigure(1);
    loglog(points_list,time_unsorted_brute,'color',plot_colors(1,:))
    hold on
    loglog(points_list,time_presorted_brute,'color',plot_colors(2,:))
    loglog(points_list,time_sort+time_presorted_search,'color',plot_colors(3,:))
    loglog(points_list,time_presorted_search,'color',plot_colors(4,:))
    loglog(points_list,(time_sort+time_presorted_search*m_val)/m_val,'color',plot_colors(5,:))
    legend('unsorted brute','presorted brute','sort then fast\_sorted\_mask','presorted fast\_sorted\_mask',...
        [sprintf('sort then 10^{%.1f}*',log10(m_val)),'fast\_sorted\_mask (scaled)'],'Location','northwest')
    hold off
    xlabel('Vector Size(n)');
    ylabel('Execution Time');
    pause(1e-6)
    last_update=ptime;

end

end
figure(1)
fprintf('\n') 

%% how many repeats?
%note this will depend on the min/max that is chosen
%from above we can see that the sort then search method
% O ~ nlog(n) +2 log(n)
%vs the brute compare
% O ~ 2n
%can never compensate for the sort time, however if we are sorting once and then doing m operations
%we can win back a speedup for large enough m
% O ~ nlog(n) + 2 m log(n)
% O ~ 2n m
%solving
%time_sort + m*time_presorted_search=m*time_unsorted_brute
%time_sort =m*time_unsorted_brute-m*time_presorted_search
%time_sort/(time_unsorted_brute-time_presorted_search) =m
figure(2)
clf
set(gcf,'color','w')
subplot(2,2,1)
loglog(points_list,time_presorted_search./time_presorted_brute,'color','k')
title('Rel. Time(Exc. Sort) fast\_sorted\_mask/brute compare')
ylabel('Relative Execution Time')
xlabel('Vector Size(n)');


subplot(2,2,2)
set(gcf,'color','w')
m_list=10.^(1:.5:4);
rel_times=arrayfun(@(m) (time_sort + m.*time_presorted_search)./ (m.*time_unsorted_brute),m_list,'UniformOutput',0);
plot_colors=viridis(numel(m_list));
if numel(rel_times)
    for ii=1:numel(rel_times)
        loglog(points_list,rel_times{ii},'color',plot_colors(ii,:))
        if ii==1, hold on ,end
    end 
end

lgd=legend(arrayfun(@(x) sprintf('m=10^{%.1f}',x),log10(m_list),'UniformOutput',0),...
    'Location','East');
lgd.Position=lgd.Position+[+0.08,0,0,0];
xlabel('Vector Size (n)');
ylabel('Relative Execution Time')
title('Rel. Time sort+m*fast\_sorted\_mask/ m*brute compare')
yl=ylim;
ylim([yl(1),2])
hold off

%now lets calculate the win factor for some values of m repeated masks 
% relative time=(time_sort + m*time_presorted_search) / m*time_unsorted_brute
subplot(2,2,3)
m_crossover=time_sort./(time_unsorted_brute-time_presorted_search);
semilogx(points_list,m_crossover,'color','k')
yl=ylim;
ylim([0,yl(2)]);
xlabel('Vector Size (n)');
ylabel('Repeated Uses (m)')
title('Repeated uses (m) for sort+fast\_sorted\_mask to win over brute compare');

subplot(2,2,4)
set(gcf,'color','w')
m_list=logspace(0,6,100);
n_desired=10^(4.0);
[~,idx]=min(abs(points_list-n_desired));
n_actual=points_list(idx);
rel_times=(time_sort(idx) + m_list.*time_presorted_search(idx))./ (m_list.*time_unsorted_brute(idx));
plot_colors=viridis(numel(m_list));
loglog(m_list,rel_times,'color','k')
xlabel('Repeated Uses (m)');
ylabel('Relative Execution Time');
title(sprintf('Rel. Time Inc. Sort for n=10^{%.1f}',log10(n_actual)))
hold off

