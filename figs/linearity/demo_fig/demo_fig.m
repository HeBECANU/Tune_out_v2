%% graphics config 
f_units='normalized';
f_pos=[0.2,0.2,0.2,0.3];
f_pos_wide=[0.2,0.2,0.25,0.3];
f_ren='painters';
[c,cl,cd] = palette(4);     % make unique colors (normal, light, dark) tones

line_sty={'-','--',':','-.'};
mark_typ={'o','s','^','d'};
str_data={'1','2'};

%% DATA
x = linspace(0,4*pi,1e3);
y{1} = sin(x);
y{2} = cos(x);

X = linspace(0,4*pi,10);
Y{1} = sin(X);
Y{2} = cos(X);
dY = 0.1*ones(size(X));

%% VIS: LINE + DATA POINTS WITH ERROR BARS
figname='waves';
h=figure('Name',figname,'Units',f_units,'Position',f_pos,'Renderer',f_ren);

hold on;

pleg=NaN(2,1);
for ii=1:2
    % NORMAL LINE PLOT
    plot(x,y{ii},'Color',cl(ii,:));
    
    % ERROR BARS
    tp=ploterr(X,Y{ii},[],dY,'o','hhxy',0);
    set(tp(1),'Marker',mark_typ{ii},...
        'MarkerFaceColor',cl(ii,:),'MarkerEdgeColor',c(ii,:),...
        'DisplayName',str_data{ii});
    set(tp(2),'Color',c(ii,:));
    
    pleg(ii)=tp(1);     % graphical object to disp in legend
end


%%% ANNNOTATION
box on;
ax=gca;
set(ax,'Layer','Top');
ax.FontSize=fontsize;
% ax.LineWidth=ax_lwidth;       % sometimes messy

% axis tight;

xlabel('x');
ylabel('y');

lgd=legend(pleg);