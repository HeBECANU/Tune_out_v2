
%% figure out the path to the 
mpath = fileparts(mfilename('fullpath'));


%% set up env

addpath(fullfile(mpath,'../../lib/Core_BEC_Analysis/lib/')) %add the path to set_up_project_path, this will change if Core_BEC_Analysis is included as a submodule
                  % in this case it should be 
%set_up_project_path('../..',{'+','figs'})

set_up_project_path('../..',{'+','figs','lib'})

hebec_constants %call the constants function that makes some globals

%% convert the bib file to a struct
bibs_dir=fullfile(mpath,'papers')
ref_st=bib_to_st(bibs_dir)


%% sort cronoligicaly 
% add a frac year field
for ii=1:numel(ref_st)
    ref_st{ii}.frac_year=ref_st{ii}.year+month_to_num(ref_st{ii}.month)/12;
end

frac_years=cellfun(@(x) x.frac_year,ref_st);
[~,sort_order]=sort(frac_years);
ref_st=ref_st(sort_order);

%%
title='Tune-out wavelengths for metastable helium';
st_idx=find(cellfun(@(x) contains(lower(x.title),lower(title)), ref_st));
tune_out=[];
tune_out.val=413.02;
tune_out.unc.tot=0.09;
tune_out.units='nm';
tune_out.type='theory';
ref_st{st_idx}.tune_out=tune_out;


title='Precision Measurement for Metastable Helium Atoms';
st_idx=find(cellfun(@(x) contains(lower(x.title),lower(title)), ref_st));
tune_out=[];
tune_out.val=413.0938;
tune_out.unc.stat=0.0009;
tune_out.unc.syst=0.002;
tune_out.unc.tot=rssq([tune_out.unc.stat,tune_out.unc.syst]);
tune_out.units='nm';
tune_out.type='exp';
ref_st{st_idx}.tune_out=tune_out;

title='Tune-out wavelength around 413 nm for the helium ';
st_idx=find(cellfun(@(x) contains(lower(x.title),lower(title)), ref_st));
%from page 5, non relativistic finite mass
tune_out_nrci=[];
tune_out_nrci.val=  413.13919;
tune_out_nrci.unc.tot=0.00002;
tune_out_nrci.units='nm';
tune_out_nrci.type='theory';
tune_out_nrci.note='nrfm';

tune_out_total=[];
tune_out_total.val=413.0859;
tune_out_total.unc.tot=0.0004;
tune_out_total.units='nm';
tune_out_total.type='theory';
tune_out_total.note='total';

ref_st{st_idx}.tune_out={tune_out_nrci,tune_out_total};

title='The variational calculation of the 413 nm 4He tune-out wavelength';
st_idx=find(cellfun(@(x) contains(lower(x.title),lower(title)), ref_st));
tune_out_nrci=[];
tune_out_nrci.val=413.038304399+0.100917093; %413.038304399(3)+0.100917093(7)
tune_out_nrci.unc.tot=rssq([0.000000003,0.000000007]);

tune_out_nrci.units='nm';
tune_out_nrci.type='theory';
tune_out_nrci.note='nrfm';

% value from p54
tune_out_total=[];
tune_out_total.val=  413.0858252;
tune_out_total.unc.tot=0.0000004;
tune_out_total.units='nm';
tune_out_total.type='theory';
tune_out_total.note='total';
ref_st{st_idx}.tune_out={tune_out_nrci,tune_out_total};
% 
title='Helium tune-out wavelength: Gauge';
st_idx=find(cellfun(@(x) contains(lower(x.title),lower(title)), ref_st));
tune_out=[];
tune_out.val=413.1392214630;
tune_out.unc.tot=0.00000000001;
tune_out.units='nm';
tune_out.type='theory';
ref_st{st_idx}.tune_out=tune_out;


title='QED and relativistic nuclear recoil corrections ';
st_idx=find(cellfun(@(x) contains(lower(x.title),lower(title)), ref_st));
tune_out=[];
% godamit, the arxiv paper gives a different value
tune_out.val=  413.09071; %including retardation correction
tune_out.unc.tot=0.00004;
tune_out.units='nm';
tune_out.type='theory';
ref_st{st_idx}.tune_out=tune_out;


%%
%%
%atomic unit conversion
atom_u=[];
% atom_u.energy=(const.hb)^2/(const.me*(const.a0^2)); 
% atom_u.time=const.hb/atom_u.energy;
atom_u.energy=4.3597447222071e-18 ;%J
atom_u.time=2.4188843265857e-17 ;%s
atom_u.polz=1.64877727436e-41;% C2⋅m2⋅J−1 
%from LYT response 190730
to=[];
to.scalar.wl.val=(const.h*const.c/(0.11030046*atom_u.energy));
to.scalar.wl.unc=0.00000002*(const.h*const.c/((0.11030046.^2)*atom_u.energy));
% not sure why the unc in E_to propagated to nm is different to the value in wavelength
[to.scalar.f.val,to.scalar.f.unc]=f2wl(to.scalar.wl.val,to.scalar.wl.unc);
fprintf(' to scalar %s nm\n %s MHz \n',...
    string_value_with_unc(to.scalar.wl.val*1e9,to.scalar.wl.unc*1e9,'b'),...
    string_value_with_unc(to.scalar.f.val*1e-6,to.scalar.f.unc*1e-6,'b'))

to.tosmht.wl.val=(const.h*const.c/(0.11030073*atom_u.energy));
to.tosmht.wl.unc=0.00000001*(const.h*const.c/((0.11030046.^2)*atom_u.energy));
% not sure why the unc in E_to propagated to nm is different to the value in wavelength
[to.tosmht.f.val,to.tosmht.f.unc]=f2wl(to.tosmht.wl.val,to.tosmht.wl.unc);
fprintf(' to tosmht %s nm\n %s MHz \n',...
    string_value_with_unc(to.tosmht.wl.val*1e9,to.tosmht.wl.unc*1e9,'b'),...
    string_value_with_unc(to.tosmht.f.val*1e-6,to.tosmht.f.unc*1e-6,'b'))

to.tens_shift.f.val=to.tosmht.f.val-to.scalar.f.val;
to.tens_shift.f.unc=sqrt((to.scalar.f.unc.^2)+(to.tosmht.f.unc.^2));

to.tens_shift.wl.val=to.tosmht.wl.val-to.scalar.wl.val;
to.tens_shift.wl.unc=sqrt((to.scalar.wl.unc.^2)+(to.tosmht.wl.unc.^2));

fprintf('shift %s nm\n %s MHz \n',...
    string_value_with_unc(to.tens_shift.wl.val*1e9,to.tens_shift.wl.unc*1e9,'b'),...
    string_value_with_unc(to.tens_shift.f.val*1e-6,to.tens_shift.f.unc*1e-6,'b'))

% We can also calculate a shift using the gradient values calculated in response 190730
polz=[];
polz.scalar_grad.au.val=7134.37; % polarizability gradient in a.u. polz per a.u. energy U
polz.scalar_grad.au.unc=0.01;
% lets convert into SI polz per hertz
polz.scalar_grad_si.val=polz.scalar_grad.au.val* ( (const.a0^3)*atom_u.polz  ) ... %
                    * ( const.h/atom_u.energy  );
polz.scalar_grad_si.unc= polz.scalar_grad.au.unc* ( (const.a0^3)*atom_u.polz  ) ... %
                    * ( const.h/atom_u.energy  );

polz.tensor_const.au.val=3.71e-3;
polz.tensor_const.au.unc=0.01e-3;
polz.tensor_const.si.val=polz.tensor_const.au.val*(const.a0^3)*atom_u.polz;
polz.tensor_const.si.unc=polz.tensor_const.au.unc*(const.a0^3)*atom_u.polz;

to.tens_shift_alt.f.val=(1/2)*polz.tensor_const.si.val/polz.scalar_grad_si.val;
to.tens_shift_alt.f.unc=to.tens_shift_alt.f.val*sqrt((polz.tensor_const.si.unc/polz.tensor_const.si.val)^2 ...
                        +(polz.scalar_grad_si.unc/polz.scalar_grad_si.val)^2) ;

[~,to.tens_shift_alt.wl.val]=f2wl(to.tosmht.f.val,to.tens_shift_alt.f.val);
[~,to.tens_shift_alt.wl.unc]=f2wl(to.tosmht.f.val,to.tens_shift_alt.f.unc)
fprintf('shift alt method %s nm\n %s MHz \n',...
    string_value_with_unc(to.tens_shift_alt.wl.val*1e9,to.tens_shift_alt.wl.unc*1e9,'b'),...
    string_value_with_unc(to.tens_shift_alt.f.val*1e-6,to.tens_shift_alt.f.unc*1e-6,'b'))

%%
% scale our values into the scalar TO
title='this work (exp.)';
st_idx=find(cellfun(@(x) contains(lower(x.title),lower(title)), ref_st));
tune_out=[];
tune_out.val=725736811-to.tens_shift_alt.f.val*1e-6;
tune_out.unc.stat=45;
tune_out.unc.syst=rssq([120,to.tens_shift_alt.f.unc*1e-6]);
tune_out.unc.tot=rssq([tune_out.unc.stat,tune_out.unc.syst]);
tune_out.units='MHz';
tune_out.type='exp';
ref_st{st_idx}.tune_out=tune_out;


title='this work (theory)';
st_idx=find(cellfun(@(x) contains(lower(x.title),lower(title)), ref_st));
tune_out=[];
tune_out.val=725736430-to.tens_shift_alt.f.val*1e-6;
tune_out.unc.tot=rssq([50,to.tens_shift_alt.f.unc*1e-6]);
tune_out.units='MHz';
tune_out.type='theory';
ref_st{st_idx}.tune_out=tune_out;



%%

to_st=[];
to_st.frac_year=[];
to_st.freq.val=[];
to_st.freq.unc=[];
to_st.is_th=[];
to_st.bib_keys={};
to_st.title={};
ii=1; % bibtex entry
jj=1; % tune out sub entry
kk=1; %out array index
while ii<= numel(ref_st)
    this_to=ref_st{ii}.tune_out;
    to_st.frac_year(kk)=ref_st{ii}.frac_year;
    to_st.bib_keys{kk}=ref_st{ii}.bib_key;
    to_st.ref{kk}=ref_st{ii};
    to_st.title{kk}=ref_st{ii}.title;

    if iscell(this_to)
        this_to_single=this_to{jj};
        jj=jj+1;
        if jj>numel(this_to)
            ii=ii+1;
            jj=1;
        end
    else
        this_to_single=this_to;
        ii=ii+1;
    end
    if strcmp(this_to_single.units,'nm')
        [to_st.freq.val(kk),to_st.freq.unc(kk)]=f2wl(this_to_single.val*1e-9,this_to_single.unc.tot*1e-9);
    elseif strcmp(this_to_single.units,'MHz')
        to_st.freq.val(kk)=this_to_single.val*1e6;
        to_st.freq.unc(kk)=this_to_single.unc.tot*1e6;
    end
    to_st.is_th(kk)= (strcmp(this_to_single.type,'theory'));
    
    to_st.type_str{kk}=this_to_single.type;
    if isfield(this_to_single,'note')
        to_st.note{kk}=this_to_single.note;
    else
        to_st.note{kk}='';
    end
    
    kk=kk+1;
end
th_mask=logical(to_st.is_th);
frac_year=to_st.frac_year;


%%
%colors_main=[[53,126,220];[33,188,44];[0,0,0]]./255;

colors_main=[[233,87,0];[33,188,44];[0,165,166]];
colors_main=colors_main./255;
lch=colorspace('RGB->LCH',colors_main(:,:));
lch(:,1)=lch(:,1)+20;
colors_detail=colorspace('LCH->RGB',lch);
%would prefer to use srgb_2_Jab here
color_shaded=colorspace('RGB->LCH',colors_main(3,:));
color_shaded(1)=125;
color_shaded=colorspace('LCH->RGB',color_shaded);

%%

stfig('convergence')
clf



plot_y_factor=1e-9;
%plot_y_shift=725700000*1e6;% Hz 
plot_y_shift=725735000*1e6;% Hz
to_st.freq.shift_scaled=(to_st.freq.val-plot_y_shift)*plot_y_factor;
errorbar(frac_year(th_mask),to_st.freq.shift_scaled(th_mask),to_st.freq.unc(th_mask)*plot_y_factor,...
        'o','CapSize',0,'MarkerSize',5,...
        'Color',colors_main(1,:),'MarkerFaceColor',colors_detail(1,:),'LineWidth',1.8)
    hold on
errorbar(frac_year(~th_mask),to_st.freq.shift_scaled(~th_mask),to_st.freq.unc(~th_mask)*plot_y_factor,...
        '^','CapSize',0,'MarkerSize',5,...
        'Color',colors_main(2,:),'MarkerFaceColor',colors_detail(2,:),'LineWidth',1.8) 

% placeholder so that the NR shows up on the legend
errorbar(2017,-100,1,...
        's','CapSize',0,'MarkerSize',5,...
        'Color',colors_main(3,:),'MarkerFaceColor',colors_detail(3,:),'LineWidth',1.8)    
    
    
    hold off
    
    
text_info=[];    
text_info.val=(to_st.freq.val-plot_y_shift);
text_info.unc=to_st.freq.unc;
text_info.frac_year=frac_year;
text_info.align=nan;
text_info.align(1:numel(text_info.val))=0;
text_info.ref_num=zeros(1,numel(text_info.val));

% get these numbers from the latex document by putting the cite comand at the end of the document
text_info.ref_num=[16,15,17,17,18,18,19,20,nan,nan];
text_info.delta_str=cell(1,numel(text_info.val));

text_info.ref_text=arrayfun(@(x) sprintf('[%u]',x),text_info.ref_num,'UniformOutput',0);
text_info.ref_text{end-1}='(theory)';
text_info.ref_text{end}='(exp.)';

for ii=1:numel(text_info.val)
    text_info.delta_str{ii}=sprintf('%s GHz\n %s',...
        string_value_with_unc(text_info.val(ii)*plot_y_factor,text_info.unc(ii)*plot_y_factor,'b'),...
        text_info.ref_text{ii});
end
text_info.frac_year=text_info.frac_year+0.2;
	
ii=1;
text_info.frac_year(ii)=text_info.frac_year(ii)+0.0;
text_info.val(ii)=to_st.freq.shift_scaled(end)*1e9+1e9;
% add an arrow to point up at the first data pt which is off screen
ha = annotation('arrow');  % store the arrow information in ha
ha.Parent = gca;           % associate the arrow the the current axes
ha.X = [0,0]+2014.3;          % the location in data units
ha.Y = [0,4]+3;   
ha.LineWidth  = 1;          % make the arrow bolder for the picture
ha.HeadWidth  = 8;
ha.HeadLength = 10;



ii=4;
text_info.frac_year(ii)=text_info.frac_year(ii)+0.1;
text_info.val(ii)=text_info.val(ii)+1.5*1e9;
text_info.align(ii)=1;

text_info.frac_year(end-2)=text_info.frac_year(end-2)-0.4;
text_info.val(end-2)=text_info.val(end-2)-0.5*1e9;
text_info.align(end-2)=2;

text_info.frac_year(end-1)=text_info.frac_year(end-1)-0.1;
text_info.val(end-1)=text_info.val(end-1)-3*1e9;
text_info.align(end-1)=1;

text_info.frac_year(end)=text_info.frac_year(end)-0.2;
text_info.val(end)=text_info.val(end)+3*1e9;
text_info.align(end)=1;


%plot out each seprately
mask=text_info.align==0;
text(text_info.frac_year(mask),text_info.val(mask)*1e-9,text_info.delta_str(mask),'HorizontalAlignment','left' )
mask=text_info.align==1;
text(text_info.frac_year(mask),text_info.val(mask)*1e-9,text_info.delta_str(mask),'HorizontalAlignment','center' )
mask=text_info.align==2;
text(text_info.frac_year(mask),text_info.val(mask)*1e-9,text_info.delta_str(mask),'HorizontalAlignment','right' )



    
numstr=sprintf('%.0f',plot_y_shift*1e-9)
ylabel(sprintf('$f_{\\mathrm{TO}}$ - %s (GHZ)',numstr))
ylim([-18,9]+(to_st.freq.val(end)-plot_y_shift)*1e-9)

%yticks(-10:5:20)
yticks((floor(min(to_st.freq.shift_scaled(2:end))/10)*10):5:(ceil(max(to_st.freq.shift_scaled(2:end))/10)*10))
xlim([2013,2021.5])

lg=legend('NR Theory','Full Theory','Experiment','Position',[0.75 0.4 0.1 0.1])
%legend_pos=lg.Position
%set(lg,'Position',legend_pos-[0.00000001,0,0,0])

axis_main=gca;
set(axis_main,'Units','normalized');
set(axis_main,'Color','none');
set(axis_main,'Position',[0.15,0.3,0.8,0.65])
set(axis_main,'XTickLabel','')
set(axis_main,'TickLength',[0.02,0.02])
set(axis_main,'LineWidth',1.2)
set(axis_main,'FontSize',12);

main_ax_pos=axis_main.Position;
plot_y_space=0.00;
axis_subplot=axes('Position',[main_ax_pos(1),0.1,main_ax_pos(3),main_ax_pos(2)-0.1-plot_y_space])
set(axis_subplot,'Color','none');
xticks(2013:2020)
box on

plot_y_factor=1e-9;
%plot_y_shift=725735000*1e6;% Hz
errorbar(frac_year(th_mask),(to_st.freq.val(th_mask)-plot_y_shift)*plot_y_factor,to_st.freq.unc(th_mask)*plot_y_factor,...
        's','CapSize',0,'MarkerSize',5,...
        'Color',colors_main(3,:),'MarkerFaceColor',colors_detail(3,:),'LineWidth',1.8)
    hold on
% errorbar(frac_year(~th_mask),(to_st.freq.val(~th_mask)-plot_y_shift)*plot_y_factor,to_st.freq.unc(~th_mask)*plot_y_factor,...
%         '^','CapSize',0,'MarkerSize',5,...
%         'Color',colors_main(2,:),'MarkerFaceColor',colors_detail(2,:),'LineWidth',1.8) 
%     hold off
ylim([-90.02,-89.88]+(to_st.freq.val(end)-plot_y_shift)*1e-9) 
set(axis_subplot,'Xlim',axis_main.XLim);
set(axis_subplot,'XTick',axis_main.XTick);
set(axis_subplot,'FontSize',axis_main.FontSize);
set(axis_subplot,'TickLength',[0.02,0.02])
set(axis_subplot,'LineWidth',1.2)



ylabels=arrayfun(@(x) sprintf('%.2f',x),axis_subplot.YTick,'UniformOutput',false)
set(axis_subplot,'YTickLabel',ylabels)
xlabel('Year')
axis_main.YLabel.Units='normalized';
% caluclate the y position so that the label is centered
axis_main.YLabel.Extent
axis_main.YLabel.Position=axis_main.YLabel.Position+[-0.05,-0.17,0];
pause(0.1)


ii=3;
text_info.frac_year(ii)=text_info.frac_year(ii)-0.4;
text_info.val(ii)=text_info.val(ii)-0.000*1e9;
text_info.align(ii)=2;

ii=5;
text_info.frac_year(ii)=text_info.frac_year(ii)+0.2;
text_info.val(ii)=text_info.val(ii)+0.03*1e9;
text_info.align(ii)=1;

ii=7;
text_info.frac_year(ii)=text_info.frac_year(ii)+0.6;
text_info.val(ii)=text_info.val(ii)+0.01*1e9;
text_info.align(ii)=1;


%plot out each seprately
mask=text_info.align==0;
text(text_info.frac_year(mask),text_info.val(mask)*1e-9,text_info.delta_str(mask),'HorizontalAlignment','left' )
mask=text_info.align==1;
text(text_info.frac_year(mask),text_info.val(mask)*1e-9,text_info.delta_str(mask),'HorizontalAlignment','center' )
mask=text_info.align==2;
text(text_info.frac_year(mask),text_info.val(mask)*1e-9,text_info.delta_str(mask),'HorizontalAlignment','right' )

%4/3 aspect ratio



fprintf('put this at the end of the document to get the red numbers to put back in this code\n %s\n',...
    sprintf('\\cite{%s},',to_st.bib_keys{:}))
fprintf('titles \n %s }\n',...
    sprintf('%s\n',to_st.title{:}))

th_diff=[];
th_diff.val=to_st.freq.val(end)-to_st.freq.val(end-1);
th_diff.unc=rssq(to_st.freq.unc(end-1:end));

sigma_th_exp=th_diff.val/th_diff.unc;

fprintf('diff between th and exp %s MHZ sigma %f \n',string_value_with_unc(th_diff.val*1e-6,th_diff.unc*1e-6,'b'),sigma_th_exp)

%% make a string output to make referenceing easy
iimax=numel(to_st.frac_year);
fprintf('\n')
for ii=1:iimax
    this_ref=to_st.ref{ii};
    type_note_str=to_st.type_str{ii};
    if ~isempty(to_st.note{ii})
        type_note_str=cat(2,type_note_str,' ',to_st.note{ii});
    end
    fprintf('place: %2u, year: %.1f, freq: %+7.2f GHz, val type: %s, bibkey: "%s", title: "%s", first auth.: %s\n',...
            ii,to_st.frac_year(ii),to_st.freq.shift_scaled(ii),...
             type_note_str,...
            to_st.bib_keys{ii},to_st.title{ii},...
            sprintf('%s %s', this_ref.author{1}.first, this_ref.author{1}.last))
end
fprintf('\n')


%%
set(gcf,'position',[493 318 600 450])
export_fig(gcf,fullfile(mpath,'to_val_convergence_history.svg'))
export_fig(gcf,fullfile(mpath,'to_val_convergence_history.eps'))
export_fig(gcf,fullfile(mpath,'to_val_convergence_history.png'))
%
%

%%

% stfig('convergence')
% clf
% 
% 
% 
% plot_y_factor=1e-9;
% %plot_y_shift=725700000*1e6;% Hz 
% plot_y_shift=725735000*1e6;% Hz
% to_st.freq.shift_scaled=(to_st.freq.val-plot_y_shift)*plot_y_factor;
% errorbar(frac_year(th_mask),to_st.freq.shift_scaled(th_mask),to_st.freq.unc(th_mask)*plot_y_factor,...
%         'o','CapSize',0,'MarkerSize',5,...
%         'Color',colors_main(1,:),'MarkerFaceColor',colors_detail(1,:),'LineWidth',1.8)
%     hold on
% errorbar(frac_year(~th_mask),to_st.freq.shift_scaled(~th_mask),to_st.freq.unc(~th_mask)*plot_y_factor,...
%         '^','CapSize',0,'MarkerSize',5,...
%         'Color',colors_main(2,:),'MarkerFaceColor',colors_detail(2,:),'LineWidth',1.8) 
%     hold off
%     
%     
% text_info=[];    
% text_info.val=(to_st.freq.val-plot_y_shift);
% text_info.unc=to_st.freq.unc;
% text_info.frac_year=frac_year;
% text_info.align=nan;
% text_info.align(1:numel(text_info.val))=0;
% text_info.ref_num=zeros(1,numel(text_info.val));
% text_info.ref_num=[14,13,20,20,21,21,12,11,nan,nan];
% text_info.delta_str=cell(1,numel(text_info.val));
% 
% text_info.ref_text=arrayfun(@(x) sprintf('[%u]',x),text_info.ref_num,'UniformOutput',0);
% text_info.ref_text{end-1}='(theory)';
% text_info.ref_text{end}='(exp.)';
% 
% for ii=1:numel(text_info.val)
%     text_info.delta_str{ii}=sprintf('%s GHz\n %s',...
%         string_value_with_unc(text_info.val(ii)*plot_y_factor,text_info.unc(ii)*plot_y_factor,'b'),...
%         text_info.ref_text{ii});
% end
% text_info.frac_year=text_info.frac_year+0.2;
% 	
% 
% text_info.frac_year(1)=text_info.frac_year(1)+0.3;
% text_info.val(1)=to_st.freq.shift_scaled(end)*1e9+13e9;
% 
% text_info.frac_year(end-2)=text_info.frac_year(end-2)-0.4;
% text_info.val(end-2)=text_info.val(end-2)-0.5*1e9;
% text_info.align(end-2)=2;
% 
% text_info.frac_year(end-1)=text_info.frac_year(end-1)-0.1;
% text_info.val(end-1)=text_info.val(end-1)-3.5*1e9;
% text_info.align(end-1)=1;
% 
% text_info.frac_year(end)=text_info.frac_year(end)-0.2;
% text_info.val(end)=text_info.val(end)+3.7*1e9;
% text_info.align(end)=1;
% 
% 
% %plot out each seprately
% mask=text_info.align==0;
% text(text_info.frac_year(mask),text_info.val(mask)*1e-9,text_info.delta_str(mask),'HorizontalAlignment','left' )
% mask=text_info.align==1;
% text(text_info.frac_year(mask),text_info.val(mask)*1e-9,text_info.delta_str(mask),'HorizontalAlignment','center' )
% mask=text_info.align==2;
% text(text_info.frac_year(mask),text_info.val(mask)*1e-9,text_info.delta_str(mask),'HorizontalAlignment','right' )
% 
% % add an arrow to point up at the first data pt which is off screen
% ha = annotation('arrow');  % store the arrow information in ha
% ha.Parent = gca;           % associate the arrow the the current axes
% ha.X = [0,0]+2014.2;          % the location in data units
% ha.Y = [0,4]+11;   
% ha.LineWidth  = 1;          % make the arrow bolder for the picture
% ha.HeadWidth  = 8;
% ha.HeadLength = 10;
% 
%     
% numstr=sprintf('%.0f',plot_y_shift*1e-9)
% ylabel(sprintf('$(\\omega_{\\mathrm{TO}}/2\\pi)$ - %s (GHZ)',numstr))
% ylim([-14.5,18]+(to_st.freq.val(end)-plot_y_shift)*1e-9)
% 
% %yticks(-10:5:20)
% yticks((floor(min(to_st.freq.shift_scaled(2:end))/10)*10):5:(ceil(max(to_st.freq.shift_scaled(2:end))/10)*10))
% xlim([2013,2021.5])
% 
% 
% axis_main=gca;
% set(axis_main,'Units','normalized');
% set(axis_main,'Color','none');
% set(axis_main,'Position',[0.15,0.3,0.8,0.65])
% set(axis_main,'XTickLabel','')
% set(axis_main,'TickLength',[0.02,0.02])
% set(axis_main,'LineWidth',1.2)
% set(axis_main,'FontSize',12);
% 
% main_ax_pos=axis_main.Position;
% plot_y_space=0.01;
% axis_subplot=axes('Position',[main_ax_pos(1),0.1,main_ax_pos(3),main_ax_pos(2)-0.1-plot_y_space])
% set(axis_subplot,'Color','none');
% xticks(2013:2020)
% box on
% 
% plot_y_factor=1e-9;
% %plot_y_shift=725735000*1e6;% Hz
% errorbar(frac_year(th_mask),(to_st.freq.val(th_mask)-plot_y_shift)*plot_y_factor,to_st.freq.unc(th_mask)*plot_y_factor,...
%         'o','CapSize',0,'MarkerSize',5,...
%         'Color',colors_main(1,:),'MarkerFaceColor',colors_detail(1,:),'LineWidth',1.8)
%     hold on
% errorbar(frac_year(~th_mask),(to_st.freq.val(~th_mask)-plot_y_shift)*plot_y_factor,to_st.freq.unc(~th_mask)*plot_y_factor,...
%         '^','CapSize',0,'MarkerSize',5,...
%         'Color',colors_main(2,:),'MarkerFaceColor',colors_detail(2,:),'LineWidth',1.8) 
%     hold off
% ylim([-90,-89.79]+(to_st.freq.val(end)-plot_y_shift)*1e-9) 
% set(axis_subplot,'Xlim',axis_main.XLim);
% set(axis_subplot,'XTick',axis_main.XTick);
% set(axis_subplot,'FontSize',axis_main.FontSize);
% set(axis_subplot,'TickLength',[0.02,0.02])
% set(axis_subplot,'LineWidth',1.2)
% ylabels=arrayfun(@(x) sprintf('%.2f',x),axis_subplot.YTick,'UniformOutput',false)
% set(axis_subplot,'YTickLabel',ylabels)
% xlabel('Year')
% axis_main.YLabel.Units='normalized';
% % caluclate the y position so that the label is centered
% axis_main.YLabel.Extent
% axis_main.YLabel.Position=axis_main.YLabel.Position+[-0.05,-0.17,0];
% pause(0.1)
% 
% 
% ii=3;
% text_info.frac_year(ii)=text_info.frac_year(ii)-0.5;
% text_info.val(ii)=text_info.val(ii)-0.00*1e9;
% text_info.align(ii)=2;
% 
% ii=5;
% text_info.frac_year(ii)=text_info.frac_year(ii)-0.4;
% text_info.val(ii)=text_info.val(ii)+0*1e9;
% text_info.align(ii)=2;
% 
% ii=7;
% text_info.frac_year(ii)=text_info.frac_year(ii)+0.2;
% text_info.val(ii)=text_info.val(ii)+0.02*1e9;
% text_info.align(ii)=1;
% 
% 
% %plot out each seprately
% mask=text_info.align==0;
% text(text_info.frac_year(mask),text_info.val(mask)*1e-9,text_info.delta_str(mask),'HorizontalAlignment','left' )
% mask=text_info.align==1;
% text(text_info.frac_year(mask),text_info.val(mask)*1e-9,text_info.delta_str(mask),'HorizontalAlignment','center' )
% mask=text_info.align==2;
% text(text_info.frac_year(mask),text_info.val(mask)*1e-9,text_info.delta_str(mask),'HorizontalAlignment','right' )
% 
% %4/3 aspect ratio
% set(gcf,'position',[493 318 600 450])
% 
% fprintf('for figure caption\n %s\n',...
%     sprintf('\\cite{%s},',to_st.bib_keys{:}))
% fprintf('titles \n %s }\n',...
%     sprintf('%s\n',to_st.title{:}))
% 
% th_diff=[];
% th_diff.val=to_st.freq.val(end)-to_st.freq.val(end-1);
% th_diff.unc=rssq(to_st.freq.unc(end-1:end));
% 
% sigma_th_exp=th_diff.val/th_diff.unc;
% 
% fprintf('diff between th and exp %s MHZ sigma %f \n',string_value_with_unc(th_diff.val*1e-6,th_diff.unc*1e-6,'b'),sigma_th_exp)
% 
% export_fig(gcf,fullfile(mpath,'to_val_convergence_history.svg'))
% export_fig(gcf,fullfile(mpath,'to_val_convergence_history.eps'))
% export_fig(gcf,fullfile(mpath,'to_val_convergence_history.png'))