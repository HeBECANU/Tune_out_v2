%test sfigure_text

%close('all')
%sfigure_text('plot3')

%plot(sin(linspace(0,1,1e3))+rand(1,1e3))


%sfigure_text(repmat('a',1,75))
%%
tic
sfigure_text('ab',0,0);
toc

%%
tic
sfigure_text(4)
toc
%%


%test_sfigure_mock_call_fun