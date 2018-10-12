%hash test
%developing ideas for a nice function cache based on hashing of the input options

Opt.Format = 'base64'; 
Opt.Method = 'MD2'; 
DataHash(anal_opts, Opt) 


%%
cache_opts=[]
cache_opts.dir='.';
function_cache(cache_opts,@(x) mean(x),1:10)