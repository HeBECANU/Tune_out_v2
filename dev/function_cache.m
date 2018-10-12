function function_cache(varargin)


%cache_opts,function_handle,var_opts
%mandatory
%cache_opts.dir
%optional
%cache_opts.no_cache


cache_opts=varargin{1};
function_handle=varargin{2};
function_opts=varargin{3:end};

hash_opt.Format = 'base64';  
hash_opt.Method = 'MD2'; 
hash_fun_args=DataHash(function_opts, hash_opt) 

fun_out=function_handle(function_opts);
log=[];
nowdt=datetime('now');
posix_time=posixtime(nowdt);

%log.iso_time=datestr(nowdt,'yyyy-mm-ddTHH:MM:SS.FFF');
                    
cache_file=fullfile(cache_opts.dir,[hash_opt.Method,'_',hash_fun_args,'.mat'])
save('import_mcp_tdc_save.mat','function_opts','fun_out','cache_opts','-v7.3');

end