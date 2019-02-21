function I_p_psi_chi=stokes_to_sph(stokes_in)

if iscell(stokes_in)
    I_p_psi_chi=vec_cellfun(@(x) stokes_to_sph_core(x{:}),stokes_in');
    %bit faster than doing
    %a=arrayfun(@(x) stokes_to_sph_core(x{:}),stokes_in,'UniformOutput',0);
    %out=cat(1,a{:});
else
    I_p_psi_chi=stokes_to_sph_core(stokes_in);
end
end


function I_p_psi_chi=stokes_to_sph_core(stokes_in)
%get the spehreical cord rep for the stokes vector
I_p_psi_chi=zeros(1,4);
I_p_psi_chi(1)=stokes_in(1);
I_p_psi_chi(2)=sqrt(sum(stokes_in(2:4).^2))/stokes_in(1);
I_p_psi_chi(3)=1/2*atan(stokes_in(3)/stokes_in(2));
I_p_psi_chi(4) = 1/2*atan(stokes_in(4)./sqrt(stokes_in(2).^2+stokes_in(3).^2));
end