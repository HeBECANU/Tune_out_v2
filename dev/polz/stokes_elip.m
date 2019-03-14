function stokes_vec = stokes_elip(phi,alpha)
    % Define a pure stokes vector by angles phi and alpha
    stokes_vec= [1,cos(phi).*cos(alpha),cos(phi).*sin(alpha),sin(phi)]'; 
end