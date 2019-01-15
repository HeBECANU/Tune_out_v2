function com_sig = phasor_tranform(bin_edges,data,setpoints,f_mod)
%Safely computes the complex fourier transform of signals 'data', and
%returns the complex phasor from the frequency bin closest to a specified
%diagnostic frequency 'f_mod' along with the value pair from the list
%'setpoints' . Needs bin_edges (1D array) to compute FFT.
% INPUTS
%   'data' (cell array of Nx3 TXY data)
%   'f_mod' (int)
%   'setpoints' (1D array)
%   'bin_edges'(1D array)
% OUTPUTS
%   Array of the form [setpoint,phasor] (float, complex)

    num_shots = length(data);
    data_range = [1,num_shots];
    com_sig = zeros(diff(data_range),2);
    bin_centres = bin_edges(1:end-1) + 0.5*mean(diff(bin_edges));
    ctr = 0;
    for i=data_range(1):data_range(2)
        ctr = ctr + 1;
        temp_t= data{i}(:,1);
        temp_set = setpoints(i);
        flux_t = histcounts(temp_t,bin_edges);
        [ft,c]= fft_tx(bin_centres,flux_t');
        [~,f_idx] = min(abs(ft(:,1)-f_mod));
        phasor = c(f_idx);
        com_sig(ctr,:) = [temp_set,phasor];
    end
    
end