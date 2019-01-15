function phasor_plot(sig,ccode)
    mean_freq = mean(sig(:,1));
    f_axis = sig(:,1);%-mean_freq;
    if ccode
        colourcode = linspace(0,1,size(sig,1));
    else
        colourcode = zeros(1,size(sig,1));
    end
    figure()
    subplot(2,1,1)
    scatter(f_axis,abs(sig(:,2)),[],colourcode,'.')
    colormap(plasma)
    title('Response magnitude')
    xlabel('Frequency (MHz)')
    ylabel('|Modulation (kHz)|')
    subplot(2,1,2)
    scatter(f_axis,angle(sig(:,2)),[],colourcode,'.')
    colormap(viridis)
    title('Response phase')
    xlabel('Frequency (MHz)')
    ylabel('\theta_{response}')

end