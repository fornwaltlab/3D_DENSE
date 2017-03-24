function res = sliceCURE(data)
    % Break the data up into slices and then show it
    res = [CURE(data(1:6,:)), CURE(data(7:12,:)), CURE(data(13:16,:))].';
end
