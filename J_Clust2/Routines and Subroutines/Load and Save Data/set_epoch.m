function [epochs, break_set_epoch, N] = set_epoch(epochs, start_time, end_time, epoch_val, unit_pts, features, ts, overlaps, waveforms)

%Description: This .m file allows the user to set "epochs" within the loaded recording session. E.g. a user may want to divide the recording
%sessions (and unit clustering) into three sections: "pre-run", "run" and "post-run"
%
%Input: 'epochs' = a 2x10 array where the first row is the start time, the second row is the end time, and each column represents an epoch,
%'start_time' = user-defined start time of epoch, 'end_time' = user-defined end time of epoch, 'epoch_val' = epoch selection from GUI drop-down menu,
%'unit_pts' = points in current clusters, 'features' = cell array of spike features, 'ts' = spike timestamps, 'overlaps' = number of temporal spike 
%waveform overlaps (within 1.5 ms), 'waveforms' = spike waveforms 
%
%Output: 'epochs', 'break_set_epoch' = variable which allows breaking out of this function, and subsequent breaking out of callback within main GUI
%'N' = cell array - empty if we are trying to view a previously saved
%epoch; otherwise, within current set time, N{1} = unit_pts, N{2} =
%features, N{3} = ts, N{4} = overlaps, N{5} = waveforms.
%

break_set_epoch = 0;
N = [];

%option to delete epoch

if epoch_val == 12
    break_set_epoch = 1;
    epoch_to_del = (inputdlg('Enter Number of Epoch to Delete (1-10):', 'Delete Epoch', [1 50]));
    if isempty(epoch_to_del)
        return;
    else
        epoch_to_del = str2num(cell2mat(epoch_to_del));
    end
    if ~epochs(2,epoch_to_del)
        disp(['No data for epoch', num2str(epoch_to_del)]);
    else
        epochs(:,epoch_to_del) = 0;
        break_set_epoch = 1;
        disp(['Epoch', num2str(epoch_to_del), ' deleted']);
    end
    return
end

%option to set a new epoch (if not, we end function b/c we have old epoch in memory)

if epochs(:,epoch_val) == 0 %if we are setting a new epoch
    k = questdlg('Set New Epoch as Currently Displayed Time?', 'Set Epoch');
    if strcmp(k, 'Yes')
        epochs(1, epoch_val) = start_time;
        epochs(2, epoch_val) = end_time;
        N{1} = unit_pts;
        N{2} = features;
        N{3} = ts;
        N{4} = overlaps;
        N{5} = waveforms;
        disp(['Epoch' num2str(epoch_val),' set from: ', num2str(start_time), '-', num2str(end_time), ' seconds']);
    else
        break_set_epoch = 1;
        return;
    end
end

end

