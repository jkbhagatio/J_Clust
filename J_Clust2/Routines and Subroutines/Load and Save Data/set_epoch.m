function [epochs, break_set_epoch] = set_epoch(epochs, start_time, end_time, epoch_val)

%Description: This .m file allows the user to set "epochs" within the loaded recording session. E.g. a user may want to divide the recording
%sessions (and unit clustering) into three sections: "pre-run", "run" and "post-run"
%
%Input: 'epochs' = a 2x10 array where the first row is the start time, the second row is the end time, and each column represents an epoch,
%'start_time' = user-defined start time of epoch, 'end_time' = user-defined end time of epoch, 'epoch_val' = epoch selection from GUI drop-down menu
%
%Output: 'epochs', 'break_set_epoch' = variable which allows breaking out of this function, and subsequent breaking out of callback within main GUI
%

break_set_epoch = 0;

%option to delete epoch
if epoch_val == 12
    epoch_to_del = (inputdlg('Enter Number of Epoch to Delete (1-10):', 'Delete Epoch', [1 50]));
    if isempty(epoch_to_del)
        break_set_epoch = 1;
        return;
    else
        epoch_to_del = str2num(cell2mat(epoch_to_del));
    end
    epochs(:,epoch_to_del) = 0;
    break_set_epoch = 1;
    return
end

%option to set a new epoch or return an old, previously defined epoch
if epochs(:,epoch_val) == 0
    k = questdlg('Set New Epoch as Currently Displayed Time?', 'Set Epoch');
    if strcmp(k, 'Yes')
        epochs(1, epoch_val) = start_time;
        epochs(2, epoch_val) = end_time;
    else
        break_set_epoch = 1;
        return;
    end
else
    epochs = epochs;
end

end

