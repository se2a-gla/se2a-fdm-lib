function [] = baseWorkspaceRestore(base_workspace_backup)
% baseWorkspaceRestore restore base workspace from backup, see
% baseWorkspaceBackup
evalin('base','clear all');
for i = 1:length(base_workspace_backup)
    assignin('base',base_workspace_backup{1,i},base_workspace_backup{2,i});
end
end