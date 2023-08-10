function engines = enginesCreate( folder_name, structure_red, axis_reversed )
% enginesCreate create preliminary engines struct for the application of
% local forces and moments in flexible aircraft models

engines.pos_c = axis_reversed.*nastranGetEngPos( folder_name );
engines.num = size(engines.pos_c,2);
T1 = zeros(2*numel(structure_red.xyz),engines.num*3);
engines.direction = repmat([1;0;0],1,engines.num);
engines.T_es = zeros(6*engines.num,length(structure_red.modal.omega_red)+6);
for i = 1:engines.num
    [~,engines.node_idx(i)] = min(vecnorm(structure_red.xyz-engines.pos_c(:,i),2,1));
    T1(6*(engines.node_idx(i)-1)+(1:3),3*(i-1)+(1:3)) = eye(3);
    engines.T_es(6*(i-1)+(1:6),:) = structure_red.modal.T(6*(engines.node_idx(i)-1)+(1:6),:);
end
engines.T_se = structure_red.modal.T' * T1;

end