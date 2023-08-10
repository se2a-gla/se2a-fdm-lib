function trim_point_out = tpConvert(trim_point_in,model_description)
% tpConvert convert trim point for use in other Simulink model

state_names	= mdlGetTypeNames(model_description,'state');
input_names	= mdlGetTypeNames(model_description,'input');
output_names = mdlGetTypeNames(model_description,'output');
num_states 	= size(state_names,1);
num_inputs	= size(input_names,1);
num_outputs	= size(output_names,1);

disp(['tpConvert: Trying to assign ',num2str(num_states), ...
    ' state values, ',num2str(num_inputs),' input values', ...
    ' and ',num2str(num_outputs),' output values.']);

trim_point_out.State 	= zeros( num_states, 1 );
trim_point_out.Input	= zeros( num_inputs, 1 );
trim_point_out.Output   = zeros( num_outputs, 1 );

trim_point_out.mdl_names = mdlGetMdlNames(model_description);

for i = 1:num_states
    val = tpGetState(trim_point_in,state_names(i));
    trim_point_out = tpSetState(trim_point_out,state_names(i),val);
end

for i = 1:num_inputs
    val = tpGetInput(trim_point_in,input_names(i));
    trim_point_out = tpSetInput(trim_point_out,input_names(i),val);
end

for i = 1:num_outputs
    val = tpGetOutput(trim_point_in,output_names(i));
    trim_point_out = tpSetOutput(trim_point_out,output_names(i),val);
end

if isfield(trim_point_in, 'Task')
    trim_point_out.Task = trim_point_in.Task;
end

disp('tpConvert: Finished!');

end

