function info = load_RNN_info(filename)
% LOAD_RNN_INFO  Load and organize trained RNN data and parameters.
%
%   INFO = LOAD_RNN_INFO(FILENAME) loads a MAT-file containing a trained
%   recurrent neural network (RNN) and associated training/testing data.
%   The function reformats inputs/outputs, collects network parameters,
%   and organizes condition indices into a convenient structure.
%
%
%   INPUT
%   -----
%   FILENAME : char or string
%       Name of the .mat file containing the trained RNN variables.
%       This function assumes that the filename contains the name of the
%       family that the RNN belongs to. Family name can be 'separate' or 'continuous'
%
%   OUTPUT
%   ------
%   INFO : struct
%       Structure containing reformatted data, indices, and network
%       parameters with the following fields:
%
%       Data:
%       -----
%       Input            - Training inputs [time x features x trials]
%       Output           - Training outputs [time x outputs x trials]
%       Test_input       - Test inputs [time x features x trials]
%       Test_Outputs     - Test outputs [time x outputs x trials]
%
%       Indices:
%       --------
%       idx_pos_train    - Position indices (training)
%       idx_dir_train    - Direction indices (training)
%       idx_pos_test     - Position indices (testing)
%       idx_dir_test     - Direction indices (testing)
%       idx_cycles_train - Cycle indices (training)
%       idx_cycles_test  - Cycle indices (testing)
%
%       idx_conditions_train - [cycle, position, direction] (training)
%       idx_conditions_test  - [cycle, position, direction] (testing)
%       idx_conds_all         - Combined training and testing conditions
%
%       Network parameters:
%       -------------------
%       net_params.B   - Concatenated bias matrix
%                        [Bipos, B, Bdir, (Btask if separate)]
%       net_params.W   - Recurrent weight matrix
%       net_params.O   - Output weight matrix
%       net_params.Ob  - Output bias
%       net_params.S0  - Initial hidden state
%
%
% Andrea Colins
% 17/12/2025

if contains(filename,'separate')

    load(filename,'inputs','B','W','O','Ob','Output','Test_input','Test_Outputs','idx_cycles_train','idx_pos_train','idx_cycles_test','idx_pos_test','Bipos','Bdir','Initial_state','idx_dir_test','idx_dir_train','Btask')
else
    load(filename,'inputs','B','W','O','Ob','Output','Test_input','Test_Outputs','idx_cycles_train','idx_pos_train','idx_cycles_test','idx_pos_test','Bipos','Bdir','Initial_state','idx_dir_test','idx_dir_train')

end

info.Input = permute(inputs,[2 1 3]);
info.Output = permute(Output,[2 1 3]);
info.Test_input = permute(Test_input,[2 1 3]);
info.Test_Outputs = permute(Test_Outputs,[2 1 3]);

info.idx_pos_test=double(idx_pos_test);
info.idx_dir_test=double(idx_dir_test);
info.idx_pos_train=double(idx_pos_train);
info.idx_dir_train=double(idx_dir_train);

%% if separate RNNs, then Btask should be there

if contains(filename,'separate')
    net_params.B=[Bipos',B',Bdir',Btask'];
else
    net_params.B=[Bipos',B',Bdir'];
end

net_params.S0=Initial_state';
net_params.W=W;
net_params.O=O;
net_params.Ob=Ob;
info.net_params=net_params;

info.idx_cycles_train = idx_cycles_train;
info.idx_cycles_test = idx_cycles_test;
info.idx_conditions_train = [idx_cycles_train,info.idx_pos_train,info.idx_dir_train];
info.idx_conditions_test = [idx_cycles_test,info.idx_pos_test,info.idx_dir_test];
info.idx_conds_all = [info.idx_conditions_train;info.idx_conditions_test];

end