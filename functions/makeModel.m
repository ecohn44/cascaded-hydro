%% ========================================================================
% GUROBI MODEL ASSEMBLY
% ========================================================================

function model = makeModel(Aeq, Aub, beq, bub, obj, lb, ub)
% Assemble a Gurobi model struct from its constraint blocks
%
% Pure: no state, no solve. Promoted out of realtimeSimulation.m so that
% realtimeDispatch.m can build a model too -- MATLAB local functions are not
% visible outside the file that defines them.
% Names are omitted: passing ~2000 name strings to Gurobi every rolling
% solve costs time and they are only used for IIS printing, which indexes
% rowNames/varNames directly.
model.A=[Aeq;Aub]; model.obj=obj; model.rhs=[beq;bub];
model.sense=[repmat('=',size(Aeq,1),1);repmat('<',size(Aub,1),1)];
model.lb=lb; model.ub=ub; model.modelsense='min';
end
