function tp_indi = tpInitController(tp_indi,gla_indi)
% init controller steady state
if length(gla_indi.mcidx)>1
    for i = 1:length(gla_indi.mcidx)
        tp_indi = tpSetState(tp_indi,['eta_d__',num2str(i)],tpGetState(tp_indi,['eta__',num2str(gla_indi.mcidx(i))]));
        tp_indi = tpSetState(tp_indi,['eta_f__',num2str(i)],tpGetState(tp_indi,['eta__',num2str(gla_indi.mcidx(i))]));
    end
else
    tp_indi = tpSetState(tp_indi,'eta_d',tpGetState(tp_indi,'eta__1'));
    tp_indi = tpSetState(tp_indi,'eta_f',tpGetState(tp_indi,'eta__1'));
end
V_Kb = [tpGetState(tp_indi,'u_Kb');tpGetState(tp_indi,'v_Kb');tpGetState(tp_indi,'w_Kb')];
V_K = norm(V_Kb,2);
tp_indi = tpSetState(tp_indi,'v_f',V_K);
tp_indi = tpSetState(tp_indi,'a_z_Kb_f',0);
tp_indi = tpSetState(tp_indi,'q_f',tpGetState(tp_indi,'q_Kb'));
end