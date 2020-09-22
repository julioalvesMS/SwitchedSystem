function printErrorFrequencyTable(data)
%PRINTOTIMIZATIONRESULT Summary of this function goes here
%   Detailed explanation goes here
    
    aux = {
    "	\begin{table}[!t]"
    "		\centering"
    "		\caption{Average Steady-State Error}"
    "		\begin{tabular}{l c c c}"
    "			"
    "			\hline"
    "			& \multicolumn{3}{c}{Voltage Error (\%)} \\"
    "			\cline{2-4}"
    "			Controller & $1MHz$ & $200kHz$ & $40kHz$\\"
    "			\hline"
    "			QNS Controller 	& $"+sprintf('%.1f',data.c2_1M)+"$ & $"+sprintf('%.1f',data.c2_200)+"$ & $"+sprintf('%.1f',data.c2_40)+"$ \\"
    "			RNS Controller  & $"+sprintf('%.1f',data.c1_1M)+"$ & $"+sprintf('%.1f',data.c1_200)+"$ & $"+sprintf('%.1f',data.c1_40)+"$ \\"
    "			RS Controller	& $"+sprintf('%.1f',data.d1_1M)+"$ & $"+sprintf('%.1f',data.d1_200)+"$ & $"+sprintf('%.1f',data.d1_40)+"$ \\"
    "			\hline"
    "		\end{tabular}"
    "		\label{tab:steady_state_error}"
    "	\end{table}"
    ""
    };
    
    fprintf('\n\n\n======= printErrorFrequencyTable =======\n\n')
    fprintf('%s\n',aux{:});
end