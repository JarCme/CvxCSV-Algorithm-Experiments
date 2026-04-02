function n_table(SNR_picked, SDR_picked, params_to_get_vals, methods)

data_met_sir = SNR_picked;
data_met_sdr = SDR_picked;

methods_to_plot = length(methods);

sep = " & ";
sep_line = " \\";
hline = " \hline";
string_tab = "";

string_tab(1,1)="\begin{table}[]";
string_tab(2,1)="\resizebox{\textwidth}{!}{%";
string_tab(end+1)="{\renewcommand{\arraystretch}{1.5}";
string_tab(end+1)="\begin{tabular}{c|c||c|c||c|c||c|c||c|c||c|c||c|c||c|c|}";
string_tab(end+1) = "\cline{2-16}";
string_tab(end+1) = " & \multicolumn{1}{ c||}{\textbf{\begin{tabular}[c]{@{}c@{}}\# frames ($N$) \\ (in seconds)\end{tabular}}}";
string_tab(end+1) = "   \multicolumn{1}{|c|}{\textbf{method}} & \multicolumn{1}{c||}{\textbf{\begin{tabular}[c]{@{}c@{}}\#CSV \\ blocks \\ ($T$)\end{tabular}}} ";

maxes_sir = zeros(length(params_to_get_vals),1);
maxes_sdr = zeros(length(params_to_get_vals),1);
for N_idx = 1:length(params_to_get_vals)
    ini = params_to_get_vals(N_idx);
    if N_idx == length(params_to_get_vals)
        % string_tab(end-1) = string_tab(end-1)+ sep + "\multicolumn{2}{c|}{$"+num2str(ini)+"$}";
        string_tab(end-1) = string_tab(end-1)+ sep + "\multicolumn{2}{c|}{\begin{tabular}[c]{@{}c@{}}"+"$"+num2str(ini)+"$ \\ ($"+(ini*1*256/4)/16000+" s$)\end{tabular}}";
    else
        % string_tab(end-1) = string_tab(end-1)+ sep + "\multicolumn{2}{c||}{$"+num2str(ini)+"$}";
        string_tab(end-1) = string_tab(end-1)+ sep + "\multicolumn{2}{c||}{\begin{tabular}[c]{@{}c@{}}"+"$"+num2str(ini)+"$ \\ ($"+(ini*1*256/4)/16000+" s$)\end{tabular}}";
    end
    string_tab(end) = string_tab(end)+sep + "\textbf{SIR}" + sep + "\textbf{SDR}";
    [~,idx_max_sir] = max(data_met_sir(:,:,N_idx));
    [~,idx_max_sdr] = max(data_met_sdr(:,:,N_idx));
    maxes_sir(N_idx) = idx_max_sir;
    maxes_sdr(N_idx) = idx_max_sdr;
end
string_tab(end-1) = string_tab(end-1) + sep_line + hline + hline;
string_tab(end) = string_tab(end) + sep_line + hline + hline;

for m_idx = 1:methods_to_plot
    m = methods(m_idx);
    string_tab(end+1) = "";
    for N_idx = 1:length(params_to_get_vals)
        ini = params_to_get_vals(N_idx);
        T = "$10$";
        if N_idx == 1
            name = char(m);
            switch char(m)
                case 'ini'
                    name = "init. LCMP";
                    T = " - ";
                case 'cvxcsv'
                    name = "CvxCSV";
                case 'bogice'
                    name = "BOGICE";
                case 'bogive'
                    name = "BOGIVE";
                case 'fastdiva'
                    name = "FastDIVA";
                case 'fastdiva_T1'
                    name = "FastDIVA";
                    T = "$1$";
                case 'fastdiva_T10'
                    name = "FastDIVA";
                case 'overiva'
                    name = "OverIVA";
                    T = " - ";
            end

           
           string_tab(end) =  string_tab(end) + "\multicolumn{1}{|c|}{"+name+"}" + sep + ...
                                "\multicolumn{1}{c||}{"+ T +"}";
        end
        
        sir_cell = data_met_sir(:,m_idx,N_idx);
        sdr_cell = data_met_sdr(:,m_idx,N_idx);

        if(m_idx == maxes_sir(N_idx))
            sir_cell = "\cellcolor{blue!15}"+num2str(sir_cell, '%.1f')+"";
        else
            sir_cell = ""+num2str(sir_cell, '%.1f')+"";
        end

        if(m_idx == maxes_sdr(N_idx))
            sdr_cell = "\cellcolor{red!15}"+num2str(sdr_cell, '%.1f')+"";   
        else
            sdr_cell = num2str(sdr_cell, '%.1f');
        end

        if N_idx == length(params_to_get_vals)
            string_tab(end) =  string_tab(end) + sep + "\multicolumn{1}{S[table-format=3.2]|}{" + sir_cell + "}" + sep + "\multicolumn{1}{S[table-format=3.2]|}{" + sdr_cell + "}";
        else
            string_tab(end) =  string_tab(end) + sep + "\multicolumn{1}{S[table-format=3.2]|}{" + sir_cell + "}" + sep + "\multicolumn{1}{S[table-format=3.2]||}{" + sdr_cell + "}";
        end
        
        
    end
    string_tab(end) = string_tab(end) + sep_line + hline;
end

string_tab(end+1,1)="\end{tabular}%";
string_tab(end+1) = "}}";
string_tab(end+1,1)="\end{table}";


for l = 1:length(string_tab) disp(string_tab(l));end

end