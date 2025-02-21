% Clear workspace and close all figures
clc;
clear;
close all;
% Indicator function: Selects which population group to analyze
indicador_funcion = 2;  % 1: Susceptible, 2: Exposed, 3: Infected, 4: Asymptomatic, 
                        % 5: Recovered, 6: Virus, 7: Vaccination, 8: Treatment
% Line styles for plotting
LT = {'-', '--mo', ':bs', '-.r*'};
data = [];
for ind = 1:4
    % Call the function that solves the epidemic model
    [A, Xx] = Solve(ind, indicador_funcion);
    
    % Plot results using a semi-logarithmic scale
    hold on;
    semilogy(Xx, A', LT{ind}, 'LineWidth', 1);
    
    % Store results for file output
    data(:, ind) = A';
end
title('Epidemic Dynamics')
xlabel('Time')
ylabel('Population Proportion')
grid on;
% Save results to a file
fileID = fopen('Treatment.dat', 'w');
% Write the header
fprintf(fileID, 'Time\tWithoutControl\tOnlyTreatment\tOnlyVaccination\tBothControls\n');
% Save data in tabular format
fprintf(fileID, '%.4f\t%.6f\t%.6f\t%.6f\t%.6f\n', [Xx; data']);
% Close the file
fclose(fileID);
% Define legends based on the selected population group
population_labels = { 
    'Susceptible', 'Exposed', 'Infected', 'Asymptomatic', 
    'Recovered', 'Virus', 'Vaccination', 'Treatment'
};
legend({
    [population_labels{indicador_funcion} ' without control'],
    [population_labels{indicador_funcion} ' only with treatment'],
    [population_labels{indicador_funcion} ' only with vaccination'],
    [population_labels{indicador_funcion} ' with both controls']
});
