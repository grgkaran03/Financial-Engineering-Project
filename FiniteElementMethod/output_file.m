function output_file()
    Ques1(-2, log(4), 100, 80, 1, 1);
    Ques1(-2, log(4), 100, 80, 1, 2);

    N_array = 50:10:500;
    M_array = 70*ones(length(N_array));

    calculate_errors(-2, 4, N_array, M_array, 1);
    calculate_errors(-2, 4, N_array, M_array, 2);
end